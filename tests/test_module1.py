#!/usr/bin/env python3
"""QA tests for TaxInd module1 - Stage 1 (validate).

All tests import the repository example table
(``data/examples/asv_abundance_table.csv``) through the Stage 1 code and
check each detection step, the end-to-end run, the CLI (exit codes,
determinism) and the compatibility of the output with the downstream
alpha-diversity module (``modules/1.1-div.py``).

Run (uv provides the interpreter and dependencies):

    uv run --with pytest --with pandas --with numpy --with openpyxl pytest tests -v
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
MODULES_DIR = REPO_ROOT / "modules"
sys.path.insert(0, str(MODULES_DIR))

import module1  # noqa: E402

EXAMPLE_TABLE = REPO_ROOT / "data" / "examples" / "asv_abundance_table.csv"
MODULE1_SCRIPT = MODULES_DIR / "module1.py"
ALPHA_SCRIPT = MODULES_DIR / "1.1-div.py"


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def raw_example() -> pd.DataFrame:
    """The example table as read by the Stage 1 reader."""
    return module1.read_table(str(EXAMPLE_TABLE))


@pytest.fixture(scope="module")
def wide_example(raw_example: pd.DataFrame) -> pd.DataFrame:
    """The example table with the first column set as the ID axis."""
    return raw_example.set_index(raw_example.columns[0])


@pytest.fixture(scope="module")
def expected_counts() -> pd.DataFrame:
    """Independent ground truth: plain pandas re-read of the example table."""
    return pd.read_csv(EXAMPLE_TABLE, index_col=0)


def make_args(outdir: Path, **overrides) -> argparse.Namespace:
    """Build a `validate`-style argparse.Namespace for run_validate()."""
    options = {
        "input": str(EXAMPLE_TABLE),
        "outdir": str(outdir),
        "sep": None,
        "sheet": None,
        "assume_structure": "auto",
        "sample_col": None,
        "taxon_col": None,
        "count_col": None,
        "transpose": False,
        "unit": "auto",
        "metadata": None,
        "metadata_id_col": 0,
    }
    options.update(overrides)
    return argparse.Namespace(**options)


def run_cli(*argv: str) -> subprocess.CompletedProcess:
    """Run modules/module1.py as a subprocess."""
    return subprocess.run(
        [sys.executable, str(MODULE1_SCRIPT), *argv],
        capture_output=True,
        text=True,
    )


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def test_sniff_separator_detects_comma():
    text = EXAMPLE_TABLE.read_text(encoding="utf-8")
    assert module1.sniff_separator(text) == ","


def test_read_table_reads_example(raw_example):
    assert raw_example.shape == (10, 101)
    assert raw_example.columns[0] == "sample_id"
    # All abundance columns are numeric.
    assert all(
        pd.api.types.is_numeric_dtype(raw_example[col])
        for col in raw_example.columns[1:]
    )
    # Feature labels are nucleotide sequences (ASVs).
    labels = [str(col) for col in raw_example.columns[1:]]
    assert all(set(label) <= set("ACGT") and len(label) >= 50 for label in labels)


def test_read_header_labels_recovers_original_header():
    labels = module1.read_header_labels(str(EXAMPLE_TABLE))
    assert labels is not None
    assert len(labels) == 101
    assert labels[0] == "sample_id"
    assert len(set(labels)) == len(labels)  # no duplicates in the example


# ---------------------------------------------------------------------------
# Detection steps
# ---------------------------------------------------------------------------


def test_detect_structure_reports_wide(raw_example):
    structure, confidence, evidence = module1.detect_structure(raw_example)
    assert structure == "wide"
    assert confidence >= 0.8
    assert evidence["first_column_unique_fraction"] == 1.0
    assert evidence["numeric_fraction_of_remaining_columns"] == 1.0


def test_detect_orientation_reports_samples_rows(wide_example):
    orientation, confidence, evidence = module1.detect_orientation(wide_example)
    assert orientation == "samples_rows"
    assert confidence >= 0.8
    assert evidence["columns"]["sequence"] == 1.0  # ASV sequences as columns
    assert evidence["rows"]["sample_like"] == 1.0   # sample-like row IDs


def test_detect_data_type_reports_counts(wide_example):
    num, _ = module1.coerce_numeric(wide_example)
    data_type, evidence = module1.detect_data_type(num)
    assert data_type == "counts"
    assert evidence["values"] == "non-negative integers"
    assert evidence["min_value"] == 0.0


def test_detect_unit_reports_asv(wide_example):
    unit, evidence = module1.detect_unit(list(wide_example.columns))
    assert unit == "ASV"
    assert evidence["sequence"] == 1.0


def test_coerce_numeric_matches_ground_truth(wide_example, expected_counts):
    num, null_mask = module1.coerce_numeric(wide_example)
    assert not null_mask.to_numpy().any()       # the example has no nulls
    assert np.isfinite(num.to_numpy()).all()    # and no non-numeric cells
    np.testing.assert_allclose(
        num.to_numpy(dtype=float), expected_counts.to_numpy(dtype=float)
    )


# ---------------------------------------------------------------------------
# End-to-end run (in process)
# ---------------------------------------------------------------------------


def test_run_validate_end_to_end(tmp_path, expected_counts):
    assert module1.run_validate(make_args(tmp_path)) == 0

    for name in ("matrix_step1.csv", "validation_report.json",
                 "validation_report.md"):
        assert (tmp_path / name).is_file(), name

    report = json.loads(
        (tmp_path / "validation_report.json").read_text(encoding="utf-8")
    )

    # Detected properties on the example table.
    assert report["properties"]["structure"]["value"] == "wide"
    assert report["properties"]["structure"]["reshaped_from_long"] is False
    assert report["properties"]["orientation"]["detected"] == "samples_rows"
    assert report["properties"]["orientation"]["transposed"] is False
    assert report["properties"]["data_type"]["value"] == "counts"
    assert report["properties"]["unit"]["value"] == "ASV"

    # Matrix stats agree with the ground truth.
    sums = expected_counts.sum(axis=1)
    assert report["matrix"]["n_samples"] == 10
    assert report["matrix"]["n_features"] == 100
    assert report["matrix"]["depth_range"] == [
        float(sums.min()), float(sums.max())
    ]
    expected_sparsity = float((expected_counts == 0).sum().sum()) / expected_counts.size
    assert report["matrix"]["sparsity"] == pytest.approx(expected_sparsity)

    # The example table is already rarefied: exactly one INFO issue.
    summary = report["issue_summary"]
    assert summary["by_severity"]["ERROR"] == 0
    assert summary["by_severity"]["WARNING"] == 0
    assert summary["by_code"] == {"ALREADY_RAREFIED": 1}
    assert report["issues"][0]["code"] == "ALREADY_RAREFIED"
    assert report["issues"][0]["severity"] == "INFO"
    assert report["metadata_check"] is None

    # The report records the input checksum.
    assert report["input"]["checksum"] == module1.sha256_of_file(EXAMPLE_TABLE)


def test_matrix_step1_roundtrip(tmp_path, expected_counts):
    module1.run_validate(make_args(tmp_path))
    matrix = pd.read_csv(tmp_path / "matrix_step1.csv", index_col=0)

    assert matrix.shape == (10, 100)
    assert list(matrix.index) == list(expected_counts.index)
    assert list(matrix.columns) == list(expected_counts.columns)
    np.testing.assert_allclose(
        matrix.to_numpy(dtype=float), expected_counts.to_numpy(dtype=float)
    )
    # Counts stay integers after the normalization.
    assert np.allclose(matrix.to_numpy(), np.round(matrix.to_numpy()))


def test_json_report_schema(tmp_path):
    module1.run_validate(make_args(tmp_path))
    report = json.loads(
        (tmp_path / "validation_report.json").read_text(encoding="utf-8")
    )

    assert set(report) == {
        "module", "stage", "generated", "input", "properties",
        "matrix", "metadata_check", "issue_summary", "issues",
    }
    assert report["module"] == "module1"
    assert report["stage"] == "validate"

    for issue in report["issues"]:
        assert set(issue) == {
            "code", "severity", "axis", "target", "detail", "value", "suggestion"
        }
        assert issue["severity"] in ("ERROR", "WARNING", "INFO")
        assert issue["code"] in module1.ISSUE_SEVERITY
        assert issue["severity"] == module1.ISSUE_SEVERITY[issue["code"]]


def test_markdown_report_contents(tmp_path):
    module1.run_validate(make_args(tmp_path))
    markdown = (tmp_path / "validation_report.md").read_text(encoding="utf-8")
    assert "Validation Report (Stage 1)" in markdown
    assert "| Structure | wide |" in markdown
    assert "| Data type | counts |" in markdown
    assert "| Unit | ASV |" in markdown
    assert "| Samples | 10 |" in markdown
    assert "| Features | 100 |" in markdown
    assert "ALREADY_RAREFIED" in markdown


# ---------------------------------------------------------------------------
# CLI (subprocess)
# ---------------------------------------------------------------------------


def test_cli_success_exit_zero(tmp_path):
    result = run_cli("validate", "-i", str(EXAMPLE_TABLE), "-o", str(tmp_path))
    assert result.returncode == 0
    assert (tmp_path / "matrix_step1.csv").is_file()
    assert (tmp_path / "validation_report.json").is_file()
    assert (tmp_path / "validation_report.md").is_file()
    assert "Stage 1 (validate) summary" in result.stdout
    assert "Data type:  counts" in result.stdout


def test_cli_missing_input_exits_one(tmp_path):
    result = run_cli(
        "validate", "-i", str(tmp_path / "missing.csv"), "-o", str(tmp_path)
    )
    assert result.returncode == 1
    assert "Error" in result.stderr


def test_cli_missing_command_exits_two():
    result = run_cli()
    assert result.returncode == 2


def test_cli_output_is_deterministic(tmp_path):
    outputs = []
    for run in range(2):
        outdir = tmp_path / f"run{run}"
        result = run_cli("validate", "-i", str(EXAMPLE_TABLE), "-o", str(outdir))
        assert result.returncode == 0
        outputs.append((outdir / "matrix_step1.csv").read_bytes())
    assert outputs[0] == outputs[1]


def test_output_feeds_alpha_diversity_module(tmp_path):
    """Stage 1 output must be directly consumable by modules/1.1-div.py."""
    module1.run_validate(make_args(tmp_path))
    result = subprocess.run(
        [sys.executable, str(ALPHA_SCRIPT), "-i", str(tmp_path / "matrix_step1.csv")],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    # 1.1-div.py writes TSV with Windows line endings; keep non-empty lines.
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    assert len(lines) == 11  # header + 10 samples
    assert lines[0].split("\t") == [
        "sample", "observed_richness", "chao1", "shannon_diversity",
        "simpson_diversity", "pielou_evenness",
    ]
    assert lines[1].startswith("S01\t")
