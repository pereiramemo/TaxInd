#!/usr/bin/env python3
"""
module1.py - Stage 1 (validate): import, validation and diagnostics.

TaxInd pipeline, module 1 (preprocessing). This stage reads a raw biological
diversity matrix (ASV/OTU/taxa abundance table in CSV, TSV or Excel format),
diagnoses its structure and contents, and writes a normalized matrix plus a
validation report that feeds Stage 2 (`clean`).

Detected properties:
  - structure   : wide (samples x taxa) vs long (tidy) format
  - orientation : samples as rows (standard) or as columns (auto-transposed)
  - data type   : counts, binary (presence/absence) or relative abundance
  - unit        : ASV (sequence identifiers) vs OTU (cluster identifiers)

Validation checks (typed, severity-classified issues):
  null values, negative values, non-numeric values, duplicated row/column
  identifiers, empty samples/features, inconsistent types, already-rarefied
  matrices, low-confidence orientation and, optionally, sample-name
  mismatches against a metadata table.

Outputs (written to --outdir, default `out/`):
  matrix_step1.csv         normalized matrix (wide, samples as rows, numeric)
  validation_report.json   machine-readable report (properties + issues)
  validation_report.md     human-readable report

Assumptions:
  - The input table has a header row.
  - In wide format, the first column is the ID axis (samples or taxa).

Exit codes:
  0  stage completed (issues, if any, are listed in the report)
  1  operational failure (missing file, unreadable/empty input, ...)
  2  command line usage error (argparse default)

Usage:
  python modules/module1.py validate -i INPUT [-o OUTDIR] [options]

Stages 2 (`clean`), 3 (`rarefy`) and 4 (`run`) are planned - see
docs/plan_module1.md.
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import hashlib
import io
import json
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import openpyxl  # noqa: F401  (only needed for Excel input)
    HAS_OPENPYXL = True
except ImportError:
    HAS_OPENPYXL = False

EXCEL_SUFFIXES = {".xlsx", ".xls"}
ENCODING_FALLBACK = ("utf-8", "utf-8-sig", "latin-1")
CANDIDATE_SEPS = (",", "\t", ";")

# Strings (lowercase, after stripping) treated as null values, besides NaN.
NULL_LIKE_STRINGS = {"", "na", "n/a", "nan", "null", "none", "-"}
INT_TOLERANCE = 1e-9          # a value within this of an integer counts as one
SUM_TOLERANCE = 1e-3          # relative-abundance rows must sum to ~1
MIN_ASV_LENGTH = 50           # minimum length of a sequence-like label
MAX_ISSUES_SHOWN_PER_CODE = 20

ERROR, WARNING, INFO = "ERROR", "WARNING", "INFO"
SEVERITY_RANK = {ERROR: 0, WARNING: 1, INFO: 2}

ISSUE_SEVERITY = {
    "NULL_VALUE": WARNING,
    "NEGATIVE_VALUE": ERROR,
    "NON_NUMERIC": ERROR,
    "DUPLICATE_ROW_ID": ERROR,
    "DUPLICATE_COL_ID": ERROR,
    "EMPTY_SAMPLE": WARNING,
    "EMPTY_FEATURE": WARNING,
    "INCONSISTENT_TYPE": WARNING,
    "ALREADY_RAREFIED": INFO,
    "LONG_FORMAT_RESOLVED": INFO,
    "ORIENTATION_LOW_CONFIDENCE": WARNING,
    "SAMPLE_NOT_IN_METADATA": WARNING,
    "METADATA_ORPHAN_SAMPLE": WARNING,
    "DUPLICATE_METADATA_ID": ERROR,
}

SEQUENCE_RE = re.compile(r"^[ACGTacgt]+$")
SAMPLE_LIKE_RE = re.compile(r"^(?:s|sam|sample|specimen)[_-]?\d+$", re.IGNORECASE)
CLUSTER_ID_RE = re.compile(r"^(?:otu|asv|zotu|esv)[_-]?\d+$", re.IGNORECASE)
TAXONOMY_RANK_RE = re.compile(r"[dkpcofgs]__")


def _fail(message: str):
    """Abort with an operational error (exit code 1)."""
    sys.exit(f"Error: {message}")


def _json_safe(value: Any) -> Any:
    """Convert numpy/NaN values into JSON-serializable Python scalars."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return None if (np.isnan(number) or np.isinf(number)) else number
    return str(value)


###############################################################################
# 2. Data model
###############################################################################


@dataclass
class Issue:
    """A single validation finding.

    axis:  "matrix" | "row" | "column" | "cell" | "metadata"
    target: row/column identifiers for cell issues, a label otherwise.
    """

    code: str
    severity: str
    axis: str
    target_id: str = ""
    row: Optional[str] = None
    column: Optional[str] = None
    detail: str = ""
    value: Any = None
    suggestion: str = ""

    def to_dict(self) -> Dict[str, Any]:
        if self.row is not None or self.column is not None:
            target: Any = {"row": self.row, "column": self.column}
        else:
            target = self.target_id
        return {
            "code": self.code,
            "severity": self.severity,
            "axis": self.axis,
            "target": target,
            "detail": self.detail,
            "value": _json_safe(self.value),
            "suggestion": self.suggestion,
        }


@dataclass
class DiversityMatrix:
    """Normalized matrix (samples as rows) plus its detected properties.

    Shared data model for the module stages; Stage 1 produces it in memory
    and persists it as `matrix_step1.csv` + `validation_report.json`, which
    Stages 2-3 load back.
    """

    df: pd.DataFrame
    structure: str
    orientation: str
    transposed: bool
    data_type: str
    unit: str
    source: str
    checksum: str


class ValidationReport:
    """Collects detected properties, matrix stats and issues."""

    def __init__(self, input_path: str, checksum: str):
        self.generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
        self.input_path = input_path
        self.input_checksum = checksum
        self.properties: Dict[str, Any] = {}
        self.stats: Dict[str, Any] = {}
        self.metadata_check: Optional[Dict[str, Any]] = None
        self.issues: List[Issue] = []

    def add(self, issue: Issue) -> None:
        self.issues.append(issue)

    def issue_summary(self) -> Dict[str, Any]:
        by_severity = {severity: 0 for severity in SEVERITY_RANK}
        by_code: Dict[str, int] = {}
        for issue in self.issues:
            by_severity[issue.severity] = by_severity.get(issue.severity, 0) + 1
            by_code[issue.code] = by_code.get(issue.code, 0) + 1
        return {"by_severity": by_severity, "by_code": by_code}

    def sorted_issues(self) -> List[Issue]:
        return sorted(
            self.issues,
            key=lambda i: (
                SEVERITY_RANK.get(i.severity, 9),
                i.code,
                str(i.target_id),
                i.row or "",
                i.column or "",
            ),
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "module": "module1",
            "stage": "validate",
            "generated": self.generated,
            "input": {"path": self.input_path, "checksum": self.input_checksum},
            "properties": self.properties,
            "matrix": self.stats,
            "metadata_check": self.metadata_check,
            "issue_summary": self.issue_summary(),
            "issues": [issue.to_dict() for issue in self.sorted_issues()],
        }

    def write_json(self, path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2)
            fh.write("\n")

    def write_markdown(self, path: Path) -> None:
        def short(value: Any, limit: int = 24) -> str:
            text = str(value)
            return text if len(text) <= limit else text[: limit - 1] + "…"

        lines: List[str] = []
        add = lines.append

        add("# TaxInd module1 - Validation Report (Stage 1)")
        add("")
        add(f"- **Input:** `{self.input_path}`")
        add(f"- **Checksum:** `{self.input_checksum}`")
        add(f"- **Generated:** {self.generated}")
        add("")

        add("## Detected properties")
        add("")
        add("| Property | Value | Notes |")
        add("|---|---|---|")
        structure = self.properties.get("structure", {})
        notes = []
        if structure.get("reshaped_from_long"):
            notes.append("reshaped from long format")
        notes.append(f"confidence {structure.get('confidence')}")
        add(f"| Structure | {structure.get('value')} | {', '.join(notes)} |")
        orientation = self.properties.get("orientation", {})
        onotes = []
        if orientation.get("forced"):
            onotes.append("forced by --transpose")
        if orientation.get("transposed"):
            onotes.append("input was transposed to samples as rows")
        onotes.append(f"confidence {orientation.get('confidence')}")
        add(f"| Orientation | {orientation.get('detected')} | {', '.join(onotes)} |")
        data_type = self.properties.get("data_type", {})
        evidence = data_type.get("evidence", {})
        add(f"| Data type | {data_type.get('value')} | "
            f"{evidence.get('reason') or evidence.get('values', '')} |")
        unit = self.properties.get("unit", {})
        uevidence = unit.get("evidence", {})
        if "forced_by" in uevidence:
            unotes = str(uevidence["forced_by"])
        else:
            unotes = "; ".join(
                f"{key.replace('_', ' ')} labels: {int(round(value * 100))}%"
                for key, value in uevidence.items()
                if isinstance(value, (int, float)) and value > 0
            ) or "no label evidence"
        add(f"| Unit | {unit.get('value')} | {unotes} |")
        add("")

        add("## Matrix summary")
        add("")
        add("| Metric | Value |")
        add("|---|---|")
        add(f"| Samples | {self.stats.get('n_samples')} |")
        add(f"| Features | {self.stats.get('n_features')} |")
        depth_range = self.stats.get("depth_range")
        if depth_range:
            add(f"| Abundance range per sample | {depth_range[0]:g} - {depth_range[1]:g} |")
        sparsity = self.stats.get("sparsity")
        if sparsity is not None:
            add(f"| Sparsity (zero cells) | {sparsity:.1%} |")
        add("")

        if self.metadata_check:
            check = self.metadata_check
            add("## Metadata cross-check")
            add("")
            add(f"- Metadata file: `{check.get('metadata_file')}` "
                 f"(ID column: `{check.get('id_column')}`)")
            add(f"- Matched samples: {check.get('n_matched')} / {self.stats.get('n_samples')}")
            add(f"- Samples missing from metadata: {check.get('n_missing_from_metadata')}")
            add(f"- Metadata samples not in matrix: {check.get('n_orphan_metadata_samples')}")
            add("")

        add("## Issues")
        add("")
        if not self.issues:
            add("No issues detected - the matrix is ready for Stage 2 (`clean`).")
        else:
            summary = self.issue_summary()
            counts = summary["by_severity"]
            add(f"**{counts.get(ERROR, 0)} error(s), "
                f"{counts.get(WARNING, 0)} warning(s), "
                f"{counts.get(INFO, 0)} info**")
            add("")
            add("| Severity | Code | Count |")
            add("|---|---|---|")
            by_code = summary["by_code"]
            for code in sorted(by_code, key=lambda c: (SEVERITY_RANK[ISSUE_SEVERITY[c]], c)):
                add(f"| {ISSUE_SEVERITY[code]} | `{code}` | {by_code[code]} |")
            add("")
            add("### Details")
            add("")
            grouped: Dict[str, List[Issue]] = {}
            for issue in self.sorted_issues():
                grouped.setdefault(issue.code, []).append(issue)
            for code, issues in grouped.items():
                add(f"#### `{code}`")
                add("")
                for issue in issues[:MAX_ISSUES_SHOWN_PER_CODE]:
                    if issue.row is not None or issue.column is not None:
                        target = f"row `{short(issue.row)}`, column `{short(issue.column)}`"
                    else:
                        target = f"`{short(issue.target_id)}`"
                    value = _json_safe(issue.value)
                    value_note = f" (value: `{short(value)}`)" if value is not None else ""
                    suggestion = f" - *{issue.suggestion}*" if issue.suggestion else ""
                    add(f"- **[{issue.severity}]** {target}{value_note}: "
                        f"{issue.detail}{suggestion}")
                hidden = len(issues) - MAX_ISSUES_SHOWN_PER_CODE
                if hidden > 0:
                    add(f"- ... and {hidden} more `{code}` issue(s)")
                add("")

        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")


###############################################################################
# 3. Input reading
###############################################################################


def sha256_of_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return f"sha256:{digest.hexdigest()}"


def read_text_with_fallback(path: Path) -> str:
    """Read a text file trying UTF-8 (BOM-aware) and Latin-1 as fallback."""
    raw = path.read_bytes()
    if raw.startswith(b"\xff\xfe") or raw.startswith(b"\xfe\xff"):
        return raw.decode("utf-16")
    for encoding in ENCODING_FALLBACK:
        try:
            text = raw.decode(encoding)
        except UnicodeDecodeError:
            continue
        if encoding == "utf-8" and text.startswith("\ufeff"):
            continue  # BOM: retry with utf-8-sig
        return text
    return raw.decode("latin-1")  # cannot fail


def sniff_separator(text: str) -> Optional[str]:
    """Guess the separator from the first non-empty line."""
    first_line = ""
    for line in text.splitlines():
        if line.strip():
            first_line = line
            break
    if not first_line:
        return None
    counts = {sep: first_line.count(sep) for sep in CANDIDATE_SEPS}
    best = max(counts, key=lambda sep: counts[sep])
    return best if counts[best] > 0 else None


def sep_from_suffix(suffix: str) -> str:
    return "\t" if suffix in (".tsv", ".txt") else ","


def resolve_sheet(sheet: Optional[str]):
    """Accept a sheet name or a 0-based index string (None -> first sheet)."""
    if sheet is None:
        return 0
    text = str(sheet).strip()
    if re.fullmatch(r"\d+", text):
        return int(text)
    return text


def read_table(path_str: str, sep: Optional[str] = None,
               sheet: Optional[str] = None) -> pd.DataFrame:
    """Read a CSV, TSV or Excel table into a DataFrame."""
    path = Path(path_str)
    if not path.is_file():
        _fail(f"file not found: {path_str}")
    suffix = path.suffix.lower()
    try:
        if suffix in EXCEL_SUFFIXES:
            if not HAS_OPENPYXL:
                _fail(f"reading Excel files requires the 'openpyxl' package: {path_str}")
            return pd.read_excel(path, sheet_name=resolve_sheet(sheet))
        text = read_text_with_fallback(path)
        effective_sep = sep or sniff_separator(text) or sep_from_suffix(suffix)
        return pd.read_csv(io.StringIO(text), sep=effective_sep)
    except Exception as exc:
        _fail(f"could not read table '{path_str}': {exc}")


def read_header_labels(path_str: str, sep: Optional[str] = None,
                       sheet: Optional[str] = None) -> Optional[List[str]]:
    """Best-effort recovery of the original (unmangled) header labels.

    pandas silently renames duplicated column names on read (e.g. a second
    'ASV-1' becomes 'ASV-1.1'); this reads the header row as data so the
    original labels can be checked for duplicates.
    """
    try:
        path = Path(path_str)
        if not path.is_file():
            return None
        suffix = path.suffix.lower()
        if suffix in EXCEL_SUFFIXES:
            if not HAS_OPENPYXL:
                return None
            raw = pd.read_excel(path, sheet_name=resolve_sheet(sheet), header=None)
            if raw.empty:
                return None
            return [str(label) for label in raw.iloc[0].tolist()]
        text = read_text_with_fallback(path)
        effective_sep = sep or sniff_separator(text) or sep_from_suffix(suffix)
        raw = pd.read_csv(io.StringIO(text), sep=effective_sep,
                          header=None, nrows=1, dtype=str)
        if raw.empty:
            return None
        return [str(label) for label in raw.iloc[0].tolist()]
    except Exception:
        return None


###############################################################################
# 4. Structure detection (wide vs long)
###############################################################################


def _is_numeric_series(series: pd.Series) -> bool:
    """True if the column is numeric, or numeric for >= 90% of its values."""
    if pd.api.types.is_numeric_dtype(series):
        return True
    filled = series.notna()
    if not filled.any():
        return False
    coerced = pd.to_numeric(series, errors="coerce")
    return coerced.notna().sum() >= 0.9 * filled.sum()


def detect_structure(df: pd.DataFrame) -> Tuple[str, float, Dict[str, Any]]:
    """Detect whether the table is in wide or long format.

    Returns (structure, confidence, evidence).
    """
    nrows, ncols = df.shape
    first = df.iloc[:, 0]
    rest = df.iloc[:, 1:]
    first_is_numeric = _is_numeric_series(first)
    rest_numeric_frac = (
        sum(1 for col in rest.columns if _is_numeric_series(rest[col])) / rest.shape[1]
        if rest.shape[1]
        else 0.0
    )
    first_unique_frac = (first.nunique() / nrows) if nrows else 1.0
    evidence = {
        "n_rows": nrows,
        "n_columns": ncols,
        "first_column_numeric": bool(first_is_numeric),
        "numeric_fraction_of_remaining_columns": round(rest_numeric_frac, 2),
        "first_column_unique_fraction": round(first_unique_frac, 2),
    }

    # Long: few columns, two or more identifier columns with repeated values
    # plus a numeric count column.
    if 3 <= ncols <= 4 and not first_is_numeric:
        numeric_cols = [c for c in df.columns if _is_numeric_series(df[c])]
        categorical = [c for c in df.columns if c not in numeric_cols]
        has_repeats = any(df[c].nunique() < nrows for c in categorical)
        if numeric_cols and len(categorical) >= 2 and has_repeats:
            classic = (
                ncols == 3
                and len(numeric_cols) == 1
                and _is_numeric_series(df.iloc[:, -1])
            )
            return "long", (0.9 if classic else 0.75), {
                **evidence,
                "numeric_columns": [str(c) for c in numeric_cols],
                "identifier_columns": [str(c) for c in categorical],
            }

    # Wide: ID-like first column plus mostly numeric remaining columns.
    if (not first_is_numeric and rest.shape[1] >= 1
            and rest_numeric_frac >= 0.8 and first_unique_frac >= 0.8):
        return "wide", min(0.99, 0.6 + 0.35 * rest_numeric_frac), evidence

    # Fallback: assume wide with low confidence.
    return "wide", (0.6 if ncols >= 5 else 0.3), evidence


def _resolve_column(df: pd.DataFrame, ref: str):
    """Resolve a column by name or 0-based position string."""
    if ref in df.columns:
        return ref
    text = str(ref)
    if re.fullmatch(r"-?\d+", text):
        index = int(text)
        if 0 <= index < df.shape[1]:
            return df.columns[index]
    _fail(f"column not found: {ref}")


def reshape_long_to_wide(df: pd.DataFrame, sample_col: Optional[str] = None,
                         taxon_col: Optional[str] = None,
                         count_col: Optional[str] = None
                         ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Reshape a long (tidy) table to wide (samples x taxa).

    Duplicate (sample, taxon) pairs keep their FIRST value (no automatic
    correction - the user decides in Stage 2). Rows with missing identifiers
    are excluded and reported.
    """
    numeric_cols = [c for c in df.columns if _is_numeric_series(df[c])]
    categorical = [c for c in df.columns if c not in numeric_cols]

    count_ref = _resolve_column(df, count_col) if count_col is not None else None
    if count_ref is None:
        if not numeric_cols:
            _fail("long format detected but no numeric count column was found")
        count_ref = numeric_cols[-1]
    sample_ref = _resolve_column(df, sample_col) if sample_col is not None else None
    if sample_ref is None:
        if not categorical:
            _fail("long format detected but no identifier column for samples was found")
        sample_ref = categorical[0]
    taxon_ref = _resolve_column(df, taxon_col) if taxon_col is not None else None
    if taxon_ref is None:
        remaining = [c for c in categorical if c != sample_ref]
        if not remaining:
            _fail("long format detected but a second identifier column "
                  "(taxa) was not found")
        taxon_ref = remaining[0]

    subset = df[[sample_ref, taxon_ref, count_ref]]
    null_keys = subset[sample_ref].isna() | subset[taxon_ref].isna()
    n_null_keys = int(null_keys.sum())
    usable = subset.loc[~null_keys]
    n_duplicate_rows = int(
        usable.duplicated(subset=[sample_ref, taxon_ref], keep=False).sum()
    )

    grouped = (
        usable.groupby([sample_ref, taxon_ref], sort=False)[count_ref]
        .agg(lambda values: values.iloc[0])
    )
    wide = grouped.unstack(fill_value=0)

    info = {
        "sample_column": str(sample_ref),
        "taxon_column": str(taxon_ref),
        "count_column": str(count_ref),
        "rows_belonging_to_duplicate_pairs": n_duplicate_rows,
        "rows_with_missing_identifiers": n_null_keys,
    }
    return wide, info


###############################################################################
# 5. Orientation detection (samples vs taxa)
###############################################################################


def _looks_taxonomic(label: str) -> bool:
    label = label.strip()
    if re.match(r"^[dkpcofgs]__", label):
        return True
    return ";" in label and bool(TAXONOMY_RANK_RE.search(label))


def _label_evidence(labels: Sequence[Any]) -> Dict[str, float]:
    strings = [str(label) for label in labels]
    n = len(strings) or 1
    return {
        "sequence": round(sum(
            1 for x in strings
            if SEQUENCE_RE.match(x) and len(x) >= MIN_ASV_LENGTH
        ) / n, 2),
        "cluster_id": round(sum(
            1 for x in strings if CLUSTER_ID_RE.match(x.strip())
        ) / n, 2),
        "taxonomy": round(sum(
            1 for x in strings if _looks_taxonomic(x)
        ) / n, 2),
        "sample_like": round(sum(
            1 for x in strings if SAMPLE_LIKE_RE.match(x.strip())
        ) / n, 2),
    }


def detect_orientation(df: pd.DataFrame) -> Tuple[str, float, Dict[str, Any]]:
    """Decide which axis holds the samples.

    Returns (orientation, confidence, evidence) where orientation is
    "samples_rows" or "samples_cols" (taxa as rows).
    """
    row_evidence = _label_evidence(list(df.index))
    col_evidence = _label_evidence(list(df.columns))
    row_taxa_score = (
        2.0 * row_evidence["sequence"]
        + row_evidence["cluster_id"] + row_evidence["taxonomy"]
    )
    col_taxa_score = (
        2.0 * col_evidence["sequence"]
        + col_evidence["cluster_id"] + col_evidence["taxonomy"]
    )

    if abs(row_taxa_score - col_taxa_score) >= 0.15:
        # Stronger taxa evidence wins: taxa in rows -> samples are columns.
        orientation = "samples_cols" if row_taxa_score > col_taxa_score else "samples_rows"
        confidence = min(0.99, 0.55 + 0.4 * abs(row_taxa_score - col_taxa_score))
    else:
        # No (or conflicting) label evidence: the taxa axis usually has
        # more labels than the sample axis.
        orientation = "samples_cols" if df.shape[0] > df.shape[1] else "samples_rows"
        confidence = 0.35
    return orientation, confidence, {
        "rows": row_evidence, "columns": col_evidence,
    }


###############################################################################
# 6. Type and unit detection
###############################################################################


def coerce_numeric(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Coerce the matrix to numeric values.

    Returns (numeric_df, null_mask):
      - null_mask marks empty cells, NaN and null-like strings ("NA", "-", ...).
      - numeric_df holds floats; cells that are neither null-like nor numeric
        become NaN (reported as NON_NUMERIC by the checks).

    Column access is positional so that duplicated column labels (reported
    as DUPLICATE_COL_ID) do not break the coercion.
    """
    n_rows, n_cols = df.shape
    numeric_arr = np.full((n_rows, n_cols), np.nan)
    null_arr = np.zeros((n_rows, n_cols), dtype=bool)
    for pos in range(n_cols):
        series = df.iloc[:, pos]
        if pd.api.types.is_numeric_dtype(series):
            null_arr[:, pos] = series.isna().to_numpy()
            numeric_arr[:, pos] = series.astype(float).to_numpy()
        else:
            lowered = series.astype(str).str.strip().str.lower()
            is_null = (series.isna() | lowered.isin(NULL_LIKE_STRINGS)).to_numpy()
            null_arr[:, pos] = is_null
            converted = pd.to_numeric(series.where(~is_null), errors="coerce")
            numeric_arr[:, pos] = converted.to_numpy(dtype=float)
    numeric = pd.DataFrame(numeric_arr, index=df.index, columns=df.columns)
    null_mask = pd.DataFrame(null_arr, index=df.index, columns=df.columns)
    return numeric, null_mask


def detect_data_type(num: pd.DataFrame) -> Tuple[str, Dict[str, Any]]:
    """Classify the numeric matrix as counts / binary / relative."""
    values = num.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    evidence: Dict[str, Any] = {"n_numeric_cells": int(finite.size)}
    if finite.size == 0:
        return "inconsistent", {**evidence, "reason": "no numeric values in the matrix"}

    evidence["min_value"] = float(finite.min())
    evidence["max_value"] = float(finite.max())

    if np.all((finite == 0) | (finite == 1)):
        return "binary", {**evidence, "values": "all cells are 0 or 1"}

    integer_like = bool(np.all(np.abs(finite - np.round(finite)) <= INT_TOLERANCE))
    non_negative = bool(np.all(finite >= 0))
    if integer_like and non_negative:
        return "counts", {**evidence, "values": "non-negative integers"}

    if non_negative and bool(np.all(finite <= 1)) and not integer_like:
        row_sums = np.nansum(values, axis=1)
        col_sums = np.nansum(values, axis=0)
        rows_ok = float(np.mean(np.abs(row_sums - 1.0) <= SUM_TOLERANCE))
        cols_ok = float(np.mean(np.abs(col_sums - 1.0) <= SUM_TOLERANCE))
        evidence["row_sums_equal_one"] = round(rows_ok, 3)
        evidence["column_sums_equal_one"] = round(cols_ok, 3)
        if rows_ok >= 0.9 or cols_ok >= 0.9:
            return "relative", {
                **evidence,
                "values": "proportions in [0, 1]; rows (or columns) sum to 1",
            }
        return "inconsistent", {
            **evidence,
            "reason": "values in [0, 1] but rows/columns do not sum to 1",
        }

    reasons = []
    if not non_negative:
        reasons.append(f"{int((finite < 0).sum())} negative value(s)")
    if not integer_like:
        reasons.append(
            f"{int((np.abs(finite - np.round(finite)) > INT_TOLERANCE).sum())} "
            "non-integer value(s)"
        )
    return "inconsistent", {**evidence, "reason": "; ".join(reasons) or "mixed values"}


def detect_unit(labels: Sequence[Any]) -> Tuple[str, Dict[str, Any]]:
    """Detect whether feature labels are ASVs (sequences) or OTUs (IDs)."""
    evidence = _label_evidence(labels)
    if evidence["sequence"] >= 0.5:
        return "ASV", evidence
    if evidence["cluster_id"] + evidence["taxonomy"] >= 0.5:
        return "OTU", evidence
    return "unknown", evidence


###############################################################################
# 7. Validation checks
###############################################################################


def check_cell_issues(raw_matrix: pd.DataFrame, num: pd.DataFrame,
                      null_mask: pd.DataFrame, report: ValidationReport) -> None:
    """Cell-level checks: null, non-numeric and negative values."""
    index = raw_matrix.index
    columns = raw_matrix.columns

    nulls = null_mask.to_numpy()
    rows, cols = np.nonzero(nulls)
    for r, c in zip(rows, cols):
        report.add(Issue("NULL_VALUE", WARNING, "cell", f"{index[r]} | {columns[c]}",
                         row=str(index[r]), column=str(columns[c]),
                         value=raw_matrix.iat[r, c],
                         detail="empty or null-like value",
                         suggestion="set a value (e.g. 0) or drop the row/column"))

    non_numeric = num.isna().to_numpy() & ~nulls
    rows, cols = np.nonzero(non_numeric)
    for r, c in zip(rows, cols):
        report.add(Issue("NON_NUMERIC", ERROR, "cell", f"{index[r]} | {columns[c]}",
                         row=str(index[r]), column=str(columns[c]),
                         value=raw_matrix.iat[r, c],
                         detail="value is not numeric",
                         suggestion="fix the value or drop the row/column"))

    negatives = num.to_numpy() < 0
    rows, cols = np.nonzero(negatives)
    for r, c in zip(rows, cols):
        report.add(Issue("NEGATIVE_VALUE", ERROR, "cell", f"{index[r]} | {columns[c]}",
                         row=str(index[r]), column=str(columns[c]),
                         value=num.iat[r, c],
                         detail="negative abundance value",
                         suggestion="set a non-negative value or drop"))


def check_duplicates(matrix: pd.DataFrame, report: ValidationReport,
                     reshaped: bool, raw_header: Optional[List[str]],
                     samples_in_columns: bool) -> None:
    """Duplicated identifiers on both axes (including pandas-mangled names)."""
    row_counts = matrix.index.value_counts()
    for sample_id, count in row_counts[row_counts > 1].items():
        positions = [i for i, label in enumerate(matrix.index) if label == sample_id]
        report.add(Issue("DUPLICATE_ROW_ID", ERROR, "row", str(sample_id),
                         detail=f"sample ID occurs {count} times (rows {positions})",
                         suggestion="rename or drop the duplicated sample"))

    col_counts = matrix.columns.value_counts()
    for feature_id, count in col_counts[col_counts > 1].items():
        positions = [i for i, label in enumerate(matrix.columns) if label == feature_id]
        report.add(Issue("DUPLICATE_COL_ID", ERROR, "column", str(feature_id),
                         detail=f"feature ID occurs {count} times (columns {positions})",
                         suggestion="rename or drop the duplicated feature"))

    # pandas renames duplicated header names on read ("a" -> "a.1"), so they
    # are only visible in the original header labels.
    if not reshaped and raw_header and len(raw_header) > 1:
        counts: Dict[str, int] = {}
        for label in raw_header[1:]:
            counts[label] = counts.get(label, 0) + 1
        for label, count in counts.items():
            if count < 2:
                continue
            if samples_in_columns:
                code, axis, human = "DUPLICATE_ROW_ID", "row", "sample"
            else:
                code, axis, human = "DUPLICATE_COL_ID", "column", "feature"
            report.add(Issue(code, ERROR, axis, label,
                             detail=(f"{human} ID '{label}' appears {count} times in "
                                     "the input header; pandas renamed the extra "
                                     "occurrence(s) (e.g. 'name.1')"),
                             suggestion="rename or drop the duplicated column"))


def check_empty(num: pd.DataFrame, report: ValidationReport) -> None:
    """Samples (rows) and features (columns) with no positive values."""
    row_positive = (num > 0).any(axis=1)
    for pos in np.nonzero(~row_positive.to_numpy())[0]:
        sample_id = str(num.index[pos])
        n_null = int(num.iloc[pos].isna().sum())
        detail = "sample has no positive values"
        if n_null:
            detail += f" ({n_null} null cell(s))"
        report.add(Issue("EMPTY_SAMPLE", WARNING, "row", sample_id, detail=detail,
                         suggestion="drop the sample or keep it"))

    col_positive = (num > 0).any(axis=0)
    for pos in np.nonzero(~col_positive.to_numpy())[0]:
        feature_id = str(num.columns[pos])
        n_null = int(num.iloc[:, pos].isna().sum())
        detail = "feature has no positive values"
        if n_null:
            detail += f" ({n_null} null cell(s))"
        report.add(Issue("EMPTY_FEATURE", WARNING, "column", feature_id, detail=detail,
                         suggestion="drop the feature or keep it"))


def check_already_rarefied(num: pd.DataFrame, report: ValidationReport) -> None:
    """INFO when every non-empty sample has the same total abundance."""
    sums = num.sum(axis=1).to_numpy()
    valid = sums[np.isfinite(sums) & (sums > 0)]
    if valid.size >= 2 and np.allclose(valid, valid[0], rtol=1e-6, atol=0):
        report.add(Issue("ALREADY_RAREFIED", INFO, "matrix", "matrix",
                         detail=(f"all {valid.size} non-empty samples have the same "
                                 f"total abundance ({valid[0]:g})"),
                         suggestion="informational - matrix appears to be rarefied already"))


def check_inconsistent_type(num: pd.DataFrame, report: ValidationReport,
                            data_type: str, evidence: Dict[str, Any]) -> None:
    """Aggregated warnings when the data type is inconsistent or imperfect."""
    if data_type == "inconsistent":
        report.add(Issue("INCONSISTENT_TYPE", WARNING, "matrix", "matrix",
                         detail=(f"data type could not be determined: "
                                 f"{evidence.get('reason', '')}"),
                         suggestion=("fix or coerce values so the matrix is "
                                     "counts, binary or relative")))
    elif data_type == "relative":
        sums = num.sum(axis=1)
        bad = sums[sums.notna() & ((sums - 1.0).abs() > SUM_TOLERANCE)]
        if len(bad):
            examples = ", ".join(str(x) for x in list(bad.index)[:5])
            report.add(Issue("INCONSISTENT_TYPE", WARNING, "matrix", "matrix",
                             detail=(f"{len(bad)} sample(s) do not sum to 1 "
                                     f"(e.g. {examples})"),
                             suggestion="fix values or re-classify the data type"))


def crosscheck_metadata(matrix: pd.DataFrame, metadata_path: str, id_col: int,
                        report: ValidationReport) -> Dict[str, Any]:
    """Cross-check sample names between the matrix and a metadata table."""
    metadata = read_table(metadata_path)
    if metadata.empty:
        _fail(f"metadata table is empty: {metadata_path}")
    if id_col < 0 or id_col >= metadata.shape[1]:
        _fail(f"--metadata-id-col {id_col} is out of range "
              f"(metadata has {metadata.shape[1]} column(s))")

    raw_ids = metadata.iloc[:, id_col]
    for pos in np.nonzero(raw_ids.isna().to_numpy())[0]:
        report.add(Issue("NULL_VALUE", WARNING, "metadata", f"metadata row {pos}",
                         detail="metadata row with a missing sample identifier",
                         suggestion="fix the metadata table"))

    ids = raw_ids.fillna("").astype(str).str.strip()
    ids = ids[ids != ""]
    unique_ids = set(ids)
    id_counts = ids.value_counts()
    for sample_id, count in id_counts[id_counts > 1].items():
        report.add(Issue("DUPLICATE_METADATA_ID", ERROR, "metadata", str(sample_id),
                         detail=f"sample ID occurs {count} times in the metadata table",
                         suggestion="rename or drop the duplicated metadata row"))

    matrix_ids = [str(label) for label in matrix.index]
    missing = [s for s in matrix_ids if s not in unique_ids]
    for sample_id in dict.fromkeys(missing):
        report.add(Issue("SAMPLE_NOT_IN_METADATA", WARNING, "row", sample_id,
                         detail="sample present in the matrix but absent from metadata",
                         suggestion="add the sample to the metadata or keep"))
    orphans = sorted(unique_ids - set(matrix_ids))
    for sample_id in orphans:
        report.add(Issue("METADATA_ORPHAN_SAMPLE", WARNING, "metadata", sample_id,
                         detail="metadata sample not present in the matrix",
                         suggestion="check the sample ID or keep"))

    return {
        "metadata_file": str(metadata_path),
        "id_column": str(metadata.columns[id_col]),
        "n_metadata_rows": int(metadata.shape[0]),
        "n_matched": len([s for s in set(matrix_ids) if s in unique_ids]),
        "n_missing_from_metadata": len(dict.fromkeys(missing)),
        "n_orphan_metadata_samples": len(orphans),
    }


###############################################################################
# 8. CLI
###############################################################################


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="module1.py",
        description=("TaxInd module 1: import, validation, cleaning and rarefaction "
                     "of diversity matrices. Stage 1 (validate) is implemented; "
                     "stages 2-4 are planned (see docs/plan_module1.md)."),
    )
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)

    validate = subparsers.add_parser(
        "validate",
        help="Stage 1: import a matrix, detect its structure, run validation checks.",
    )
    validate.add_argument(
        "-i", "--input", required=True,
        help="Input abundance table (CSV, TSV or Excel).",
    )
    validate.add_argument(
        "-o", "--outdir", default="out",
        help="Output directory (default: out).",
    )
    validate.add_argument(
        "--sep", default=None,
        help="Column separator. Auto-detected if not provided.",
    )
    validate.add_argument(
        "--sheet", default=None,
        help="Excel sheet name or 0-based index (default: first sheet).",
    )
    validate.add_argument(
        "--assume-structure", choices=("auto", "wide", "long"), default="auto",
        help="Skip structure detection and assume the given layout (default: auto).",
    )
    validate.add_argument(
        "--sample-col", default=None,
        help="[long format] Column with sample IDs (name or 0-based position).",
    )
    validate.add_argument(
        "--taxon-col", default=None,
        help="[long format] Column with taxon/ASV IDs (name or 0-based position).",
    )
    validate.add_argument(
        "--count-col", default=None,
        help="[long format] Column with abundance values (name or 0-based position).",
    )
    validate.add_argument(
        "--transpose", action="store_true",
        help="Force swapping rows and columns (overrides orientation detection).",
    )
    validate.add_argument(
        "--unit", choices=("auto", "OTU", "ASV"), default="auto",
        help="Override unit detection (default: auto).",
    )
    validate.add_argument(
        "--metadata", default=None,
        help="Metadata table (CSV, TSV or Excel) to cross-check sample names.",
    )
    validate.add_argument(
        "--metadata-id-col", type=int, default=0,
        help="0-based column of the metadata table holding sample IDs (default: 0).",
    )
    return parser


###############################################################################
# 9. Stage 1 driver
###############################################################################


def run_validate(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    if not input_path.is_file():
        _fail(f"input file not found: {args.input}")

    report = ValidationReport(str(input_path), sha256_of_file(input_path))

    # ---- read ---------------------------------------------------------------
    raw = read_table(args.input, sep=args.sep, sheet=args.sheet)
    if raw.empty:
        _fail("input table has no data rows")
    if raw.shape[1] < 2:
        _fail("input table needs an ID column plus at least one abundance column")
    raw_header = read_header_labels(args.input, sep=args.sep, sheet=args.sheet)

    # ---- structure -----------------------------------------------------------
    if args.assume_structure == "auto":
        structure, structure_conf, structure_ev = detect_structure(raw)
    else:
        structure, structure_conf = args.assume_structure, 1.0
        structure_ev = {"forced_by": "--assume-structure"}

    reshape_info: Optional[Dict[str, Any]] = None
    if structure == "long":
        wide, reshape_info = reshape_long_to_wide(
            raw, args.sample_col, args.taxon_col, args.count_col
        )
        detail = (f"reshaped with sample column '{reshape_info['sample_column']}', "
                  f"taxon column '{reshape_info['taxon_column']}', "
                  f"count column '{reshape_info['count_column']}'")
        if reshape_info["rows_belonging_to_duplicate_pairs"]:
            detail += (f"; {reshape_info['rows_belonging_to_duplicate_pairs']} row(s) "
                       "belong to duplicate (sample, taxon) pairs - the first value "
                       "was kept")
        report.add(Issue("LONG_FORMAT_RESOLVED", INFO, "matrix", "matrix",
                         detail=detail, suggestion="informational"))
        if reshape_info["rows_with_missing_identifiers"]:
            report.add(Issue(
                "NULL_VALUE", WARNING, "matrix", "long format rows",
                detail=(f"{reshape_info['rows_with_missing_identifiers']} long-format "
                        "row(s) had missing sample/taxon identifiers and were "
                        "excluded from the matrix"),
                suggestion="fix the identifiers in the input file"))
    else:
        wide = raw.set_index(raw.columns[0])

    if wide.shape[1] == 0:
        _fail("input table has no abundance columns")

    # ---- orientation -----------------------------------------------------------
    if args.transpose:
        detected_orientation = "forced by --transpose"
        orientation_conf, orientation_ev = 1.0, {"forced_by": "--transpose"}
        wide = wide.T
        transposed, forced = True, True
    else:
        detected_orientation, orientation_conf, orientation_ev = detect_orientation(wide)
        transposed = detected_orientation == "samples_cols"
        forced = False
        if transposed:
            wide = wide.T
        if orientation_conf < 0.5:
            report.add(Issue(
                "ORIENTATION_LOW_CONFIDENCE", WARNING, "matrix", "matrix",
                detail=(f"orientation '{detected_orientation}' decided with low "
                        f"confidence ({orientation_conf:.2f}); label evidence was "
                        "inconclusive and dimensions were used as tie-breaker"),
                suggestion="rerun with --transpose if samples and taxa are swapped"))

    matrix = wide
    reshaped = reshape_info is not None

    # ---- label normalization -------------------------------------------------------
    def _is_missing_label(label: Any) -> bool:
        return label is None or (isinstance(label, float) and label != label)

    for labels, axis, human in ((matrix.index, "row", "sample"),
                                 (matrix.columns, "column", "feature")):
        for pos, label in enumerate(labels):
            if _is_missing_label(label):
                report.add(Issue("NULL_VALUE", WARNING, axis,
                                 f"(missing {human} label at position {pos})",
                                 detail=f"missing {human} identifier",
                                 suggestion="provide an identifier or drop"))
    matrix.index = pd.Index(
        ["" if _is_missing_label(x) else str(x) for x in matrix.index],
        name=matrix.index.name,
    )
    matrix.columns = pd.Index(
        ["" if _is_missing_label(x) else str(x) for x in matrix.columns],
        name=matrix.columns.name,
    )

    # ---- numeric coercion ---------------------------------------------------------
    raw_matrix = matrix.copy()
    num, null_mask = coerce_numeric(matrix)

    # ---- type & unit ---------------------------------------------------------------
    data_type, data_type_ev = detect_data_type(num)
    if args.unit != "auto":
        unit, unit_ev = args.unit, {"forced_by": "--unit"}
    else:
        unit, unit_ev = detect_unit(list(matrix.columns))

    report.properties = {
        "structure": {
            "value": structure,
            "confidence": round(float(structure_conf), 2),
            "reshaped_from_long": reshaped,
            "reshape": reshape_info,
            "evidence": structure_ev,
        },
        "orientation": {
            "detected": detected_orientation,
            "confidence": round(float(orientation_conf), 2),
            "transposed": bool(transposed),
            "forced": bool(forced),
            "evidence": orientation_ev,
        },
        "data_type": {"value": data_type, "evidence": data_type_ev},
        "unit": {"value": unit, "evidence": unit_ev},
    }

    # ---- checks ---------------------------------------------------------------------
    check_cell_issues(raw_matrix, num, null_mask, report)
    check_duplicates(matrix, report, reshaped, raw_header, samples_in_columns=transposed)
    check_empty(num, report)
    check_already_rarefied(num, report)
    check_inconsistent_type(num, report, data_type, data_type_ev)

    # ---- metadata cross-check ----------------------------------------------------------
    if args.metadata:
        report.metadata_check = crosscheck_metadata(
            matrix, args.metadata, args.metadata_id_col, report
        )

    # ---- matrix stats --------------------------------------------------------------------
    row_sums = num.sum(axis=1)
    finite_sums = row_sums[np.isfinite(row_sums)]
    depth_range = ([float(finite_sums.min()), float(finite_sums.max())]
                   if len(finite_sums) else None)
    values = num.to_numpy(dtype=float)
    finite_cells = int(np.isfinite(values).sum())
    sparsity = float((values == 0).sum()) / finite_cells if finite_cells else None
    report.stats = {
        "n_samples": int(matrix.shape[0]),
        "n_features": int(matrix.shape[1]),
        "depth_range": depth_range,
        "sparsity": sparsity,
    }

    # ---- outputs ---------------------------------------------------------------------------
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    out = num.copy()
    out.index.name = out.index.name or "sample_id"
    matrix_path = outdir / "matrix_step1.csv"
    out.to_csv(matrix_path, index=True, float_format="%.10g")

    json_path = outdir / "validation_report.json"
    md_path = outdir / "validation_report.md"
    report.write_json(json_path)
    report.write_markdown(md_path)

    # ---- console summary --------------------------------------------------------------------
    counts = report.issue_summary()["by_severity"]
    print("Stage 1 (validate) summary")
    print(f"  Input:       {args.input}")
    print(f"  Structure:   {structure}" + (" (reshaped from long)" if reshaped else ""))
    print("  Orientation: samples as rows"
          + (" (input was transposed)" if transposed else ""))
    print(f"  Data type:  {data_type}")
    print(f"  Unit:       {unit}")
    print(f"  Matrix:     {report.stats['n_samples']} samples x "
          f"{report.stats['n_features']} features")
    print(f"  Issues:     {counts[ERROR]} error(s), {counts[WARNING]} warning(s), "
          f"{counts[INFO]} info")
    if counts[ERROR] or counts[WARNING]:
        print("  Next step:  run Stage 2 ('clean') to review and resolve the issues.")
    print("Outputs:")
    print(f"  Matrix:     {matrix_path}")
    print(f"  Report:     {json_path}")
    print(f"  Report:     {md_path}")
    return 0


###############################################################################
# 10. Run Main
###############################################################################


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "validate":
        sys.exit(run_validate(args))
    # Future stages: clean (Stage 2), rarefy (Stage 3), run (Stage 4).
    build_parser().error(f"unknown command: {args.command}")


if __name__ == "__main__":
    main()
