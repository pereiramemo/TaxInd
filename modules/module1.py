#!/usr/bin/env python3
"""
module1.py - TaxInd module 1: import, validation and cleaning of diversity
matrices (preprocessing stage of the TaxInd pipeline).

Implemented stages:

  validate (Stage 1)
      Reads a raw biological diversity matrix (ASV/OTU/taxa abundance table
      in CSV, TSV or Excel format), diagnoses its structure and contents,
      and writes a normalized matrix plus a validation report that feeds
      Stage 2.

        Detected properties:
          - structure   : wide (samples x taxa) vs long (tidy) format
          - orientation : samples as rows (standard) or columns (transposed)
          - data type   : counts, binary (presence/absence), relative
          - unit        : ASV (sequences) vs OTU (cluster identifiers)

        Validation checks (typed, severity-classified issues): null values,
        negative values, non-numeric values, duplicated identifiers, empty
        samples/features, inconsistent types, already-rarefied matrices,
        low-confidence orientation and, optionally, sample-name mismatches
        against a metadata table.

  clean (Stage 2)
      Reviews the issues detected by Stage 1 and corrects the matrix.
      Interactive terminal session (default) or replay of a saved decisions
      file (--decisions). Actions: set, drop-row, drop-col, rename, coerce,
      keep (accept as is), ignore (postpone). Every decision is recorded
      (decisions.json + change_log.json) and the cleaned matrix is
      re-validated (postclean report).

Outputs (written to --outdir):
  Stage 1: matrix_step1.csv, validation_report.{json,md}
  Stage 2: matrix_clean.csv, decisions.json, change_log.json,
           postclean_report.md

Assumptions:
  - The input table has a header row.
  - In wide format, the first column is the ID axis (samples or taxa).

Exit codes:
  0  stage completed (issues, if any, are listed in the report)
  1  operational failure (missing file, unreadable/empty input, ...)
  2  command line usage error (argparse default)

Usage:
  python modules/module1.py validate -i INPUT [-o OUTDIR] [options]
  python modules/module1.py clean -i MATRIX [--decisions FILE] [options]

Stages 3 (`rarefy`) and 4 (`run`) are planned - see docs/plan_module1.md.
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import hashlib
import io
import json
import re
import shutil
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
# 8. Cleaning (Stage 2)
###############################################################################

ACTIONS = {"set", "drop-row", "drop-col", "rename", "coerce", "keep", "ignore"}
BULK_ACTIONS = {"set", "drop-row", "drop-col", "coerce", "keep", "ignore"}
ACTION_ALIASES = {
    "s": "set", "dr": "drop-row", "dc": "drop-col",
    "rn": "rename", "c": "coerce", "k": "keep", "i": "ignore",
}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _issue_target_ids(issue: Dict[str, Any]) -> Tuple[Optional[str], Optional[str]]:
    """(row, column) identifiers an issue refers to; (None, None) if n/a."""
    axis = issue.get("axis")
    target = issue.get("target")
    if isinstance(target, dict):
        row = target.get("row")
        column = target.get("column", target.get("col"))
        return (str(row) if row is not None else None,
                str(column) if column is not None else None)
    if axis == "row":
        return (str(target) if target is not None else None, None)
    if axis == "column":
        return None, (str(target) if target is not None else None)
    return None, None


def _decision_target(issue: Dict[str, Any]) -> Dict[str, str]:
    """Structured decision target derived from a report issue."""
    row, column = _issue_target_ids(issue)
    target: Dict[str, str] = {}
    if row is not None:
        target["row"] = row
    if column is not None:
        target["column"] = column
    if (not target and issue.get("axis") == "metadata"
            and issue.get("target") is not None):
        target["id"] = str(issue["target"])
    return target


def _decision_key(decision: Dict[str, Any]) -> Tuple[str, str]:
    """Stable key matching a decision to the issue it resolves."""
    return (str(decision.get("issue")),
            json.dumps(decision.get("target") or {}, sort_keys=True))


def _issue_key(issue: Dict[str, Any]) -> Tuple[str, str]:
    """Stable key of a (post-clean) issue for decision matching."""
    return (str(issue.get("code")),
            json.dumps(_decision_target(issue), sort_keys=True))


def _available_actions(issue: Dict[str, Any]) -> List[str]:
    """Actions that make sense for a given issue."""
    axis = issue.get("axis")
    if axis == "cell":
        return ["set", "drop-row", "drop-col", "keep", "ignore"]
    if axis == "row":
        return ["drop-row", "rename", "keep", "ignore"]
    if axis == "column":
        return ["drop-col", "rename", "keep", "ignore"]
    if axis == "metadata":
        return ["keep", "ignore"]
    if issue.get("code") == "INCONSISTENT_TYPE":
        return ["coerce", "keep", "ignore"]
    return ["keep", "ignore"]


def _suggested_action(issue: Dict[str, Any],
                      matrix: pd.DataFrame) -> Optional[str]:
    """Default action shown for an issue (applied when the user presses enter)."""
    code = issue.get("code")
    if code == "NULL_VALUE":
        row, column = _issue_target_ids(issue)
        if row is not None and column is not None:
            return "set 0"
        return None
    if code == "EMPTY_SAMPLE":
        return "drop-row"
    if code == "EMPTY_FEATURE":
        return "drop-col"
    if code == "INCONSISTENT_TYPE":
        row, column = _issue_target_ids(issue)
        if row is None and column is None:
            num, _ = coerce_numeric(matrix)
            finite = num.to_numpy(dtype=float)
            finite = finite[np.isfinite(finite)]
            if finite.size and bool(np.any(
                    np.abs(finite - np.round(finite)) > INT_TOLERANCE)):
                return "coerce"
        return None
    return None


def _parse_action_tokens(tokens: List[str]) -> Tuple[Optional[str], Optional[Any], Optional[str]]:
    """Parse action tokens -> (action, value, error)."""
    if not tokens:
        return None, None, "no action given"
    head = ACTION_ALIASES.get(tokens[0].lower(), tokens[0].lower())
    if head == "set":
        if len(tokens) < 2:
            return None, None, "usage: set <number>"
        try:
            number = float(tokens[1])
        except ValueError:
            return None, None, f"'{tokens[1]}' is not a number"
        if not np.isfinite(number):
            return None, None, "value must be a finite number"
        return "set", number, None
    if head == "rename":
        if len(tokens) < 2:
            return None, None, "usage: rename <new_id>"
        return "rename", tokens[1], None
    if head in ("drop-row", "drop-col", "coerce", "keep", "ignore"):
        return head, None, None
    return None, None, f"unknown action '{tokens[0]}' - type ? for help"


def _print_clean_help() -> None:
    print("  Actions:")
    print("    set <number>   set a new value for the cell")
    print("    drop-row       drop the sample row (every occurrence of the ID)")
    print("    drop-col       drop the feature column (renamed variants too)")
    print("    rename <id>    rename the 2nd..Nth occurrence of a duplicated ID,")
    print("                    or the row/column itself when it is not duplicated")
    print("    coerce         [INCONSISTENT_TYPE] round every non-integer cell")
    print("    keep           accept the data as is (reject the invalidation)")
    print("    ignore         postpone (the issue stays in the report)")
    print("    all <action>   apply <action> to every remaining issue of this code")
    print("    q              finish early (remaining issues are ignored)")
    print("    <enter>        apply the suggested action, if any")


def _new_decision(issue: Dict[str, Any], action: str,
                  value: Any = None) -> Dict[str, Any]:
    decision = {
        "issue": issue.get("code"),
        "target": _decision_target(issue),
        "action": action,
    }
    if value is not None:
        decision["value"] = value
    return decision


def interactive_clean_session(matrix: pd.DataFrame,
                               pending: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Interactive session: the user decides what to do with each issue.

    Returns the decision list (one decision per reviewed issue; issues not
    reached because of an early quit are recorded as ignored).
    """
    decisions: List[Dict[str, Any]] = []
    resolved = [False] * len(pending)
    quit_requested = False
    aborted_by_eof = False
    answers = 0
    index = 0

    print("Type ? at any prompt for the list of actions.")
    while index < len(pending) and not quit_requested:
        if resolved[index]:
            index += 1
            continue
        issue = pending[index]
        code = issue.get("code")
        row, column = _issue_target_ids(issue)
        suggestion = _suggested_action(issue, matrix)

        print("-" * 72)
        header = f"({index + 1}/{len(pending)}) [{issue.get('severity')}] {code}"
        if row is not None:
            header += f" - row '{row}'"
        if column is not None:
            header += f" column '{column}'"
        print(header)
        if issue.get("value") is not None:
            print(f"        value: {issue.get('value')!r}")
        if issue.get("detail"):
            print(f"        {issue['detail']}")
        if suggestion:
            print(f"        suggested: {suggestion}")

        while True:
            try:
                line = input("  action> ").strip()
            except EOFError:
                aborted_by_eof = True
                quit_requested = True
                break
            if not line:
                if suggestion:
                    line = suggestion
                else:
                    print("  (no suggested action for this issue - type ? for help)")
                    continue
            tokens = line.split()
            first = tokens[0].lower()
            if first in ("?", "h", "help"):
                _print_clean_help()
                continue
            if first in ("q", "quit", "exit"):
                quit_requested = True
                break
            if first in ("all", "a"):
                action, value, error = _parse_action_tokens(tokens[1:])
                if error:
                    print(f"  {error}")
                    continue
                if action not in BULK_ACTIONS:
                    print(f"  '{action}' cannot be applied in bulk")
                    continue
                applied = 0
                for pos in range(index, len(pending)):
                    if resolved[pos] or pending[pos].get("code") != code:
                        continue
                    if action not in _available_actions(pending[pos]):
                        continue
                    decisions.append(_new_decision(pending[pos], action, value))
                    resolved[pos] = True
                    applied += 1
                answers += applied
                print(f"  applied '{action}' to {applied} issue(s) of code {code}")
                if resolved[index]:
                    break
                continue
            action, value, error = _parse_action_tokens(tokens)
            if error:
                print(f"  {error}")
                continue
            if action not in _available_actions(issue):
                print(f"  action '{action}' is not applicable here "
                      f"(available: {', '.join(_available_actions(issue))})")
                continue
            if action == "set" and (row is None or column is None):
                print("  'set' applies to cell issues only")
                continue
            if action == "rename":
                new_id = str(value)
                if new_id in matrix.index or new_id in matrix.columns:
                    print(f"  identifier '{new_id}' already exists")
                    continue
            decisions.append(_new_decision(issue, action, value))
            resolved[index] = True
            answers += 1
            break

        if not quit_requested:
            index += 1

    if quit_requested:
        for pos in range(len(pending)):
            if not resolved[pos]:
                decisions.append(_new_decision(pending[pos], "ignore"))
                resolved[pos] = True

    if aborted_by_eof and answers == 0:
        _fail("interactive session received no input - "
               "pass --decisions for non-interactive cleaning")

    kept = sum(1 for d in decisions if d["action"] == "keep")
    ignored = sum(1 for d in decisions if d["action"] == "ignore")
    print("-" * 72)
    print(f"Session finished: {len(decisions)} decision(s) "
          f"({kept} kept, {ignored} ignored, "
          f"{len(decisions) - kept - ignored} action(s)).")
    return decisions


def write_decisions_file(path: Path, matrix_path: Path,
                         decisions: List[Dict[str, Any]]) -> None:
    """Write the decisions file (replayable; records the input checksum)."""
    document = {
        "input": matrix_path.name,
        "input_checksum": sha256_of_file(matrix_path),
        "generated": _now(),
        "decisions": decisions,
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(document, fh, indent=2)
        fh.write("\n")


def load_decisions_file(path: Path, matrix_path: Path) -> List[Dict[str, Any]]:
    """Load a decisions file; hard error on input checksum mismatch."""
    try:
        with open(path, encoding="utf-8") as fh:
            document = json.load(fh)
    except json.JSONDecodeError as exc:
        _fail(f"could not parse decisions file '{path}': {exc}")
    if not isinstance(document, dict) or not isinstance(document.get("decisions"), list):
        _fail(f"malformed decisions file '{path}'")
    checksum = document.get("input_checksum")
    if not checksum:
        _fail(f"decisions file '{path}' has no input checksum")
    if checksum != sha256_of_file(matrix_path):
        _fail(f"decisions file '{path}' was generated for a different matrix "
              f"(file expects {checksum}, input is {sha256_of_file(matrix_path)})")
    return document["decisions"]


def _validate_decision_list(decisions: List[Dict[str, Any]]) -> None:
    """Hard validation of a decisions file (must replay exactly)."""
    for position, decision in enumerate(decisions, start=1):
        action = decision.get("action")
        target = decision.get("target") or {}
        row = target.get("row")
        column = target.get("column", target.get("col"))
        value = decision.get("value")
        if action not in ACTIONS:
            _fail(f"decision #{position}: unknown action '{action}'")
        if action == "set":
            if row is None or column is None:
                _fail(f"decision #{position}: 'set' needs a row and a column")
            if (isinstance(value, bool) or not isinstance(value, (int, float))
                    or not np.isfinite(value)):
                _fail(f"decision #{position}: 'set' needs a finite numeric value")
        elif action == "rename":
            if row is None and column is None:
                _fail(f"decision #{position}: 'rename' needs a row or a column")
            if not value:
                _fail(f"decision #{position}: 'rename' needs a new identifier")
        elif action == "drop-row" and row is None:
            _fail(f"decision #{position}: 'drop-row' needs a row")
        elif action == "drop-col" and column is None:
            _fail(f"decision #{position}: 'drop-col' needs a column")


def _column_candidates(columns: pd.Index, label: str) -> List[str]:
    """Column labels matching an ID, including pandas-renamed variants."""
    exact: List[str] = []
    mangled: List[str] = []
    for column in columns:
        column = str(column)
        if column == label:
            exact.append(column)
        elif column.startswith(label + ".") and column[len(label) + 1:].isdigit():
            mangled.append(column)
    return exact + mangled


def apply_decisions(matrix: pd.DataFrame, decisions: List[Dict[str, Any]]
                    ) -> Tuple[pd.DataFrame, List[Dict[str, Any]], List[str]]:
    """Apply decisions to a matrix.

    Returns (clean_matrix, change_log, notes):
      - change_log: one entry per matrix modification.
      - notes: decisions that could not be applied (e.g. targets removed by
        an earlier decision), reported in the post-clean report.
    """
    clean = matrix.copy()
    change_log: List[Dict[str, Any]] = []
    notes: List[str] = []

    for decision in decisions:
        action = decision.get("action")
        issue = str(decision.get("issue", ""))
        target = decision.get("target") or {}
        row_id = target.get("row")
        column_id = target.get("column", target.get("col"))
        value = decision.get("value")

        if action in ("keep", "ignore"):
            continue  # recorded in decisions.json; no matrix change

        elif action == "set":
            row_positions = [p for p, label in enumerate(clean.index) if label == row_id]
            col_positions = [p for p, label in enumerate(clean.columns) if label == column_id]
            if len(row_positions) != 1 or len(col_positions) != 1:
                notes.append(f"set skipped: '{row_id}'/'{column_id}' "
                             "not present or duplicated")
                continue
            before = clean.iat[row_positions[0], col_positions[0]]
            clean.iat[row_positions[0], col_positions[0]] = float(value)
            change_log.append({
                "timestamp": _now(), "issue": issue,
                "target": {"row": row_id, "column": column_id},
                "action": "set",
                "before": _json_safe(before), "after": float(value),
            })

        elif action == "drop-row":
            occurrences = [label for label in clean.index if label == row_id]
            if not occurrences:
                notes.append(f"drop-row skipped: '{row_id}' not present")
                continue
            n_dropped = len(occurrences)
            clean = clean.drop(index=sorted(set(occurrences)))
            change_log.append({
                "timestamp": _now(), "issue": issue,
                "target": {"row": row_id}, "action": "drop-row",
                "before": f"{n_dropped} row(s) '{row_id}'",
                "after": f"{clean.shape[0]} sample(s) remain",
            })

        elif action == "drop-col":
            candidates = _column_candidates(clean.columns, column_id)
            if not candidates:
                notes.append(f"drop-col skipped: '{column_id}' not present")
                continue
            clean = clean.drop(columns=candidates)
            change_log.append({
                "timestamp": _now(), "issue": issue,
                "target": {"column": column_id}, "action": "drop-col",
                "before": f"column(s) {candidates}",
                "after": f"{clean.shape[1]} feature(s) remain",
            })

        elif action == "rename":
            new_id = str(value)
            if new_id in clean.index or new_id in clean.columns:
                notes.append(f"rename skipped: '{new_id}' already exists")
                continue
            if row_id is not None:
                positions = [p for p, label in enumerate(clean.index)
                             if label == row_id]
                if not positions:
                    notes.append(f"rename skipped: '{row_id}' not present")
                    continue
                labels = list(clean.index)
                renamed: List[str] = []
                if len(positions) == 1:
                    labels[positions[0]] = new_id
                    renamed.append(new_id)
                else:
                    for occurrence, pos in enumerate(positions[1:], start=1):
                        replacement = (new_id if occurrence == 1
                                       else f"{new_id}.{occurrence}")
                        labels[pos] = replacement
                        renamed.append(replacement)
                clean.index = pd.Index(labels, name=clean.index.name)
                change_log.append({
                    "timestamp": _now(), "issue": issue,
                    "target": {"row": row_id}, "action": "rename",
                    "before": row_id, "after": renamed,
                })
            elif column_id is not None:
                positions = _column_candidates(clean.columns, column_id)
                if not positions:
                    notes.append(f"rename skipped: '{column_id}' not present")
                    continue
                order = [p for p, label in enumerate(clean.columns)
                         if label in positions]
                labels = list(clean.columns)
                renamed = []
                if len(order) == 1:
                    labels[order[0]] = new_id
                    renamed.append(new_id)
                else:
                    for occurrence, pos in enumerate(order[1:], start=1):
                        replacement = (new_id if occurrence == 1
                                       else f"{new_id}.{occurrence}")
                        labels[pos] = replacement
                        renamed.append(replacement)
                clean.columns = pd.Index(labels, name=clean.columns.name)
                change_log.append({
                    "timestamp": _now(), "issue": issue,
                    "target": {"column": column_id}, "action": "rename",
                    "before": column_id, "after": renamed,
                })
            else:
                notes.append("rename: missing row/column in decision")

        elif action == "coerce":
            if row_id is not None and column_id is not None:
                row_positions = [p for p, label in enumerate(clean.index) if label == row_id]
                col_positions = [p for p, label in enumerate(clean.columns)
                                 if label == column_id]
                if len(row_positions) != 1 or len(col_positions) != 1:
                    notes.append(f"coerce skipped: '{row_id}'/'{column_id}' "
                                  "not present or duplicated")
                    continue
                cell = clean.iat[row_positions[0], col_positions[0]]
                if cell is None or (isinstance(cell, float) and np.isnan(cell)):
                    clean.iat[row_positions[0], col_positions[0]] = 0.0
                    change_log.append({
                        "timestamp": _now(), "issue": issue,
                        "target": {"row": row_id, "column": column_id},
                        "action": "coerce", "before": None, "after": 0.0,
                    })
                elif isinstance(cell, float) and abs(cell - round(cell)) > INT_TOLERANCE:
                    rounded = float(round(cell))
                    clean.iat[row_positions[0], col_positions[0]] = rounded
                    change_log.append({
                        "timestamp": _now(), "issue": issue,
                        "target": {"row": row_id, "column": column_id},
                        "action": "coerce", "before": cell, "after": rounded,
                    })
                else:
                    notes.append(f"coerce skipped at '{row_id}'/'{column_id}': "
                                  "nothing to convert")
            else:
                values, _ = coerce_numeric(clean)
                numeric = values.to_numpy(dtype=float)
                mask = np.isfinite(numeric) & (
                    np.abs(numeric - np.round(numeric)) > INT_TOLERANCE)
                positions = np.nonzero(mask)
                count = int(positions[0].size)
                if count == 0:
                    notes.append("coerce: no non-integer cells to round")
                    continue
                examples = [f"{values.index[r]} | {values.columns[c]}"
                            for r, c in zip(positions[0][:5], positions[1][:5])]
                clean = values.round()
                change_log.append({
                    "timestamp": _now(), "issue": issue,
                    "target": {}, "action": "coerce",
                    "before": f"{count} non-integer cell(s)",
                    "after": "rounded to integers",
                    "examples": examples,
                })
        else:
            notes.append(f"unknown action '{action}'")

    return clean, change_log, notes


def _postclean_validation(matrix: pd.DataFrame, matrix_clean_path: Path,
                          metadata: Optional[str] = None,
                          metadata_id_col: int = 0) -> ValidationReport:
    """Re-run detection and checks on the cleaned matrix."""
    report = ValidationReport(str(matrix_clean_path),
                               sha256_of_file(matrix_clean_path))
    raw_matrix = matrix.copy()
    num, null_mask = coerce_numeric(matrix)
    data_type, data_type_ev = detect_data_type(num)
    unit, unit_ev = detect_unit(list(matrix.columns))
    report.properties = {
        "structure": {"value": "wide", "confidence": 1.0,
                      "reshaped_from_long": False, "reshape": None,
                      "evidence": {"cleaned_matrix": True}},
        "orientation": {"detected": "samples_rows", "confidence": 1.0,
                        "transposed": False, "forced": False,
                        "evidence": {"normalized_in_stage": 1}},
        "data_type": {"value": data_type, "evidence": data_type_ev},
        "unit": {"value": unit, "evidence": unit_ev},
    }
    check_cell_issues(raw_matrix, num, null_mask, report)
    check_duplicates(matrix, report, reshaped=True, raw_header=None,
                     samples_in_columns=False)
    check_empty(num, report)
    check_already_rarefied(num, report)
    check_inconsistent_type(num, report, data_type, data_type_ev)
    if metadata:
        report.metadata_check = crosscheck_metadata(
            matrix, metadata, metadata_id_col, report)
    values = num.to_numpy(dtype=float)
    finite_cells = int(np.isfinite(values).sum())
    sums = num.sum(axis=1)
    finite_sums = sums[np.isfinite(sums)]
    report.stats = {
        "n_samples": int(matrix.shape[0]),
        "n_features": int(matrix.shape[1]),
        "depth_range": ([float(finite_sums.min()), float(finite_sums.max())]
                         if len(finite_sums) else None),
        "sparsity": (float((values == 0).sum()) / finite_cells
                     if finite_cells else None),
    }
    return report


def write_postclean_markdown(path: Path, report: ValidationReport,
                             statuses: List[str],
                             decision_counts: Dict[str, int],
                             changes_count: int, notes: List[str],
                             before_stats: Dict[str, Any]) -> None:
    """Human-readable post-clean report (Stage 2)."""

    def short(value: Any, limit: int = 24) -> str:
        text = str(value)
        return text if len(text) <= limit else text[: limit - 1] + "…"

    lines: List[str] = []
    add = lines.append
    add("# TaxInd module1 - Post-clean Report (Stage 2)")
    add("")
    add(f"- **Cleaned matrix:** `{report.input_path}`")
    add(f"- **Checksum:** `{report.input_checksum}`")
    add(f"- **Generated:** {report.generated}")
    add("")
    add("## Cleaning summary")
    add("")
    add("| Item | Count |")
    add("|---|---|")
    add(f"| Decisions | {decision_counts['total']} |")
    add(f"| Kept (accepted as is) | {decision_counts['kept']} |")
    add(f"| Ignored (pending) | {decision_counts['ignored']} |")
    add(f"| Matrix changes | {changes_count} |")
    if notes:
        add(f"| Skipped decisions | {len(notes)} |")
    add("")
    add("## Matrix dimensions")
    add("")
    add("| | Before | After |")
    add("|---|---|---|")
    add(f"| Samples | {before_stats.get('n_samples')} | "
        f"{report.stats.get('n_samples')} |")
    add(f"| Features | {before_stats.get('n_features')} | "
        f"{report.stats.get('n_features')} |")
    add("")
    add("## Re-detected properties")
    add("")
    add("| Property | Value | Notes |")
    add("|---|---|---|")
    data_type = report.properties.get("data_type", {})
    evidence = data_type.get("evidence", {})
    add(f"| Data type | {data_type.get('value')} | "
        f"{evidence.get('reason') or evidence.get('values', '')} |")
    unit = report.properties.get("unit", {})
    add(f"| Unit | {unit.get('value')} | |")
    add("")
    add("## Issues after cleaning")
    add("")
    if not report.issues:
        add("No issues detected in the cleaned matrix - ready for Stage 3 (`rarefy`).")
    else:
        add("| Severity | Code | Target | Status |")
        add("|---|---|---|---|")
        pairs = sorted(
            zip(report.issues, statuses),
            key=lambda pair: (SEVERITY_RANK.get(pair[0].severity, 9),
                              pair[0].code, str(pair[0].target_id)),
        )
        for issue, status in pairs:
            row, column = _issue_target_ids(issue.to_dict())
            if row is not None and column is not None:
                target = f"{short(row)} | {short(column)}"
            elif row is not None:
                target = short(row)
            elif column is not None:
                target = short(column)
            else:
                target = short(issue.target_id) or "-"
            add(f"| {issue.severity} | `{issue.code}` | {target} | {status} |")
    if notes:
        add("")
        add("## Skipped decisions")
        add("")
        for note in notes:
            add(f"- {note}")
    add("")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def run_clean(args: argparse.Namespace) -> int:
    matrix_path = Path(args.input)
    if not matrix_path.is_file():
        _fail(f"input file not found: {args.input}")
    report_path = (Path(args.report) if args.report
                   else matrix_path.parent / "validation_report.json")
    if not report_path.is_file():
        _fail(f"validation report not found: {report_path} "
              "(run 'validate' first, or pass --report)")

    with open(report_path, encoding="utf-8") as fh:
        stage1 = json.load(fh)
    issues = stage1.get("issues", [])
    pending = [i for i in issues if i.get("severity") in (ERROR, WARNING)]
    info_count = len(issues) - len(pending)

    matrix = pd.read_csv(matrix_path, index_col=0)
    if matrix.empty:
        _fail("matrix is empty")
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)

    outdir = Path(args.outdir) if args.outdir else matrix_path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    if args.decisions:
        decisions_path = Path(args.decisions)
        if not decisions_path.is_file():
            _fail(f"decisions file not found: {args.decisions}")
        decisions = load_decisions_file(decisions_path, matrix_path)
        _validate_decision_list(decisions)
        saved_decisions_path = outdir / "decisions.json"
        if decisions_path.resolve() != saved_decisions_path.resolve():
            shutil.copyfile(decisions_path, saved_decisions_path)
    else:
        print("Stage 2 (clean) - interactive session")
        print(f"  Matrix: {matrix_path} "
              f"({matrix.shape[0]} samples x {matrix.shape[1]} features)")
        print(f"  Issues: {len(pending)} to review"
              + (f" (+{info_count} informational, no action needed)"
                 if info_count else ""))
        print()
        try:
            decisions = interactive_clean_session(matrix, pending)
        except KeyboardInterrupt:
            print("\nAborted - no changes were applied.")
            return 1
        saved_decisions_path = outdir / "decisions.json"
        write_decisions_file(saved_decisions_path, matrix_path, decisions)
        print(f"\nDecisions saved to: {saved_decisions_path} (saved before applying)")

    clean, changes, notes = apply_decisions(matrix, decisions)

    matrix_clean_path = outdir / "matrix_clean.csv"
    out = clean.copy()
    out.index.name = out.index.name or "sample_id"
    out.to_csv(matrix_clean_path, index=True, float_format="%.10g")

    change_log_path = outdir / "change_log.json"
    with open(change_log_path, "w", encoding="utf-8") as fh:
        json.dump(changes, fh, indent=2)
        fh.write("\n")

    if clean.shape[0] == 0 or clean.shape[1] == 0:
        print("Warning: cleaning removed every sample or feature; "
              "the cleaned matrix is empty.")

    postclean_path = outdir / "postclean_report.md"
    postclean = _postclean_validation(clean, matrix_clean_path,
                                      metadata=args.metadata,
                                      metadata_id_col=args.metadata_id_col)
    kept_keys = {_decision_key(d) for d in decisions if d.get("action") == "keep"}
    ignored_keys = {_decision_key(d) for d in decisions if d.get("action") == "ignore"}
    statuses = []
    for issue in postclean.issues:
        key = _issue_key(issue.to_dict())
        if key in kept_keys:
            statuses.append("kept (accepted)")
        elif key in ignored_keys:
            statuses.append("ignored (pending)")
        else:
            statuses.append("remaining")

    decision_counts = {
        "total": len(decisions),
        "kept": sum(1 for d in decisions if d.get("action") == "keep"),
        "ignored": sum(1 for d in decisions if d.get("action") == "ignore"),
    }
    write_postclean_markdown(
        postclean_path, postclean, statuses, decision_counts,
        len(changes), notes, stage1.get("matrix", {}),
    )

    counts = postclean.issue_summary()["by_severity"]
    print("Stage 2 (clean) summary")
    print(f"  Matrix:      {matrix_path.name} -> {matrix_clean_path.name}")
    print(f"  Decisions:   {decision_counts['total']} "
          f"({decision_counts['kept']} kept, {decision_counts['ignored']} ignored, "
          f"{decision_counts['total'] - decision_counts['kept'] - decision_counts['ignored']} action(s))")
    if notes:
        print(f"  Skipped:     {len(notes)} decision(s) could not be applied (see report)")
    print(f"  Postclean:   {counts[ERROR]} error(s), {counts[WARNING]} warning(s), "
          f"{counts[INFO]} info")
    print("Outputs:")
    print(f"  Matrix:     {matrix_clean_path}")
    print(f"  Decisions:  {saved_decisions_path}")
    print(f"  Change log: {change_log_path}")
    print(f"  Report:     {postclean_path}")
    print("Next step: run Stage 3 ('rarefy') on the cleaned matrix.")
    return 0


###############################################################################
# 9. CLI
###############################################################################


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="module1.py",
        description=("TaxInd module 1: import, validation, cleaning and rarefaction "
                     "of diversity matrices. Stages 1 (validate) and 2 (clean) are "
                     "implemented; stages 3-4 are planned (see docs/plan_module1.md)."),
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

    clean = subparsers.add_parser(
        "clean",
        help="Stage 2: review issues and correct the matrix "
             "(interactive, or replay a decisions file).",
    )
    clean.add_argument(
        "-i", "--input", required=True,
        help="Matrix to clean (matrix_step1.csv from Stage 1).",
    )
    clean.add_argument(
        "-o", "--outdir", default=None,
        help="Output directory (default: directory of the input matrix).",
    )
    clean.add_argument(
        "--report", default=None,
        help="Stage 1 validation report "
             "(default: validation_report.json next to the input matrix).",
    )
    clean.add_argument(
        "--decisions", default=None,
        help="Decisions file to replay non-interactively (generated by a "
             "previous interactive session).",
    )
    clean.add_argument(
        "--metadata", default=None,
        help="Metadata table to repeat the sample-name cross-check.",
    )
    clean.add_argument(
        "--metadata-id-col", type=int, default=0,
        help="0-based column of the metadata table holding sample IDs (default: 0).",
    )
    return parser


###############################################################################
# 10. Stage 1 driver
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
# 11. Run Main
###############################################################################


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "validate":
        sys.exit(run_validate(args))
    if args.command == "clean":
        sys.exit(run_clean(args))
    # Future stages: rarefy (Stage 3), run (Stage 4).
    build_parser().error(f"unknown command: {args.command}")


if __name__ == "__main__":
    main()
