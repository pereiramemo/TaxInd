# module1.py — Implementation Plan

**Status:** approved — validated against the TaxInd project plan
**Target file:** `modules/module1.py` (English: code, docstrings, CLI help, reports)
**Position in pipeline:** `module1.py` (import + validation + cleaning + rarefaction) → `1.1-div.py` / `1.2-div.py` / `1.3-div.py`

---

## 1. Overview

`module1.py` is the preprocessing stage of TaxInd. It takes a raw biological
diversity matrix (CSV, TSV or Excel), diagnoses its structure and quality,
lets the user resolve the detected problems, rarefies the result and exports
a clean matrix that is directly consumable by the downstream diversity
modules.

**Responsibilities:** import → validation → cleaning → rarefaction → export.

**Confirmed design decisions**

| Decision | Choice |
|---|---|
| Stage 2 interaction | Hybrid: interactive terminal session that generates a **replayable decisions file** (reproducible, automatable via `--decisions`) |
| Stage 3 default | Depth = minimum sample sum; **single random draw** per sample (classic rarefaction); fixed `--seed` for reproducibility |
| Stage 4 content | Export, consolidated pipeline report, integration with modules 1.1/1.2/1.3, tests, README update |
| File name | `modules/module1.py` |
| Metadata validation | Yes — optional `--metadata` cross-check in Stage 1 (sample-name matching). Taxonomy cross-validation deferred to the module that consumes taxonomy |
| Transformations | Stage 3 emits `relative` and `binary` variants of the rarefied matrix, in addition to counts |
| Containerization | Deferred — pinned-deps YAML + Dockerfile per module will be added later, project-wide |

### Alignment with the TaxInd pipeline plan

Validated against the project-level plan (Nextflow-orchestrated, containerized
modules, each runnable as an independent script). `module1.py` covers workflow
steps 1–2 (*data validation* → *rarefaction and data transformations*) and
stays callable as a standalone script (Nextflow compatibility: §2). Gamma
diversity and environmental metrics are future modules, not part of module1
(see §8).

---

## 2. Architecture

### CLI (subcommands, one per stage + full pipeline)

```bash
# Stage 1 — validate & diagnose
python modules/module1.py validate -i INPUT [-o OUTDIR] [--metadata FILE] [options]

# Stage 2 — clean (interactive or replayed from a decisions file)
python modules/module1.py clean -i MATRIX [-o OUTDIR] [--decisions FILE] [options]

# Stage 3 — rarefy (+ optional transformations)
python modules/module1.py rarefy -i MATRIX [-o OUTDIR] [--transform relative binary] [options]

# Full pipeline (stage 1 + non-interactive clean + stage 3)
python modules/module1.py run -i INPUT [-o OUTDIR] [--metadata FILE] [--transform ...] [options]
```

### Output directory layout

```
out/
├── matrix_step1.csv             # normalized matrix after stage 1 (wide, samples as rows)
├── validation_report.json        # machine-readable (issues + metadata)
├── validation_report.md          # human-readable
├── matrix_clean.csv              # after stage 2
├── decisions.json                # user decisions (replayable)
├── change_log.json               # record of every modification
├── postclean_report.md           # re-validation after cleaning
├── rarefied_matrix.csv           # final output (stage 3, counts)
├── rarefied_matrix_relative.csv  # optional --transform relative
├── rarefied_matrix_binary.csv    # optional --transform binary
├── rarefaction_report.md         # depth, seed, excluded samples, before/after
└── pipeline_report.md            # consolidated (stage 4)
```

### Data model (in-memory)

| Object | Contents |
|---|---|
| `DiversityMatrix` | `DataFrame` + metadata: `structure` (wide/long), `orientation`, `data_type`, `unit` (OTU/ASV/unknown), source file |
| `Issue` (dataclass) | `code`, `severity`, `axis`, `target_id`, `detail`, `value`, `suggestion` |
| `ValidationReport` | list of issues + summary stats; `to_json()` / `to_markdown()` |
| `ChangeLog` | `timestamp`, `issue_code`, `target`, `action`, `before`, `after` |

### Severity taxonomy

| Severity | Meaning | Examples |
|---|---|---|
| `ERROR` | blocks analysis | non-numeric values, negative values, duplicated IDs |
| `WARNING` | probably fixable | nulls, empty samples/features, inconsistent type, sample-name mismatches |
| `INFO` | observation | detected structure/orientation, matrix already rarefied |

### Dependencies

`numpy`, `pandas` (already used by the project) + `openpyxl` (Excel input).
No other new dependencies.

### Nextflow compatibility (design constraints)

Although orchestration is deferred, the module stays pipeline-ready:

- `clean` runs non-interactively whenever `--decisions` is provided
  (interactive mode is never required for automation).
- Documented exit codes: `0` success, `1` validation failure / user error,
  `2` usage error (argparse default).
- Machine-readable JSON reports with stable schemas
  (`validation_report.json`, `decisions.json`).
- All outputs written inside `-o OUTDIR`; nothing depends on stdout parsing.

---

## 3. Stage 1 — Import & Validation

### 3.1 Reading

`read_table(path, sep=None, sheet=None) -> DataFrame`

- Formats: `.csv`, `.tsv`, `.txt`, `.xlsx`, `.xls`.
- Separator: sniffed from the first lines (comma / tab / semicolon); extension
  used only as fallback. `--sep` overrides.
- Encoding fallback chain: `utf-8` → `utf-8-sig` → `latin-1`.

### 3.2 Structure detection — wide vs long

`detect_structure(df) -> ("wide" | "long", confidence)`

- **long:** ~3 columns with pattern (categorical, categorical, numeric) and
  repeated values in the categorical columns.
- **wide:** first column is ID-like, remaining columns numeric.
- If long → `reshape_long_to_wide(df, sample_col, taxon_col, count_col)`.
  Column names auto-detected, overridable via `--sample-col/--taxon-col/--count-col`.
- Override: `--assume-structure {auto,wide,long}`.

### 3.3 Orientation detection — samples vs OTU/ASV

`detect_orientation(df) -> ("samples_rows" | "samples_cols", confidence)`

Heuristics (combined into a confidence score, reported in the validation
report):

1. Axis whose labels are pure `ACGT` and length ≥ 50 → taxa axis (ASVs).
2. Axis matching sample-like patterns (`S\d+`, `SAM\d+`, `sample_*`, ...).
4. Tie-breaker: taxa axis usually has more unique labels than the sample axis.
4. If confidence is low → emit `ORIENTATION_LOW_CONFIDENCE` (WARNING) and
   let the user decide with `--transpose`.

### 3.4 Data type detection

`detect_data_type(df) -> ("counts" | "binary" | "relative" | "inconsistent", evidence)`

| Detected | Rule |
|---|---|
| `binary` | all values ∈ {0, 1} |
| `counts` | numeric, ≥ 0, integer within tolerance (`abs(x − round(x)) ≤ 1e-9`) |
| `relative` | floats in [0, 1], non-integer, row sums ≈ 1 (or column sums ≈ 1) |
| `inconsistent` | none of the above → one issue per offending cell/value class |

### 3.5 Unit detection — OTU vs ASV

`detect_unit(df) -> ("ASV" | "OTU" | "unknown", evidence)`

- `ASV`: column (or row) labels are nucleotide sequences (`^[ACGT]+$`, length ≥ 50) —
  matches the repo's example table.
- `OTU`: labels like `OTU\d+`, `ASV\d+` (short IDs), or taxonomy strings
  (`k__...;p__...;c__...`).
- `unknown` otherwise (user override: `--unit {auto,OTU,ASV}`).

### 3.6 Validation checks → issue codes

| Code | Severity | Trigger | Suggested fix |
|---|---|---|---|
| `NULL_VALUE` | WARNING | empty cell, `NaN`, `"NA"`, `"N/A"`, `"-"`, whitespace-only | `set 0` or drop |
| `NEGATIVE_VALUE` | ERROR | value < 0 | fix value or drop |
| `NON_NUMERIC` | ERROR | cell fails numeric coercion | fix or drop |
| `DUPLICATE_ROW_ID` | ERROR | repeated sample ID | rename or drop |
| `DUPLICATE_COL_ID` | ERROR | repeated feature ID | rename or drop |
| `EMPTY_SAMPLE` | WARNING | all-zero row | drop or keep |
| `EMPTY_FEATURE` | WARNING | all-zero column | drop or keep |
| `INCONSISTENT_TYPE` | WARNING | decimals in counts / values > 1 in binary / sums ≠ 1 in relative | coerce or reclassify |
| `ALREADY_RAREFIED` | INFO | all row sums equal (within tolerance) | informational |
| `LONG_FORMAT_RESOLVED` | INFO | reshaped long → wide | — |
| `ORIENTATION_LOW_CONFIDENCE` | WARNING | heuristics inconclusive | `--transpose` |
| `SAMPLE_NOT_IN_METADATA` | WARNING | matrix sample missing from metadata table | add to metadata or `keep` |
| `METADATA_ORPHAN_SAMPLE` | WARNING | metadata sample missing from matrix | check ID or `keep` |
| `DUPLICATE_METADATA_ID` | ERROR | repeated sample ID in the metadata table | rename or drop |

### 3.7 Metadata cross-validation (optional)

`--metadata FILE` (CSV/TSV, read with the same `read_table()` helper;
`--metadata-id-col N`, default 0) validates that sample names match between
the abundance matrix and the metadata table:

- `SAMPLE_NOT_IN_METADATA` — sample present in the matrix but absent in metadata.
- `METADATA_ORPHAN_SAMPLE` — metadata sample not present in the matrix.
- `DUPLICATE_METADATA_ID` — repeated sample IDs in the metadata table.

The metadata table is only **validated**, never modified; module1 does not
use it beyond the cross-check (environmental metrics are future modules).
If provided, the cross-check is repeated in the Stage 2 post-clean report
(dropping samples in Stage 2 may surface new `METADATA_ORPHAN_SAMPLE` issues).

### 3.8 Stage 1 outputs

- `matrix_step1.csv` — normalized matrix (wide, samples as rows, numeric).
- `validation_report.json` + `validation_report.md`:
  detected structure/orientation/type/unit, summary stats
  (n_samples, n_taxa, depth range, sparsity), issues grouped by severity,
  and (if `--metadata`) the sample-name matching summary.

---

## 4. Stage 2 — Cleaning & Manual Correction (hybrid)

### 4.1 Flow

1. Load matrix + `validation_report.json` from stage 1.
2. Show invalidation list grouped by severity (`ERROR` → `WARNING` → `INFO`),
   each with code, location, current value and suggested fix.
3. Interactive session: the user decides per issue (see actions below).
   Bulk shortcuts available (e.g. *apply suggested fix to all `NULL_VALUE`*).
4. Decisions are saved to `decisions.json` (with input checksum) **before**
   being applied — the session can be replayed non-interactively.
5. `apply_decisions()` → cleaned matrix + `change_log.json`.
6. Re-validate the cleaned matrix → `postclean_report.md`
   (only user-`keep`ed issues should remain; if `--metadata` was provided,
   the sample-name cross-check is repeated — dropped samples may surface
   new `METADATA_ORPHAN_SAMPLE` issues).

Non-interactive mode: `--decisions decisions.json` replays the exact same
cleaning without prompts (checksum mismatch on the input → hard error).

### 4.2 Actions per issue

| Action | Meaning |
|---|---|
| `set <value>` | manually modify a value |
| `drop-row` | remove a sample (drops every occurrence of a duplicated ID) |
| `drop-col` | remove a feature (drops pandas-renamed variants too) |
| `rename <new_id>` | resolve a duplicated ID |
| `coerce` | convert type (e.g. `"3"` → `3`) |
| `keep` | **reject** the invalidation — accept the data as is |
| `ignore` | postpone (issue stays in the report) |

### 4.3 Decisions file schema (`decisions.json`)

```json
{
  "input": "matrix_step1.csv",
  "input_checksum": "sha256:...",
  "generated": "2026-11-11T12:00:00",
  "decisions": [
    {"issue": "NEGATIVE_VALUE", "target": {"row": "S03", "column": "ASV-12"}, "action": "set", "value": 1},
    {"issue": "DUPLICATE_ROW_ID", "target": {"row": "S07"}, "action": "drop-row"},
    {"issue": "EMPTY_FEATURE", "target": {"column": "ACGTTA..."}, "action": "keep"}
  ]
}
```

### 4.4 Change log schema (`change_log.json`)

```json
[
  {
    "timestamp": "2026-11-11T12:00:01",
    "issue": "NEGATIVE_VALUE",
    "target": {"row": "S03", "col": "ASV-12"},
    "action": "set",
    "before": -4,
    "after": 1
  }
]
```

### 4.5 Stage 2 outputs

`matrix_clean.csv`, `decisions.json`, `change_log.json`, `postclean_report.md`.

---

## 5. Stage 3 — Rarefaction & transformations

### 5.1 Algorithm (default: classic single-draw rarefaction)

```python
rng = np.random.default_rng(seed)
for each sample:
    pool = np.repeat(np.arange(n_taxa), counts)   # one entry per individual
    draw = rng.permutation(pool)[:depth]           # subsample WITHOUT replacement
    rarefied = np.bincount(draw, minlength=n_taxa)
```

- **Depth (default `auto`):** minimum positive sample sum → keeps all samples.
  Alternatives: `min`, `median`, or an explicit integer.
- **Single draw per sample** (confirmed default); `--iterations N` averages N
  draws and rounds (documented as non-classical).
- **Reproducibility:** `--seed` with a fixed documented default.
- **Underflow** (sample sum < depth): `--on-underflow {exclude, keep}`.
  A sample whose sum equals the depth passes through unchanged.

### 5.2 Input requirements

- Rarefaction is only defined for `counts`. `binary` is treated as 0/1 counts
  (depth must be ≤ number of positive features). `relative` → hard error with
  a message suggesting to provide absolute counts.

### 5.3 Transformations (post-rarefaction)

`--transform {relative,binary}` — repeatable flag; each value emits an
additional output file derived from the **rarefied** matrix:

| Transform | Rule | Output |
|---|---|---|
| `relative` | each row divided by its (post-rarefaction) sum | `rarefied_matrix_relative.csv` |
| `binary` | `x > 0 → 1` (presence/absence) | `rarefied_matrix_binary.csv` |

Default: no transformation (counts only). The counts matrix
(`rarefied_matrix.csv`) is always written. The `binary` output feeds the
*presence/absence alpha/beta/gamma* analyses of the pipeline; the `relative`
output feeds abundance-based analyses on compositional data.

### 5.4 Stage 3 outputs

- `rarefied_matrix.csv` — integer counts, samples as rows, first column = sample IDs
  (exact input format of `1.1-div.py`, `1.2-div.py`, `1.3-div.py`).
- `rarefied_matrix_relative.csv` / `rarefied_matrix_binary.csv` — only for
  the transforms requested via `--transform`.
- `rarefaction_report.md` — chosen depth, seed, excluded samples,
  per-sample before/after totals, transforms generated.

---

## 6. Stage 4 — Export, Reports & Integration

1. **`pipeline_report.md`** — consolidated view: input file + checksum,
   parameters and results of every stage, final dimensions.
2. **Compatibility guarantee** — the rarefied matrices feed directly into the
   existing modules:
   ```bash
   python modules/1.1-div.py -i out/rarefied_matrix.csv
   python modules/1.2-div.py -i out/rarefied_matrix.csv
   python modules/1.3-div.py -i out/rarefied_matrix.csv
   # presence/absence variant for binary analyses:
   python modules/1.3-div.py -i out/rarefied_matrix_binary.csv
   ```
3. **Tests** — `tests/test_module1.py`:
   - detection: wide/long, orientation, data type, unit (uses
     `data/examples/asv_abundance_table.csv` + synthetic fixtures);
   - validation: nulls, negatives, non-numeric, duplicates, empty rows/cols,
     metadata sample-name matching (matches, orphans, duplicates);
   - rarefaction: row sums == depth, no negatives, values ≤ original,
     same seed → identical output, underflow exclusion;
   - transformations: relative row sums == 1, binary values ∈ {0, 1};
   - cleaning: decisions file replay produces identical matrix.
4. **README** — add module 1 section; fix existing path references
   (`examples/...` → `data/examples/asv_abundance_table.csv`).

---

## 7. Deliverables checklist

- [x] **Stage 1:** `validate` subcommand — reading, structure/orientation/
      type/unit detection, issue checks, optional metadata cross-check,
      reports. *(implemented; QA suite `tests/test_module1.py` — 17 tests
      driven by `data/examples/asv_abundance_table.csv`, all passing)*
- [x] **Stage 2:** `clean` subcommand — interactive session, decisions file,
      actions, change log, post-clean report. *(implemented; QA suite
      extended to 24 tests, all passing)*
- [ ] **Stage 3:** `rarefy` subcommand — rarefaction + relative/binary
      transformations + report.
- [ ] **Stage 4:** `run` subcommand, consolidated report, tests, README update.

Each stage is implemented and reviewed before moving to the next one.

---

## 8. Non-goals / future work

- No automatic imputation of missing values (every change is a user decision).
- No taxonomy-aware features (merge by taxonomy, etc.); taxonomy
  cross-validation is deferred to the module that consumes taxonomic
  annotation.
- Interactive session is plain terminal prompts (no TUI framework).
- `module1.py` does not rarefy presence/absence tables into counts.
- Gamma diversity and environmental metrics — future modules (module1's
  output already supports them).
- Containerization (pinned-deps YAML + Dockerfile per module) and Nextflow
  orchestration — deferred to a later, project-wide stage covering all
  modules.
