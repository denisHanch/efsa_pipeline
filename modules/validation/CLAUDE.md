# Validation Module — CLAUDE.md

Comprehensive reference for the `modules/validation/` directory. This module
validates and standardises genomic input files (FASTA/GenBank genomes, FASTQ/BAM
reads, GFF/GTF/BED features) before they enter the Nextflow pipeline. At the end
it emits `validated_params.json`, which Nextflow consumes via `-params-file`.

---

## Directory layout

```
modules/validation/
├── main.py                        # CLI entry point — orchestrates full workflow
├── validation.sh                  # Bash wrapper with default paths and safety clears
├── nextflow_params_handler.py     # Builds and serialises validated_params.json
├── test_nextflow_params.py        # Tests for nextflow_params_handler
├── README.md                      # User-facing how-to and config guide pointer
├── TODO.md                        # Open work items
└── validation-pkg/
    ├── setup.py
    ├── requirements.txt           # biopython, numpy, pysam, structlog, typing_extensions
    ├── requirements-dev.txt       # pytest
    ├── validation_pkg/
    │   ├── __init__.py            # Functional API (validate_genome, validate_reads, …); __all__ exports ValidationReport and __version__ individually
    │   ├── config_manager.py      # JSON config loader → Config / *Config dataclasses
    │   ├── exceptions.py          # Exception hierarchy rooted at ValidationError
    │   ├── report.py              # ValidationReport — collects results, writes report.txt
    │   ├── utils/
    │   │   ├── base_settings.py   # BaseSettings, BaseOutputMetadata, BaseValidatorSettings
    │   │   ├── base_validator.py  # BaseValidator abstract class
    │   │   ├── file_handler.py    # Compression, format detection, file I/O utilities
    │   │   ├── formats.py         # Enums: CodingType, GenomeFormat, ReadFormat, FeatureFormat
    │   │   ├── path_utils.py      # Path resolution + path-traversal security
    │   │   └── sequence_stats.py  # N50 calculation
    │   └── validators/
    │       ├── genome_validator.py        # GenomeValidator + GenomeOutputMetadata
    │       ├── read_validator.py          # ReadValidator + ReadOutputMetadata
    │       ├── feature_validator.py       # FeatureValidator + FeatureOutputMetadata
    │       ├── interfile_genome.py        # genomexgenome_validation + GenomeXGenomeSettings
    │       └── interfile_read.py          # readxread_validation + ReadXReadSettings
    ├── validation_utils/
    │   ├── logger.py              # ValidationLogger (structlog-based); setup_logging()
    │   ├── ref_defragment.py      # defragment_reference() — unsupported workaround
    │   └── tests/
    │       └── test_ref_defragment.py
    └── tests/
        ├── test_config_manager.py
        ├── test_genome_validator.py
        ├── test_read_validator.py
        ├── test_feature_validator.py
        ├── test_interfile_genome.py
        ├── test_interfile_genome_characterization.py
        ├── test_interfile_read.py
        ├── test_logger.py
        ├── test_report.py
        ├── test_exceptions.py
        ├── test_formats.py
        ├── test_file_handler.py
        ├── test_path_utils.py
        ├── test_path_sanitization.py
        └── test_settings.py
```

---

## Running the module

### Bash wrapper (normal use)
```bash
# from the repo root (default config)
./validation.sh

# custom config path
./validation.sh --config path/to/config.json
```
`validation.sh` runs with `set -euo pipefail`. All outputs go to `data/valid/`
directly; previous runs are preserved via auto-incremented filenames
(`validation_001.log`, `report_1.txt`, etc.). The script also creates a
timestamped sub-directory (`data/valid/run_YYYYMMDD_HHMMSS/`) and exports its
path as `$EFSA_VALIDATION_RUN_DIR`. `main.py` reads this env var and writes the
log, report, and all validated files there. When called directly (without the
shell wrapper), `main.py` creates its own timestamped sub-directory. Only
`validated_params.json` is always written to the static `data/valid/` path so
Nextflow can load it via a stable `-params-file` path.

### Directly
```bash
python3 ./modules/validation/main.py data/inputs/config.json
```

### Tests
```bash
cd modules/validation/validation-pkg/
pytest tests/
```

---

## Outputs

Per-run files go to `data/valid/run_YYYYMMDD_HHMMSS/`; only `validated_params.json` lands directly in `data/valid/`.

| File | Description |
|---|---|
| `run_*/validation.log` | Structured JSON log (structlog). Auto-incremented if exists. |
| `run_*/report.txt` | Human-readable validation report with per-file statistics. |
| `run_*/` validated genome/reads/feature files | Standardised copies of every input file. |
| `validated_params.json` | Nextflow `-params-file`; always at static `data/valid/` path. |

---

## High-level workflow (`main.py`)

```
1.  Parse sys.argv[1] as config_path
2.  setup_logging() → data/valid/validation.log
3.  ConfigManager.load(config_path) → Config          # returns 1 on failure
4.  Instantiate per-validator Settings objects
5.  Initialise fatal_errors list + register_required_failure / register_missing_output helpers
6.  validate_genome(ref_genome_config, ref_settings)              # required
      → register_missing_output on success, register_required_failure on ValidationError
7.  validate_genome(mod_genome_config, mod_settings)              # optional
8.  genomexgenome_validation(ref_res, mod_res, gxg_settings)      # only if both present
9.  validate_genome(ref_plasmid_config, plasmid_settings)         # optional
10. validate_genome(mod_plasmid_config, plasmid_settings)         # optional
11. validate_reads(reads_configs, reads_settings)                  # required
      → register_missing_output per read, register_required_failure on ValidationError
12. readxread_validation(reads_res, readxread_settings)            # only if reads validated
13. validate_feature(ref_feature_config, ref_feature_settings)    # optional
14. validate_feature(mod_feature_config, mod_feature_settings)    # optional
15. if fatal_errors: report.add_fatal_errors(fatal_errors)
16. report.flush(format='text')
17. nf_params.build_params(validation_results) → NextflowParams
18. nf_params.write_params(params, data/valid/validated_params.json)
19. sys.exit(0)   # always 0; errors surfaced through report and log
```

**Exit code is always 0.** Fatal errors (required inputs that failed) are
accumulated in `fatal_errors`, forwarded to the report via `add_fatal_errors()`,
and cause `overall_status` to show `✗ FAILED`. `main.py` returns 1 **only** if
config loading itself fails.

---

## Config file (`data/inputs/config.json`)

Full specification is in `validation-pkg/docs/CONFIG_GUIDE.md`.

### Minimal example
```json
{
  "ref_genome_filename":  {"filename": "ref.fasta.gz"},
  "reads": [
    {"filename": "sample_R1.fastq.gz", "ngs_type": "illumina"},
    {"filename": "sample_R2.fastq.gz", "ngs_type": "illumina"}
  ]
}
```

### Full example
```json
{
  "ref_genome_filename":  {"filename": "ref.fasta.gz"},
  "mod_genome_filename":  {"filename": "mod.fasta"},
  "ref_plasmid_filename": {"filename": "plasmid_ref.fasta"},
  "mod_plasmid_filename": {"filename": "plasmid_mod.fasta"},
  "reads": [
    {"filename": "illumina_R1.fastq.gz", "ngs_type": "illumina"},
    {"filename": "illumina_R2.fastq.gz", "ngs_type": "illumina"},
    {"filename": "ont_reads.fastq.gz",   "ngs_type": "ont"},
    {"filename": "pacbio_reads.bam",     "ngs_type": "pacbio"}
  ],
  "ref_feature_filename": {"filename": "ref.gff3"},
  "mod_feature_filename": {"filename": "mod.gff3"},
  "options": {
    "validation_level": "strict",
    "threads": 4,
    "type": "prokaryote"
  }
}
```

### `options` keys
| Key | Values | Default | Notes |
|---|---|---|---|
| `validation_level` | `"strict"`, `"trust"`, `"minimal"` | `"trust"` | Controls depth of validation. **Case-insensitive** — stored as uppercase internally. |
| `threads` | positive int or `null` | `null` (auto-detect) | Warns if > CPU cores or > 16 |
| `logging_level` | `"DEBUG"`, `"INFO"`, `"WARNING"`, `"ERROR"` | `"INFO"` | Console log level. **Case-insensitive.** |
| `type` | `"prokaryote"`, `"eukaryote"` | `"prokaryote"` | Organism type hint. **Case-insensitive** — stored as uppercase internally. |

---

## Configuration loading (`config_manager.py`)

### Key constants
```python
ALLOWED_GLOBAL_OPTIONS  = {'threads', 'validation_level', 'logging_level', 'type'}
ALLOWED_FILE_OPTIONS    = {'threads', 'validation_level'}
MAX_RECOMMENDED_THREADS = 16
DEFAULT_THREADS         = 8
```

### Config dataclasses

**`BaseValidatorConfig`**
```
filename    : str
filepath    : Path    # resolved absolute path
basename    : str     # filename without all extensions (e.g. "genome" from "genome.fasta.gz")
coding_type : CodingType   # GZIP / BZIP2 / NONE
output_dir  : Path    # always data/valid/
global_options : Dict[str, Any]
```

**`GenomeConfig`** (extends `BaseValidatorConfig`)
```
detected_format  : GenomeFormat   # FASTA or GENBANK
n_sequence_limit : Optional[int]  # default 5; forbidden for plasmid configs
```

**`ReadConfig`** (extends `BaseValidatorConfig`)
```
detected_format : ReadFormat    # FASTQ or BAM
ngs_type        : str           # "ILLUMINA", "ONT", "PACBIO" — uppercased and validated in __post_init__
```

**`FeatureConfig`** (extends `BaseValidatorConfig`)
```
detected_format : FeatureFormat   # GFF, GTF, or BED
```

**`Config`** — top-level container
```
ref_genome   : GenomeConfig          # required
reads        : List[ReadConfig]      # required, non-empty
mod_genome   : Optional[GenomeConfig]
ref_plasmid  : Optional[GenomeConfig]
mod_plasmid  : Optional[GenomeConfig]
ref_feature  : Optional[FeatureConfig]
mod_feature  : Optional[FeatureConfig]

Properties (from options dict):
  .threads         → int | None
  .logging_level   → str
  .validation_level → str
  .type            → str | None
```

### Option flow

```
config.json "options"
      ↓
_parse_options() → stores only explicitly set keys in config.options
      ↓
config.options passed as global_options to each file config
      ↓
BaseValidator.__init__() reads global_options.get("validation_level") → None if not set
      ↓
BaseValidator._initialize_defaults() uppercases validation_level; fills None → 'TRUST', None threads → DEFAULT_THREADS (8)
```

**Important:** defaults are not stored in `config.options` when an option is omitted from the JSON.
The `Config` properties (`config.validation_level`, `config.type`, etc.) return defaults via
`.get(key, default)`, but `config.options` itself may be an incomplete dict.

### File-level override

Individual file entries in `config.json` can override `threads` and `validation_level` only:

```json
{"filename": "ref.fasta", "validation_level": "trust", "threads": 4}
```

`_merge_options()` merges these into the file's `global_options`. Any other keys in a file entry
are logged as warnings and ignored. `ngs_type` in read entries is extracted before calling
`_parse_file_config` and passed via `extra_fields` to avoid this warning.

### `ConfigManager.load(config_path, cli_options=None, logger=None) → Config`
Static method. `logger` is optional — a `ValidationLogger` instance can be injected
(useful in tests to assert on logged calls); defaults to a fresh `ValidationLogger()`.

Full validation pipeline:
1. Read and JSON-parse the file
2. `_validate_required_fields()` — requires `ref_genome_filename` and non-empty `reads`
3. `_setup_output_directory()` — creates `config_dir.parent / "valid"`
4. `_parse_genome_configs()` — ref (required), mod/plasmids (optional)
5. `_parse_reads_configs()` — supports both `filename` and `directory` keys;
   **raises `ValueError` if a `directory` entry for `ngs_type="ont"` / `"ONT"` or
   `ngs_type="pacbio"` / `"PACBIO"` contains more than one file** — merge reads first
6. `_parse_feature_configs()` — ref_feature, mod_feature (optional)
7. `_parse_options()` — validates and normalises global options
8. Returns fully populated `Config`

Raises: `ValidationFileNotFoundError`, `ConfigurationError`, `ValueError`

---

## Validation modes

| Mode | Scope | Stats collected | Use case |
|---|---|---|---|
| `minimal` | Extension/compression check only | Input/output paths, elapsed time | Trust existing data, maximum speed |
| `trust` | Sample first 10 sequences/reads | Counts, basic per-file stats | Fast QC on reasonably clean data |
| `strict` (default) | Every sequence/read | Full stats: GC%, N50, per-read lengths | Default; comprehensive QC |

In `minimal` mode the input file is copied as-is if format and compression already
match expectations; otherwise only an extension/coding check is performed.

---

## Exception hierarchy (`exceptions.py`)

```
ValidationError                    # base
├── ConfigurationError             # bad config.json
├── FileNotFoundError              # missing input file
├── FileFormatError                # unrecognised or malformed format
│   ├── FastaFormatError
│   ├── GenBankFormatError
│   ├── BedFormatError
│   ├── GffFormatError
│   ├── FastqFormatError
│   └── BamFormatError
├── CompressionError               # decompression/compression failure
├── GenomeValidationError          # genome-specific validation failure
├── FeatureValidationError         # feature-specific validation failure
├── ReadValidationError            # read-specific validation failure
└── InterFileValidationError       # inter-file consistency failure
```

All exceptions live in `validation_pkg.exceptions` and are importable from
`validation_pkg`.

---

## Format & compression enums (`utils/formats.py`)

### `CodingType`
```
GZIP, BZIP2, NONE

.to_extension()  → ".gz" / ".bz2" / ""
.normalize(val)  → CodingType (accepts string, enum, or None)
```

### `GenomeFormat`
```
FASTA, GENBANK

.to_biopython()  → "fasta" / "genbank"
.to_extension()  → ".fasta" / ".genbank"

Aliases:  fa, fna, faa → FASTA
          gb, gbk, genbank → GENBANK
```

### `ReadFormat`
```
FASTQ, BAM

.to_biopython()  → "fastq" / "bam"
.to_extension()  → ".fastq" / ".bam"

Aliases:  fq → FASTQ
```

### `FeatureFormat`
```
GFF, GTF, BED

.to_biopython()  → "gff" / "gtf" / "bed"
.to_extension()  → ".gff" / ".gtf" / ".bed"

Aliases:  gff3 → GFF
          gff2 → GTF
```

---

## Base classes (`utils/base_settings.py`, `utils/base_validator.py`)

### `BaseSettings` (abstract base for all Settings and Metadata)
```python
.copy() → BaseSettings
.update(**kwargs) → BaseSettings       # returns new instance, validates field names
.to_dict() → Dict[str, Any]
.from_dict(data) → BaseSettings        # classmethod
```

### `BaseOutputMetadata` (dataclass, extends `BaseSettings`)
```
input_file         : str
output_file        : str
output_filename    : str
validation_level   : str
elapsed_time       : float

.format_common_fields(indent) → List[str]   # for report generation
```

### `BaseValidatorSettings` (dataclass, extends `BaseSettings`)
```
coding_type           : CodingType     # normalised on __post_init__
output_filename_suffix : Optional[str]
output_subdir_name    : Optional[str]
```

### `BaseValidator` (abstract)

**Constructor receives:** `config: Any`, `settings: Optional[BaseSettings]`, `logger: Optional[ValidationLogger]`
Extracts from config: `output_dir`, `validation_level`, `threads`, `input_path`.
Default threads = 8 if not specified.
Logger defaults to `ValidationLogger(name=self._validator_type)` — so genome/read/feature
validators log under `"genome"`, `"read"`, `"feature"` respectively.

**`run() → OutputMetadata`** — main entry point
```
build output path (_build_output_path)
start timer
if validation_level == "minimal":
    _handle_minimal_mode()       # copy if format/compression match
else:
    _run_validation()            # subclass implementation
stop timer, fill metadata, return
```

**Abstract interface** (must be implemented by every subclass):
```python
@property _validator_type → str           # "genome" / "read" / "feature"
@property OutputMetadata → Type
@property _output_format → str           # "fasta" / "fastq" / "gff"
@property _expected_format → Any
@property _get_validator_exception() → Type[Exception]
def _run_validation() → Path
def _write_output() → Path
def _fill_output_metadata(output_path: Path) → None
```

---

## GenomeValidator (`validators/genome_validator.py`)

### `GenomeValidator.Settings` (dataclass)
```
# Validation thresholds
allow_empty_sequences : bool = False
allow_empty_id        : bool = False

# Genome-level editing
is_plasmid                   : bool = False   # marks file as plasmid in report
plasmid_split                : bool = False   # split multi-plasmid FASTA to separate files
plasmids_to_one              : bool = False   # merge plasmids into one FASTA (ref default: True)
main_longest                 : bool = False   # sort sequences longest-first (ref default: True)
main_first                   : bool = False   # keep first sequence at top after sorting
replace_id_with              : Optional[str]  # replace all IDs with this prefix
replace_id_with_incremental  : Optional[str]  # replace with "chr1", "chr2", … (both default: "chr")
min_sequence_length          : int = 100      # filter sequences shorter than this

# From BaseValidatorSettings
coding_type            : CodingType
output_filename_suffix : Optional[str]   # ref / mod / plasmid
output_subdir_name     : Optional[str]
```

Raises `GenomeValidationError` if:
- `plasmid_split` and `plasmids_to_one` are both True
- `main_longest` and `main_first` are both True

### `GenomeOutputMetadata` (dataclass)
```
# Always present
input_file, output_file, validation_level, elapsed_time   (from BaseOutputMetadata)
num_sequences              : int
longest_sequence_length    : int
longest_sequence_id        : str
sequence_ids               : List[str]
sequence_lengths           : List[int]
num_sequences_filtered     : int          # count removed by min_sequence_length
plasmid_count              : int
plasmid_filenames          : List[str]
fragmented                 : bool         # True when sequence count > n_sequence_limit

# Strict mode only
total_genome_size  : int
gc_content         : float    # percentage
n50                : int

.format_statistics(indent, input_settings) → List[str]   # used by ValidationReport
```

### Internal processing order
```
_parse_file()        # BioPython parse → self.sequences
_validate_sequences() # duplicates, empty IDs, min_length filter
_apply_edits()       # plasmid split/merge, reorder, replace IDs
_write_output()      # FASTA + compression
_fill_output_metadata()
```

---

## ReadValidator (`validators/read_validator.py`)

### Constants
```python
TRUST_MODE_SAMPLE_SIZE      = 10
PARALLEL_CHUNK_MULTIPLIER   = 4
MIN_CHUNK_SIZE               = 1000
ILLUMINA_PAIRED_END_PATTERNS = [...]  # regexes to detect R1/R2 suffixes
```

### `ReadValidator.Settings` (dataclass)
```
# Validation
check_invalid_chars  : bool = False   # reject non-ATCGN bases
allow_empty_id       : bool = False
allow_duplicate_ids  : bool = True

# BAM handling
keep_bam   : bool = False    # keep original BAM alongside FASTQ output
ignore_bam : bool = False    # copy BAM to output dir without FASTQ conversion
                             # (downstream tools must handle BAM directly)
                             # NOTE: does NOT return None — returns real BAM copy path

# Output organisation
outdir_by_ngs_type : bool = False   # write to illumina/, ont/, pacbio/ subdirs

# From BaseValidatorSettings
coding_type : CodingType    # default gz (reads always stored compressed)
output_filename_suffix, output_subdir_name
```

### `ReadOutputMetadata` (dataclass)
```
# Always present
input_file, output_file, validation_level, elapsed_time
base_name                  : str     # filename without R1/R2 suffix
read_number                : int     # 1 or 2 (0 if not detected)
ngs_type                   : str     # "ILLUMINA", "ONT", "PACBIO" (uppercased from input)
illumina_pairing_detected  : str     # "illumina" if R1/R2 pattern matched; else ""
num_reads                  : int

# Strict mode only
n50                    : int
total_bases            : int
mean_read_length       : float
longest_read_length    : int
shortest_read_length   : int
```

### Paired-end detection
`_detect_illumina_pairing(filename)` iterates `ILLUMINA_PAIRED_END_PATTERNS` looking
for R1/R2 markers. Returns `(base_name, read_number)`. Stored in metadata so
`readxread_validation` can pair files.

---

## FeatureValidator (`validators/feature_validator.py`)

### `FeatureValidator.Settings` (dataclass)
```
sort_by_position  : bool = True       # sort features by seqname + start
check_coordinates : bool = True       # validate start ≤ end, non-negative
replace_id_with   : Optional[str]     # replace seqname prefix ("chr")

# From BaseValidatorSettings
coding_type, output_filename_suffix, output_subdir_name
```

### `Feature` (dataclass — internal)
```
seqname, feature_type, score, strand, source, frame, attributes : str
start, end : int

.length → int   (end - start)
```

### `FeatureOutputMetadata` (dataclass)
```
input_file, output_file, validation_level, elapsed_time
num_features  : int
feature_types : List[str]    # unique feature types (gene, CDS, exon, …)
sequence_ids  : List[str]    # unique seqnames referenced
```

### Trust vs strict coordinate validation
- **Trust:** sample-validates first 10 features only
- **Strict:** validates every feature

### gffread fallback
`_parse_input()` first tries `gffread` for format normalisation. If `gffread` is
unavailable or exits non-zero, the validator **falls back to `_parse_gff()`**
(direct GFF3 line parser) on the already-decompressed input — it no longer
silently produces an empty output file. A `WARNING` validation issue is recorded.
If the file is genuinely unparseable, both paths return 0 features.

---

## Inter-file validators

### `genomexgenome_validation` (`validators/interfile_genome.py`)

**Settings:**
```python
@dataclass
class GenomeXGenomeSettings:
    same_number_of_sequences : bool = True
    same_sequence_ids        : bool = False
    same_sequence_lengths    : bool = False   # requires same_sequence_ids=True
    characterize             : bool = True    # run minimap2 alignment
```

**Function signature:**
```python
genomexgenome_validation(
    ref_genome_result  : GenomeOutputMetadata,
    mod_genome_result  : GenomeOutputMetadata,
    settings           : GenomeXGenomeSettings
) → Dict[str, Any]
```

**Return value:**
```python
{
    'passed': bool,
    'warnings': List[str],
    'errors': List[str],      # if non-empty, also raises InterFileValidationError
    'metadata': {
        'ref_num_sequences'  : int,
        'mod_num_sequences'  : int,
        'common_sequence_ids': List[str],
        'ref_only_ids'       : List[str],
        'mod_only_ids'       : List[str],
        'length_mismatches'  : {
            seq_id: {'ref_length': int, 'mod_length': int, 'difference': int}
        },
        'contigs_found'  : bool,
        'plasmids_found' : bool,
        'contig_files'   : List[str],    # minimap2 alignment output paths
        'plasmid_file'   : Optional[str]
    }
}
```

**Called by `main.py` only when:**
- Both `ref_genome_res` and `mod_genome_res` are not None
- `mod_genome_res.fragmented` is False

### `readxread_validation` (`validators/interfile_read.py`)

**Settings:**
```python
@dataclass
class ReadXReadSettings:
    pair_end_basename : bool = True    # check R1/R2 pairing
    allow_missing_r1  : bool = False   # R2 without R1 → warning instead of error
```

**Function signature:**
```python
readxread_validation(
    reads_results : List[ReadOutputMetadata],
    settings      : ReadXReadSettings
) → Dict[str, Any]
```

**Return value:**
```python
{
    'passed': bool,
    'warnings': List[str],
    'errors': List[str],      # if non-empty, also raises ReadValidationError
    'metadata': {
        'pairs_checked'   : int,
        'complete_pairs'  : List[str],   # base names with R1+R2 present
        'missing_r1'      : List[str],   # R2 with no matching R1
        'duplicate_r1'    : List[str],   # multiple R1s for same base name
        'duplicate_r2'    : List[str]    # multiple R2s for same base name
    }
}
```

**Called by `main.py` only when:** `reads_res is not None`

---

## Logger (`validation_utils/logger.py`)

Not a singleton — each call to `ValidationLogger()` creates a new instance, but all
instances share the same underlying structlog/stdlib configuration set by `setup_logging`.

### Setup
```python
from validation_utils.logger import ValidationLogger, setup_logging

# Called once in main.py to configure handlers:
logger = setup_logging(console_level='DEBUG', log_file=Path('data/valid/validation.log'))

# Anywhere else — creates a new instance bound to the shared structlog config:
logger = ValidationLogger()              # name defaults to "validation"
logger = ValidationLogger(name="genome") # used by BaseValidator subclasses
```

### Logger name
`ValidationLogger(name=...)` sets the structlog logger name (appears in log records).
The stdlib routing logger is always `"validation"` regardless of the instance name.
`BaseValidator` passes `self._validator_type` as the name automatically.

### Log methods
```python
logger.debug(message, **context)
logger.info(message, **context)
logger.warning(message, **context)   # also appends to validation_issues
logger.error(message, **context)     # also appends to validation_issues
logger.critical(message, **context)  # also appends to validation_issues
logger.add_validation_issue(level, category, message, details)
```

### Console vs file output
- **File** (`validation.log`): JSON, includes tracebacks (`format_exc_info`)
- **Console**: human-readable coloured output, tracebacks suppressed

### Timers
```python
logger.start_timer('genome_ref')
elapsed = logger.stop_timer('genome_ref')   # returns float seconds
logger.add_file_timing(input_file, validator_type, elapsed)
logger.display_file_timings_summary()
```

### Thread safety
All shared state (`validation_issues`, `_timers`, `file_timings`) is protected by
`threading.Lock`.

### Log file auto-increment
If `validation.log` already exists it creates `validation_001.log`, `validation_002.log`,
and so on — the old log is never overwritten.

---

## ValidationReport (`report.py`)

### Usage in `main.py`
```python
report = ValidationReport(output_dir / "report.txt")
report.write(ref_genome_res, file_type="genome")
report.write(reads_res, file_type="read")           # accepts single or list
report.write(genomexgenome_res, file_type="genomexgenome")
report.write(readxread_res, file_type="readxread")
if fatal_errors:
    report.add_fatal_errors(fatal_errors)           # propagate required-input failures
report.flush(format='text')   # writes report.txt; also accepts 'json'
```

### `add_fatal_errors(errors: List[str]) → None`
Appends file-level fatal error messages (required validators that raised or
produced no output). These cause `overall_status` to be `✗ FAILED` and are
printed as a **Fatal File-Level Errors** block in the text report and included
as `"fatal_errors"` in the JSON summary.

### `overall_status` logic
```python
overall_status = "✓ PASSED" if interfile_failed == 0 and not self.fatal_errors else "✗ FAILED"
```

### Internal record types
- **`FileValidationRecord`** — wraps a single validator's `OutputMetadata`
  - `validator_type`: "genome", "read", or "feature"
  - Renders: title, common fields (input/output/time), statistics, non-default settings
- **`InterFileValidationRecord`** — wraps an inter-validator result dict
  - `validation_type`: "genomexgenome" or "readxread"
  - Renders: status (PASSED/FAILED), errors, warnings, metadata

### Report auto-increment
Like the log file, `report_1.txt`, `report_2.txt` are created if `report.txt` exists.

---

## Utility functions

### Path security (`utils/path_utils.py`)
All path operations go through these functions to prevent directory traversal.

```python
resolve_filepath(base_dir, filename) → Path
    # Resolves relative path; raises ConfigurationError if outside base_dir

sanitize_path_component(component, allow_slashes=False) → str
    # Blocks: "..", absolute paths, null bytes, Windows reserved names
    # Raises ValueError on violation

build_safe_output_dir(base_dir, subdir_name=None, create=False) → Path
    # sanitize + resolve + optional mkdir

strip_all_extensions(filename, path=None) → str
    # "genome.fasta.gz" → "genome"

get_incremented_path(path, separator="_") → Path
    # "report.txt" → "report_1.txt" if "report.txt" exists
```

### File handling (`utils/file_handler.py`)
```python
detect_compression_type(filepath) → CodingType
    # from extension: .gz, .bz2, or NONE

detect_file_format(filepath, format_enum) → Enum
    # from extension, using enum aliases

open_file_with_coding_type(filepath, coding_type, mode) → TextIO
    # transparent decompression on open

open_compressed_writer(filepath, coding_type, threads) → TextIO
    # transparent compression on write

check_tool_available(tool_name) → bool
    # thread-safe, cached check via shutil.which

get_compression_command(coding_type, mode, threads) → Tuple[str, List[str]]
    # prefers pigz/pbzip2 over gzip/bzip2 when available

convert_file_compression(src, src_coding, dst, dst_coding, threads)
    # convenience: gz_to_bz2, bz2_to_gz, none_to_gz, gz_to_none, …

parse_config_file_value(value, field_name) → Tuple[str, Dict]
    # normalises config entry (plain string or {"filename": …, …} dict)
```

### Sequence statistics (`utils/sequence_stats.py`)
```python
calculate_n50(lengths: List[int]) → int
    # Sort descending, cumulative sum until ≥ 50% of total length.
    # Returns the length at which cumulative sum crosses the 50% threshold.
```

---

## `nextflow_params_handler.py`

### `NextflowParams` dataclass
```
# general_options (always written to JSON)
run_ref_x_mod     : bool = False   # both genomes validated + not fragmented + gxg passed
run_truvari       : bool = False   # reserved; always False
run_illumina      : bool = False   # illumina reads validated
run_nanopore      : bool = False   # ont reads validated
run_pacbio        : bool = False   # pacbio reads validated
contig_file_size  : int  = 0       # len(gxg_metadata['contig_files'])
run_vcf_annotation: bool = False   # ref_feature present

# input_output_options
nanopore_fastq      : Optional[str] = None   # always written (may be null)
ref_fasta_validated : Optional[str] = None   # omitted from JSON when None
mod_fasta_validated : Optional[str] = None   # omitted from JSON when None
pacbio_fastq        : Optional[str] = None   # omitted from JSON when None
gff                 : Optional[str] = None   # omitted from JSON when None

.to_dict() → dict
```

### `_path(meta) → Optional[str]` (internal)
Safe path extractor — **never returns the string `"None"`**:
```python
def _path(meta):
    if meta is None:
        return None
    value = getattr(meta, "output_file", None)
    if value is None:
        return None
    value = str(value).strip()
    return value if value and value != "None" else None
```

### `build_params(validation_results: dict) → NextflowParams`
Expected keys in `validation_results`:
```
ref_genome    : GenomeOutputMetadata | None
mod_genome    : GenomeOutputMetadata | None
genomexgenome : dict (from genomexgenome_validation) | None
reads         : List[ReadOutputMetadata] | None
ref_feature   : FeatureOutputMetadata | None
```

Before building flags, reads are pre-filtered:
```python
reads = [r for r in (validation_results.get("reads") or [])
         if r is not None and _path(r) is not None]
```

`run_ref_x_mod` is True only if:
- `ref_path` and `mod_path` are both non-None (and non-`"None"`)
- Neither genome has `fragmented=True`
- `gxg.get("passed", False)` is True

### `write_params(params, path: Path) → None`
Writes JSON with 2-space indent, UTF-8, `ensure_ascii=False`.

---

## Settings used in `main.py` (reference)

```python
# Reference genome — required
ref_genome_settings = GenomeValidator.Settings(
    plasmids_to_one=True,
    main_longest=True,
    coding_type=None,
    output_filename_suffix='ref',
    replace_id_with_incremental='chr',
    min_sequence_length=100
)

# Modified genome — optional
mod_genome_settings = GenomeValidator.Settings(
    plasmids_to_one=False,
    coding_type=None,
    output_filename_suffix='mod',
    replace_id_with_incremental='chr',
    min_sequence_length=100
)

# Plasmid genomes
plasmid_settings = GenomeValidator.Settings(
    is_plasmid=True,
    plasmids_to_one=True,
    coding_type=None,
    output_filename_suffix='plasmid'
)

# Reads
reads_settings = ReadValidator.Settings(
    coding_type='gz',
    outdir_by_ngs_type=True
)

# Reference feature annotations
ref_feature_settings = FeatureValidator.Settings(
    sort_by_position=False,
    check_coordinates=False,
    replace_id_with='chr',
    coding_type=None,
    output_filename_suffix='ref'
)

# Modified feature annotations
mod_feature_settings = FeatureValidator.Settings(
    sort_by_position=False,
    check_coordinates=False,
    replace_id_with='chr',
    coding_type=None,
    output_filename_suffix='mod'
)

# Inter-genome
genomexgenome_settings = GenomeXGenomeSettings(
    characterize=True,
    same_sequence_ids=False,
    same_number_of_sequences=False
)

# Inter-read
readxread_settings = ReadXReadSettings()   # all defaults
```

---

## Functional API (`validation_pkg/__init__.py`)

Thin wrappers — prefer these over instantiating validators directly.

```python
validate_genome(genome_config, settings=None)   → GenomeOutputMetadata
validate_genomes(genome_configs, settings=None) → List[GenomeOutputMetadata]

validate_read(read_config, settings=None)       → ReadOutputMetadata
validate_reads(read_configs, settings=None)     → List[ReadOutputMetadata]

validate_feature(feature_config, settings=None) → FeatureOutputMetadata
validate_features(feature_configs, settings=None) → List[FeatureOutputMetadata]
```

---

## Adding a new validator — checklist

1. Create `validators/new_validator.py`
2. Define `NewOutputMetadata(BaseOutputMetadata)` with validator-specific fields
3. Define `NewValidator.Settings(BaseValidatorSettings)` — only fields that differ from defaults
4. Subclass `BaseValidator`; implement all abstract methods
5. Add format support in `utils/formats.py` if needed
6. Register new exception subclass in `exceptions.py` if needed
7. Export from `validation_pkg/__init__.py`
8. Add `write(result, file_type="new_type")` handling in `report.py`
9. Write tests in `validation-pkg/tests/test_new_validator.py`

---

## Package metadata

```
name    : validation_pkg
version : 0.1.0
author  : Dominika Bohuslavova
license : EUPL-1.2
python  : ≥ 3.10

dependencies:
  biopython        == 1.85
  numpy            == 2.2.6
  pysam            == 0.23.3
  structlog        == 25.4.0
  typing_extensions== 4.15.0

dev:
  pytest == 8.4.2
```

---

## Input scenarios

The supported input scenarios (single contig, fragmented assembly below/above limit,
multi-sequence reference) are documented in `docs/validation/OVERVIEW.md`.

---

## Known gaps / TODO

- See `TODO.md` for outstanding items
