# Validation Module — CLAUDE.md

Comprehensive reference for the `modules/validation/` directory. This module
validates and standardises genomic input files (FASTA/GenBank genomes, FASTQ/BAM
reads) before they enter the Nextflow pipeline. At the end
it emits `validated_params.json`, which Nextflow consumes via `-params-file`.

---

## Directory layout

```
modules/validation/
├── main.py                        # CLI entry point — orchestrates full workflow
├── validation.sh                  # Bash wrapper with default paths and safety clears
├── README.md                      # User-facing how-to and config guide pointer
├── utils/
│   ├── nextflow_params_handler.py # Builds and serialises validated_params.json
│   └── ref_defragment.py          # Unsupported workaround: merge fragmented reference
└── validation-pkg/
    ├── setup.py
    ├── requirements.txt           # biopython, numpy, pysam, structlog, typing_extensions
    ├── requirements-dev.txt       # pytest
    ├── validation_pkg/
    │   ├── __init__.py            # Public exports (validator classes, inter-file functions, logging); no functional-API wrappers
    │   ├── config_manager.py      # JSON config loader → Config / *Config dataclasses
    │   ├── exceptions.py          # Exception hierarchy rooted at ValidationError
    │   ├── utils/
    │   │   ├── base_settings.py   # BaseSettings, BaseOutputMetadata, BaseValidatorSettings
    │   │   ├── base_validator.py  # BaseValidator abstract class
    │   │   ├── file_handler.py    # Compression, format detection, file I/O utilities
    │   │   ├── formats.py         # Enums: CodingType, GenomeFormat, ReadFormat, OrganismType, ValidationLevel, LoggingLevel, NgsType
    │   │   ├── logger.py          # ValidationLogger singleton (structlog-based)
    │   │   └── path_utils.py      # Path resolution + path-traversal security
    │   └── validators/
    │       ├── genome_validator.py        # GenomeValidator + GenomeOutputMetadata
    │       ├── read_validator.py          # ReadValidator + ReadOutputMetadata
    │       ├── interfile_genome.py        # genomexgenome_validation + GenomeXGenomeSettings
    │       └── interfile_read.py          # readxread_validation + ReadXReadSettings
    └── tests/                     # not tracked in git (untracked locally)
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
`validation.sh` runs with `set -euo pipefail`. Each run writes to its own
timestamped sub-directory (`data/valid/run_YYYYMMDD_HHMMSS/`); previous runs
are **preserved**. The sub-directory path is exported as
`$VALIDATION_RUN_DIR` for downstream steps.

All options (threads, validation level, logging level, organism type) are set
exclusively via `config.json` — there are no CLI flags for them.

### Directly
```bash
python3 ./modules/validation/main.py data/inputs/config.json
```

### Tests
Test files live under `validation-pkg/tests/` and `utils/tests/` locally but are **not tracked in git**.
```bash
cd modules/validation/validation-pkg/
pytest tests/
```

---

## Outputs

| File | Location | Description |
|---|---|---|
| `validation_<run_id>.log` | `data/outputs/logs/` | Structured JSON log (structlog). Falls back to `validation.log` when no run ID. Auto-incremented if exists. |
| Validated genome/reads/feature files | `data/valid/run_YYYYMMDD_HHMMSS/` | Standardised copies of every input file. |
| `validated_params.json` | `data/valid/` | Nextflow `-params-file`; consumed by the main pipeline. |

---

## High-level workflow (`main.py`)

```
1.  Parse sys.argv → config_path
2.  setup_logging() → data/outputs/logs/validation_<run_id>.log
3.  ConfigManager.load(config_path) → Config   # returns 1 on failure
4.  Instantiate per-validator Settings objects
5.  GenomeValidator(ref_genome_config, ref_settings).run()        # required
      → logger.error on ValidationError or missing output_file
6.  GenomeValidator(mod_genome_config, mod_settings).run()        # optional
7.  genomexgenome_validation(ref_res, mod_res, gxg_settings)      # only if both present AND neither fragmented
8.  GenomeValidator(ref_plasmid_config, plasmid_settings).run()   # optional
9.  GenomeValidator(mod_plasmid_config, plasmid_settings).run()   # optional
10. [ReadValidator(rc, reads_settings).run() for rc in reads]     # required
      → logger.error on ValidationError or missing output_file
11. readxread_validation(reads_res, readxread_settings)            # only if reads validated
12. nf_params.build_params(validation_results) → NextflowParams
13. nf_params.write_params(params, data/valid/validated_params.json)
14. sys.exit(0)   # always 0; errors surfaced through log
```

**Exit code is always 0.** Validation failures are logged via `logger.error`.
`main.py` returns 1 **only** if config loading fails.

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
| `validation_level` | `"strict"`, `"trust"`, `"minimal"` | `"trust"` | Controls depth of validation |
| `threads` | positive int or `null` | `null` (auto-detect) | Warns if > CPU cores or > 16 |
| `logging_level` | `"DEBUG"`, `"INFO"`, `"WARNING"`, `"ERROR"` | `"INFO"` | Console log level |
| `type` | `"prokaryote"`, `"eukaryote"` | `"prokaryote"` | Organism type hint |

All options are set exclusively via `config.json` — there are no CLI equivalents.

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
ngs_type        : str           # "illumina", "ont", "pacbio" — validated in __post_init__
```

**`Config`** — top-level container
```
ref_genome   : GenomeConfig          # required
reads        : List[ReadConfig]      # required, non-empty
mod_genome   : Optional[GenomeConfig]
ref_plasmid  : Optional[GenomeConfig]
mod_plasmid  : Optional[GenomeConfig]

Properties (from options dict):
  .threads              → int | None
  .logging_level        → str         # default "INFO"
  .validation_level     → str         # default "trust"
  .type                 → str         # default "prokaryote"
```

### `ConfigManager.load(config_path: str, cli_options: Dict[str, Any] = None) → Config`
Static method. `cli_options` dict (from `_parse_cli_args`) is applied as fallback
for keys not set in `config.json`. Full validation pipeline:
1. Read and JSON-parse the file
2. `_parse_options(data, config, cli_options=cli_options)` — validates and
   normalises global options; config.json values override CLI values
3. `_validate_required_fields()` — requires `ref_genome_filename` and non-empty `reads`
4. `_setup_output_directory()` — creates `config_dir.parent / "outputs" / "valid"`
5. `_parse_genome_configs()` — ref (required), mod/plasmids (optional)
6. `_parse_reads_configs()` — supports both `filename` and `directory` keys;
   each file in the directory becomes a separate `ReadConfig` entry (all ngs types
   support multiple files per directory)
7. Returns fully populated `Config`

Raises: `ValidationFileNotFoundError`, `ConfigurationError`, `ValueError`

---

## Validation modes

| Mode | Scope | Use case |
|---|---|---|
| `minimal` | Extension/compression check only | Trust existing data, maximum speed |
| `trust` | Sample first 10 sequences/reads | Fast validation on reasonably clean data |
| `strict` (default) | Every sequence/read | Default; full validation + total genome size |

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
├── ReadValidationError            # read-specific validation failure
└── InterFileValidationError       # inter-file consistency failure
```

All exceptions live in `validation_pkg.exceptions` and are importable from
`validation_pkg`.

---

## Format & compression enums (`utils/formats.py`)

All enums expose a `normalize(value)` classmethod that accepts strings (case-insensitive), the enum itself, or `None` (returns the default).

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

Aliases:  fa, fna → FASTA
          gb, gbk, genbank → GENBANK
```

### `ReadFormat`
```
FASTQ, BAM

.to_biopython()  → "fastq" / "bam"
.to_extension()  → ".fastq" / ".bam"

Aliases:  fq → FASTQ
```

### `OrganismType`
```
PROKARYOTE, EUKARYOTE

.normalize(val)  → OrganismType (default: PROKARYOTE)
.value           → "prokaryote" / "eukaryote"
```

### `ValidationLevel`
```
STRICT, TRUST, MINIMAL

.normalize(val)  → ValidationLevel (default: TRUST)
.value           → "strict" / "trust" / "minimal"
```

### `LoggingLevel`
```
DEBUG, INFO, WARNING, ERROR

.normalize(val)  → LoggingLevel (default: INFO)
.value           → "DEBUG" / "INFO" / "WARNING" / "ERROR"
```

### `NgsType`
```
ILLUMINA, ONT, PACBIO

.normalize(val)  → NgsType (no default — raises ValueError if None)
.value           → "illumina" / "ont" / "pacbio"
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

**Constructor receives:** `config: Any`, `settings: Optional[BaseSettings]`
Extracts from config: `output_dir`, `validation_level`, `threads`, `input_path`.
Default threads = 8 if not specified.

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
sequence_ids               : List[str]
sequence_lengths           : List[int]   # critical: summed by build_params() for genome size
num_sequences_filtered     : int          # count removed by min_sequence_length
plasmid_count              : int
plasmid_filenames          : List[str]
plasmid_output_paths       : List[str]
fragmented                 : bool         # True when sequence count >= n_sequence_limit

# Strict mode only
total_genome_size  : int

```

### Internal processing order
```
_parse_file()           # BioPython parse → self.sequences
_validate_sequences()   # duplicates, empty IDs, min_length filter;
                        # if len(seqs) >= n_sequence_limit: copy file as-is,
                        #   call _deduplicate_copied_ids(), set fragmented=True, return early
_apply_edits()          # plasmid split/merge, reorder, replace IDs
_write_output()         # FASTA + compression
_fill_output_metadata()
```

### `_deduplicate_copied_ids()`
Called immediately after `copy_file()` in the fragmented early-exit path.
Renames duplicate sequence IDs in-place in the output file: second occurrence of
`id` becomes `id_1`, third becomes `id_2`, etc. Logs a WARNING for each renamed
pair. Updates `self.sequences` to the deduplicated list.

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
allow_empty_id       : bool = False

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
ngs_type                   : str     # "illumina", "ont", "pacbio"
input_format               : str     # "fastq" or "bam"
illumina_pairing_detected  : str     # "illumina" if R1/R2 pattern matched; else ""
```

### Paired-end detection
`_detect_illumina_pairing(filename)` iterates `ILLUMINA_PAIRED_END_PATTERNS` looking
for R1/R2 markers. Returns `(base_name, read_number)`. Stored in metadata so
`readxread_validation` can pair files.

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
- `ref_genome_res.fragmented` is False **and** `mod_genome_res.fragmented` is False

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

## Logger (`utils/logger.py`)

Singleton pattern. There is exactly one `ValidationLogger` instance throughout a run.

### Setup
```python
from validation_pkg import setup_logging, get_logger

logger = setup_logging(console_level='DEBUG', log_file=Path('data/outputs/logs/validation_<run_id>.log'))
logger = get_logger()   # retrieve singleton anywhere
```

### Log methods
```python
logger.debug(message, **context)
logger.info(message, **context)
logger.warning(message, **context)   # also appends to validation_issues
logger.error(message, **context)     # also appends to validation_issues
logger.critical(message, **context)  # also appends to validation_issues
logger.add_validation_issue(level, category, message, details)
```

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
If `validation.log` already exists it creates `validation_1.log`, `validation_2.log`,
and so on — the old log is never overwritten.

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

---

## `utils/nextflow_params_handler.py`

### `NextflowParams` dataclass
```
# general_options (always written to JSON)
run_ref_x_mod     : bool = False   # both genomes validated + not fragmented + gxg passed
run_illumina      : bool = False   # illumina reads validated
run_nanopore      : bool = False   # ont reads validated
run_pacbio        : bool = False   # pacbio reads validated
contig_file_size  : int  = 0       # len(gxg_metadata['contig_files'])
validation_timestamp : str = ""    # YYYYMMDD_HHMMSS
ref_genome_size_bp : Optional[int] = None  # validated reference genome size, bp
mod_genome_size_bp : Optional[int] = None  # validated modified genome size, bp

# input_output_options (omitted from JSON when None)
ref_fasta_validated : Optional[str] = None
mod_fasta_validated : Optional[str] = None
ref_plasmid_fasta   : Optional[str] = None
mod_plasmid_fasta   : Optional[str] = None

# input_output_options — lists (always present, may be empty)
illumina_fastqs : List[str]
ont_fastqs      : List[str]
ont_bams        : List[str]
pacbio_fastqs   : List[str]
pacbio_bams     : List[str]
contig_files    : List[str]

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

Genome-size params:
- `ref_genome_size_bp` and `mod_genome_size_bp` are derived from validator metadata.
- The builder prefers `total_genome_size` when available, otherwise it sums `sequence_lengths`.
- These values are consumed by `restructure_sv_tbl` to populate `pct_of_ref_genome` and `pct_of_mod_genome` without asking users to pass genome sizes manually.

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

---

## Adding a new validator — checklist

1. Create `validators/new_validator.py`
2. Define `NewOutputMetadata(BaseOutputMetadata)` with validator-specific fields
3. Define `NewValidator.Settings(BaseValidatorSettings)` — only fields that differ from defaults
4. Subclass `BaseValidator`; implement all abstract methods
5. Add format support in `utils/formats.py` if needed
6. Register new exception subclass in `exceptions.py` if needed
7. Export from `validation_pkg/__init__.py`
8. Wire result into `build_params()` in `nextflow_params_handler.py` if the output path is needed by Nextflow
9. Write tests in `validation-pkg/tests/test_new_validator.py`

---

## Bug fixes (changelog)

### `n_sequence_limit` off-by-one (`genome_validator.py`)
**Symptom:** ref.fa or mod.fa with exactly `n_sequence_limit` sequences (e.g. 5 with
default limit=5) was treated as NOT fragmented. `main_longest` ran on a multi-sequence
reference, splitting it into a chromosome + plasmid file instead of copying as-is
(scenario 4 per OVERVIEW.md).

**Fix:** changed `len(self.sequences) > n_sequence_limit` → `>= n_sequence_limit`.

**Semantics after fix:** sequences `< limit` → normal processing; sequences `>= limit`
→ `fragmented=True`, file copied as-is, `_deduplicate_copied_ids()` called.
Default `n_sequence_limit=5` means 1–4 sequences → normal, 5+ → fragmented.

---

### GXG called when ref was fragmented (`main.py`)
**Symptom:** `genomexgenome_validation` (minimap2 alignment) was called even when
`ref_genome_res.fragmented=True` — only `mod.fragmented` was checked.

**Fix:** condition now requires BOTH to be not fragmented:
```python
if (mod_genome_res is not None and ref_genome_res is not None
        and not getattr(mod_genome_res, 'fragmented', False)
        and not getattr(ref_genome_res, 'fragmented', False)):
```

---

### Duplicate sequence IDs in fragmented output (`genome_validator.py`)
**Symptom:** when a genome was copied as-is (fragmented path), duplicate IDs were
preserved in the output file, causing `samtools sort` to fail with "duplicate entry
in sam header".

**Fix:** `_deduplicate_copied_ids()` is called immediately after `copy_file()` in
both fragmented early-exit paths. Also applied in `interfile_genome.py`
(`_deduplicate_plasmid_ids()`) for plasmid output files.

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
