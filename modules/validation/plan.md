# Plan: Refactor nextflow_params_handler.py to use a Dataclass

## Context
`build_params` currently returns a plain `dict`, making params opaque and hard to
introspect. The user wants a typed `NextflowParams` dataclass for IDE support,
explicit field definitions, and `repr`, while `write_params` (and JSON output)
must continue to emit a plain dict.

---

## Files to Modify
- `modules/validation/nextflow_params_handler.py` — main change
- `modules/validation/test_nextflow_params.py` — update assertions from dict to attribute access

---

## Implementation

### 1. `nextflow_params_handler.py`

Add imports:
```python
from dataclasses import dataclass, field
from typing import Optional
```

Define the dataclass **above** `build_params`:
```python
@dataclass
class NextflowParams:
    # validated_inputs
    ref_fasta_avail: bool = False
    mod_fasta_avail: bool = False
    # general_options
    run_ref_x_mod: bool = False
    run_syri: bool = False
    run_truvari: bool = False
    run_illumina: bool = False
    run_nanopore: bool = False
    run_pacbio: bool = False
    contig_file_size: int = 0
    run_vcf_annotation: bool = False
    # input_output_options — always present (may be None)
    nanopore_fastq: Optional[str] = None
    # input_output_options — omitted from JSON when None
    ref_fasta_validated: Optional[str] = None
    mod_fasta_validated: Optional[str] = None
    pacbio_fastq: Optional[str] = None
    gff: Optional[str] = None

    def to_dict(self) -> dict:
        """Produce a dict for Nextflow -params-file JSON.
        Optional path fields are excluded when None; nanopore_fastq is always included."""
        _OPTIONAL_PATHS = {"ref_fasta_validated", "mod_fasta_validated", "pacbio_fastq", "gff"}
        result = {
            "ref_fasta_avail": self.ref_fasta_avail,
            "mod_fasta_avail": self.mod_fasta_avail,
            "run_ref_x_mod": self.run_ref_x_mod,
            "run_syri": self.run_syri,
            "run_truvari": self.run_truvari,
            "run_illumina": self.run_illumina,
            "run_nanopore": self.run_nanopore,
            "run_pacbio": self.run_pacbio,
            "contig_file_size": self.contig_file_size,
            "run_vcf_annotation": self.run_vcf_annotation,
            "nanopore_fastq": self.nanopore_fastq,
        }
        for k in _OPTIONAL_PATHS:
            v = getattr(self, k)
            if v is not None:
                result[k] = v
        return result
```

Update `build_params` signature: `-> NextflowParams` and return a dataclass
instance using keyword arguments (identical logic, just constructing
`NextflowParams(...)` instead of a dict).

Update `write_params` signature: `params: NextflowParams` → call
`params.to_dict()` before `json.dump`.

### 2. `test_nextflow_params.py`

All assertions currently use `p["key"]` — replace with `p.key`.
Membership tests (`"key" in p`, `"key" not in p`) become `p.key is not None` /
`p.key is None`.

Specific replacements:
| Old | New |
|-----|-----|
| `p["run_syri"]` | `p.run_syri` |
| `p["run_ref_x_mod"]` | `p.run_ref_x_mod` |
| `p["ref_fasta_avail"]` | `p.ref_fasta_avail` |
| `p["ref_fasta_validated"]` | `p.ref_fasta_validated` |
| `"ref_fasta_validated" in p` | `p.ref_fasta_validated is not None` |
| `"mod_fasta_validated" not in p` | `p.mod_fasta_validated is None` |
| `"pacbio_fastq" not in p` | `p.pacbio_fastq is None` |
| `"gff" not in p` | `p.gff is None` |
| (all other `p["key"]`) | `p.key` |

For `TestWriteParams`, `write_params` now takes a `NextflowParams` — pass the
dataclass directly (no change to call sites). `data["run_syri"]` checks on the
loaded JSON remain valid. The `params["run_syri"]` reference in the test body
becomes `params.run_syri`.

---

## Verification
1. Run existing tests: `python -m pytest test_nextflow_params.py -v`
2. Confirm all tests pass without modification to business logic.
