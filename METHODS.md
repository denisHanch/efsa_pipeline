# METHODS

This file is intentionally non-duplicative.

## Canonical Parameter Reference

Bioinformatics tool parameters, scientific rationale, and citations are documented in:

- `docs/nextflow/tool-parameters.md`

Use that file as the single source of truth for analysis thresholds and tool options.

## What This File Covers

This file only documents runtime resource policy and small-machine debugging.

### Resource Limit Policy

Resource governance has two layers:

1. Nextflow requests in `nextflow.config`:
   - `cpus`, `memory`, `time`
2. Tool thread usage in process scripts:
   - `task.cpus`

Rules:

- Keep finite `time` for all heavy processes to prevent indefinite hangs.
- Cap process CPUs via `Math.min(params.max_cpu as int, <cap>)`.
- Size memory/time by process class (light, mapping, sorting, SV calling).

### Small Machine Debugging (resource request too high)

If the host has less RAM than requested (for example 32 GB host vs 64/96 GB request), tasks usually fail to start or are OOM-killed (`exit 137`), and pipeline execution stops because `errorStrategy = 'terminate'`.

Quick checks:

1. `data/outputs/logs/trace.tsv`:
   - look at `status`, `exit`, `memory`, `peak_rss`, `realtime`
2. failing task logs:
   - `.command.err`, `.command.log` in that task `work/` directory

Typical error signatures:

- `Killed`
- `Out of memory`
- `Cannot allocate memory`
- exit code `137`

Recovery:

1. Run with lower global CPU pressure: `--max_cpu <small_value>`
2. Lower per-process `memory`/`time` requests in `nextflow.config` for the target machine
3. Resume: `nextflow run main.nf ... -resume`

Recommended operational setup:

- maintain separate profiles for low-resource and high-resource environments
  (for example `small`, `prokaryote`, `eukaryote`).
