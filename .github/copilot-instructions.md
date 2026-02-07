## EFSA Pipeline — Copilot instructions (concise)

- Big picture
  - The project is a Nextflow DSL2 pipeline orchestrating three main domains: short-read processing, long-read processing (PacBio/ONT) and reference-vs-modified genome comparison. Entry point: `main.nf` (DSL2). Workflows live under `workflows/` and reusable processes live under `modules/`.
  - Validation and preprocessing are done by the Python validation package in `modules/validation/` which writes cleaned inputs to `data/valid/`. The Nextflow pipeline reads from `data/valid/` and writes to `data/outputs/`.
  - Containerization is central: developers run the pipeline inside the project's Docker image (see `Dockerfile`, `docker-compose.yml`, `run_container.sh`). Nextflow container images are configured in `nextflow.config` under `params.containers`.

- Quick commands (real examples used in repo)
  - Start developer container (recommended): `./run_container.sh`
  - Validate inputs: `python3 ./modules/validation/main.py ./data/inputs/config.json` (creates `data/valid/`)
  - Run pipeline (resumable):
    `nextflow run main.nf -process.containerOptions "-u $(id -u):$(id -g)" --max_cpu $(nproc) -resume`
  - Run single workflow: `nextflow run workflows/long_read.nf -resume` (see `README.md` for examples)

- Important file references (where to look/change behavior)
  - `main.nf` — orchestrates which workflows run, input detection, and final aggregation.
  - `workflows/` — contains workflow entrypoints (e.g. `long_read.nf`, `short_read.nf`, `fasta_ref_x_mod.nf`). Use these to modify high-level data flow.
  - `modules/` — building blocks (mapping.nf, qc.nf, sv_calling.nf, logs.nf). Edit here to change process commands and channels.
  - `nextflow.config` — global params (containers map, default dirs, docker enable flag). Update `params.containers` to change tool images.
  - `modules/validation/` — Python validation package and tests; input config format in `data/inputs/config_template.json` and docs in `validation/CONFIG_GUIDE.md`.
  - `modules/utils/create_sv_output.py` — final SV table assembler called from `main.nf`.

- Project conventions and patterns (project-specific)
  - Data layout: `data/inputs/` (raw), `data/valid/` (validated), `data/outputs/` (results). Readers and workflows expect this structure.
  - Filenames: the validation config expects keys such as `ref_genome_filename`, `mod_genome_filename`, `ref_plasmid_filename`, and read subdirectories `illumina/`, `pacbio/`, `ont/` under `data/inputs`.
  - Centralized container map: tool images live in `nextflow.config` under `params.containers` — update versions there rather than inside many process definitions.
  - Nextflow DSL2 idioms: workflows use `take:`, `main:`, `emit:`, and `Channel` construction; prefer returning channels (`emit`) from workflows for downstream composition.
  - Error handling: `process.errorStrategy = 'ignore'` is set globally in `nextflow.config` — be careful when changing this as it affects pipeline robustness and retries.

- Integration points & external dependencies
  - Many processes run in Docker containers listed in `nextflow.config`. Ensure images exist and are accessible (registry or local build). Some images are `ecomolegmo/*` and `staphb/*`.
- Do not edit anything in the validation package (`modules/validation/`) — it is managed by a separate CI process.

- Debugging & development tips
  - Always run validation first (creates `data/valid/`). Skipping it will usually make Nextflow find no inputs.
  - Use `-resume` to continue after fixing tasks. Use `-with-report` and `-with-timeline` for diagnostics and performance analysis.
  - Run inside the provided container (`./run_container.sh`) so all required binaries and Python packages are available (the host may differ).
  - Check logs under `data/outputs/logs/` and validation logs under `logs/` for detailed messages.

- When changing behavior
  - To add/upgrade a tool, update `params.containers` in `nextflow.config` and, if necessary, update `Dockerfile` or rebuild the specific tool image.
  - To add a new high-level pipeline, create a new `workflows/*.nf` using DSL2 patterns and include it from `main.nf` (see how `long_read.nf` is included).

- Please review
  - Tell me if you'd like more detail about any module, or want me to merge additional content from the README or validation docs into this file.
