# Nextflow Configuration

The file `nextflow.config` is the central configuration for the pipeline. It defines default parameters, per-process resource limits and container images, execution profiles, and reporting options.

## Parameters (`params`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `config_json` | `data/inputs/config.json` | Path to the input configuration JSON |
| `out_dir` | `data/outputs` | Output directory for results |
| `log_dir` | `data/outputs/logs` | Directory for pipeline logs and reports |
| `max_cpu` | `1` | Maximum CPUs available per process (override with `--max_cpu`) |
| `clean_work` | `true` | Remove Nextflow `work/` directory after a successful run |
| `run_truvari` | `false` | Enable cross-caller VCF comparison with Truvari |
| `ref_plasmid_fasta` | `null` | Optional reference plasmid FASTA |
| `mod_plasmid_fasta` | `null` | Optional modified plasmid FASTA |
| `gff` | `null` | GFF3 feature file (overridden by validated params at runtime) |
| `help` | `false` | Print help message and exit |

## Per-Process Configuration

Each bioinformatics tool has a `withName` block that specifies its **container image** and **resource limits** together:

```groovy
withName: minimap2 {
    container = 'ecomolegmo/minimap2:v2.30@sha256:...'
    cpus = Math.min(params.max_cpu as int, 48)
    memory = '64 GB'
    time = '24h'
}
```

**Key conventions:**

- **Container pinning**: All images use `@sha256:` digests for reproducibility — the exact image bytes are locked, not just the tag.
- **CPU scaling**: Uses `Math.min(params.max_cpu, N)` so processes scale to available hardware but never exceed a tool-specific cap.
- **Centralized updates**: To upgrade a tool version, change it in `nextflow.config` only — individual process definitions do not specify containers.

### Default Resource Limits

Processes without explicit resource overrides use:

| Resource | Default |
|----------|---------|
| CPUs | 1 |
| Memory | 8 GB |
| Time | 8h |

## Error Strategy

```groovy
errorStrategy = 'terminate'
```

The pipeline terminates on the first process failure. Change to `'finish'` to let running tasks complete before stopping, or `'ignore'` to skip failures (use with caution).

## Docker

Docker execution is enabled by default:

```groovy
docker {
    enabled = true
    runOptions = '--user $(id -u):$(id -g)'
}
```

Containers run as the current user to avoid root-owned output files.

## Profiles

| Profile | Description |
|---------|-------------|
| `standard` | Default — uses parameters as defined above |
| `test` | Sets `config_json` to `data/inputs/test/config.json` for the bundled test dataset |

Run with a profile:

```bash
nextflow run main.nf -profile test --max_cpu $(nproc)
```

## Reports and Tracing

Nextflow can generate execution diagnostics (report, timeline, trace) via command-line flags. See [Running the Pipeline — Nextflow Options](running-pipeline.md#nextflow-options) for the flags and [Logs](../outputs/logs.md) for descriptions of the generated files.

## Manifest

The `manifest` block declares pipeline metadata (name, version, license, Nextflow version requirement) used by Nextflow and registries.
