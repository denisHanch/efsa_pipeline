# Pipeline Runtime Messages & Mapping Summary

During execution, the pipeline prints progress messages indicating which workflow is currently running and what type of reads are being processed.

## Runtime Status Messages

When the pipeline is running, you will see real-time messages like:

```text
ℹ️  Running pipeline: processing long-pacbio reads → mapping to the reference & modified fasta.

ℹ️  Running pipeline: processing long-ont reads → mapping to the reference & modified fasta.

ℹ️  Running pipeline: processing short reads → mapping to the reference & modified fasta.

```

These messages help track the execution order and confirm that all three pipelines are being executed as expected.

## Unmapped Reads Statistics

After mapping, the pipeline reports the number and percentage of **unmapped reads** for each analysis.
This is useful for assessing mapping efficiency and data quality.

### Example Output

```text
📊 short-mod mapping:
    Unmapped reads: 19,880 (2.06%)
    Total input reads: 963,427

📊 short-ref mapping:
    Unmapped reads: 25,360 (2.63%)
    Total input reads: 963,362

📊 short-ref-plasmid mapping against plasmid:
    Unmapped reads: 4,677 (0.49%)
    Total input reads: 963,362

📊 ont/long-ref mapping:
    Unmapped reads: 52,745 (3.04%)
    Total input reads: 1,732,734

📊 ont/long-mod mapping:
    Unmapped reads: 48,145 (2.64%)
    Total input reads: 1,825,876

📊 pacbio/long-mod mapping:
    Unmapped reads: 41,596 (2.28%)
    Total input reads: 1,826,736

📊 pacbio/long-ref mapping:
    Unmapped reads: 47,472 (2.74%)
    Total input reads: 1,733,973
```

### Interpretation

**Unmapped reads** represent sequences that did not align to the provided reference or modified FASTA files.

A low percentage of unmapped reads indicates:

- High mapping quality
- Good reference/assembly quality
- Low contamination or sequencing errors

If the percentage of unmapped reads is unusually high, this may indicate:

- Poor read quality
- Inadequate or incomplete reference
- Contamination
- Incorrect input file selection

## Pipeline Execution Summary

The Nextflow pipelines ran successfully and produced the expected outputs. Each step completed without errors:

```text
✅ The ref_mod processing pipeline completed successfully.

✅ The long-read processing pipeline completed successfully.

✅ The short-read processing pipeline completed successfully.

✅ The execution of main.nf processing pipeline completed successfully.

📋 Process execution manifest: data/outputs/logs/process_manifest.txt
```

At the end of every run (success or failure), a **process execution manifest** is generated at `data/outputs/logs/process_manifest.txt`. This file lists every process that ran, its status (COMPLETED/FAILED/CACHED), and its exit code.

## Error Handling

When a failure occurs, the pipeline stops immediately, copies all process logs to `data/outputs/logs/`, and generates a `process_manifest.txt` with exit codes. See [Logs](../outputs/logs.md) for details on log files and the process manifest format.

```text
❌ The execution of main.nf processing pipeline failed: ...
Check the process execution manifest in data/outputs/logs/process_manifest.txt for details on which processes failed.
```

## Removal of Nextflow Work Directory

When the pipeline is executed with the option `--clean_work true` or when the parameter clean_work is set to true in the nextflow.config:

```text
params {
    clean_work = true
}
```

Nextflow automatically removes the temporary `work/` directory after successful completion and logs a message to confirm this.

```text
ℹ️ Nextflow `work/` directory was removed.
```

**Notes:**

- The `work/` directory contains intermediate files and temporary outputs generated during pipeline execution.
- Removing it saves disk space while retaining all final results in the `out_dir`.
- If you want to keep intermediate files for debugging or inspection, set: `clean_work = false` in the params section of the nextflow.config or use `--clean_work false` when running the pipeline.

## See Also

- [Running the Pipeline](running-pipeline.md) - Pipeline execution commands
- [Output Directory Structures](directory-structures.md) - Where to find results