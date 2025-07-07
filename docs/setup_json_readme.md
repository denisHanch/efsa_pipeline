# Documentation of `setup.json`

The `setup.json` file is a configuration file that is generated and updated during the input parsing process. It is used by the main analysis pipeline (based on Nextflow).

The primary purpose of this file is to store information about:
- Job status and identification
- File locations and extensions
- Used compression method
- Preprocessing steps applied to the input files  so they are suitable for entry into the main pipeline.

This file should be read-only to provide brief information about files and the data.

---

## Structure

The structure of the `setup.json` file is as follows:
``` json
{
    "status": "SUCCESS/ERROR_XY/WARNING", // Indicates whether the pipeline can proceed; detailed info is stored in the log file.
    "jobID": "timestamp", // Unique identifier to trace the job in the log files.
    "files": {
        "reference_fasta": {
            "fasta_filepath": "/absolute/path/to/Rtimestamp.fa.??", // Absolute path to the reformatted FASTA file, ready for use in the pipeline.
            "original_filepath": "/absolute/path/to/reference_genome.??", // Absolute path to the original reference file.
            "reformat_bool": true, // Indicates whether the original file was reformatted.
            "reformat_from": "gbk", // Format of the original file before reformatting.
            "compressed_bool": true, // Indicates whether the file is compressed.
            "compression_type": "gz", // Compression method used.
            "recompressed_bool": true, // Indicates whether recompression was necessary.
            "recompressed_from": "bz2/gz/xz/zip/zst/??" // Original compression format before recompression.
        },
        "reference_features": { // Null if no feature file is found.
            "gff_filepath": "/absolute/path/to/Rtimestamp.gff.??", // Absolute path to the GFF file (features).
            "original_filepath": "/absolute/path/to/reference_genome.??", // Absolute path to the original feature file.
            "reformat_bool": true, // Indicates whether the original file was reformatted.
            "reformat_from": "gtf/bed", // Format of the original file before reformatting.
            "info": {
                // Important metadata extracted from the GFF file.
            },
            "compressed_bool": true, // Indicates whether the file is compressed.
            "compression_type": "gz", // Compression method used.
            "recompressed_bool": true, // Indicates whether recompression was necessary.
            "recompressed_from": "bz2/gz/xz/zip/zst/??" // Original compression format before recompression.
        },
        "modification_fasta": {
            "fasta_filepath": "/absolute/path/to/Mtimestamp.fa.??", // Absolute path to the reformatted FASTA file, ready for use in the pipeline.
            "original_filepath": "/absolute/path/to/mod_genome.??", // Absolute path to the original modified genome file.
            "reformat_bool": true, // Indicates whether the original file was reformatted.
            "reformat_from": "gbk", // Format of the original file before reformatting.
            "compressed_bool": true, // Indicates whether the file is compressed.
            "compression_type": "gz", // Compression method used.
            "recompressed_bool": true, // Indicates whether recompression was necessary.
            "recompressed_from": "bz2/gz/xz/zip/zst/??" // Original compression format before recompression.
        },
        "modification_features": { // Null if no feature file is found; may be generated automatically.
            "gff_filepath": "/absolute/path/to/Mtimestamp.gff.??", // Absolute path to the GFF file (features).
            "original_filepath": "/absolute/path/to/mod_genome.??", // Absolute path to the original feature file.
            "reformat_bool": true, // Indicates whether the original file was reformatted.
            "reformat_from": "gtf/bed", // Format of the original file before reformatting.
            "info": {
                // Important metadata extracted from the GFF file.
            },
            "compressed_bool": true, // Indicates whether the file is compressed.
            "compression_type": "gz", // Compression method used.
            "recompressed_bool": true, // Indicates whether recompression was necessary.
            "recompressed_from": "bz2/gz/xz/zip/zst/??" // Original compression format before recompression.
        },
        "modification_reads": {
            "fastq_filepath": "/absolute/path/to/Mtimestamp.fastq.??", // Absolute path to the reformatted FASTQ file.
            "original_filepath": "/absolute/path/to/reads.??", // Absolute path to the original reads file.
            "reformat_bool": true, // Indicates whether the original file was reformatted.
            "reformat_from": "bam", // Format of the original file before reformatting.
            "compressed_bool": true, // Indicates whether the file is compressed.
            "compression_type": "gz", // Compression method used.
            "recompressed_bool": true, // Indicates whether recompression was necessary.
            "recompressed_from": "bz2/gz/xz/zip/zst/??" // Original compression format before recompression.
        }
    }
}
```
## Explanation

### Status
Serves as a quick indicator of whether preprocessing was successful. It can take one of the following values: **SUCCESS**, **WARNING**, or **ERROR_XY**, where `XY` is replaced by the corresponding error code from the list below.

If everything completes without issues, the status will be **SUCCESS**.

A **WARNING** indicates that one or more non-critical issues occurred, such as:
  - Feature extraction was incomplete, or the feature file was not found.
  - TODO

Warnings are informative only and do not prevent the pipeline from continuing. They serve as a notification to the user but are not considered fatal errors.

An **ERROR_XY** represents a critical failure during input parsing or processing. If an status of this type occurs, the pipeline should not proceed until the issue is reviewed and resolved.

The table below summarizes the error types:

| ERROR | Description |
|-------|-------------|
| 01 | Reference genome file was not recognized |
| 02 | Modified genome file was not recognized |
| 03 | Reads file was not recognized |
| 04 | Reformatting from format X to Y failed. Details: ... |
| 05 | Unknown conversion from format X to format Y |
| 06 | Conversion from format X to format Y was not successful. Details: ... |
| ?? | TODO |


---

### JobID

A generated timestamp used as a unique identifier throughout the pipeline. It is especially useful for referencing jobs in the log files.

---

## Compression

Compression must be validated and possibly reformatted depending on the tools used in the pipeline, as not all tools support all compression types. 

For this reason:
- Genome FASTA files are preferably left uncompressed.
- Read files should be standardized using the gzip compression method with the `.gz` extension.
