# Tool version pins

This file documents explicit version pins used in custom tool images to improve reproducibility.

## Python packages (pinned)

| Tool image | Package | Version |
|---|---|---|
| tools/cutesv/Dockerfile | cuteSV | 2.1.3 |
| tools/sniffles/Dockerfile | sniffles | 2.7.3 |
| tools/debreak/Dockerfile | pysam | 0.23.3 |
| tools/pbzip2/Dockerfile | pbzip2 | 1.1.13 |
| tools/minimap2/Dockerfile | minimap2 | 2.30 |
| tools/mosdepth/Dockerfile | mosdepth | 0.3.11 |
| tools/trimgalore/Dockerfile | trimgalore | 0.6.10 |
| tools/validation/Dockerfile | validation_pkg | 0.1.0 (local) |
| tools/validation/Dockerfile | structlog | 25.5.0 |
| tools/validation/Dockerfile | py3-pandas (Alpine) | 2.3.3-r0 |

## Git sources (pinned)

| File | Repository | Pinned ref | Resolved commit |
|---|---|---|---|
| tools/debreak/Dockerfile | https://github.com/Maggi-Chen/DeBreak.git | v1.2 | 21d9b4294089ccd9df0feec976239551f4359096 |
| tools/syri/Dockerfile | https://github.com/schneebergerlab/syri.git | v1.7.1 | 61e5aecbb506553c22f88649015013fad28fb05e |
| tools/survivor/Dockerfile | https://github.com/fritzsedlazeck/SURVIVOR.git | 1.0.7 | 07404b74d3fe42b20f1362f4d8625f284b9d536c |
| Dockerfile | https://github.com/gpertea/gffread.git | v0.12.7 | 5647f076f616c583ea32fd19629d549d70fc43ea |

## Notes

- Tag pins are significantly more reproducible than floating branch HEAD installs.
- For maximum immutability, consider replacing tag pins with commit SHA checkout in Docker builds.
- For Nextflow container references, prefer image digests (for example image@sha256:...) over mutable tags when publishing production runs.
