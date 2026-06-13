# Custom tool build version pins

This file documents explicit version pins used while building custom tool images. It is not a license inventory and does not duplicate the container digest table in `docs/project/third-party-licenses.md`.

## Package and binary pins

| Build file | Component | Version / ref | Pin source |
|---|---|---|
| tools/cutesv/Dockerfile | cuteSV | 2.1.3 | `pip wheel cuteSV==2.1.3` |
| tools/debreak/Dockerfile | DeBreak | v1.2 | `ARG DEBREAK_REF=v1.2` |
| tools/debreak/Dockerfile | minimap2 | 2.30 | `ARG MINIMAP2_VERSION=2.30` |
| tools/debreak/Dockerfile | pysam | 0.23.3 | `ARG PYSAM_VERSION=0.23.3` |
| tools/gffread/Dockerfile | gffread | v0.12.7 | `ARG GFFREAD_REF=v0.12.7` |
| tools/minimap2/Dockerfile | minimap2 | 2.30 | `ARG MINIMAP2_VERSION=2.30` |
| tools/mosdepth/Dockerfile | mosdepth | 0.3.11 | `ARG MOSDEPTH_VERSION=0.3.11` |
| tools/pbzip2/Dockerfile | pbzip2 | 1.1.13 | source tarball URL |
| tools/sniffles/Dockerfile | sniffles | 2.7.3 | `ARG SNIFFLES_VERSION=2.7.3` |
| tools/survivor/Dockerfile | SURVIVOR | 07404b74d3fe42b20f1362f4d8625f284b9d536c | `ARG SURVIVOR_REF=...` |
| tools/syri/Dockerfile | SyRI | v1.7.1 | `ARG SYRI_REF=v1.7.1` |
| tools/trimgalore/Dockerfile | cutadapt | 4.9 | `ARG CUTADAPT_VERSION=4.9` |
| tools/trimgalore/Dockerfile | TrimGalore | 0.6.10 | `ARG TRIMGALORE_VERSION=0.6.10` |
| tools/validation/Dockerfile | validation package | local source | `COPY modules/validation/ /tmp/validation/` |
| tools/validation/Dockerfile | structlog | 25.5.0 | scripts venv install |
| tools/validation/Dockerfile | py3-pandas (Alpine) | 2.3.3-r0 | `apk add py3-pandas=2.3.3-r0` |

## Validation image copied binaries

| Consuming build file | Copied image | Digest prefix | Copied binary |
|---|---|---|---|
| tools/validation/Dockerfile | ecomolegmo/pbzip2:v1.1.13 | `sha256:4a308661...` | `/usr/bin/pbzip2` |
| tools/validation/Dockerfile | ecomolegmo/minimap2:v2.30 | `sha256:50d38b71...` | `/usr/local/bin/minimap2` |
| tools/validation/Dockerfile | ecomolegmo/gffread:v0.12.7 | `sha256:dad98757...` | `/usr/local/bin/gffread` |

## Notes

- `tools/cutesv/Dockerfile` pins cuteSV itself, but builds `pysam` without an explicit version pin.
- `tools/syri/Dockerfile` installs SyRI dependencies (`pandas`, `numpy`, `pysam`, `scipy`, `cython`, `igraph`, `matplotlib`, `psutil`) without explicit version pins.
- Git tag refs are more reproducible than floating branch installs, but commit SHA checkouts are stronger when an upstream tag is not immutable.
- Runtime Nextflow container image tags and digests are tracked in `nextflow.config` and summarized in `docs/project/third-party-licenses.md`.
