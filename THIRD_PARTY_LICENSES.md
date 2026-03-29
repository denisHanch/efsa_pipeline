# Third-Party Licenses

This document lists all third-party tools, libraries, and software components used in the EFSA Pipeline, along with their respective licenses and source locations.

The EFSA Pipeline itself is licensed under the **European Union Public Licence v1.2 (EUPL-1.2)**.

---

## License Compatibility Summary

EUPL-1.2 is a copyleft license. Its Appendix explicitly lists compatible licenses for combined/derivative work distribution (GPL-2.0, GPL-3.0, AGPL-3.0, LGPL-2.1, LGPL-3.0, MPL-2.0, EPL-1.0, CeCILL-2.1, and others). Key notes:

| Compatibility Class | Licenses | Status |
|---|---|---|
| Permissive (no restrictions) | MIT, BSD-2-Clause, BSD-3-Clause, Apache-2.0, Artistic-2.0, PSF-2.0 | Compatible |
| Weak copyleft | LGPL-2.1, LGPL-3.0 | Compatible (EUPL-1.2 Appendix) |
| Strong copyleft | GPL-2.0, GPL-3.0 | Compatible (EUPL-1.2 Appendix) |
| Database / data licenses | Varies by SNPEff database source | Verify per database |

> **Important:** The pipeline invokes all bioinformatics tools as separate processes inside isolated Docker containers. No tool code is statically or dynamically linked into the pipeline source. This limits copyleft obligations to distribution of the containers themselves, not the pipeline scripts.

> **Governmental / Commercial Use:** All tools listed below use OSI-approved open-source licenses. None impose academic-only or non-commercial restrictions at the license level. However, some tools' **databases** (e.g., ClinVar, dbSNP used by SNPEff) may carry additional terms — verify before distribution.

---

## 1. Bioinformatics Tools (Containerized)

### Sequence Alignment & Mapping

| Tool | Version | License | Source URL |
|---|---|---|---|
| minimap2 | 2.30 | MIT | https://github.com/lh3/minimap2 |
| BWA | 0.7.17 | GPL-3.0-or-later | https://github.com/lh3/bwa |
| MUMmer (nucmer, delta-filter, show-coords) | 4.0.1 | Artistic-2.0 | https://github.com/mummer4/mummer |

### Structural Variant Calling

| Tool | Version | License | Source URL |
|---|---|---|---|
| CuteSV | 2.1.3 | MIT | https://github.com/tjiangHIT/cuteSV |
| Sniffles | 2.7.3 | MIT | https://github.com/fritzsedlazeck/Sniffles |
| DeBreak | 1.2 (commit-pinned) | MIT | https://github.com/Maggi-Chen/DeBreak |
| DELLY | 1.7.3 | BSD-3-Clause | https://github.com/dellytools/delly |
| SURVIVOR | commit 07404b74 | MIT | https://github.com/fritzsedlazeck/SURVIVOR |
| FreeBayes | 1.2.0 | MIT | https://github.com/freebayes/freebayes |

### Variant Processing & Annotation

| Tool | Version | License | Source URL |
|---|---|---|---|
| bcftools | 1.23 | MIT | https://github.com/samtools/bcftools |
| SAMtools | 1.23 | MIT | https://github.com/samtools/samtools |
| HTSlib (bgzip, tabix) | 1.23 | MIT | https://github.com/samtools/htslib |
| Picard | 3.4.0 | MIT | https://github.com/broadinstitute/picard |
| SnpEff | 4.1k | LGPL-3.0-only | https://github.com/pcingola/SnpEff |
| Truvari | 5.4.0 | MIT | https://github.com/ACEnglish/truvari |

### Genome Comparison

| Tool | Version | License | Source URL |
|---|---|---|---|
| SyRI | 1.7.1 | MIT | https://github.com/schneebergerlab/syri |
| gffread | 0.12.7 | MIT | https://github.com/gpertea/gffread |

### Quality Control

| Tool | Version | License | Source URL |
|---|---|---|---|
| FastQC | 0.11.9 | GPL-3.0-or-later | https://github.com/s-andrews/FastQC |
| MultiQC | 1.33 | GPL-3.0-or-later | https://github.com/MultiQC/MultiQC |
| NanoPlot | 1.46.2 | MIT | https://github.com/wdecoster/NanoPlot |

### Read Processing

| Tool | Version | License | Source URL |
|---|---|---|---|
| TrimGalore | 0.6.10 | GPL-3.0-or-later | https://github.com/FelixKrueger/TrimGalore |
| cutadapt | 4.9 | MIT | https://github.com/marcelm/cutadapt |

### Compression Utilities

| Tool | Version | License | Source URL |
|---|---|---|---|
| pbzip2 | 1.1.13 | BSD-3-Clause | http://compression.ca/pbzip2/ |
| pigz | system | zlib License | https://zlib.net/pigz/ |

---

## 2. Python Packages

### Runtime Dependencies

| Package | Version | License | Source URL |
|---|---|---|---|
| biopython | 1.85 | BSD-3-Clause AND Biopython License | https://biopython.org |
| numpy | 2.2.6 | BSD-3-Clause | https://numpy.org |
| pysam | 0.23.3 | MIT | https://github.com/pysam-developers/pysam |
| structlog | 25.4.0 | MIT OR Apache-2.0 | https://github.com/hynek/structlog |
| typing_extensions | 4.15.0 | PSF-2.0 | https://github.com/python/typing_extensions |

### Development / Test Dependencies

| Package | Version | License | Source URL |
|---|---|---|---|
| pytest | 8.4.2 | MIT | https://github.com/pytest-dev/pytest |

---

## 3. Workflow Runtime

| Component | Version | License | Source URL |
|---|---|---|---|
| Nextflow | 25.10.4 | Apache-2.0 | https://github.com/nextflow-io/nextflow |
| nf-schema (plugin) | 2.4.2 | Apache-2.0 | https://github.com/nextflow-io/nf-schema |
| Docker Engine | 29.3.1 | Apache-2.0 | https://github.com/moby/moby |
| Alpine Linux | 3.23 | Mixed (see note) | https://alpinelinux.org |
| OpenJDK | 17 (JRE) | GPL-2.0 with Classpath Exception | https://openjdk.org |

> **Note on Alpine Linux:** Alpine itself is MIT-licensed, but it packages third-party software under their respective licenses.

---

## 4. Documentation Tools

| Component | Version | License | Source URL |
|---|---|---|---|
| MkDocs Material | latest | MIT | https://github.com/squidfunk/mkdocs-material |
| mkdocs-mermaid2-plugin | 10.9.0 | MIT | https://github.com/fralau/mkdocs-mermaid2-plugin |

---

## 5. Container Base Images

| Image | Version / SHA256 | Source |
|---|---|---|
| staphb/bcftools | 1.23@sha256:22f4ddfd... | https://github.com/StaPH-B/docker-builds |
| staphb/samtools | 1.23@sha256:ed378537... | https://github.com/StaPH-B/docker-builds |
| staphb/htslib | 1.23@sha256:400e7c4e... | https://github.com/StaPH-B/docker-builds |
| staphb/mummer | 4.0.1@sha256:f4106644... | https://github.com/StaPH-B/docker-builds |
| staphb/nanoplot | 1.46.2@sha256:824dbbe1... | https://github.com/StaPH-B/docker-builds |
| staphb/multiqc | 1.33@sha256:4dbb26ba... | https://github.com/StaPH-B/docker-builds |
| biocontainers/bwa | v0.7.17_cv1@sha256:9479b73e... | https://biocontainers.pro |
| biocontainers/fastqc | v0.11.9_cv8@sha256:82e5fa41... | https://biocontainers.pro |
| biocontainers/freebayes | v1.2.0-2-deb_cv1@sha256:1ce0e08e... | https://biocontainers.pro |
| biocontainers/snpeff | v4.1k_cv3@sha256:200e8122... | https://biocontainers.pro |
| broadinstitute/picard | 3.4.0@sha256:dda3ad26... | https://hub.docker.com/r/broadinstitute/picard |
| dellytools/delly | v1.7.3@sha256:045aba10... | https://github.com/dellytools/delly |
| ecomolegmo/cutesv | v1.0.4@sha256:9c29e67b... | Internal build (based on cuteSV MIT) |
| ecomolegmo/debreak | v1.0.3@sha256:f2d8f9b5... | Internal build (based on DeBreak MIT) |
| ecomolegmo/minimap2 | v2.30@sha256:50d38b71... | Internal build (based on minimap2 MIT) |
| ecomolegmo/sniffles | v1.0.3@sha256:d3875e4e... | Internal build (based on Sniffles MIT) |
| ecomolegmo/survivor | v1.0.3@sha256:90263d6b... | Internal build (based on SURVIVOR MIT) |
| ecomolegmo/syri | v1.0.4@sha256:789d5cdf... | Internal build (based on SyRI MIT) |
| ecomolegmo/trimgalore | v1.0.1@sha256:802cc917... | Internal build (based on TrimGalore GPL-3.0) |
| ecomolegmo/truvari | v1.0.3@sha256:ae346a3b... | Internal build (based on Truvari MIT) |

---

## License Texts

Full license texts for each component are available at their respective source URLs. The licenses referenced in this document are:

- **MIT** — https://opensource.org/licenses/MIT
- **Apache-2.0** — https://www.apache.org/licenses/LICENSE-2.0
- **GPL-3.0** — https://www.gnu.org/licenses/gpl-3.0.html
- **GPL-2.0** — https://www.gnu.org/licenses/gpl-2.0.html
- **LGPL-3.0** — https://www.gnu.org/licenses/lgpl-3.0.html
- **BSD-3-Clause** — https://opensource.org/licenses/BSD-3-Clause
- **Artistic-2.0** — https://opensource.org/licenses/Artistic-2.0
- **PSF-2.0** — https://docs.python.org/3/license.html
- **EUPL-1.2** — https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12

---

*Last updated: 2026-03-29*
*Generated from Nextflow config, Dockerfiles, and requirements.txt in the efsa_pipeline repository.*
