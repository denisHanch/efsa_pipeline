#   Tool LICENSEs summary

##   Mummer4
-   version used: 4.0.1
-   [latest release](https://github.com/mummer4/mummer/releases)
-   LICENSE: [artistic license 2.0](https://opensource.org/license/artistic-2-0), open source project free for all organizations, both for- and non-profit
-   last paper [here](https://doi.org/10.1371/journal.pcbi.1005944)
-   requirements:   fig2dev (3.2.3), gnuplot (4.0), xfig (3.2)


##  trimgalore
-   version used:   ??
-   [latest release](https://github.com/FelixKrueger/TrimGalore/releases/tag/0.6.10)
-   LICENSE: [GNU General Public License (GPL) Version 3](https://opensource.org/license/gpl-3-0)
-   requirements:   cutadapt, fastqc
-   no origin paper

##  FastQC (part of trimgalore)
-   version used:   ??
-   LICENSE: [GNU General Public License (GPL) Version 3](https://opensource.org/license/gpl-3-0)
-   requirements:   Picard tools, java > 1.8.x
-   no origin paper

##  Cutadapt (part of trimgalore)
-   version used: ??
-   [code](https://github.com/marcelm/cutadapt/)
-   LICENSE: [MIT](https://opensource.org/license/mit)
-   paper [here](DOI:10.14806/ej.17.1.200 )

##  picard tools
-   version used: ??
-   [latest release](https://github.com/broadinstitute/picard/releases/tag/3.4.0)
-   LICENSE: [MIT](https://opensource.org/license/mit)
-   not published in paper
-   requirements:   java > 1.8.x

``` bibtex
@misc{Picard2019toolkit,
  title = {Picard toolkit},
  year = {2019},
  publisher = {Broad Institute},
  journal = {Broad Institute, GitHub repository},
  howpublished = {\url{https://broadinstitute.github.io/picard/}}
}
```

##  bwa
-   version used: v0.7.17-r1188 **NEWER AVAILABLE**
-   [latest release](https://github.com/lh3/bwa/releases/tag/v0.7.19)
-   LICENSE: [GNU General Public License (GPL) Version 3](https://opensource.org/license/gpl-3-0)
-   paper on bwa-short [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC2705234/), on bwa package [here](https://arxiv.org/abs/1303.3997) **ONLY ON ARCHIVE** 
-   requirements:   BWT-SW under GPL, other are under MIT license

##  bwa-mem
-   version used: v0.7.17-r1188 **NEWER AVAILABLE**
-   [latest release](https://github.com/lh3/bwa/releases/tag/v0.7.19)
-   LICENSE: [GNU General Public License (GPL) Version 3](https://opensource.org/license/gpl-3-0)
-   paper [here](https://arxiv.org/abs/1303.3997) **ONLY ON ARCHIVE** - maybe we can think about this [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) with [paper](https://ieeexplore.ieee.org/document/8820962). 
-   requirements:   same as bwa

##  bcftools
-   version: 1.21
-   [latest release](https://github.com/samtools/bcftools/releases/tag/1.22)
-   LICENSE: [MIT](https://opensource.org/license/mit)
-   [paper](https://doi.org/10.1093/gigascience/giab008)
-   requirements:   htslib 

##  samtools
-   version: 1.19.2
-   [latest release](https://github.com/samtools/samtools/releases/tag/1.22)
-   LICENSE: [MIT](https://opensource.org/license/mit)
-   [paper](https://doi.org/10.1093/gigascience/giab008)
-   requirements:   htslib 

##   htslib
-   version: 1.19.1
-   [latest release](https://github.com/samtools/htslib/releases/tag/1.22)
-   LICENSE: [MIT](https://opensource.org/license/mit)
-   [paper](https://doi.org/10.1093/gigascience/giab007)
-   requirements:   zlib 

##  SnpEff
-   version: 1.19.1
-   [latest release](https://github.com/samtools/htslib/releases/tag/1.22)
-   LICENSE: [MIT](https://opensource.org/license/mit)
-   [paper](https://www.tandfonline.com/doi/full/10.4161/fly.19695)
-   requirements:   Java

``` bibtex
@article{cingolani2012program,
title={A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3},
author={Cingolani, P. and Platts, A. and Coon, M. and Nguyen, T. and Wang, L. and Land, S.J. and Lu, X. and Ruden, D.M.},
    journal={Fly},
    volume={6},
    number={2},
    pages={80-92},
    year={2012}
}
```

##  Freebayes
- version: ??
- [latest release](https://github.com/freebayes/freebayes/releases/tag/v1.3.10)
- LICENSE: [MIT](https://opensource.org/license/mit)
- [paper](https://arxiv.org/abs/1207.3907) **ONLY ARCHIVE**
- requirements: bc, samtools, parallel, meson, ninja-build, libvcflib-tools, vcftools (or vcflib) (can installed with all dependencies)

##  multiqc
-   version:  ??
-   [latest release](https://github.com/MultiQC/MultiQC/releases/tag/v1.30)
-   LICENSE: [GNU General Public License (GPL) Version 3](https://opensource.org/license/gpl-3-0)
-   [paper](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)

##  Delly
-   version:  ??
-   [latest release](https://github.com/dellytools/delly/releases/tag/v1.3.3)
-   LICENSE: [BSD 3-Clause License](https://opensource.org/license/bsd-3-clause)
-   [paper](https://academic.oup.com/bioinformatics/article/28/18/i333/245403?login=true)

##  plotsr
-   version:  ??
-   [latest release](https://github.com/schneebergerlab/plotsr/releases/tag/v1.1.0)
-   LICENSE: [MIT](https://opensource.org/license/mit)
-   [paper](https://academic.oup.com/bioinformatics/article/38/10/2922/6569079?login=true)

