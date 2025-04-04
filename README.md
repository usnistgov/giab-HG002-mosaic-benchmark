# TODO

- clean up code to use mosaic database (or subset of it)
- document lod analysis and add code for figure
- document mosaic database generation
- revise figures: figure 3 - drop panel B? larger font?, figure 2 - remove/ improve annotations cleaner axis labels, fig 4 larger text and clearer boxplot; fig 1 simplify?

# HG002 Mosaic Benchmark Set

This repository contains the methodology used to characterize baseline mosaic variants in the [Genome In A Bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle) reference cell line (HG002) using 
a trio-based approach and a somatic caller to create a new HG002 mosaic benchmark. 

GIAB reference WGS data from HG002 (son) was passed to Strelka2 as the tumor data and the combined parents’ data (HG003-father and HG004-mother) was passed as normal.
Potential mosaics identified from the Illumina 300X tumor-normal callset were evaluated using high coverage BGI, Element, and Pacbio HiFi Revio orthogonal datasets.

## Methods Overview

1. Establish limit of detection (LOD) using GIAB reference mixtures
2. Generate potential mosaics list with [Illumina 300X tumor-normal callset](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/mosaic_v1.10/GRCh38/SNV/SupplementaryFiles/HG002_GRCh38_Strelka2-Ill300X.vcf.gz), [GIABv4.2.1 small variant germline benchmark](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/), and [GRCh38 stratifications](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh38/)
3. Database creation that includes orthogonal datasets and subsequent filtering to identify mosaic candidates
4. Manual curation of candidates
5. External validation to refine draft benchmark

<!--[Trio-based methods](https://github.com/usnistgov/giab-HG002-mosaic-benchmark/blob/main/figures/Figure_1_trio_based_methodology.png) -->

<img src="https://github.com/usnistgov/giab-HG002-mosaic-benchmark/blob/main/figures/Figure_1_trio_based_methodology.png" alt="Trio-based methods" width="700">
  
## Benchmark location and usage

HG002 mosaic benchmark v1.1 is available at the [GIAB ftp site](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/mosaic_v1.10/GRCh38/SNV/). 

Use of this GIAB somatic resource includes, but is not limited to:

  * benchmarking mosaic variant callers
  * negative controls for WGS somatic callers or targeted clinical sequencing in tumor-only mode
  * benchmarking somatic callers in tumor-normal mode using GIAB mixtures
  * dataset for germline researchers to filter low fraction somatic variants from their data
  * benchmarking some types of off-target edits

## Repository contents

READMEs, scripts, and figures from [*A robust benchmark for detecting low-frequency variants in the HG002 Genome In A Bottle NIST reference material*](https://doi.org/10.1101/2024.12.02.625685) are contained herein. 

## Software requirements

All software information is located in [Supplementary Table 8](https://www.biorxiv.org/content/10.1101/2024.12.02.625685v1.supplementary-material).

## Sharing/Access Information

Licenses/restrictions placed on the data, or limitations of reuse: Publicly released data are freely available for reuse without embargo.

## Links to publicly accessible locations of data

Info for all GIAB reference and orthogonal datasets are found at [Supplementary Table 1](https://www.biorxiv.org/content/10.1101/2024.12.02.625685v1.supplementary-material).
GIAB Ashkenzai Jewish (AJ) trio WGS data (high coverage, Illumina 300X) are at the GIAB ftp links below. 

* HG002 (son) [fastq](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/)
* HG003 (father) [fastq](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/HG003_HiSeq300x_fastq/)
* HG004 (mother) [fastq](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/HG004_HiSeq300x_fastq/)

Reference
* GRCh38 [fasta](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz), [fai](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz.fai)

## General Information

Authors

* Principal Investigator: Justin Zook, NIST, [jzook@nist.gov](mailto:jzook@nist.gov)
* Nate Olson, NIST, [nathanael.olson@nist.gov](mailto:nathanael.olson@nist.gov)
* Camille Daniels, MDIC, [cdaniels@mdic.org](mailto:cdaniels@mdic.org)
* Adetola Abdulkadir, MDIC, [aabdulkadir@mdic.org](mailto:aabdulkadir@mdic.org)


<!-- 
Information to include in the README
1. Software or Data description
   - Statements of purpose and maturity
   - Description of the repository contents
   - Technical installation instructions, including operating
     system or software dependencies
2. Contact information
   - PI name, NIST OU, Division, and Group names
   - Contact email address at NIST
   - Details of mailing lists, chatrooms, and discussion forums,
     where applicable
3. Related Material
   - URL for associated project on the NIST website or other Department
     of Commerce page, if available
   - References to user guides if stored outside of GitHub
4. Directions on appropriate citation with example text
5. References to any included non-public domain software modules,
   and additional license language if needed, *e.g.* [BSD][li-bsd],
   [GPL][li-gpl], or [MIT][li-mit]

The more detailed your README, the more likely our colleagues
around the world are to find it through a Web search. For general
advice on writing a helpful README, please review
[*Making Readmes Readable*][18f-guide] from 18F and Cornell's
[*Guide to Writing README-style Metadata*][cornell-meta].

[18f-guide]: https://github.com/18F/open-source-guide/blob/18f-pages/pages/making-readmes-readable.md
[cornell-meta]: https://data.research.cornell.edu/content/readme
[gh-cdo]: https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-code-owners
[gh-mdn]: https://github.github.com/gfm/
[gh-nst]: https://github.com/usnistgov
[gh-odi]: https://odiwiki.nist.gov/ODI/GitHub.html
[gh-osr]: https://github.com/usnistgov/opensource-repo/
[gh-ost]: https://github.com/orgs/usnistgov/teams/opensource-team
[gh-rob]: https://odiwiki.nist.gov/pub/ODI/GitHub/GHROB.pdf
[gh-tpl]: https://github.com/usnistgov/carpentries-development/discussions/3
[li-bsd]: https://opensource.org/licenses/bsd-license
[li-gpl]: https://opensource.org/licenses/gpl-license
[li-mit]: https://opensource.org/licenses/mit-license
[nist-code]: https://code.nist.gov
[nist-disclaimer]: https://www.nist.gov/open/license
[nist-s-1801-02]: https://inet.nist.gov/adlp/directives/review-data-intended-publication
[nist-open]: https://www.nist.gov/open/license#software
[wk-rdm]: https://en.wikipedia.org/wiki/README -->


## Other

- Installing and compiling Rust readcount parser following 
https://doc.rust-lang.org/stable/book/ch01-01-installation.html?highlight=MacOS

```bash
brew install rustup  
rustup default stable
## Adding to path
echo 'export PATH="/opt/homebrew/opt/rustup/bin:$PATH"' >> /Users/nolson/.zshrc
```

Created variant tsv using `cut -f1-4 data/igv_curated_variants_RM.tsv > curated_variants.tsv`
Modify rust to include usage, ref_count, alt_count, depth

Different bam readcount pipeline run input tsvs
Dec 20  2023 ortho_bam_NIST_RM_libs.tsv
giab_id	lib	path	ref
RM	HG002	illumina	/projects/camille_mosaic/strelka2_rerun-GRCh38/HG002.GRCh38.300x.bam	GRCh38
RM	HG002	element_standard	/projects/camille_mosaic/element/new_NIST_RM_libs/GRCh38-GIABv3_GAT-APP-C138_element.bam	GRCh38
RM	HG002	element_long	/projects/camille_mosaic/element/new_NIST_RM_libs/GRCh38-GIABv3_GAT-APP-C148_element.bam	GRCh38
RM	HG002	bgiseq	/projects/data/bgiseq-high-coverage/hg38_result/upload/alignment/HG002/HG002.bam	GRCh38
RM	HG002	pacbio_revio_hifi	/projects/camille_mosaic/pacbio_revio_hifi/new_NIST_RM_lib/HG002_PacBio-HiFi-Revio_GRCh38-GIABv3_20231031.bam	GRCh38
RM	HG002	pacbio_sequel_ccsA_10kb	/projects/camille_mosaic/pacbio_sequel_ccs/HG002.10kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam	GRCh38
RM	HG002	pacbio_sequel_ccsB_15kb	/projects/camille_mosaic/pacbio_sequel_ccs/HG002.15kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam	GRCh38

Jun 14  2023 ortho_bam_plus_illumina.tsv
giab_id	lib	path	ref
RM	HG002	illumina	/projects/camille_mosaic/strelka2_rerun-GRCh38/HG002.GRCh38.300x.bam	GRCh38
CO	HG002	element	/projects/camille_mosaic/element/HG002_PCR-free_2x150_100X.bam	GRCh38
RM	HG002	bgiseq	/projects/data/bgiseq-high-coverage/hg38_result/upload/alignment/HG002/HG002.bam	GRCh38
CO	HG002	pacbio_revio_hifiA	/projects/data/pacbio_revio-hifi_ajtrio/HG002.m84005_220827_014912_s1.GRCh38.bam	GRCh38
CO	HG002	pacbio_revio_hifiB	/projects/data/pacbio_revio-hifi_ajtrio/HG002.m84005_220919_232112_s2.GRCh38.bam	GRCh38
CO	HG002	pacbio_revio_hifiC	/projects/data/pacbio_revio-hifi_ajtrio/HG002.m84011_220902_175841_s1.GRCh38.bam	GRCh38
CO	HG002	pacbio_revio_hifiD	/projects/camille_mosaic/hprc_pacbio_revio-hifi/sorted_aligned_GRCh38_m84039_230117_233243_s1.hifi_reads.default.bam	GRCh38
RM	HG002	pacbio_sequel_ccsA_10kb	/projects/camille_mosaic/pacbio_sequel_ccs/HG002.10kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam	GRCh38
CO	HG002	pacbio_sequel_ccsB_11kb	/projects/camille_mosaic/pacbio_sequel_ccs/HG002_GRCh38.haplotag.10x.bam	GRCh38
CO	HG002	pacbio_sequel_ccsC_15kb	/projects/camille_mosaic/pacbio_sequel_ccs/HG002.15kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam	GRCh38
CO	HG002	pacbio_sequel_ccsD_15kb_20kb	/projects/camille_mosaic/pacbio_sequel_ccs/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam	GRCh38