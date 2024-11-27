# HG002 Mosaic Benchmark Set

This repository contains the methodology used to characterize baseline mosaic variants in the [Genome In A Bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle) reference cell line (HG002) using 
a trio-based approach and Strelka2 somatic caller to create a new HG002 mosaic benchmark.

## Methods Overview

1. Establish limit of detection (LOD) using GIAB reference mixtures
2. Generate potential mosaic list with tumor-normal callset, GIABv4.2.1 small variant germline benchmark, and GRCh38 stratifications
3. Database creation and filtering to identify mosaic candidates
4. Manual curation of candidates
5. External validation to refine draft benchmark
  
## Benchmark details and usage

HG002 mosaic benchmark v1.1 was generated using GRCh38 as reference and is available at the [GIAB ftp site](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/mosaic_v1.10/GRCh38/SNV/). 

Use of this GIAB somatic resource includes, but is not limited to:

  * benchmarking mosaic variant callers
  * negative controls for WGS somatic callers or targeted clinical sequencing in tumor-only mode
  * benchmarking somatic callers in tumor-normal mode using GIAB mixtures
  * dataset for germline researchers to filter low fraction somatic variants from their data
  * benchmarking some types of off-target edits

## Repository contents

READMEs, scripts, and figures from [*A robust benchmark for detecting low-frequency variants in the HG002 Genome In A Bottle NIST reference material*](URL) are contained herein. 

## Software requirements

All software information is located in [Supplementary Table 8](URL).

##Sharing/Access Information

Licenses/restrictions placed on the data, or limitations of reuse: Publicly released data are freely available for reuse without embargo.

##Links to publicly accessible locations of data

Info for all GIAB reference and orthogonal datasets are found at [Supplementary Table 1](URL).
GIAB Ashkenzai Jewish (AJ) trio WGS data (high coverage, Illumina 300X) are at the GIAB ftp links below. 

* HG002-son [fastq](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/)
* HG003-father [fastq](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/HG003_HiSeq300x_fastq/)
* HG004-mother [fastq](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/HG004_HiSeq300x_fastq/)

Reference
* GRCh38 [fasta](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz), [fai](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz.fai)

## General Information

Author Information

* Principal Investigator: Justin Zook, NIST, [jzook@nist.gov](mailto:jzook@nist.gov)
* Nate Olson, NIST,[nathanael.olson@nist.gov](mailto:nathanael.olson@nist.gov)
* Camille Daniels, MDIC,[cdaniels@mdic.org](mailto:cdaniels@mdic.org)
* Adetola Abdulkadir, MDIC, [cdaniels@mdic.org](mailto:aabdulkadir@mdic.org)


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
