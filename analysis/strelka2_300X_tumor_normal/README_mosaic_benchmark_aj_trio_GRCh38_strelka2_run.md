#  GIAB AJ trio Illumina 300X Strelka2 run - HG002 mosaic benchmark  

Purpose: Generate a somatic callset with AJ trio Illumina 300X bams, GRCh38 genome reference, and Strelka2 caller to determine potential mosaics for HG002 mosaic benchmark generation

## Pre-processing

AJ trio Illumina 300X NovoAlign GRCh38 bams and indexes were downloaded to NIST Linux from [GIAB ftp site](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/). Bam files are already sorted.

HG002 (son) [bam](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam) [bai](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam.bai)

HG003 (father) [bam](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.300x.bam) [bai](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.300x.bam.bai)

HG004 (mother) [bam](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.300x.bam) [bai](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.300x.bam.bai)

## Merge and index parental bams to create a normal bam 

HG003 and HG004 bams were merged via samtools 1.15

```
samtools merge -@ 14 -o samtools_merged_HG003_HG004_GRCh38.300X.bam HG003_GRCh38.300X.bam HG004_GRCh38.300X.bam

samtools index -@ 14 samtools_merged_HG003_HG004_GRCh38.300X.bam
```

## Get GRCh38 genome 

GRCh38 files were downloaded from [GIAB ftp site](https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/)
to NIST Linux

GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta.gz
GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta.fai
GRCh38.bed

#Split bed file by chromosome (1-22, X, and Y)

#samtools view -b GRCh38.bed chr# > chr#.bed

```
for chr in {1..22} X Y; do
    samtools view -h input.bed "chr${chr}" > "chr${chr}.bed"
done
```

GRCh38 fasta, index, and bed files were uploaded to DNANexus using mamba env DNAnexus_CLI (dpxy client)


## Strelka2 variant calling

Strelka2 app was run by Camille on DNAnexus web GUI. The workflow is designed to take tumor/normal bam, bai, and reference genome files as input and performs somatic variant calling. 

Inputs 
tumor bam: HG002.GRCh38.300x.bam (son)
normal bam: samtools_merged_HG003_HG004.300x.bam (combined parents)
reference: GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta.gz
targets bed: beds/GRCh38_chr_#.bed 

Output is sent to a DNAnexus directory and downloaded locally.

To create a single vcf, all snv and indel vcf.gz outputs for each chr were concatenated using bcftools v1.16. 

```
bcftools concat -f vcf_list.txt -o HG002_GRCh38.vcf.gz

tabix -p vcf HG002_GRCh38.vcf.gz

bcftools view -H HG002_GRCh38.vcf.gz | wc -l
```

1,273,474 total variants called


## Strelka2 callset filtering to remove normal variants

bcftools was used to filter the Strelka2 vcf to create a vcf with normal variants removed

Create vcf that just has the normals flagged

```
bcftools filter -s VARIANT_DETECTED_IN_NORMAL -e "FMT/TAR[0:0] > 5 & FILTER != 'PASS'" HG002_GRCh38.vcf.gz > HG002_GRCh38_normals_flagged.vcf.gz
```

Create vcf with the normals removed

```
bcftools filter -e "FILTER == 'VARIANT_DETECTED_IN_NORMAL'" HG002_GRCh38_normals_flagged.vcf.gz > HG002_GRCh38_normals_removed.vcf.gz

bcftools view -H HG002_GRCh38_normals_removed.vcf.gz | wc -l
```

1,120,703 variants 


## Sort, index, and reformat normals removed vcf for vcf eval run

```
bcftools sort -o HG002_GRCh38_normals_removed_sorted.vcf.gz HG002_GRCh38_normals_removed.vcf.gz

tabix HG002_GRCh38_normals_removed_sorted.vcf.gz

sh reformat_strelka2.sh HG002_GRCh38_normals_removed_sorted.vcf.gz <prefix(i.e. HG002) > HG002_reformatted.vcf.gz
```
cp 

## Compare strelka2 normals removed vcf to GIAB GRCh38 4.2.1 benchmark
## Run using which mamba env - ADD this info and amend directories to match this repo

```
rtg vcfeval -b ~/git_repos/mosaic-benchmarking/resources/giab_benchmark_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
-c ~/git_repos/mosaic-benchmarking/data/strelka2_wgs_parents-vs-child/HG002_reformatted.vcf.gz \
--evaluation-regions ~/git_repos/mosaic-benchmarking/resources/giab_benchmark_v4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
--squash-ploidy \
-o  ~/git_repos/mosaic-benchmarking/data/strelka2_wgs_parents-vs-child/vcfeval_comparison \
-t ~/git_repos/mosaic-benchmarking/resources/GRCh38_hs38d1.sdf \
--decompose \
--all-records \
--ref-overlap \
--threads=3    
```

tp.vcf.gz - 389,494 true positive variants
fp.vcf.gz - 425,679 false positive variants

The false positive vcf.gz was filtered further to identify potential mosaics.

## Merge AJ trio V4.2.1 benchmark vcfs (this was already completed)

```
zsh merge_benchmark_vcf.zsh
```
Output in data/panel_design

GRCh38_aj-trio_v4.2.1_merged-variants.vcf.gz

## Create list of potential mosaics

The make_panel.zsh script was run to generate 
-AJ trio v4.2.1 benchmark bed intersection 
-AJ trio v4.2.1 benchmark complex union bed
-AJ trio v4.2.1 mosaic target bed (creating by subtracting complex union bed from benchmark bed intersection)
-Intersects vcfeval output (HG002 fp vcf.gz) with mosaic target bed to generate vcf with variants 
that occur in mosaic target regions

From mosaic-benchmarking/, run

```
zsh ~/git_repos/mosaic-benchmarking/scratch/make-draft-benchmark/make_panel.zsh
```

Outputs
aj_trio_complex_union.bed
GRCh38_aj_trio_v.4.2.1_mosaic-target.bed
GRCh38_HG002_v0.0.1_mosaic-benchmark.vcf - 366,728 variants (fall within mosaic target regions; file used to create potential mosaics database)

