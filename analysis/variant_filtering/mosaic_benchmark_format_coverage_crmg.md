

#### Final benchmark VCF formatting

```{r format_vcf_for_igv_python}

#Convert new benchmark VCF from Strelka2 to IGV compatible format and vcf4.2 version
#Move GT field to start of FORMAT column and corresponding value to start of HG002 column) - use strelka2_to_igv_format_4.2.py

python ./strelka2vcf_to_IGV_format_4.2.py <input.vcf> <output.vcf>

```

##### Add new VCF header (v.4.2)

```{r reheader_vcf}

#Get input vcf header in benchmark vcf creation folder

bcftools view -h input.vcf > header.txt

#Open header.txt in text editor and amend, reheader original vcf

bcftools reheader -h header.txt -o <output.vcf> <input.vcf>

```

#### VII. hap.py comparison of new mosaic benchmark to a callset

```{r add_mosaics_to_giabv4.2.1_benchmark}

#Add new mosaic benchmark variants to giabv4.2.1 benchmark

bcftools concat -a -D -o HG002_GRCh38_mosaics_RM_v1.0_plus_giabv4.2.1_benchmark.vcf.gz new_mosaic_benchmark_igv_new_header.vcf.gz HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz

```

##### hap.py command

```{r happy_run}

hap.py HG002_GRCh38_mosaics_RM_v1.0_plus_giabv4.2.1_benchmark.vcf.gz HG002_GRCh38_mosaic_benchmark_RM_v1.0.vcf.gz -f HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -r GRCh38.fa --engine=vcfeval --engine-vcfeval-template GRCh38.sdf -o happy_run_HG002_GRCh38_mosaic_benchmark_RM_v1.0_plus_GIABv4.2.1_GIAB_stratifications --stratification HG002_stratifications.tsv
```

##### Get % of GRCh38 genome covered by benchmark variants

```{r mosaic_genome_coverage}

# Determine mosaic benchmark variant coverage in GRCh38 genome (includes gaps) 

bedtools genomecov -i HG002_GRCh38_mosaic_benchmark_RM_v1.0.bed -g human.b38.genome > HG002_mosaic_variant_GRCh38_coverage.txt

# Bottom of output file

#genome	1	2457537662	3088286401	0.795761

2.45Gbp and 79.5% of genome covered (includes gaps)

# Determine mosaic benchmark variant coverage in GRCh38 - USING NON-GAPPED GENOME 
# See non-gap-ref-size.md to generate a GRCh38 genome non-gapped bed to only retain autosomes; code also generates gaps.bed

# Sum the bases in the benchmark and sum the bases in the non-gap-autosomes, and divide the former by the latter
awk '{sum+=$3-$2} END {print sum}' your.bed

awk '{sum+=$3-$2} END {print sum}' GRCh38_nongap_autosomes_Nate.bed

# 2745187818 Gbp

# 2457537662/2745187818 = 89.5% of GRCh 38 genome covered (excludes gaps, autosomes only)

bedtools genomecov -i HG002_GRCh38_mosaic_benchmark_RM_v1.0.vcf -g GRCh38_nongap_autosomes.bed > HG002_mosaic_variant_GRCh38_no_gaps_coverage.txt

```

##### Get GRCh38 coding regions overlapped by benchmark variants

```{r mosaic_genome_coverage}

# Determine GRCh38 coding regions that benchmark variants are found

#Get GRCh38 ref seq CDS BED

wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38@all/Funct#ional/GRCh38_refseq_cds.bed.gz

gzip -d cds.file.gz

bedtools intersect -a HG002_GRCh38_mosaic_benchmark_RM_v1.0.vcf -b GRCh38_refseq_cds.bed  > HG002_GRCh38_mosaic_variants_RM_GRCh38_cds.txt

```

#### VIII. Determine which variants hit medically relevant genes (CMRGs) 

```{r cmrg_analysis}

# Get list of 5,026 medically relevant genes used for CMRG benchmark set --> 
wget https://github.com/usnistgov/giab-cmrg-benchmarkset/blob/master/data/manually_created_files/GRCh38_mrg_full_gene.bed

#Sort mrg bed 

sort -k1,1V -k2,2n GRCh38_mrg_full_gene.bed > GRCh38_mrg_full_gene_sorted.bed 

#Intersect mosaic benchmark vcf with mrg bed
#Determine which mosaic benchmark variants overlap which MRG genes

bedtools intersect -a HG002_GRCh38_mosaic_benchmark_RM_v1.0.vcf -b GRCh38_mrg_full_gene.bed -wo > HG002_GRCh38_mosaic_benchmark_RM_v1.0_mrg_intersect.bed.txt

#Determine coverage for mosaic benchmark regions that overlap MRG regions

bedtools coverage -a GRCh38_mrg_full_gene.bed -b HG002_GRCh38_mosaic_benchmark_RM_v1.0.bed > HG002_GRCh38_mosaic_benchmark_RM_v1.0_mrg_coverage.txt

# Get separate bed of MRGs with >90% coverage by benchmark regions
# Add | awk '$NF>=0.9' to end of command above and wc -l new bed to get # of mrgs

#Determine # of bases in benchmark regions that are in MRGs
#Add | awk -F'\t' '{sum += $6} END {print "Total bases in medically relevant genes in mosaic benchmark regions: " sum}' to end of command above

#Calculate fraction of MRG bases covered by benchmark
#Add  | awk -F'\t' '{sum7 += $7; sum6 +=$6} END {print sum6/sum7}' to end of command above

```

