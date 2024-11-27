## Making Tumor/ Normal Mixtures - for GIAB Mosaic Variant Power Analysis

#wkdir="/Volumes/giab/analysis/giab-mosaic-variants/data"
wkdir="data"
#inbam_dir="/Volumes/giab/data/alignment/AshkenazimTrio/Illumina/NHGRI_Illumina300X_AJtrio_novoalign_bams"
inbam_dir="data/bams"

## Defining Variables and Wildcards ###################################################
HG=["HG002","HG003"]
TN=["normal", "tumor_af00", "tumor_af01", "tumor_af05", 
    "tumor_af10", "tumor_af25", "tumor_af50"]

## Wildcard Constraints
wildcard_constraints:
    hg="HG002|HG003"

rule all:
    input: 
        expand(wkdir + "/subsampled_bams/chrom20/300X/HG002_517.{prop}.cov", 
               prop = ["01", "05", "10", "25", "50"]),
        expand(wkdir + "/subsampled_bams/chrom20/300X/HG003_517.{prop}.cov", 
               prop = ["25", "40", "45", "49", "50"]),
        wkdir + "/subsampled_bams/chrom20/300X/HG003_531.50.cov",
        expand(wkdir + "/mixtures/chrom20/300X/{tn}_sorted.bam.bai", tn = TN),
        expand(wkdir + "/mixtures/chrom20/300X/{tn}_sorted.cov", tn = TN)

## Adding Read Groups to input bams ###################################################
rule tumor_read_group:
    input: inbam_dir + "/HG002.hs37d5.300X.chr20.bam"
    output:
        bam=wkdir + "/chrom20/300X/HG002_wrg.bam",
        bai=wkdir + "/chrom20/300X/HG002_wrg.bam.bai"
    params: sm="HG002", id="HG002-300X"
    conda: "mosaic.yaml"
    shell: """
            picard AddOrReplaceReadGroups \
                I={input} \
                O={output.bam} \
                RGID={params.id} \
                RGLB=all \
                RGPL=ILLUMINA \
                RGPU=all \
                RGSM={params.sm}

            samtools index {output.bam}
    """

rule normal_read_group:
    input: inbam_dir + "/HG003.hs37d5.300X.chr20.bam"
    output: 
        bam=wkdir + "/chrom20/300X/HG003_wrg.bam",
        bai=wkdir + "/chrom20/300X/HG003_wrg.bam.bai"
    params: sm="HG003", id="HG003-300X"
    conda: "mosaic.yaml"
    shell: """
           picard AddOrReplaceReadGroups \
                I={input} \
                O={output.bam} \
                RGID={params.id} \
                RGLB=all \
                RGPL=ILLUMINA \
                RGPU=all \
                RGSM={params.sm}

            samtools index {output.bam}
    """

### Subsetting Bams for Mixtures ########################################################
rule subset_bam:
    input: 
        bam=wkdir + "/chrom20/300X/{hg}_wrg.bam", 
        bai=wkdir + "/chrom20/300X/{hg}_wrg.bam.bai"
    output: 
        bam=wkdir + "/subsampled_bams/chrom20/300X/{hg}_{frac}.bam",
        bai=wkdir + "/subsampled_bams/chrom20/300X/{hg}_{frac}.bam.bai"	
    params: frac = "{frac}"
    conda: "mosaic.yaml"
    shell: """
	samtools view -bh -s {params.frac} {input.bam} > {output.bam}
	samtools index {output.bam} 
    """

### Subset  QC/ Sanity Check ########################################################
rule index_subset_bam:
    input: wkdir + "/subsampled_bams/chrom20/300X/{hg}_{frac}.bam"
    output: wkdir + "/subsampled_bams/chrom20/300X/{hg}_{frac}.bam.bai"
    conda: "mosaic.yaml"
    shell: "samtools index {input}"
        
rule index_stats:
    input: 
        bam=wkdir + "/subsampled_bams/chrom20/300X/{hg}_{frac}.bam",
        bai=wkdir + "/subsampled_bams/chrom20/300X/{hg}_{frac}.bam"
    output: wkdir + "/subsampled_bams/chrom20/300X/{hg}_{frac}_idxstats.txt"
    conda: "mosaic.yaml"
    shell: """
        samtools idxstats {input.bam} > {output}
    """

rule depth_stats:
    input: wkdir + "/subsampled_bams/chrom20/300X/{hg}_{frac}.bam"
    output: wkdir + "/subsampled_bams/chrom20/300X/{hg}_{frac}.cov"
    conda: "mosaic.yaml"
    shell: "samtools depth {input} > {output}"

### Make Mixtures #######################################################################
rule normal_bam:
    input: wkdir + "/chrom20/300X/HG003_wrg.bam"
    output: wkdir + "/mixtures/chrom20/300X/normal.bam"
    shell: "cp {input} {output}"

rule af50_bam:
    input: wkdir + "/chrom20/300X/HG002_wrg.bam"
    output: wkdir + "/mixtures/chrom20/300X/tumor_af50.bam"
    shell: "cp {input} {output}"

rule af25_bam:
    input: 
        wkdir + "/subsampled_bams/chrom20/300X/HG002_517.50.bam", 
        wkdir + "/subsampled_bams/chrom20/300X/HG003_517.50.bam"
    output: wkdir + "/mixtures/chrom20/300X/tumor_af25.bam"
    params: "-f -O bam"
    threads: 3
    wrapper: "0.38.0/bio/samtools/merge"

rule af10_bam:
    input: 
        wkdir + "/subsampled_bams/chrom20/300X/HG002_517.20.bam", 
        wkdir + "/subsampled_bams/chrom20/300X/HG003_517.80.bam"
    output: wkdir + "/mixtures/chrom20/300X/tumor_af10.bam"
    params: "-f -O bam"
    threads: 3
    wrapper: "0.38.0/bio/samtools/merge"

rule af05_bam:
    input: 
        wkdir + "/subsampled_bams/chrom20/300X/HG002_517.10.bam", 
        wkdir + "/subsampled_bams/chrom20/300X/HG003_517.90.bam"
    output: wkdir + "/mixtures/chrom20/300X/tumor_af05.bam"
    params: "-f -O bam"
    threads: 3
    wrapper: "0.38.0/bio/samtools/merge"

rule af01_bam:
    input: 
        wkdir + "/subsampled_bams/chrom20/300X/HG002_517.02.bam", 
        wkdir + "/subsampled_bams/chrom20/300X/HG003_517.98.bam"
    output: wkdir + "/mixtures/chrom20/300X/tumor_af01.bam"
    params: "-f -O bam"
    threads: 3
    wrapper: "0.38.0/bio/samtools/merge"

rule af00_bam:
    input: wkdir + "/chrom20/300X/HG003_wrg.bam"
    output: wkdir + "/mixtures/chrom20/300X/tumor_af00.bam"
    params: "-f -O bam"
    threads: 3
    wrapper: "0.38.0/bio/samtools/merge"

## Sort and index t-n mixtures ############################################################
rule sort_tn_bam:
    input: wkdir + "/mixtures/chrom20/300X/{tn}.bam"
    output: wkdir + "/mixtures/chrom20/300X/{tn}_sorted.bam"
    conda: "mosaic.yaml"
    shell: "samtools sort -m 4G -@8 {input} > {output}"

rule index_sorted_tn_bam:
    input: wkdir + "/mixtures/chrom20/300X/{tn}_sorted.bam"
    output: wkdir + "/mixtures/chrom20/300X/{tn}_sorted.bam.bai"
    conda: "mosaic.yaml"
    shell: "samtools index {input}"

### Mix Depth ############################################################################
rule mix_depth:
    input: wkdir + "/mixtures/chrom20/300X/{tn}_sorted.bam"
    output: wkdir + "/mixtures/chrom20/300X/{tn}_sorted.cov"
    conda: "mosaic.yaml"
    shell: "samtools depth {input} > {output}"
