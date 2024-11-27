## Reformatting strelka2 and loFreq chrom 20 150X mixture VCFs
TN=["tumor_af00", "tumor_af01", "tumor_af05", 
    "tumor_af10", "tumor_af25", "tumor_af50"]


vcf_dir="data/vcfs/mixtures/"

rule all:
    input: 
        expand("processed/loFreq_mix_chr20-300X_{tn}.vcf.gz", tn = TN),
        expand("processed/strelka_mix_chr20-300X_{tn}.vcf.gz", tn = TN)


rule reformat_lofreq:
    input: 
        snvs=vcf_dir + "loFreq_chr20_300X_results/{tn}_sorted.somatic-snvs.vcf.gz"
    output: "processed/loFreq_mix_chr20-300X_{tn}.vcf.gz"
    params: prefix="loFreq_mix_chr20-300X_{tn}"
    conda: "mosaic.yaml"
    shell: """
            bash scripts/reformat_loFreq.sh \
                {input.snvs} {params.prefix}
    """

rule reformat_strelka:
    input: 
        snvs= vcf_dir + "strelka_chr20_300X_results/{tn}_sorted.snvs.vcf.gz",
        indels= vcf_dir + "strelka_chr20_300X_results/{tn}_sorted.indels.vcf.gz"
    output: "processed/strelka_mix_chr20-300X_{tn}.vcf.gz"
    params: prefix="strelka_mix_chr20-300X_{tn}"
    conda: "mosaic.yaml"
    shell: """
            bash scripts/reformat_strelka2.sh \
                {input.snvs} {input.indels} {params.prefix}
    """
