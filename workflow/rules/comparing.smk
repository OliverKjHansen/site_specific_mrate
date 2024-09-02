configfile: "../config/config.yaml"

snp_type = config["snp_type"]
indel_type = config["indel_type"]

rule all:
    input:
        expand(["../../MakeLogRegInput/annotated_datasets/combined6/{snp_type}_GC_repli_recomb_meth_0.002_all3_long_hg38.dat.gz", #change
                "../../MakeLogRegInput/annotated_datasets/all_possible/{snp_type}_GC_repli_recomb_meth_0_long_hg38.dat.gz",
                "../../MakeLogRegInput/annotated_datasets/combined6/{indel_type}_GC_repli_recomb_meth_0.002_long_hg38.dat.gz"
                "../output/EvenOddSplit/9mer_{snp_type}_LassoBestModel.RData",
                "../output/EvenOddSplit/genomic_{snp_type}_LassoBestModel.RData",
                "../output/EvenOddSplit/complete_{snp_type}_LassoBestModel.RData",
                "../output/EvenOddSplit/3models_{snp_type}_levels.RData",
                "../output/EvenOddSplit/3models_{snp_type}_predictions.RData"],snp_type = snp_type, indel_type = indel_type)


rule EvenOddChromosomeSplit_snps:
    input:
        data = "../MakeLogRegInput/annotated_datasets/combined6/{muttype}_GC_repli_recomb_meth_0.002_all3_long_hg38.dat.gz"
    resources:
        threads=2,
        time=450,
        mem_mb=50000 
    output:
        model_9mer = "output/EvenOddSplit/9mer_{muttype}_LassoBestModel.RData",
        model_genomic = "output/EvenOddSplit/genomic_{muttype}_LassoBestModel.RData",
        model_complete = "output/EvenOddSplit/complete_{muttype}_LassoBestModel.RData",
        levels = "output/EvenOddSplit/3models_{muttype}_levels.RData",
        predictions = "output/EvenOddSplit/3models_{muttype}_predictions.RData",
    shell:"""
    Rscript scripts/EvenOddSplit.R {input.data} 
    """