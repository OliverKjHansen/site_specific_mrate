


rule all:
    input:
        expand(["../MakeLogRegInput/annotated_datasets/combined6/{muttype}_GC_repli_recomb_meth_0.002_all3_long_hg38.dat.gz",
                "../MakeLogRegInput/annotated_datasets/all_possible/{muttype}_GC_repli_recomb_meth_0_long_hg38.dat.gz",
                "output/EvenOddSplit/9mer_{muttype}_LassoBestModel.RData",
                "output/EvenOddSplit/genomic_{muttype}_LassoBestModel.RData",
                "output/EvenOddSplit/complete_{muttype}_LassoBestModel.RData",
                "output/EvenOddSplit/3models_{muttype}_levels.RData",
                "output/EvenOddSplit/3models_{muttype}_predictions.RData"],muttype = mut_type)


rule EvenOddChromosomeSplit:
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