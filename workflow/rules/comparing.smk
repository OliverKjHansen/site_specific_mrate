configfile: "../../config/config.yaml"

type_ = config["type"]
paths = config["mutation_files"]

models = config["models"]

rule all:
    input:
        expand(["output/EvenOddSplit/{modeltype}_{mutationtype}_LassoBestModel.RData",
        "output/EvenOddSplit/{modeltype}_{muttype}_predictions.RData"],
        mutationtype = type_, modeltype = models)


rule EvenOddChromosomeSplit:
    input: 
        data = lambda wc: paths[wc.mutationtype]    
    resources:
        threads=2,
        time=450,
        mem_mb=50000
    params: 
        modeltype = "{modeltype}"
    output:
        model_type = "output/EvenOddSplit/{modeltype}_{mutationtype}_LassoBestModel.RData",
        predictions = "output/EvenOddSplit/{modeltype}_{mutationtype}_predictions.RData",
    shell:"""
    Rscript scripts/EvenOddSplit.R {input.data} {output.model} {params.modeltype} {output.predictions}
    """

# rule model_selection:
#     input:
#         trainingfile = lambda wc: paths[wc.mutationtype]    
#     resources:
#         threads=2,
#         time=450,
#         mem_mb=50000 
#     output:
#         model_9mer = "output/EvenOddSplit/9mer_{muttype}_LassoBestModel.RData",
#         model_genomic = "output/EvenOddSplit/genomic_{muttype}_LassoBestModel.RData",
#         model_complete = "output/EvenOddSplit/complete_{muttype}_LassoBestModel.RData",
#         levels = "output/EvenOddSplit/3models_{muttype}_levels.RData",
#         predictions = "output/EvenOddSplit/3models_{muttype}_predictions.RData",
#     shell:"""
#     Rscript scripts/EvenOddSplit.R {input.data} 
#     """