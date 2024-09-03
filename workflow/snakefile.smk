configfile: "../config/config.yaml"

mutationtype = config["type"]
paths = config["mutation_files"]
logmodel = config["log_model"]


#FILE = glob_wildcards("../../MakeLogRegInput/annotated_datasets/combined6/{type}")
#print(FILE)

# haploinsufficiency analysis 
#include: "rules/haploinsufficiency.smk"

# compare models analysis
#include: "rules/compare_models.smk"
#lool at the yaml file to change stufff for  a smarter solution 

rule all:
    input:
        expand(["../output/models/{mutationtype}_{logmodel}_LassoBestModel.RData",
                "../output/levels/{mutationtype}_{logmodel}_levels.RData"
                ], mutationtype = mutationtype,logmodel = logmodel)

## add a rule which finds the reference context for the logistic regression 
# should make it a part of the training script
#maybe is should add model summary plots within this rule?

rule training_models:
    input: 
        trainingfile = lambda wc: paths[wc.mutationtype]
        #log_model = lambda wc: log_model
    resources:
        threads=4,
        time=250,
        mem_mb=50000
    conda: "envs/callrv2.yaml"
    output:
        model = "../output/models/{mutationtype}_{logmodel}_LassoBestModel.RData", # add no intercept_model
        levels = "../output/levels/{mutationtype}_{logmodel}_levels.RData"
    shell:"""
    Rscript scripts/modeltraining.R {input.trainingfile} {wildcards.mutationtype} {wildcards.logmodel} {output.model}
    """

#placeholder
# rule model_check_plotting:s
#     snp_model: expand(["../output/models/snp/{snp_type}_LassoBestModel.RData"]
#     indel_model: expand(["../output/models/indel/{indel_type}_LassoBestModel.RData"]
#     output: "plots/inter/summary.txt" # maybe its best to use a dummy_file instead of every plot
#     shell:"""
#    Rscript scripts/model_analysis.R {input.snp_model}
#    Rscript scripts/model_analysis.R {input.indel_model}
#    finshed {date} >> {output}
#    """

# rule prediction:
#     input:
#         snp_model = "output/{muttype}_LassoBestModel.RData",
#         #indel_model = ???
#         pred_data = "../MakeLogRegInput/annotated_datasets/all_possible/{muttype}_GC_repli_recomb_meth_0_long_hg38.dat.gz"
#     resources:
#         threads=2,
#         time=420,
#         mem_mb=15000 
#     params: 
#         levels = "output/{muttype}_levels.RData",
#     output:
#         predictions = "output/{muttype}_predictions.RData"
#     shell:"""
#     Rscript scripts/predicting.R {input.pred_data} {params.levels} {input.model}   
#     """