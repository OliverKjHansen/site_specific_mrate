
configfile: '../config/config.yaml'

snp_type = config["snp_type"]

# haploinsufficiency analysis 
#include: "rules/haploinsufficiency.smk"

# compare models analysis
#include: "rules/compare_models.smk"


#inputfiles, be smarter about this # can i be smarter when I point to the inputfiles
#"/home/oliver/MutationAnalysis/MakeLogRegInput/annotated_datasets/combined6/
#/MakeLogRegInput/annotated_datasets/all_possible/{muttype}_GC_repli_recomb_meth_0_long_hg38.dat.gz"

rule all:
    input:
        expand(["../../MakeLogRegInput/annotated_datasets/combined6/{snp_type}_GC_repli_recomb_meth_0.002_all3_long_hg38.dat.gz", #change
                "../../MakeLogRegInput/annotated_datasets/all_possible/{snp_type}_GC_repli_recomb_meth_0_long_hg38.dat.gz", #change
                "../output/models/{snp_type}_LassoBestModel.RData",
                "../output/levels/{snp_type}_levels.RData"#,
                #"output/{snp_type}_predictions.RData"
                ],snp_type = snp_type)

## add a rule which finds the reference context for the logistic regression 
# should make it a part of the training script
#maybe is should add model summary plots within this rule?

rule snp_models:
    input:
        trainingfile = "../../MakeLogRegInput/annotated_datasets/combined6/{snp_type}_GC_repli_recomb_meth_0.002_all3_long_hg38.dat.gz"
    resources:
        threads=2,
        time=450,
        mem_mb=50000 
    conda: "envs/callrv2.yaml"
    output:
        model = "../output/models/{snp_type}_LassoBestModel.RData", # add no intercept_model
        levels = "../output/levels/{snp_type}_levels.RData"
    shell:"""
    Rscript scripts/snp_model.R {input.trainingfile} {output.model}
    """
#placeholder for when we are ready to make the indel models
# rule indel_models:
#     input:
#         trainingfile = "../MakeLogRegInput/annotated_datasets/combined6/{muttype}_GC_repli_recomb_meth_0.002_all3_long_hg38.dat.gz"
#     resources:
#         threads=2,
#         time=450,
#         mem_mb=50000 
#     output:
#         model = "output/{muttype}_LassoBestModel.RData",
#         levels = "output/{muttype}_levels.RData"
#     shell:"""
#     Rscript scripts/full_model.R {input.trainingfile} 
#     """

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

#placeholder
# rule plotting:
#     input:
#     output:
#     shell: