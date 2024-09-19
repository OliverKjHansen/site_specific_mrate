configfile: "../config/config.yaml"

mutationtypes = config["mutationtype"]
variationtypes = config["variationtype"]

paths = config["mutationfiles"] # maybe do this smarter so the inputdata is not in the configfile
logmodels = config["logmodel"] # intercept or no_intercept(nobeta)
models = config["models"]
possible_variants_path = config["possible_variants_path"]

#gencode = config["gencode"]
genome2bit = config["hg382bit"]
genomebedfile = config["hg38bedfile"]
annotation_parameters = config["annotation_parameters"]

bck_kmer= config["bck_kmer"]
mut_translations = config["mut_translations"]

include: "rules/kmerpapa_crossvalidation.smk"
penalty = config["penalty_values"]
pseudo = config["pseudo_counts"]
# haploinsufficiency analysis 
#include: "rules/haploinsufficiency.smk"
# compare models analysis
#include: "rules/compare_models.smk"

glorific_path = "/home/oliver/.cache/pypoetry/virtualenvs/glorific-JUSrNIgv-py3.12/bin/glorific" #change this later

# ["../output/KmerCount/{variationtype}_unmutated_kmers.tsv",
#                 "../output/KmerCount/{mutationtype}_mutated_kmers.tsv",

rule all:
    input:
        expand([#"../output/KmerPaPa/{mutationtype}_cv/{mutationtype}_{penalty}_{pseudo}_PaPa_cv.tsv",
                #"../output/KmerPaPa/{mutationtype}_cv/{mutationtype}_parametergridfile.txt",
                #"plots/KmerPaPa/{mutationtype}_alpha_and_pseudo_loglike.pdf",
                "../output/KmerPaPa/best_partition/{mutationtype}_best_papa.txt"], mutationtype = mutationtypes),# penalty = penalty, pseudo = pseudo), #could make this its own workflow
        expand(["../output/AnnotatedMutations/{mutationtype}_annotated.dat.gz",], mutationtype = mutationtypes),
                #"../output/AnnotatedPossibleVariants/{mutationtype}_possible_lof_annotated.dat.gz"], mutationtype = mutationtypes),
        expand(["../output/models/{mutationtype}_{logmodel}_LassoBestModel.RData"], mutationtype = mutationtypes,logmodel = logmodels),
        expand([ "../output/EvenOddSplit/{modeltype}_{mutationtype}_summary.RData"], mutationtype = mutationtypes, modeltype = models),
        expand(["../output/CodingSplit/coding_{modeltype}_{mutationtype}_summary.RData",
                "../output/CodingSplit/noncoding_{modeltype}_{mutationtype}_summary.RData"], mutationtype = mutationtypes, modeltype = models)


# need to fix the sorting and future filterong of files or do some magic with touching output if i refilter and sort everytime??
# rule SortingFilesAndFilter: ## it is easier to do as a rule because snakemake restarts everyting if nit
#     # input:
#     #     mutations = "../resources/{mutationtype}denovo.tsv", #hardcoded should change
#     #     callability = genomebedfile # maybe run some blacklist filterning on this, # A callability file that show which regions are good to filter on # might be a place holder
#     resources:
#         threads=4,
#         time=120,
#         mem_mb=5000
#     output:
#         mutations = "../resources/{mutationtype}denovo.tsv", #hardcoded should change
#         callability = genomebedfile
#     shell:"""
#     sort -k1,1 -o {output.mutations} {output.mutations}
#     sort -k1,1 -o {output.callability} {output.callability}
#     """

#rule CreateAnnotationTrack

#add the partion-flag -p kmerpap_output when i decide to run it
# does not work with sex chromosomes
rule AnnotateMutations: # we annotate existing mutations and generate a dataset of sampled from the full genome look at downsample parameter
    input: 
        ref_genome = genome2bit,
        kmerpartition = "../resources/papa_files/autosome_{mutationtype}_4.txt", #hardcoded should change # for now it is prevoius kmerpapafiles
        mutations = "../resources/{mutationtype}denovo.tsv", #hardcoded should change Raw mutations
        callability = genomebedfile, # maybe run some blacklist filterning on this, # A callability file that show which regions are good to filter on # might be a place holder
    resources:
        threads=4,
        time=120,
        mem_mb=5000
    #conda: "envs/glorific.yaml" # add if i can get glorific to be in a conda environment, will work with pip too
    params: 
        glorific = glorific_path, # when conda compatibility fixed remive this
        var_type = lambda wc: mut_translations[wc.mutationtype][0],
        downsample = "0.002",
        annotationfile= annotation_parameters #I could make this myself 
    output:
        annotated_mutations = "../output/AnnotatedMutations/{mutationtype}_annotated.dat.gz"
    shell:"""
    {params.glorific} {params.var_type} {input.ref_genome} {input.callability} {input.mutations} -a {params.annotationfile} -d {params.downsample} -p {input.kmerpartition} -l --verbose | gzip > {output.annotated_mutations}
    """

rule training_models:
    input: 
        trainingfile = lambda wc: paths[wc.mutationtype] # hardcoded, should at some point be the output from the AnnotateMutations rule
        #annotated_mutations = "../output/AnnotatedMutations/{mutationtype}_annotated.dat.gz" 
    resources:
        threads=4,
        time=250,
        mem_mb=80000
    conda: "envs/callrv2.yaml"
    output:
        model = "../output/models/{mutationtype}_{logmodel}_LassoBestModel.RData", # add no intercept_model
    shell:"""
    Rscript scripts/modeltraining.R {input.trainingfile} {wildcards.mutationtype} {wildcards.logmodel} {output.model}
    """

# rule GeneratePossibleMutations: # use Genovo and glorific
#     input: 
#         ref_genome = genome2bit,
#         kmerpartition = "../resources/papa_files/autosome_{mutationtype}_4.txt", #hardcoded should change
#         possible_lof = "../resources/unzipped_sorted/{mutationtype}.gencode_v42_LoF_SNVs_sorted.txt", #hardcoded should change
#         callability = genomebedfile, # maybe run some blacklist filterning on this, # A callability file that show which regions are good to filter on # might be a place holder
#         annotationfile= annotation_parameters #I could make this myself 
#     resources:
#         threads=4,
#         time=120,
#         mem_mb=5000
#     #conda: "envs/glorific.yaml" # add if i can get glorific to be in a conda environment, will work with pip too
#     params: 
#         glorific = glorific_path, # when conda compatibility fixed remive this
#         var_type = lambda wc: mut_translations[wc.mutationtype][0],
#         downsample = "0"
#     output:
#         possible_mutations = "../output/AnnotatedPossibleVariants/{mutationtype}_possible_lof_annotated.dat.gz"
#     shell:"""
#     {params.glorific} {params.var_type} {input.ref_genome} {input.callability} {input.possible_lof} -a {input.annotationfile} -d {params.downsample} -l --verbose | gzip > {output.annotated_mutations}
#     genovo 
#     sort -k1,1 -o {input.possible_lof} {input.possible_lof} # remember to sort the output or else the whole pipeline will run again
#     """
# touch empty
#/home/oliver/.cache/pypoetry/virtualenvs/glorific-JUSrNIgv-py3.12/bin/glorific indel files/hg38.2bit ../site_specific_mrate/resources/cds_42.bed empty -a ../MakeLogRegInput/parameter_files/hg38/GC_repli_recomb_meth.txt -d 1 -l --verbose 

#add the partion-flag -p kmerpap_output when i decide to run it
# rule AnnotatePossibleMutations:
#     input: 
#         ref_genome = genome2bit,
#         kmerpartition = "../resources/papa_files/autosome_{mutationtype}_4.txt", #hardcoded should change
#         possible_lof = lambda wc: possible_variants_path[wc.mutationtype], # hardcoded, should at some point be the output from the GeneratePossible rule
#         callability = genomebedfile, # maybe run some blacklist filterning on this, # A callability file that show which regions are good to filter on # might be a place holder
#         annotationfile= annotation_parameters #I could make this myself 
#     resources:
#         threads=4,
#         time=120,
#         mem_mb=5000
#     #conda: "envs/glorific.yaml" # add if i can get glorific to be in a conda environment, will work with pip too
#     params: 
#         glorific = glorific_path, # when conda compatibility fixed remive this
#         var_type = lambda wc: mut_translations[wc.mutationtype][0],
#         downsample = "0"
#     output:
#         annotated_mutations = "../output/AnnotatedPossibleVariants/{mutationtype}_possible_lof_annotated.dat.gz"
#     shell:"""
#     {params.glorific} {params.var_type} {input.ref_genome} {input.callability} {input.possible_lof} -a {input.annotationfile} -d {params.downsample} -l --verbose | gzip > {output.annotated_mutations}
#     """

# rule Prediction:
#     input:
#         snp_model = "output/{muttype}_nobeta_LassoBestModel.RData",
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

#rule ScalePredictions 