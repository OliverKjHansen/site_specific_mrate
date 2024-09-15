configfile: "../config/config.yaml"

mutationtypes = config["mutationtype"]
paths = config["mutationfiles"] # maybe do this smarter so the inputdata is not in the configfile
logmodels = config["logmodel"] # intercept or no_intercept(nobeta)
models = config["models"]

#gencode = config["gencode"]
genome2bit = config["hg382bit"]
genomebedfile = config["hg38bedfile"]


variationtypes = config["variationtype"]
radius_background = config["radiustype"]




# haploinsufficiency analysis 
#include: "rules/haploinsufficiency.smk"
# compare models analysis
#include: "rules/compare_models.smk"

glorific = "/home/oliver/.cache/pypoetry/virtualenvs/glorific-JUSrNIgv-py3.12/bin/glorific" #change this later

rule all:
    input:
        expand(["../output/KmerCount/{mutationtype}_unmutated_kmers.tsv",
                "../output/KmerCount/{mutationtype}_mutated_kmers.tsv",
                "../output/KmerPaPa/{mutationtype}_PaPa.tsv"], mutationtype = mutationtypes),
        expand(["../output/models/{mutationtype}_{logmodel}_LassoBestModel.RData"], mutationtype = mutationtypes,logmodel = logmodels),
        expand([ "../output/EvenOddSplit/{modeltype}_{mutationtype}_summary.RData"], mutationtype = mutationtypes, modeltype = models),
        expand(["../output/CodingSplit/coding_{modeltype}_{mutationtype}_summary.RData",
                "../output/CodingSplit/noncoding_{modeltype}_{mutationtype}_summary.RData"], mutationtype = mutationtypes, modeltype = models)

#This rule takes a list of mutations in #?# and creates the optimale KmerPaPa partition
#Where the mutation file should be vcf-like text file where the first four columns are: Chromosome, Position, Ref_Allele, Alt_Allele
###BIG PROBLEM### It seems that the mutations in the mutation doesnt match the reference genome from the 2bit file
rule KmerCount:
    input:
        mutationfile = "../resources/{mutationtype}denovo.tsv",
        regions = genomebedfile,
        genome = genome2bit
    resources:
        threads=4,
        time=250,
        mem_mb=80000
    params: 
        variationtype = lambda wc: variationtypes[wc.mutationtype][0], #this check if the mutationtype is a indel or snv
        breaktype = lambda wc: variationtypes[wc.mutationtype][1], #this determind how the breakpoit is modellen. Only importent for indels. should be empty for snv
        radius = lambda wc: radius_background[wc.mutationtype] #sets the radius from the configfile
    conda: "envs/kmercounter.yaml"
    output: 
        backgroundcount = "../output/KmerCount/{mutationtype}_unmutated_kmers.tsv",
        kmercount = "../output/KmerCount/{mutationtype}_mutated_kmers.tsv"
    shell:"""
        kmer_counter background --bed {input.regions} {params.radius} {input.genome} > {output.backgroundcount}
        kmer_counter {params.variationtype} -r 4 {input.genome} {input.mutationfile} {params.breaktype} > {output.kmercount}
    """
#--reverse_complement_method middle??

#Maybe the there could be inserted some crossvalidation here some wehere
#  incorporate superpatterns -s SUPER_PATTERN, --super_pattern SUPER_PATTERN
rule KmerPapa:
    input:
        backgroundcount = "../output/KmerCount/{mutationtype}_unmutated_kmers.tsv",
        kmercount = "../output/KmerCount/{mutationtype}_mutated_kmers.tsv"
    resources:
        threads=4,
        time=250,
        mem_mb=80000
    conda: "envs/kmerpapa.yaml"
    output: 
        kmerpartition = "../output/KmerPaPa/{mutationtype}_PaPa.tsv" 
    shell:"""
    kmerpapa --positive {input.kmercount} --background {input.backgroundcount} --penalty_values 3 5 6 --pseudo_counts 0.5 1 10 > {output.kmerpartition}
    """

##
#rule CreatingAnnotationTrack

#rule AnnotatingMutations:


rule training_models:
    input: 
        trainingfile = lambda wc: paths[wc.mutationtype]
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

#rule AnnotatingPossibleVariants: 

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

