configfile: "../config/config.yaml"

mutationtypes = config["mutationtype"]
variationtypes = config["variationtype"]

#wildcards-
paths = config["mutationfiles"] # maybe do this smarter so the inputdata is not in the configfile
logmodels = config["logmodel"] # intercept or no_intercept(nobeta)
models = config["models"]
mutationfiles = config["mutationfiles"]
#possible_variants_path = config["possible_variants_path"] ##should be removed when i can generate all the files myself
#annotated_variants_path = config["annotated_variants_path"] #should be removed when i can generate all the files myself
chromosomes = config["chromosomes"]

#gencode = config["gencode"]
genome2bit = config["hg382bit"]
genomebedfile = config["hg38bedfile_strict"]
annotation_parameters = config["annotation_parameters"]

bck_kmer= config["bck_kmer"]
mut_translations = config["mut_translations"]

#kmerpapa crossvalidation
include: "rules/kmerpapa_crossvalidation.smk"
penalty = config["penalty_values"]
alpha = config["alpha_values"]

# haploinsufficiency analysis 
#include: "rules/haploinsufficiency.smk"

# compare models analysis
include: "rules/compare_models.smk"

glorific_path = "/home/oliver/.cache/pypoetry/virtualenvs/glorific-JUSrNIgv-py3.12/bin/glorific" #change this later

# ["../output/KmerCount/{variationtype}_unmutated_kmers.tsv",
#                 "../output/KmerCount/{mutationtype}_mutated_kmers.tsv",

#chrom pos context repli GC_1k GC_10k GC_100k recomb_decode recomb_pyrho meth CpG_I CDS mut A C G T

rule all:
    input:
        expand([#"../output/KmerPaPa/{mutationtype}_cv/{mutationtype}_{penalty}_{pseudo}_PaPa_cv.tsv",
                #"../output/KmerPaPa/{mutationtype}_cv/{mutationtype}_parametergridfile.txt",
                #"plots/KmerPaPa/{mutationtype}_penalty_and_pseudo_loglike.pdf",
                "../output/KmerPaPa/best_partition/{mutationtype}_best_papa.txt"], mutationtype = mutationtypes),# penalty = penalty, pseudo = pseudo), #could make this its own workflow
        expand(["../output/AnnotatedMutations/{mutationtype}_annotated.txt.gz"], mutationtype = mutationtypes),
        expand(["../output/Models/{mutationtype}_{logmodel}_LassoBestModel.RData"], mutationtype = mutationtypes, logmodel = logmodels),
        expand(["../output/PossibleMutations/variants/{chromosome}_possible_variants.txt.gz",
                "../output/PossibleMutations/LoF/{mutationtype}/{chromosome}_possibleLoF.txt",
                "../output/PossibleMutations/annotatedLoF/{mutationtype}/{chromosome}_possibleLoF_annotated.tsv",
                "../output/Predictions/{mutationtype}/{logmodel}/{mutationtype}_{logmodel}_{chromosome}_predictions.tsv",
                "../output/Predictions/{mutationtype}/{logmodel}/{mutationtype}_{logmodel}_{chromosome}_transcript_predictions.tsv"], mutationtype = mutationtypes,logmodel = logmodels, chromosome = chromosomes),
               # "../output/Transcripts/{mutationtype}_{logmodel}_predictions_small.tsv",
               # "../output/Transcripts/expected/expected_transcript_1se_{logmodel}.tsv",
               # "../output/Transcripts/expected/expected_transcript_min_{logmodel}.tsv"], mutationtype = mutationtypes,logmodel = logmodels),
        expand([ "../output/EvenOddSplit/{modeltype}_{mutationtype}_summary.tsv"], mutationtype = mutationtypes, modeltype = models)
                #"../output/CodingSplit/coding_{modeltype}_{mutationtype}_summary.RData",
                #"../output/CodingSplit/noncoding_{modeltype}_{mutationtype}_summary.RData"], mutationtype = mutationtypes, modeltype = models)

#add the partion-flag -p kmerpap_output when i decide to run it
# does not work with sex chromosomes
rule AnnotateMutations: # we annotate existing mutations and generate a dataset of sampled from the full genome look at downsample parameter
    input: 
        refgenome = genome2bit,
        kmerpartition = "../output/KmerPaPa/best_partition/{mutationtype}_best_papa.txt",
        mutations = lambda wc: mutationfiles[wc.mutationtype],
        callability = genomebedfile,
        annotationfile = annotation_parameters
    resources:
        threads=8,
        time=440,
        mem_mb=5000
    #conda: "envs/glorific.yaml" # add if i can get glorific to be in a conda environment, will work with pip too
    params: 
        glorific = glorific_path, # when conda compatibility fixed remive this
        var_type = lambda wc: mut_translations[wc.mutationtype][0],
        downsample = "0.002"
    output:
        annotated_mutations = "../output/AnnotatedMutations/{mutationtype}_annotated.txt.gz"
    shell:"""
    {params.glorific} {params.var_type} {input.refgenome} {input.callability} {input.mutations} -a {input.annotationfile} -d {params.downsample} -p {input.kmerpartition} -l --verbose | gzip > {output.annotated_mutations}
    """

rule TrainingModels:
    input: 
        #trainingfile = lambda wc: paths[wc.mutationtype] # hardcoded, should at some point be the output from the AnnotateMutations rule
        annotated_mutations = "../output/AnnotatedMutations/{mutationtype}_annotated.txt.gz" 
    resources:
        threads=8,
        time=880,
        mem_mb=150000
    conda: "envs/callrv2.yaml"
    params: #add a list of genomic features
    output:
        model = "../output/Models/{mutationtype}_{logmodel}_LassoBestModel.RData", # add no intercept_model
    shell:"""
    Rscript scripts/modeltraining.R {input.annotated_mutations} {wildcards.mutationtype} {wildcards.logmodel} {output.model}
    """

rule PossibleVariants:
    input:
        transcript_file = "../resources/gencode.v42.annotation.gff3.gz", # change to a config variable
        refgenome = genome2bit
    resources:
        threads=4,
        time=120,
        mem_mb=5000
    output:
        transcript_pr_chromosome = temp("../output/PossibleMutations/variants/{chromosome}_transcript_regions.txt"),
        possible_variants_pr_chromosome = "../output/PossibleMutations/variants/{chromosome}_possible_variants.txt.gz"
    shell:"""
    gunzip --stdout {input.transcript_file} | awk '$1 == "{wildcards.chromosome}"' - | /home/oliver/.cargo/bin/genovo --action transform --gff3 - --genomic-regions {output.transcript_pr_chromosome}
    /home/oliver/.cargo/bin/genovo --action possible_mutations --genomic-regions {output.transcript_pr_chromosome} --genome {input.refgenome} | tail -n +3 - | gzip > {output.possible_variants_pr_chromosome}
    """
# rule PossibleMutations:
#     input:
#         possible_mutations_pr_chromosome = expand(["../output/PossibleMutations/chr/{chromosomes}_variants.txt.gz"], chromosomes = chromosomes)
#     resources:
#         threads=4,
#         time=120,
#         mem_mb=20000
#     output:
#         possible_mutations = "../output/PossibleMutations/all_chromosomes_possible.txt.gz"
#     shell:"""
#     cat {input.possible_mutations_pr_chromosome} > {output.possible_mutations}
#     """
# the possible mutations is based on transcripts. This means that if two transripts overlap then the same possible mutation will occur multiple times in the final file. i take uniqe lines to remove dup mutations
rule PossibleLoF:
    input:
        possible_variants = "../output/PossibleMutations/variants/{chromosome}_possible_variants.txt.gz" 
    resources:
        threads=4,
        time=120,
        mem_mb=10000
    output:
        unfiltered = temp("../output/PossibleMutations/LoF/{mutationtype}/{chromosome}_possibleLoF_unfiltered.txt"),
        possibleLoF = "../output/PossibleMutations/LoF/{mutationtype}/{chromosome}_possibleLoF.txt"
    shell:"""
    python scripts/splitting_to_mutationtypes.py {input.possible_variants} {wildcards.mutationtype} > {output.unfiltered}
    awk '{{print $0}}' {output.unfiltered} | sort -k1,1 -k2,2n | uniq > {output.possibleLoF}
    """

rule AnnotatePossibleVariants: # indels are alrady annotated 
    input: 
        refgenome = genome2bit,
        kmerpartition = "../output/KmerPaPa/best_partition/{mutationtype}_best_papa.txt", 
        possibleLoF = "../output/PossibleMutations/LoF/{mutationtype}/{chromosome}_possibleLoF.txt",
        callability = genomebedfile,
        annotationfile= annotation_parameters 
    resources:
        threads=4,
        time=120,
        mem_mb=10000
    #conda: "envs/glorific.yaml" # add if i can get glorific to be in a conda environment, will work with pip too
    params: 
        glorific = glorific_path, # when conda compatibility fixed remive this
        var_type = lambda wc: mut_translations[wc.mutationtype][0],
        downsample = "0"
    output:
        annotatedLoF = "../output/PossibleMutations/annotatedLoF/{mutationtype}/{chromosome}_possibleLoF_annotated.tsv"
    shell:"""
    {params.glorific} {params.var_type} {input.refgenome} {input.callability} {input.possibleLoF} -a {input.annotationfile} -d {params.downsample} -p {input.kmerpartition} -l --verbose | uniq > {output.annotatedLoF}
    """
#annotated_lof = lambda wc: annotated_variants_path[wc.mutationtype] # old
rule Prediction:
    input:
        model = "../output/Models/{mutationtype}_{logmodel}_LassoBestModel.RData",
        annotatedLoF = "../output/PossibleMutations/annotatedLoF/{mutationtype}/{chromosome}_possibleLoF_annotated.tsv",
        possibleLoF = "../output/PossibleMutations/LoF/{mutationtype}/{chromosome}_possibleLoF.txt"
    resources:
        threads=4,
        time=360,
        mem_mb=20000 
    conda: "envs/callrv2.yaml"
    params: 
        levels = "../output/Models/{mutationtype}_{logmodel}_levels.RData"
    output:
        predictions = "../output/Predictions/{mutationtype}/{logmodel}/{mutationtype}_{logmodel}_{chromosome}_predictions.tsv",
        wtranscripts = "../output/Predictions/{mutationtype}/{logmodel}/{mutationtype}_{logmodel}_{chromosome}_transcript_predictions.tsv"
    shell:"""
    Rscript scripts/predictions.R {input.annotatedLoF} {input.model} {params.levels} {wildcards.logmodel} {output.predictions}
    Rscript scripts/joiningtranscripts.R {input.possibleLoF} {output.predictions} {output.wtranscripts}
    """

#rule ScalePredictions 