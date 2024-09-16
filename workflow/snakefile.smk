configfile: "../config/config.yaml"

mutationtypes = config["mutationtype"]
variationtypes = config["variationtype"]

paths = config["mutationfiles"] # maybe do this smarter so the inputdata is not in the configfile
logmodels = config["logmodel"] # intercept or no_intercept(nobeta)
models = config["models"]

#gencode = config["gencode"]
genome2bit = config["hg382bit"]
genomebedfile = config["hg38bedfile"]

bck_kmer= config["bck_kmer"]
mut_translations = config["mut_translations"]




# haploinsufficiency analysis 
#include: "rules/haploinsufficiency.smk"
# compare models analysis
#include: "rules/compare_models.smk"

glorific = "/home/oliver/.cache/pypoetry/virtualenvs/glorific-JUSrNIgv-py3.12/bin/glorific" #change this later

rule all:
    input:
        expand(["../output/KmerCount/{variationtype}_unmutated_kmers.tsv",
                "../output/KmerCount/{mutationtype}_mutated_kmers.tsv",
                "../output/KmerPaPa/{mutationtype}_PaPa.tsv"], mutationtype = mutationtypes, variationtype = variationtypes),
        expand(["../output/models/{mutationtype}_{logmodel}_LassoBestModel.RData"], mutationtype = mutationtypes,logmodel = logmodels),
        expand([ "../output/EvenOddSplit/{modeltype}_{mutationtype}_summary.RData"], mutationtype = mutationtypes, modeltype = models),
        expand(["../output/CodingSplit/coding_{modeltype}_{mutationtype}_summary.RData",
                "../output/CodingSplit/noncoding_{modeltype}_{mutationtype}_summary.RData"], mutationtype = mutationtypes, modeltype = models)
rule BackgroundKmerCount:
    input:
        regions = genomebedfile,
        genome = genome2bit
    resources:
        threads=4,
        time=250,
        mem_mb=40000
    params:
        bck_kmer = lambda wc: bck_kmer[wc.variationtype] #sets the radius from the configfile
    conda: "envs/kmercounter.yaml"
    output: 
        backgroundcount = "../output/KmerCount/{variationtype}_unmutated_kmers.tsv"
    shell:"""
        kmer_counter background --bed {input.regions} {params.bck_kmer} {input.genome} > {output.backgroundcount}
    """

#Where the mutation file should be vcf-like text file where the first four columns are: Chromosome, Position, Ref_Allele, Alt_Allele
rule MutationsKmerCount:
    input:
        mutationfile = "../resources/{mutationtype}denovo.tsv",
        regions = genomebedfile,
        genome = genome2bit
    resources:
        threads=4,
        time=250,
        mem_mb=10000
    params: 
        vat_type = lambda wc: mut_translations[wc.mutationtype][0], #this check if the mutationtype is a indel or snv
        sample = lambda wc: mut_translations[wc.mutationtype][1], # too sample, only relevant for indels
        breaktype = lambda wc:mut_translations[wc.mutationtype][2] #this determind how the breakpoit is modellen. Only importent for indels. should be empty for snv #sets the radius from the configfile
    conda: "envs/kmercounter.yaml"
    output: 
        kmercount = "../output/KmerCount/{mutationtype}_mutated_kmers.tsv"
    shell:"""
        kmer_counter {params.vat_type} {params.sample} -r 4 {input.genome} {input.mutationfile} {params.breaktype} > {output.kmercount}
    """
#--reverse_complement_method middle??

#Maybe the there could be inserted some crossvalidation here some wehere
s_pattern = {"A2C": "--super_pattern NNNNANNNN", "A2G": "--super_pattern NNNNANNNN", "A2T": "--super_pattern NNNNANNNN", 
            "C2A": "--super_pattern NNNNCNNNN", "C2G": "--super_pattern NNNNCNNNN", "C2T": "--super_pattern NNNNCNNNN", 
            "insertion": "", "deletion": ""} ## im too tired to do this in a smart way

rule KmerPaPa:
    input:
        backgroundcount=lambda wc: "../output/KmerCount/{}_unmutated_kmers.tsv".format(mut_translations[wc.mutationtype][0]), ##i tried to be smart but it took me 20 mintutes to realise i had to put the 0-index here 
        kmercount=lambda wc: "../output/KmerCount/{mutationtype}_mutated_kmers.tsv"
    resources:
        threads=8,
        time=600,
        mem_mb=100000 #more memory
    params:
        bck_kmer = lambda wc: s_pattern[wc.mutationtype]
    conda: "envs/kmerpapa.yaml"
    output: 
        kmerpartition = "../output/KmerPaPa/{mutationtype}_PaPa.tsv" 
    shell:"""
    kmerpapa --positive {input.kmercount} --background {input.backgroundcount} {params.bck_kmer} --penalty_values 3 5 6 --pseudo_counts 0.5 1 10 > {output.kmerpartition}
    """

##
#rule CreatingAnnotationTrack

# rule AnnotatingMutations:


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

