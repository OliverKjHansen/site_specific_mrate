

#rule ScalePredictionsBySynonymous: # better to split into multiple rules

glorific_path = "/home/oliver/.cache/pypoetry/virtualenvs/glorific-JUSrNIgv-py3.12/bin/glorific" #change this later

rule PossibleSyn:
    input:
        possible_variants = "../output/PossibleMutations/variants/{chromosome}_possible_variants.txt.gz" 
    resources:
        threads=4,
        time=120,
        mem_mb=10000
    output:
        unfiltered = temp("../output/PossibleMutations/syn/{mutationtype}/{chromosome}_possiblesyn_unfiltered.txt"),
        possibleSyn = "../output/PossibleMutations/syn/{mutationtype}/{chromosome}_possibleSyn.txt",
    shell:"""
    python scripts/splitting_to_synonymous.py {input.possible_variants} {wildcards.mutationtype} > {output.unfiltered}
    awk '{{print $0}}' {output.unfiltered} | sort -k1,1 -k2,2n | uniq > {output.possibleSyn}
    """

rule AnnotatePossibleSyn:
    input: 
        refgenome = genome2bit,
        kmerpartition = "../output/KmerPaPa/best_partition/{mutationtype}_best_papa.txt", 
        possibleSyn = "../output/PossibleMutations/syn/{mutationtype}/{chromosome}_possibleSyn.txt",
        callability = genomebedfile,
        annotationfile= annotation_parameters 
    resources:
        threads=8,
        time=880,
        mem_mb=80000
    params:
        glorific = glorific_path, # when conda compatibility fixed remive this
        var_type = lambda wc: mut_translations[wc.mutationtype][0],
        downsample = "0"
    output:
        annotateSyn = "../output/PossibleMutations/annotatedsyn/{mutationtype}/{chromosome}_possiblesyn_annotated.tsv"
    shell:"""
    {params.glorific} {params.var_type} {input.refgenome} {input.callability} {input.possibleSyn} -a {input.annotationfile} -d {params.downsample} -p {input.kmerpartition} -l --verbose | uniq > {output.annotateSyn}
    """

rule PredictSyn: # make a flag that makes empty insertion and deletion files
    input:
        model = "../output/Models/{mutationtype}_{logmodel}_LassoBestModel.RData",
        annotateSyn = "../output/PossibleMutations/annotatedsyn/{mutationtype}/{chromosome}_possiblesyn_annotated.tsv",
        possibleSyn = "../output/PossibleMutations/syn/{mutationtype}/{chromosome}_possibleSyn.txt"
    resources:
        threads=8,
        time=120,
        mem_mb=40000
    conda: "../envs/callrv2.yaml"
    params: 
        levels = "../output/Models/{mutationtype}_{logmodel}_levels.RData"
    output:
        predictions = "../output/Predictions/syn/{mutationtype}/{logmodel}/{mutationtype}_{logmodel}_{chromosome}_predictions.tsv",
        wtranscripts = "../output/Predictions/syn/{mutationtype}/{logmodel}/{mutationtype}_{logmodel}_{chromosome}_transcript_predictions.tsv"
    shell:"""
    if [[ {wildcards.mutationtype} == "insertion" || {wildcards.mutationtype} == "deletion" ]]; then
        touch {output.predictions}
        touch {output.wtranscripts}
    else
        Rscript scripts/predictions.R {input.annotateSyn} {input.model} {params.levels} {wildcards.logmodel} {output.predictions}
        Rscript scripts/joiningtranscripts.R {input.possibleSyn} {output.predictions} {output.wtranscripts}
    fi
    """ 

rule ExpectedSynAllMutationtypes:
    input:
        wtranscripts = expand(["../output/Predictions/syn/{mutationtype}/{{logmodel}}/{mutationtype}_{{logmodel}}_{{chromosome}}_transcript_predictions.tsv"], mutationtype = mutationtypes)
    resources:
        threads=4,
        time=60,
        mem_mb=50000
    conda: "../envs/callrv2.yaml"
    output:
        wtranscripts = "../output/Predictions/syn/alltypes/{logmodel}_transcript_predictions_{chromosome}.tsv"
    shell:"""
    awk 'FNR>1 || NR==1' {input.wtranscripts} > {output.wtranscripts}
    """

rule ObservedSynonymous: # run again
    input:
        genovo = "../output/Genovo/genovo_results_{chromosome}.txt",
        canonical = "../resources/mart_export.txt",
        wtranscripts = "../output/Predictions/syn/alltypes/{logmodel}_transcript_predictions_{chromosome}.tsv",  
    resources:
        threads=8,
        time=120,
        mem_mb=50000
    conda: "../envs/callrv2.yaml"
    output:
        #genovo_tmp = temp("../output/Genovo/{logmodel}_all_genovo_results.txt"),
        #wtranscripts = temp("../output/Predictions/{logmodel}_transcript_predictions.tsv"),
        synonymosrate = "../output/ObservedExpected/syn/{logmodel}/{logmodel}_observed_expected_{chromosome}.tsv"
    shell:"""
    Rscript scripts/synonymousrate.R {input.genovo} {input.canonical} {input.wtranscripts} {output.synonymosrate}
    """
rule Synonymosrate:
    input:  
        synonymosrate =  expand(["../output/ObservedExpected/syn/{{logmodel}}/{{logmodel}}_observed_expected_{chromosome}.tsv"], chromosome = chromosomes)
    resources:
        threads=4,
        time=120,
        mem_mb=150000
    conda: "../envs/callrv2.yaml"
    output:
        chromosomes = "../output/ObservedExpected/syn/all/synonymous_{logmodel}_observed_expected.tsv"
        #roc = "plots/{logmodel}_haploinsufficiency.pdf"
    shell:"""
    awk 'FNR>1 || NR==1' {input.synonymosrate} > {output.chromosomes}
    """