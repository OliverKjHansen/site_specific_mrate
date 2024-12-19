configfile: "../config/config.yaml"

# mutationtypes = config["type"]
# paths = config["mutation_files"]
# models = config["models"]
#gencode = config["gencode"]
glorific_path = "/home/oliver/.cache/pypoetry/virtualenvs/glorific-JUSrNIgv-py3.12/bin/glorific" #change this later

rule EvenOddChromosomeModels:
    input: 
        #data = lambda wc: paths[wc.mutationtype]
        annotated_mutations = "../output/AnnotatedMutations/{mutationtype}_annotated.txt.gz"
    resources:
        threads=8,
        time=450,
        mem_mb=100000
    output:
        tmp_mutations = temp("../output/AnnotatedMutations/{mutationtype}_annotated_{modeltype}.txt"),
        summary = "../output/EvenOddSplit/{modeltype}_{mutationtype}_summary.tsv"
    shell:"""
    zcat {input.annotated_mutations} > {output.tmp_mutations}
    Rscript scripts/EvenOddSplit.R {output.tmp_mutations} {wildcards.modeltype} {wildcards.mutationtype} {output.summary}
    """


rule CompareBackgroundKmerCount:
    input:
        callability = evengenomebed,
        refgenome = genome2bit
    resources:
        threads=4,
        time=120,
        mem_mb=5000 
    params:
        bck_kmer = lambda wc: bck_kmer[wc.variationtype] #sets the radius from the configfile
    conda: "../envs/kmercounter.yaml"
    output: 
        backgroundcount = "../output/Compare/KmerCount/{variationtype}_unmutated_kmers.tsv"
    shell:"""
    kmer_counter background --bed {input.callability} {params.bck_kmer} {input.refgenome} > {output.backgroundcount}
    """

#Where the mutation file should be vcf-like text file where the first four columns are: Chromosome, Position, Ref_Allele, Alt_Allele
#fitler out denovo with the strict file. 
rule CompareMutationsKmerCount:
    input:
        mutationfile = lambda wc: mutationfiles[wc.mutationtype],  ##change later #placeholder
        callability = evengenomebed,
        refgenome = genome2bit
    resources:
        threads=1,
        time=60,
        mem_mb=5000 
    params: 
        var_type = lambda wc: mut_translations[wc.mutationtype][0], #this check if the mutationtype is a indel or snv
        sample = lambda wc: mut_translations[wc.mutationtype][1], # too sample, only relevant for indels
        breaktype = lambda wc:mut_translations[wc.mutationtype][2] #this determind how the breakpoint is modeled. Only importent for indels. should be empty for snv #sets the radius from the configfile
    conda: "../envs/kmercounter.yaml"
    output: 
        filtered_denovo = temp("../output/Compare/{mutationtype}denovo_filtered.tsv"),
        kmercount = "../output/Compare/KmerCount/{mutationtype}_mutated_kmers.tsv"
    shell:"""
    awk -v OFS="\\t" '{{print $1,$2-1,$2,$3,$4}}' {input.mutationfile} | bedtools intersect -a - -b {input.callability} | awk -v OFS="\\t" '{{print $1,$3,$4,$5}}' > {output.filtered_denovo}
    kmer_counter {params.var_type} {params.sample} -r 4 {input.refgenome} {output.filtered_denovo} {params.breaktype} > {output.kmercount}
    """

s_pattern = {"A2C": "--super_pattern NNNNANNNN", "A2G": "--super_pattern NNNNANNNN", "A2T": "--super_pattern NNNNANNNN", 
            "C2A": "--super_pattern NNNNCNNNN", "C2G": "--super_pattern NNNNCNNNN", "C2T": "--super_pattern NNNNCNNNN", 
            "insertion": "", "deletion": ""} ## im too tired to do this in a smart way

rule CompareKmerPaPaCrossvalidation:
    input:
        backgroundcount =lambda wc: "../output/Compare/KmerCount/{}_unmutated_kmers.tsv".format(mut_translations[wc.mutationtype][0]), ##i tried to be smart but it took me 20 mintutes to realise i had to put the 0-index here 
        kmercount = "../output/Compare/KmerCount/{mutationtype}_mutated_kmers.tsv"
    resources:
        threads=8,
        time=480,
        mem_mb=150000 #more memory
    params:
        super_pattern = lambda wc: s_pattern[wc.mutationtype]
    conda: "../envs/kmerpapa.yaml"
    output: 
        cv_values = "../output/Compare/KmerPaPa/{mutationtype}_cv/{mutationtype}_{penalty}_{alpha}_PaPa_cv.tsv" 
    shell:"""
    kmerpapa --positive {input.kmercount} --background {input.backgroundcount} {params.super_pattern} --penalty_values {wildcards.penalty} --pseudo_counts {wildcards.alpha} --CV_only --nfolds 3 --CVfile {output}
    """
#Finding the best KmerPaPa partition
#properly make a gridplot or a minimum function

rule CompareFindingBestParameter:
    input:
        crossvalidation = expand(["../output/Compare/KmerPaPa/{{mutationtype}}_cv/{{mutationtype}}_{penalty}_{alpha}_PaPa_cv.tsv"], penalty = penalty, alpha = alpha)
    resources:
        threads=1,
        time=60,
        mem_mb=5000 #more memory
    conda: "../envs/rplotting.yaml"
    output: 
        parameter_grid = "../output/Compare/KmerPaPa/{mutationtype}_cv/{mutationtype}_parametergridfile.txt",
        r_plot = "plots/Compare/KmerPaPa/{mutationtype}_penalty_and_alpha_loglike.pdf",
        best_parameters = "../output/Compare/KmerPaPa/best_parameters/{mutationtype}_penalty_and_alpha_best_parameters.txt" # just to output to BestKmerPaPa
    shell:"""
    awk 'FNR>1 || NR==1' {input.crossvalidation} > {output.parameter_grid}
    Rscript scripts/parameter_search_plot.R {wildcards.mutationtype} {output.parameter_grid} {output.r_plot} {output.best_parameters}
    """

rule CompareBestKmerPaPa:
    input:
        backgroundcount = lambda wc: "../output/Compare/KmerCount/{}_unmutated_kmers.tsv".format(mut_translations[wc.mutationtype][0]), ##i tried to be smart but it took me 20 mintutes to realise i had to put the 0-index here 
        kmercount = "../output/Compare/KmerCount/{mutationtype}_mutated_kmers.tsv",
        best_parameters = "../output/Compare/KmerPaPa/best_parameters/{mutationtype}_penalty_and_alpha_best_parameters.txt" # kinda of working as a dummy
    resources:
        threads=8,
        time=480,
        mem_mb=150000 #more memory
    params:
        super_pattern = lambda wc: s_pattern[wc.mutationtype],
        alpha = lambda wc: int(open("../output/Compare/KmerPaPa/best_parameters/{}_penalty_and_alpha_best_parameters.txt".format(wc.mutationtype)).read().split()[0]),
        penalty = lambda wc: int(open("../output/Compare/KmerPaPa/best_parameters/{}_penalty_and_alpha_best_parameters.txt".format(wc.mutationtype)).read().split()[1])
    conda: "../envs/kmerpapa.yaml"
    output: 
        kmerpapa = "../output/Compare/KmerPaPa/best_partition/{mutationtype}_best_papa.txt",
        kmerpapatypes = "../output/Compare/KmerPaPa/best_partition/{mutationtype}_rate.txt"
    shell:"""
    kmerpapa --positive {input.kmercount} --background {input.backgroundcount} {params.super_pattern} --penalty_values {params.penalty} --pseudo_counts {params.alpha} -o {output.kmerpapa}
    sed -i '1s/^/#/' {output.kmerpapa}
    awk 'FNR>1 || NR==0 {{print "{wildcards.mutationtype}",$1,$4}}' {output.kmerpapa} > {output.kmerpapatypes}
    """

rule CompareAnnotateTrainingData: # we annotate mutations from the even chromosomes as training data
    input: 
        refgenome = genome2bit,
        kmerpartition = "../output/Compare/KmerPaPa/best_partition/{mutationtype}_best_papa.txt",
        mutations = lambda wc: mutationfiles[wc.mutationtype],
        callability = evengenomebed,
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
        annotated_mutations = "../output/Compare/AnnotatedTrainingData/{mutationtype}_annotated.txt.gz"
    shell:"""
    {params.glorific} {params.var_type} {input.refgenome} {input.callability} {input.mutations} -a {input.annotationfile} -d {params.downsample} -p {input.kmerpartition} -l --verbose | gzip > {output.annotated_mutations}
    """
rule CompareTrainingModels: # train models on even chromosomes
    input: 
        kmerpartition = "../output/Compare/KmerPaPa/best_partition/{mutationtype}_best_papa.txt",
        annotated_mutations = "../output/Compare/AnnotatedTrainingData/{mutationtype}_annotated.txt.gz" 
    resources:
        threads=8,
        time=880,
        mem_mb=15000
    conda: "../envs/callrv2.yaml"
    params: #add a list of genomic features
    output:
        model = "../output/Compare/Models/{mutationtype}_{compare}_LassoBestModel.RData" # add no intercept_model
    shell:"""
    Rscript scripts/compare_modeltraining.R {input.annotated_mutations} {wildcards.mutationtype} {input.kmerpartition} {wildcards.compare} {output.model}
    """

rule CompareAnnotateTestData: # we annotate mutations from the odd genome as test data
    input: 
        refgenome = genome2bit,
        kmerpartition = "../output/Compare/KmerPaPa/best_partition/{mutationtype}_best_papa.txt",
        mutations = lambda wc: mutationfiles[wc.mutationtype],
        callability = oddgenomebed,
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
        annotated_mutations = "../output/Compare/AnnotatedTestData/{mutationtype}_annotated.txt.gz"
    shell:"""
    {params.glorific} {params.var_type} {input.refgenome} {input.callability} {input.mutations} -a {input.annotationfile} -d {params.downsample} -p {input.kmerpartition} -l --verbose | gzip > {output.annotated_mutations}
    """
rule ComparePrediction:
    input:
        model = "../output/Compare/Models/{mutationtype}_{compare}_LassoBestModel.RData",
        annotateddata = "../output/Compare/AnnotatedTestData/{mutationtype}_annotated.txt.gz",
        kmerpartition = "../output/Compare/KmerPaPa/best_partition/{mutationtype}_best_papa.txt"
    resources:
        threads=4,
        time=360,
        mem_mb=20000 
    conda: "../envs/callrv2.yaml"
    params: 
        levels = "../output/Compare/Models/{mutationtype}_{compare}_levels.RData"
    output:
        predictions = "../output/Compare/Predictions/{mutationtype}_{compare}_predictions.tsv"
    shell:"""
    Rscript scripts/compare_predictions.R {input.annotateddata} {input.kmerpartition} {input.model} {wildcards.compare} {params.levels} {output.predictions}    
    """
# rule CreatingCodingRegions:
#     input: 
#         gencode = gencode # change a bit
#     resources:
#         threads=2,
#         time=250,
#         mem_mb=20000
#     conda: "../envs/bedtools.yaml"
#     output:
#         coding = "../resources/exons_44.bed"
#     shell:"""
#     gunzip -c {input.gencode} | grep 'transcript_type "protein_coding"' |\
#     awk '($3=="exon") {{printf("%s\\t%s\\t%s\\n",$1,int($4)-1,$5);}}' |\
#     sort -T . -t $'\\t' -k1,1 -k2,2n |\
#     bedtools merge > {output.coding}
#     """

# rule CodingNonCodingSplit:
#     input:
#         data = lambda wc: paths[wc.mutationtype],
#         exon_bedfile = "../resources/exons_44.bed"
#     resources:
#         threads=2,
#         time=200,
#         mem_mb=40000
#     conda: "../envs/bedtools.yaml"
#     output:
#         mut_bedfile = temp("../output/CodingSplit/{mutationtype}.bed"),
#         coding = "../output/CodingSplit/coding_{mutationtype}.dat",
#         noncoding = "../output/CodingSplit/noncoding_{mutationtype}.dat"
#  #       summary = "../output/CodingSplit/{modeltype}_{mutationtype}_summary.RData"
#     shell:"""
#     zcat {input.data} | awk -v OFS="\\t" '{{$1=$1}}1' - | awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $2+1, $0}}' - | cut -f-3,6- -> {output.mut_bedfile}
#     intersectBed -wa -a {output.mut_bedfile} -b {input.exon_bedfile} | cut -f-2,4- | expand -t 1 - -> {output.coding}
#     intersectBed -v -a {output.mut_bedfile} -b {input.exon_bedfile} | cut -f-2,4- | expand -t 1 - -> {output.noncoding}
#     sed -i "1i $(zcat {input.data} | head -n 1)" {output.coding}
#     sed -i "1i $(zcat {input.data} | head -n 1)" {output.noncoding}
#    """

# rule CodingNonCodingModels:
#     input: 
#         coding = "../output/CodingSplit/coding_{mutationtype}.dat",
#         noncoding = "../output/CodingSplit/noncoding_{mutationtype}.dat"
#     resources:
#         threads=4,
#         time=450,
#         mem_mb=80000
#     output:
#         summary_coding = "../output/CodingSplit/coding_{modeltype}_{mutationtype}_summary.RData",
#         summary_noncoding = "../output/CodingSplit/noncoding_{modeltype}_{mutationtype}_summary.RData"
#     shell:"""
#     Rscript scripts/EvenOddSplit.R {input.coding} {wildcards.modeltype} {wildcards.mutationtype} {output.summary_coding}
#     Rscript scripts/EvenOddSplit.R {input.noncoding} {wildcards.modeltype} {wildcards.mutationtype} {output.summary_noncoding}
#     """

# rule model_selection:
#     input:
#         summary = expand(["../output/EvenOddSplit/{modeltype}_{mutationtype}_model_summary.RData"])
#     resources:
#         threads=2,
#         time=450,
#         mem_mb=50000 
#     output:
#        msAIC = "plots/modelselection_AIC.pdf",
#        msBIC = "plots/modelselection_BIC.pdf"
#     shell:"""
#     touch {output.msAIC}
#     touch {output.msBIC}
#     """
    # shell:"""
    # Rscript scripts/modelselection.R {input.summary} {output.msAIC}
    # """