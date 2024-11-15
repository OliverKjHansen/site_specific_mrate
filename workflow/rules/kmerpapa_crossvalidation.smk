#configfile: "../config/config.yaml"

# mutationtypes = config["mutationtype"]
# variationtypes = config["variationtype"]

# #gencode = config["gencode"]
# genome2bit = config["hg382bit"]
# genomebedfile = config["hg38bedfile"]

# bck_kmer= config["bck_kmer"] #change to superpattern
# mut_translations = config["mut_translations"]

penalty = config["penalty_values"]
alpha = config["alpha_values"]

rule BackgroundKmerCount:
    input:
        callability = genomebedfile,
        refgenome = genome2bit 
    resources:
        threads=4,
        time=120,
        mem_mb=5000 
    params:
        bck_kmer = lambda wc: bck_kmer[wc.variationtype] #sets the radius from the configfile
    conda: "../envs/kmercounter.yaml"
    output: 
        backgroundcount = "../output/KmerCount/{variationtype}_unmutated_kmers.tsv"
    shell:"""
    kmer_counter background --bed {input.callability} {params.bck_kmer} {input.refgenome} > {output.backgroundcount}
    """

#Where the mutation file should be vcf-like text file where the first four columns are: Chromosome, Position, Ref_Allele, Alt_Allele
#fitler out denovo with the strict file. 
rule MutationsKmerCount:
    input:
        mutationfile = lambda wc: mutationfiles[wc.mutationtype],  ##change later #placeholder
        callability = genomebedfile,
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
        filtered_denovo = temp("../output/Filtered/{mutationtype}denovo_filtered.tsv"),
        kmercount = "../output/KmerCount/{mutationtype}_mutated_kmers.tsv"
    shell:"""
    awk -v OFS="\\t" '{{print $1,$2-1,$2,$3,$4}}' {input.mutationfile} | bedtools intersect -a - -b {input.callability} | awk -v OFS="\\t" '{{print $1,$3,$4,$5}}' > {output.filtered_denovo}
    kmer_counter {params.var_type} {params.sample} -r 4 {input.refgenome} {input.mutationfile} {params.breaktype} > {output.kmercount}
    """

s_pattern = {"A2C": "--super_pattern NNNNANNNN", "A2G": "--super_pattern NNNNANNNN", "A2T": "--super_pattern NNNNANNNN", 
            "C2A": "--super_pattern NNNNCNNNN", "C2G": "--super_pattern NNNNCNNNN", "C2T": "--super_pattern NNNNCNNNN", 
            "insertion": "", "deletion": ""} ## im too tired to do this in a smart way

rule KmerPaPaCrossvalidation:
    input:
        backgroundcount =lambda wc: "../output/KmerCount/{}_unmutated_kmers.tsv".format(mut_translations[wc.mutationtype][0]), ##i tried to be smart but it took me 20 mintutes to realise i had to put the 0-index here 
        kmercount = "../output/KmerCount/{mutationtype}_mutated_kmers.tsv"
    resources:
        threads=8,
        time=480,
        mem_mb=150000 #more memory
    params:
        super_pattern = lambda wc: s_pattern[wc.mutationtype]
    conda: "../envs/kmerpapa.yaml"
    output: 
        cv_values = "../output/KmerPaPa/{mutationtype}_cv/{mutationtype}_{penalty}_{alpha}_PaPa_cv.tsv" 
    shell:"""
    kmerpapa --positive {input.kmercount} --background {input.backgroundcount} {params.super_pattern} --penalty_values {wildcards.penalty} --pseudo_counts {wildcards.alpha} --CV_only --nfolds 3 --CVfile {output}
    """
#Finding the best KmerPaPa partition
#properly make a gridplot or a minimum function

rule FindingBestParameter:
    input:
        crossvalidation = expand(["../output/KmerPaPa/{{mutationtype}}_cv/{{mutationtype}}_{penalty}_{alpha}_PaPa_cv.tsv"], penalty = penalty, alpha = alpha)
    resources:
        threads=1,
        time=60,
        mem_mb=5000 #more memory
    conda: "../envs/rplotting.yaml"
    output: 
        parameter_grid = "../output/KmerPaPa/{mutationtype}_cv/{mutationtype}_parametergridfile.txt",
        r_plot = "plots/KmerPaPa/{mutationtype}_penalty_and_alpha_loglike.pdf",
        best_parameters = "../output/KmerPaPa/best_parameters/{mutationtype}_penalty_and_alpha_best_parameters.txt" # just to output to BestKmerPaPa
    shell:"""
    awk 'FNR>1 || NR==1' {input.crossvalidation} > {output.parameter_grid}
    Rscript scripts/parameter_search_plot.R {wildcards.mutationtype} {output.parameter_grid} {output.r_plot} {output.best_parameters}
    """

rule BestKmerPaPa:
    input:
        backgroundcount =lambda wc: "../output/KmerCount/{}_unmutated_kmers.tsv".format(mut_translations[wc.mutationtype][0]), ##i tried to be smart but it took me 20 mintutes to realise i had to put the 0-index here 
        kmercount = "../output/KmerCount/{mutationtype}_mutated_kmers.tsv",
        best_parameters = "../output/KmerPaPa/best_parameters/{mutationtype}_penalty_and_alpha_best_parameters.txt" # kinda of working as a dummy
    resources:
        threads=8,
        time=480,
        mem_mb=150000 #more memory
    params:
        super_pattern = lambda wc: s_pattern[wc.mutationtype],
        penalty = lambda wc: int(open("../output/KmerPaPa/best_parameters/{}_penalty_and_alpha_best_parameters.txt".format(wc.mutationtype)).read().split()[0]),
        alpha = lambda wc: int(open("../output/KmerPaPa/best_parameters/{}_penalty_and_alpha_best_parameters.txt".format(wc.mutationtype)).read().split()[1])
    conda: "../envs/kmerpapa.yaml"
    output: 
        kmerpapa = "../output/KmerPaPa/best_partition/{mutationtype}_best_papa.txt"
    shell:"""
    kmerpapa --positive {input.kmercount} --background {input.backgroundcount} {params.super_pattern} --penalty_values {params.penalty} --pseudo_counts {params.alpha} -o {output.kmerpapa}
    sed -i '1s/^/#/' {output.kmerpapa}
    """
