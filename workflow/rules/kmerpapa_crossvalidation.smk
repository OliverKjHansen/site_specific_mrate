#configfile: "../config/config.yaml"

# mutationtypes = config["mutationtype"]
# variationtypes = config["variationtype"]

# #gencode = config["gencode"]
# genome2bit = config["hg382bit"]
# genomebedfile = config["hg38bedfile"]

# bck_kmer= config["bck_kmer"] # change to superpattern
# mut_translations = config["mut_translations"]

rule BackgroundKmerCount:
    input:
        regions = genomebedfile, # might be a placeholder # maybe add blacklist filtering
        ref_genome = genome2bit #change to ref_genome
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
        kmer_counter background --bed {input.regions} {params.bck_kmer} {input.ref_genome} > {output.backgroundcount}
    """

#Where the mutation file should be vcf-like text file where the first four columns are: Chromosome, Position, Ref_Allele, Alt_Allele
rule MutationsKmerCount:
    input:
        mutationfile = "../resources/{mutationtype}denovo.tsv", ##change later #placeholder
        regions = genomebedfile,
        ref_genome = genome2bit
    resources:
        threads=1,
        time=60,
        mem_mb=5000 
    params: 
        var_type = lambda wc: mut_translations[wc.mutationtype][0], #this check if the mutationtype is a indel or snv
        sample = lambda wc: mut_translations[wc.mutationtype][1], # too sample, only relevant for indels
        breaktype = lambda wc:mut_translations[wc.mutationtype][2] #this determind how the breakpoit is modellen. Only importent for indels. should be empty for snv #sets the radius from the configfile
    conda: "../envs/kmercounter.yaml"
    output: 
        kmercount = "../output/KmerCount/{mutationtype}_mutated_kmers.tsv"
    shell:"""
        kmer_counter {params.var_type} {params.sample} -r 4 {input.ref_genome} {input.mutationfile} {params.breaktype} > {output.kmercount}
    """
#--reverse_complement_method middle??

#Maybe the there could be inserted some crossvalidation here some wehere
s_pattern = {"A2C": "--super_pattern NNNNANNNN", "A2G": "--super_pattern NNNNANNNN", "A2T": "--super_pattern NNNNANNNN", 
            "C2A": "--super_pattern NNNNCNNNN", "C2G": "--super_pattern NNNNCNNNN", "C2T": "--super_pattern NNNNCNNNN", 
            "insertion": "", "deletion": ""} ## im too tired to do this in a smart way

#remeber to insert the values for pseudo counts and crossvalidation- Right now it is just to see if it works
#make it a loop that trains a grid for each mutation type
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
        cv_values = "../output/KmerPaPa/{mutationtype}_cv/{mutationtype}_{penalty}_{pseudo}_PaPa_cv.tsv" 
    shell:"""
    kmerpapa --positive {input.kmercount} --background {input.backgroundcount} {params.super_pattern} --penalty_values {wildcards.penalty} --pseudo_counts {wildcards.pseudo} --CV_only --nfolds 3 --CVfile {output}
    """
#Finding the best KmerPaPa partition
#properly make a gridplot or a minimum function
#Rule 



