import sys
import gzip
# """
# this script takes a file of possible variants and splits it into the 8 mutationtype: A2C, A2G, A2T, C2G, C2A, C2T, deletions, insertions
# the input should be a genovo output formatted like this
# chr1    69420   C       G       ENST00000641515.2       Synonymous      TAGCCATGG
# chr1    69420   C       T       ENST00000641515.2       Synonymous      TAGCCATGG
# chr1    69420   C       CNNN    ENST00000641515.2       InFrameIndel    TAGCCATGG
# chr1    69420   C       CN      ENST00000641515.2       FrameshiftIndel TAGCCATGG
# chr1    69421   A       C       ENST00000641515.2       Missense        AGCCATGGG
# """

#ONLY HANDLES .GZ FILES, could proberæy change quick

base_table = {"A":"T", 
              "T":"A", 
              "G":"C", 
              "C":"G",
              "insertion":"insertion",
              "deletion":"deletion"}


def splitting_to_indels(file, mutationtype):
    with gzip.open(file, "rt") as f:
        for line in f:
            chrom, pos, ref, alt, transcript, outcome, sequence = line.split("\n")[0].split("\t")[0:7]
            if len(alt) > 1:
                if outcome == "FrameshiftIndel":
                    print(line.split("\n")[0])


def splitting_to_snps(file, mutationtype):
    mut_ref = mutationtype[0]
    mut_alt = mutationtype[2]
    with gzip.open(file, "rt") as f:
        for line in f:
            chrom, pos, ref, alt, transcript, outcome, sequence = line.split("\n")[0].split("\t")[0:7]
            if (ref == mut_ref) and (alt == mut_alt):
                if (outcome == "SpliceSite") or (outcome == "Nonsense"):
                    print(line.split("\n")[0])
                #print(chrom,pos,ref,alt, sep ="\t")
            if (ref == base_table[mut_ref] and alt == base_table[mut_alt]):
                if (outcome == "SpliceSite") or (outcome == "Nonsense"):
                    print(line.split("\n")[0])
                #print(chrom,pos, ref, alt, sep ="\t")
                #print(chrom,pos,ref, alt, base_table[ref], base_table[alt])

def main_split(file, mutationtype):
    if len(mutationtype) == 3:
        splitting_to_snps(file, mutationtype)
    if len(mutationtype) > 3:
        splitting_to_indels(file, mutationtype)

if __name__ == '__main__':
    file = sys.argv[1]
    mutationtype = sys.argv[2]
    #output = sys.argv[3]
    main_split(file, mutationtype)

