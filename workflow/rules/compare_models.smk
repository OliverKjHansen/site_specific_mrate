configfile: "../config/config.yaml"

# mutationtypes = config["type"]
# paths = config["mutation_files"]
# models = config["models"]
gencode = config["gencode"]

rule EvenOddChromosomeModels:
    input: 
        #data = lambda wc: paths[wc.mutationtype]
        annotated_mutations = "../output/AnnotatedMutations/{mutationtype}_annotated.txt.gz"
    resources:
        threads=4,
        time=450,
        mem_mb=80000
    output:
        tmp_mutations = temp("../output/AnnotatedMutations/{mutationtype}_annotated.txt")
        summary = "../output/EvenOddSplit/{modeltype}_{mutationtype}_summary.tsv"
    shell:"""
    zcat {input.annotated_mutations} > {output.tmp_mutations}
    Rscript scripts/EvenOddSplit.R {output.tmp_mutations} {wildcards.modeltype} {wildcards.mutationtype} {output.summary}
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