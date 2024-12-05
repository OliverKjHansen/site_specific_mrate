rule ObservedMutations:
    input:
        mutations = "../resources/mutations/gnomad/vcfs/gnomad.genomes.v4.1.sites.{chromosome}.vcf.bgz",
        transcript_file = "../resources/gencode.v42.annotation.gff3.gz",
        callability = genomebedfile,
        refgenome = genome2bit,
        point_mutation_probabilities = "../output/KmerPaPa/best_partition/snv_mutation_rate_papa.txt",
        indel_mutation_probabilities = "../output/KmerPaPa/best_partition/indel_mutation_rate_papa.txt"
    resources:
        threads=8,
        time=880,
        mem_mb=150000
    conda: "../envs/kmercounter.yaml"
    output:
        filtered_mutations = "../resources/mutations/gnomad/derived/gnomad.filtered_{chromosome}.vcf",
        genovo = "../output/Genovo/genovo_results_{chromosome}.txt"
    shell:"""
    gunzip -c {input.mutations} | bedtools intersect -a - -b {input.callability} | awk -v OFS="\t" '{{print $1,$2,$4,$5}}' > {output.filtered_mutations}
    sed -i $"s/\\t/ /g" {output.filtered_mutations}
    gunzip --stdout {input.transcript_file} | awk '$1 == "{wildcards.chromosome}"' - | \
    /home/oliver/.cargo/bin/genovo \
	--gff3 - \
	--observed-mutations {output.filtered_mutations} \
	--genome {input.refgenome} \
	--point-mutation-probabilities {input.point_mutation_probabilities} \
	--indel-mutation-probabilities {input.indel_mutation_probabilities} \
	--significant-mutations {output.genovo}
    """

# FILTER=<ID=AC0,Description="Allele count is zero after filtering out low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)">FILTER=<ID=AS_VQSR,Description="Failed VQSR filtering thresholds of -2.502 for SNPs and -0.7156 for indels">
# FILTER=<ID=InbreedingCoeff,Description="Inbreeding coefficient < -0.3">
# FILTER=<ID=PASS,Description="Passed all variant filters">

# only take the header from the first file
rule ExpectedAllMutationtypes:
    input:
        wtranscripts = expand(["../output/Predictions/{mutationtype}/{{logmodel}}/{mutationtype}_{{logmodel}}_{{chromosome}}_transcript_predictions.tsv"], mutationtype = mutationtypes)
    resources:
        threads=4,
        time=60,
        mem_mb=50000
    conda: "../envs/callrv2.yaml"
    output:
        wtranscripts = "../output/Predictions/alltypes/{logmodel}_transcript_predictions_{chromosome}.tsv"
    shell:"""
    awk 'FNR>1 || NR==1' {input.wtranscripts} > {output.wtranscripts}
    """
rule ObservedExpectedRatio: # run again
    input:
        genovo = "../output/Genovo/genovo_results_{chromosome}.txt",
        canonical = "../resources/mart_export.txt",
        wtranscripts = "../output/Predictions/alltypes/{logmodel}_transcript_predictions_{chromosome}.tsv",  
        haploinsufficiency = "../resources/ClinGen_gene_curation_list_GRCh38.tsv"
    resources:
        threads=8,
        time=120,
        mem_mb=50000
    conda: "../envs/callrv2.yaml"
    output:
        #genovo_tmp = temp("../output/Genovo/{logmodel}_all_genovo_results.txt"),
        #wtranscripts = temp("../output/Predictions/{logmodel}_transcript_predictions.tsv"),
        oeratio = "../output/ObservedExpected/{logmodel}/{logmodel}_observed_expected_{chromosome}.tsv"
    shell:"""
    Rscript scripts/observedexpected.R {input.genovo} {input.canonical} {input.wtranscripts} {input.haploinsufficiency} {output.oeratio}
    """
rule HaploinsufficiencyROCCurve:
    input:  
        oeratio =  expand(["../output/ObservedExpected/{{logmodel}}/{{logmodel}}_observed_expected_{chromosome}.tsv"], chromosome = chromosomes)
    resources:
        threads=4,
        time=120,
        mem_mb=150000
    conda: "../envs/callrv2.yaml"
    output:
        chromosomes = "../output/ObservedExpected/all/{logmodel}_observed_expected.tsv"
        #roc = "plots/{logmodel}_haploinsufficiency.pdf"
    shell:"""
    awk 'FNR>1 || NR==1' {input.oeratio} > {output.chromosomes}
    """
#Rscript scripts/observedexpected.R {output.roc}