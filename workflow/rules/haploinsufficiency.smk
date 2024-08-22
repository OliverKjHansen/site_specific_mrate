
rule all:
    input:
        expand(["output/haploinsufficiency/merged_predictionfiles.csv",
                'output/haploinsufficiency/merged_data.csv',
                'output/haploinsufficiency/expected_n_mutations.csv'],muttype = mut_type)



rule RData2CSV:
    input:
        predictions = "output/{muttype}_predictions.RData"
    output:
        csv = "output/csvformat/{muttype}_predictions.csv"
    shell:"""
    Rscript scripts/r2csv.R {input.predictions}
    """

rule creating_expected_counts:
    input:
        transcripts = "files/gencode_v42_LoF_SNV_w_transcript_collapsed.csv",
        merged_predictionfiles = expand("output/csvformat/{muttype}_predictions.csv", muttype = mut_type)
    resources:
        threads=2,
        time=450,
        mem_mb=50000 
    output:
        merged_predictionfiles = "output/haploinsufficiency/merged_predictionfiles.csv",
        output_csv_file = 'output/haploinsufficiency/merged_data.csv',
        transcriptid_csv_file = 'output/haploinsufficiency/expected_n_mutations.csv'
    shell:"""
    awk 'FNR>1 || NR==1' {input.merged_predictionfiles} >> {output.merged_predictionfiles} 
    python scripts/merging_df.py {output.merged_predictionfiles} {input.transcripts} 
    """