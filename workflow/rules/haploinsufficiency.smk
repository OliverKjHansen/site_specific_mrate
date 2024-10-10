rule AnnotatingTranscripts:
    input:
        transcript_file = "../resources/test_cds.bed", # placeholder
        predictions = "../output/DeleteLate_Predictions/{mutationtype}_{logmodel}_predictions.tsv", #change later
    resources:
        threads=4,
        time=120,
        mem_mb=50000
    conda: "../envs/bedtools.yaml"
    output: 
        shorten_file = temp("../output/Transcripts/{mutationtype}_{logmodel}_shorten.tsv"), 
        transcripts_predictions_tmp =  temp("../output/Transcripts/{mutationtype}_{logmodel}_predictions_tmp.tsv"), # should be able to do this in on go
        transcripts_predictions =  "../output/Transcripts/{mutationtype}_{logmodel}_predictions.tsv",
        small_file = "../output/Transcripts/{mutationtype}_{logmodel}_predictions_small.tsv"
    shell:"""
    awk -v OFS="\\t" '{{print $1,$2,$3,$13}}' {input.transcript_file} > {output.shorten_file}
    awk -v OFS="\\t" '{{ $2=$2-1"\\t"$2; print }}' {input.predictions} > {output.transcripts_predictions_tmp}
    head -n 1 {output.transcripts_predictions_tmp} > {output.transcripts_predictions}
    sed -i '1 s/$/\\tchrom\\tint_start\\tint_end\\ttranscript_id/' {output.transcripts_predictions}
    bedtools intersect -wo -a {output.transcripts_predictions_tmp} -b {output.shorten_file} | awk -v OFS="\\t" '{{$NF=""; print $0}}' - >> {output.transcripts_predictions}
    awk -v OFS="\\t" '{{print $1,$3,$20,$21,$25,"{wildcards.mutationtype}"}}' {output.transcripts_predictions} > {output.small_file}
    """
rule SummingTranscripts:
    input:
        transcripts_predictions_all =  expand(["../output/Transcripts/{mutationtype}_{{logmodel}}_predictions_small.tsv"], mutationtype = mutationtypes)
    resources:
        threads=4,
        time=120,
        mem_mb=100000 
    conda: "../envs/callrv2.yaml"
    output: 
        all_types = "../output/Transcripts/expected/alltypes_{logmodel}_predictions_small.tsv",
        routput_1se = "../output/Transcripts/expected/expected_transcript_1se_{logmodel}.tsv",
        routput_min = "../output/Transcripts/expected/expected_transcript_min_{logmodel}.tsv"
    shell:"""
    awk FNR!=1 {input.transcripts_predictions_all} > {output.all_types}
    Rscript scripts/summing_transcripts.R {output.all_types} {output.routput_1se} {output.routput_min} 
    """