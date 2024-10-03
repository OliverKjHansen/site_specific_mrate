
rule AnnotatingTranscripts:
    input:
        transcript_file = "../resources/test_cds.bed", # placeholder
        predictions = "../output/Predictions/{mutationtype}_{logmodel}_predictions.tsv",
        header_file = "../resources/header_file.txt"
    resources:
        threads=4,
        time=120,
        mem_mb=5000 
    conda: "../envs/bedtools.yaml"
    output: 
        shorten_file = temp("../output/Transcripts/{mutationtype}_{logmodel}_shorten.tsv"), 
        transcripts_predictions_tmp =  temp("../output/Transcripts/{mutationtype}_{logmodel}_predictions_tmp.tsv"), # should be able to do this in on go
        transcripts_predictions =  "../output/Transcripts/{mutationtype}_{logmodel}_predictions.tsv"
    shell:"""
    awk -v OFS="\\t" '{{print $1,$2,$3,$13}}' {input.transcript_file} > {output.shorten_file}
    awk -v OFS="\\t" '{{ $2=$2-1"\\t"$2; print }}' {input.predictions} > {output.transcripts_predictions_tmp}
    head -n 1 {output.transcripts_predictions_tmp} > {output.transcripts_predictions}
    sed -i '1 s/$/\\tchrom\\tint_start\\tint_end\\ttranscript_id/' {output.transcripts_predictions}
    bedtools intersect -wo -a {output.transcripts_predictions_tmp} -b {output.shorten_file} | awk -v OFS="\\t" '{{$NF=""; print $0}}' - >> {output.transcripts_predictions}
    """