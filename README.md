# site_specific_mrate
We present a model that has  an improved mutation prediction, which we can use it to gain novel insights of the human genome

The following project is structered as suggested in the snakemake documentation: "https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html"

#tree -r -L 3
#add comment to the tree levels 

```bash
.
workflow
│   ├── snakefile.smk
│   ├── scripts
│   ├── rules
│   ├── plots
│   ├── logs
│   └── envs
├── resources
│   ├── snv_mutation_rate_model.txt
│   ├── readme_downloads.sh
│   ├── prepare_dnm
│   ├── old
│   ├── observed_mutations_DDD_2017_chr22.txt
│   ├── mutations
│   ├── mart_export.txt
│   ├── __MACOSX
│   ├── indel_mutation_rate_model.txt
│   ├── hg38.bed
│   ├── hg38.2bit
│   ├── genomeannotations
│   ├── gencode.v44.annotation.gtf.gz
│   ├── gencode.v42.annotation.gtf.gz
│   ├── gencode.v42.annotation.gff3.gz
│   ├── ClinGen_gene_curation_list_GRCh38.tsv
│   ├── annotationtracks.txt
│   ├── annotationtracks_others.txt
│   ├── 20160622.allChr.pilot_mask.bed
│   ├── 20160622.allChr.mask_noXY_strict.bed
│   ├── 20160622.allChr.mask_noXY.bed
│   └── 20160622.allChr.mask.bed
├── README.md
├── output
│   ├── Predictions
│   ├── PossibleMutations
│   ├── ObservedExpected
│   ├── Models
│   ├── KmerPaPa
│   ├── KmerCount
│   ├── Genovo
│   ├── EvenOddSplit
│   ├── CodingSplit
│   └── AnnotatedMutations
├── LICENSE
└── config
    └── config.yaml
```

| mutation type name            | numeric code |
|-------------------------------|--------------|
| unknown                       | 0            |
| synonymous                    | 1            |
| missense                      | 2            |
| nonsense                      | 3            |
| loss of stop codon            | 4            |
| loss of start codon           | 5            |
| loss of canonical splice site | 6            |
| Mutation within an intron     | 7            |
| Non-frameshift indel          | 8            |
| Frameshift indel              | 9            |
|                               |              |