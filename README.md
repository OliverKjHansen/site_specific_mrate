# site_specific_mrate
We present a model that has  an improved mutation prediction, which we can use it to gain novel insights of the human genome

The following project is structered as suggested in the snakemake documentation: "https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html"

#tree -r -L 3
#add comment to the tree levels 

```bash
.
├── workflow
│   ├── snakefile.smk 
│   ├── scripts
│   ├── rules
│   │   ├── haploinsufficiency.smk 
│   │   └── comparing.smk
│   ├── plots
│   └── envs
├── resources
├── README.md
├── output
│   └── int
├── LICENSE
└── config
    └── config.yaml
```