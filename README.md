# site_specific_mrate
We present a model that has  an improved mutation prediction, which we can use it to gain novel insights of the human genome

The following project is structered as suggested in the snakemake documentation: "https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html"

<!-- ├── .gitignore
├── README.md
├── LICENSE.md
├── workflow
│   ├── rules
|   │   ├── module1.smk
|   │   └── module2.smk
│   ├── envs
|   │   ├── tool1.yaml
|   │   └── tool2.yaml
│   ├── scripts
|   │   ├── script1.py
|   │   └── script2.R
│   ├── notebooks
|   │   ├── notebook1.py.ipynb
|   │   └── notebook2.r.ipynb
│   ├── report
|   │   ├── plot1.rst
|   │   └── plot2.rst
|   └── Snakefile
├── config
│   ├── config.yaml
│   └── some-sheet.tsv
├── results
└── resources -->
.
├── config
│   └── config.yaml
├── LICENSE
├── README.md
├── resources
├── results
└── workflow
    ├── envs
    ├── plots
    ├── rules
    ├── scripts
    └── snakefile.smk