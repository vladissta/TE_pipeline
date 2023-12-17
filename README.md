# TE pipeline

## Table of Contents

- [About](#about)
- [Getting Started](#getting_started)
- [Usage](#usage)

## About <a name = "about"></a>

Pipeline to identify transposable elements _de novo_ in inputed sequenced samples  
(specified as *input_files* variable in **variables_config.yaml** file)

## Getting Started <a name = "getting_started"></a>


Clone repository
``` bash
git clone https://github.com/vladissta/TE_pipeline
cd TE_pipeline
```

Create conda environment
``` bash
conda env create -f environment.yaml
conda activate te_venv
```

## Usage <a name = "usage"></a>

1. Edit **variables_config.yaml** file (input paths of reference genome, hieriarchy file, etc.)  
2. Think if you should use snakemake profile – **profile_snakemake directory**  
**[!]** Current profile includes settings to execute snakemake on **slurm scheduler**  

If profile is useful for you:
``` bash
snakemake --profile profile_snakemake/ 
```
If not:

``` bash
snakemake
```

It is recommended firstly test pipeline by executing **dry run**:
``` bash
snakemake -n
```

## Author <a name = "author"></a>

[@vladissta](https://github.com/vladissta)