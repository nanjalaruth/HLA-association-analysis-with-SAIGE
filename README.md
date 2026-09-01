# HLA PheWAS using SAIGE (Nextflow pipeline)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![SAIGE](https://img.shields.io/badge/SAIGE-1.1.6.3-blue.svg)](https://saigegit.github.io/SAIGE-doc/)

## Introduction

The pipeline runs a phenome-wide association study (PheWAS) of imputed HLA alleles, amino acid residues and SNPs in the MHC region against a panel of binary phenotypes using SAIGE, a mixed-model method that accounts for sample relatedness and case–control imbalance.
See the [Paper](https://www.nature.com/articles/s41588-018-0184-y) and the [GitHub](https://github.com/saigegit/SAIGE) repo.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It runs inside the SAIGE Singularity container, so SAIGE itself does not need to be installed. Each phenotype is analysed in parallel against the same sparse GRM.

The pipeline takes in a genome-wide PLINK binary fileset (`.bed/.bim/.fam`) for relatedness estimation and the null model, a chromosome 6 VCF of imputed HLA variants for association testing, and a single phenotype–covariate file.

### Pipeline steps

| Step | Process | SAIGE script | Description |
| ---- | ------- | ------------ | ----------- |
| 0 | `create_sparse_grm` | `createSparseGRM.R` | Builds a sparse genetic relationship matrix from 5 000 random genome-wide markers with a relatedness cut-off of 0.05 |
| 1 | `saige_logistic_step_1` | `step1_fitNULLGLMM.R` | Fits the null logistic mixed model for each phenotype using the sparse GRM, adjusting for sex, age, age² and age×sex, and estimates the variance ratio |
| 2 | `saige_logistic_step_2` | `step2_SPAtests.R` | Tests each variant in the chromosome 6 VCF (MAF ≥ 0.01, MAC ≥ 1) with saddlepoint approximation and Firth correction for rare or unbalanced variants (`--pCutoffforFirth=0.1`, `--is_fastTest=TRUE`) |

> **Status:** the current `main.nf` executes step 0 only; steps 1 and 2 are written and commented out in the `workflow` block. Uncomment them once the sparse GRM has been built and checked.

## Installation

1. Nextflow

```
wget -qO- https://get.nextflow.io | bash
```

2. Singularity/Apptainer and the SAIGE container

```
singularity pull saige_1.1.6.3.sif docker://wzhou88/saige:1.1.6.3
```

Then set `process.container` in the `singularity` profile of `nextflow.config` to the path of the downloaded image, and add any data directories that need mounting to `singularity.runOptions` (e.g. `-B /path/to/data`).

## Input files

### Genome-wide genotypes (`whole_plink`)

PLINK binary fileset covering the whole genome (or a LD-pruned subset of common variants). It is used for the sparse GRM and for the null model in step 1. This should **not** be restricted to the MHC.

### Chromosome 6 VCF (`chr6_vcf`)

Bgzipped VCF of imputed HLA variants with a `.csi` or `.tbi` index. Genotypes are read from the `GT` field with `ref-first` allele ordering. Variant IDs should follow the usual HLA imputation convention (`HLA_*`, `AA_*`, `SNP_*`, `rs*`) so results can be split by variant type downstream.

### Sample list (`whole_ids`)

Plain-text file with one sample ID per line, in the same order as the VCF header. Passed to SAIGE as `--sampleFile`.

### Phenotype–covariate file (`cov_pheno`)

A single tab-delimited file containing both the covariates and every phenotype, with `IID` as the sample ID column. The covariate columns are fixed in `main.nf` and must be named exactly `sex`, `age`, `age_squared` and `age_sex`. Each phenotype column must be named `<label>_pheno`, where `<label>` is the name given in the config; the pipeline appends `_pheno` automatically.

```
IID     sex  age  age_squared  age_sex  T1D_pheno  RA_pheno  asthma_pheno
1001    1    54   2916         54       0          1         NA
1002    0    61   3721         0        1          0         0
```

Binary phenotypes must be coded 0/1 with `NA` for missing.

## Running the pipeline

The pipeline does not require installation as `Nextflow` will automatically fetch it from `GitHub`.

### Configure your data

Copy `conf/test.config` and edit the paths to point to your files:

```groovy
params {
    whole_plink = [
        ['MyCohort',
         '/path/to/genomewide.bed',
         '/path/to/genomewide.bim',
         '/path/to/genomewide.fam'],
    ]

    chr6_vcf = [
        ['MyCohort',
         '/path/to/chr6_hla_imputed.vcf.gz',
         '/path/to/chr6_hla_imputed.vcf.gz.csi'],
    ]

    whole_ids = "/path/to/sample_ids.txt"

    cov_pheno = [
        ["T1D",    "/path/to/cov_pheno.tsv"],
        ["RA",     "/path/to/cov_pheno.tsv"],
        ["asthma", "/path/to/cov_pheno.tsv"],
    ]

    outdir = "./output"
}
```

`cov_pheno` is a list of `[phenotype_label, cov_pheno_file]`; the same file can be reused for every phenotype. The dataset name (first element of `whole_plink` and `chr6_vcf`) must match, as it is used to join the GRM, null model and VCF channels.

### Run on a Slurm cluster

```
nextflow run nanjalaruth/HLA-PheWAS-SAIGE -profile singularity,slurm -c <path to your edited config file> -resume
```

A submission script is provided in `saige.sh`, which launches the whole pipeline from a single high-memory job:

```
sbatch saige.sh
```

### Run locally

```
nextflow run nanjalaruth/HLA-PheWAS-SAIGE -profile singularity -c <path to your edited config file> -resume
```

## To run the updated version of this pipeline, run:

```
nextflow pull nanjalaruth/HLA-PheWAS-SAIGE
```

## Arguments

### Required Arguments

| Argument       | Usage                                              | Description |
| -------------- | -------------------------------------------------- | ----------- |
| -profile       | <singularity,slurm,standard,debug>                 | Configuration profile to use. `singularity` is required to run the SAIGE container; combine with `slurm` to submit each process to the scheduler |
| -c             | <path to config>                                   | Config file defining the parameters below |
| --whole_plink  | [[dataset, bed, bim, fam], ...]                    | Genome-wide PLINK fileset for the sparse GRM and null model |
| --chr6_vcf     | [[dataset, vcf.gz, index], ...]                    | Chromosome 6 VCF of imputed HLA variants and its index |
| --whole_ids    | <path to ids.txt>                                  | Sample ID file matching the VCF |
| --cov_pheno    | [[label, cov_pheno_file], ...]                     | Phenotype labels and the phenotype–covariate file |

### Optional Arguments

| Argument       | Default                        | Description |
| -------------- | ------------------------------ | ----------- |
| --outdir       | `./output`                     | Directory for published results |
| --tracedir     | `${outdir}/pipeline_info`      | Directory for the Nextflow timeline, report, trace and DAG |
| --max_memory   | `128.GB`                       | Upper bound for process memory |
| --max_cpus     | `16`                           | Upper bound for process CPUs |
| --max_time     | `240.h`                        | Upper bound for process wall time |

All three SAIGE processes are labelled `bigmem` (120 GB, 9 CPUs, 24 h, doubled on retry). Adjust the `withLabel` blocks in `nextflow.config` if your cluster limits differ.

## Output

```
output/
├── Regression_results/
│   ├── <dataset>_<pheno>_saige_out.rda                    # null model (step 1)
│   ├── <dataset>_<pheno>_saige_out.varianceRatio.txt      # variance ratio (step 1)
│   └── SAIGE/
│       └── <dataset>_<pheno>_binary.SAIGE.vcf.genotype.txt   # association results (step 2)
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    ├── execution_trace.txt
    └── pipeline_dag.png
```

The sparse GRM (`<dataset>_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx` and its `.sampleIDs.txt`) is kept in the Nextflow `work/` directory and reused across phenotypes; uncomment the `publishDir` line in `create_sparse_grm` to copy it to `outdir`.

The step 2 output is the standard SAIGE results table (`CHR, POS, MarkerID, Allele1, Allele2, AC_Allele2, AF_Allele2, N, BETA, SE, Tstat, p.value, ...`), one file per phenotype. Classical HLA alleles can be extracted with `grep '^6.*HLA_'` or by filtering `MarkerID` in R.

## Notes

- SAIGE step 2 is run with `--LOCO=FALSE` because the sparse GRM is built genome-wide while only chromosome 6 is tested.
- The `--minMAF=0.01` threshold in step 2 excludes rare HLA alleles; lower it in `main.nf` if these are of interest.
- The MHC is highly polymorphic; the GRM should be built from genome-wide markers, not from the HLA panel itself.

## Related pipelines

- [hla_typing_using_HLA-LA](https://github.com/nanjalaruth/hla_typing_using_HLA-LA) — HLA typing from sequencing data
- [HLA-association-analysis-with-Regenie](https://github.com/nanjalaruth/HLA-association-analysis-with-Regenie) — the same association design using Regenie instead of SAIGE

## Support

I track open tasks using github's [issues](https://github.com/nanjalaruth/HLA-PheWAS-SAIGE/issues)
