# nf-report pipeline

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-reclaim/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

**scbirlab/nf-report** is a Nextflow pipeline for compound Repurposing by Ortholgy (RepOrt).

**Table of contents**

- [Processing steps](#processing-steps)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Credit](#credit)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Further help](#further-help)

## Processing steps

**scbirlab/nf-report** carries out the following steps:

1. Fetch all high confidence targets with inhibitors from ChEMBL
2. Fetch all protein sequences for targets from UniProt
3. Fetch proteomes of query organisms from UniProt
4. BLAST target protein sequences against query organism proteomes
5. Annotate LOEUF genomic constraint scores for human targets
6. Output inhibitors for targets meeting identity, coverage, and LOEUF cutoffs from ChEMBL

In parallel, the pipeline fetches data to assess LOEUF cutoffs for toxicity:

1. Fetch all inhibitors with cell IC50 or CC50 _and_ target biochemical $K_i$ from ChEMBL
2. Fetch all targets of these inhibitors from ChEMBL
3. Filter down to inhibitors that have target biochemical $K_i$ from ChEMBL
4. Output cell line IC50 paired with target biochemical $K_i$

## Requirements

You need Nextflow and either Anaconda, Singularity, or Docker to be installed.

### First time using Nextflow?

#### Crick users

If you're at the Crick **or your shared cluster has Nextflow and Singularity already installed**, try:

```bash
module load Nextflow Singularity
```

#### Others

Otherwise, if it's your first time using Nextflow on your system, you can install it using `conda`:

```bash
conda install -c bioconda nextflow
```

You may need to set the `NXF_HOME` environment variable. For example,

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

To make this a permanent change, you can do something like the following:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

## Quick start

The easiest way to get going is by specifying parameters on the command-line:

```bash
nextflow run scbirlab/nf-report \
    --organism_id 243273  \
    --min_identity 0.3 \
    --min_coverage 0.5 \
    --min_pchembl 7.0
```

Here's what the flags mean:
- `--organism_id`: The Taxon ID of the organism, whih you can find at NCBI or UniProt
- `--min_identity` (optional): minimum amino acid identity for orthology
- `--min_coverage` (optional): minimum coverage for orthology
- `--min_pchembl` (optional): minimum reported pChEMBL (potency) for inhibitors

[Other options](#command-line-usage) are available.

### Running with Singularity, Docker, or Conda

**scbirlab/nf-reclaim** runs on a Singularity container engine by default to ensure software versions are consistent. If you have 
docker installed, you can run using `-with-docker` to use it instead, or if you have Conda you can run `-with-conda`.

### Running more than one query in parallel

Make a [sample sheet (see below)](#sample-sheet) with columns representing the flags above, and, 
optionally, a [`nextflow.config` file](#inputs) in the directory where you want the pipeline to run. 
Then simply run:

```bash 
nextflow run scbirlab/nf-reclaim
```

### Pipeline versions

If you want to run a particular tagged version of the pipeline, such as `v0.0.2`, you can do so using

```bash 
nextflow run scbirlab/nf-reclaim -r v0.0.2
```

For help, use `nextflow run scbirlab/nf-reclaim --help`.

The first time you run the pipeline on your system, the software dependencies in `environment.yml` will be installed. 
This may take several minutes.

## Inputs

### Command-line usage

The pipeline can be run with command-line arguments:

```bash
nextflow run scbirlab/nf-reclaim --organism_id <taxon ID>
```

The following parameters are **required**:

```bash
--organism_id       Taxon ID for organism
# or if using sample sheet
--sample_sheet      CSV listing Taxon ID for multiple organisms
```

The following parameters have default options, and are **optional**.

- `min_identity = 35`: minimum amino acid identity for orthology
- `min_coverage = 0.7`: minimum sequence coverage for orthology
- `min_loeuf = 0.515`: minimum genomic constraint for human targets
- `min_pchembl = 6.0`: minimum inhibitor pChEMBL
- `gnomad_version = "4.1"`: which gNOMAD version to use for LOEUF values
- `tox_cell_lines = ["HCT116","HEK293T",...,"CHO"]`: cell lines to fetch toxicity data
- `outputs = "outputs"`: output directory

### Sample sheet

You can run multiple combinations in one command using a sample sheet. The sample sheet is a **CSV** file with one row per combination of parameters to run. 
```bash
nextflow run scbirlab/nf-reclaim --sample_sheet path/to/sample-sheet.csv
```

#### Sample sheet structure

Here is an example of the sample sheet to find all the mycoplasma orthologous 
putative inhibitors:

| organism_id | proteome_name           |
| ----------- | ----------------------- |
| 243273      | "Mycoplasma genitalium" |

Further examples are in the `test` directory of this repository.

### Config-file usage (recommended)

For reproducibility, self-documentation, and to save typing, 
parameters with the same names as the command line flags above can be 
provided in a `nextflow.config` file in the working directory. For example:

```groovy
params {
    organism_id = "243273"
}
```

Or with a sample sheet:

```groovy
params {
    sample_sheet = "path/to/sample-sheet.csv"
}
```

## Outputs

Outputs are saved in the `output` folder defined above. 

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/nf-reclaim/issues).

## Further help

Here are the pages of the software and databases used by this pipeline.

Databases:

- [ChEMBL](https://www.ebi.ac.uk/chembl/) for inhibitors and targets
- [UniProt](https://www.uniprot.org/) for protein sequences
- [NCBI Genbank](https://www.ncbi.nlm.nih.gov/genbank/) for taxonomy

Software:

- [diamond](https://github.com/bbuchfink/diamond) to BLAST many-against-many protein sequences.
