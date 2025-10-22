# Data collection
This section describes how to download data from NCBI using CLI

## Downloading a dataset from NCBI
This section describes how to download genomic and transcriptome data from NCBI using the CLI

### Download the required tooling
With a version of conda (I'm using miniforge3), create an environment and install needed CLI tooling
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/

Create the conda environment

`conda create -n ncbi-cli`

Activate the conda environment

`conda activate ncbi-cli`

Install needed NCBI CLI tooling

`conda install ncbi-datasets-cli`

### Downloading data via CLI
- Navigate to NCBI Genome https://www.ncbi.nlm.nih.gov/datasets/genome/
- Search for a crop (Glycine max, Triticum aestivum, Arachis hypogaea)
- Within an assembly page: `datasets` will be at the top
    - select and copy the text

Once you have the download command, activate your ncbi-cli environment and donwload the data

`conda activate ncbi-cli`

An example of the to-copy command from NCBI

`datasets download genome accession GCF_000004515.6 --include gff3,rna,cds,protein,genome,seq-report`
