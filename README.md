# Workflow : ONT-VariantCalling

This workflow do VariantCalling on Oxford Nanopore Technologies data.

## Installation

### Clone the repo

```bash
git clone --recursive https://github.com/MobiDL/ONT-VariantCalling.git
```

### Install dependencies

__Conda is recommended:__ [see here to install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

You will need to create a new environment based on conda.

```bash
conda env create -f environment.yml
conda activate ONT-VariantCalling
```

Following specifications on [Clair](https://github.com/HKU-BAL/Clair#option-1-bioconda),
you will need to install `intervaltree` package :

```bash
pypy3 -m ensurepip
pypy3 -m pip install --no-cache-dir intervaltree==3.0.2
conda deactivate
```

## Configure your inputs

You can create your input file replacing editing the template or creating your own inputs file.

### Minimal

```json
{
	"ONT_VariantCalling.fastqPath": "String",
	"ONT_VariantCalling.refFa": "File",
	"ONT_VariantCalling.refFai": "File",
	"ONT_VariantCalling.modelPath": "String",
	"ONT_VariantCalling.outputPath": "String",
}
```

### Extended

A full option templates is provided (inputs.tpl.json).

This template is separating in 4 categories (blank line) :
1. Global pipeline inputs (i.e. minimal)
2. Global pipeline options
3. Specific tasks inputs
4. Specific tasks options

## Launch

### Local

```bash
conda activate ONT-VariantCalling
cromwell run \
	--inputs /path/to/inputs.json \
	ONT-VariantCalling.wdl
```
