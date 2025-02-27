# wf-human-sv
This workflow provides an easy way to call structural variants in human genomic data.
## Introduction

The pipeline performs the following steps:
* Maps reads using [lra](https://github.com/ChaissonLab/LRA)
* Calls variants using [cuteSV](https://github.com/tjiangHIT/cuteSV)
* Filters variants by minimum/maximum length, read support, or type (e.g. insertion, deletion, etc.)
* Optionally evaluates yours calls against a truthset using [truvari](https://github.com/spiralgenetics/truvari)
## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-human-sv --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* A sorted, indexed VCF file containing the SV calls made.
* A sorted, indexed BAM file containing the alignments used to make the calls. 
* an HTML report document detailing the primary findings of the workflow.
## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)