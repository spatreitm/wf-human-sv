# wf-human-sv

This workflow provides an easy way to call structural variants in human genomic data.

The pipeline performs the following steps:
*  Maps reads using lra
*  Calls variants using cuteSV
*  Filters variants by minimum/maximum length, read support, or type (e.g. insertion, deletion, etc.)
*  Optionally evaluates yours calls against a truthset using truvari


## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

> See the sections below for installation of these prerequisites in various scenarios.
> It is not required to clone or download the git repository in order to run the workflow.

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-human-sv --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* a simple text file providing a summary of sequencing reads,
* an HTML report document detailing the primary findings of the workflow.


### Supported installations and GridION devices

Installation of the software on a GridION can be performed using the command

`sudo apt install ont-nextflow`

This will install a java runtime, Nextflow and docker. If *docker* has not already been
configured the command below can be used to provide user access to the *docker*
services. Please logout of your computer after this command has been typed.

`sudo usermod -aG docker $USER`

### Installation on Ubuntu devices

For hardware running Ubuntu the following instructions should suffice to install
Nextflow and Docker in order to run the workflow.

1. Install a Jva runtime environment (JRE):

   ```sudo apt install default-jre```

2. Download and install Nextflow may be downloaded from https://www.nextflow.io:

   ```curl -s https://get.nextflow.io | bash```

   This will place a `nextflow` binary in the current working directory, you 
   may wish to move this to a location where it is always accessible, e.g:

   ```sudo mv nextflow /usr/local/bin```

3. Install docker and add the current user to the docker group to enable access:

   ```
   sudo apt install docker.io
   sudo usermod -aG docker $USER
   ```

## Running the workflow

The `wf-human-sv` workflow can be controlled by the following parameters. The `fastq` parameter
is the most important parameter: it is required to identify the location of the
sequence files to be analysed.

**Parameters:**

- `fastq` specifies a *directory* or *file* path to FASTQ file(s) (required)
- `reference` specifies a *file* path to your reference genome (required)
- `mode` switch between standard or benchmark modes. (default: standard)
- `out_dir` the path for the output (default: output)

There are many more detailed settings that you can access by editing the config (See below).

In benchmarking mode, once the SV calls have been made they will be compared against a 
user-supplied truthset, generating numbers for precision and recall.

To run the workflow using Docker containers supply the `-profile standard`
argument to `nextflow run`:

> The command below uses test data available from the [github repository](https://github.com/epi2me-labs/wf-human-sv/tree/master/test_data)
> It can be obtained with `git clone https://github.com/epi2me-labs/wf-human-sv`.

```
# run the pipeline with the test data
OUTPUT=output
nextflow run epi2me-labs/wf-human-sv \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --sample my_sample
    --fastq test_data/reads.fastq.gz \
    --reference test_data/chr21.fna.gz
    --out_dir ${OUTPUT}
```

The output of the pipeline will be found in `./output` for the above
example. This directory contains the nextflow working directories alongside
the two primary outputs of the pipeline: a `seqs.txt` file containing a summary
of all reads, and a `report.html` file summarising the workflows calculations.

### Running the workflow with Conda

To run the workflow using conda rather than docker, simply replace 

    -profile standard 

with

    -profile conda

in the command above.

### Configuration and tuning

> This section provides some minimal guidance for changing common options, see
> the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for further details.

The default settings for the workflow are described in the configuration file `nextflow.config`
found within the git repository. The default configuration defines an *executor* that will 
use a specified maximum CPU cores (four at the time of writing) and RAM (eight gigabytes).

If the workflow is being run on a device other than a GridION, the available memory and
number of CPUs may be adjusted to the available number of CPU cores. This can be done by
creating a file `my_config.cfg` in the working directory with the following contents:

```
executor {
    $local {
        cpus = 4
        memory = "8 GB"
    }
}
```

and running the workflow providing the `-c` (config) option, e.g.:

```
# run the pipeline with custom configuration
nextflow run epi2me-labs/wf-human-sv \
    -c my_config.cfg \
    ...
```

The contents of the `my_config.cfg` file will override the contents of the default
configuration file. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html)
for more information concerning customized configuration.

**Using a fixed conda environment**

By default, Nextflow will attempt to create a fresh conda environment for any new
analysis (for reasons of reproducibility). This may be undesirable if many analyses
are being run. To avoid the situation a fixed conda environment can be used for all
analyses by creating a custom config with the following stanza:

```
profiles {
    // profile using conda environments rather than docker
    // containers
    fixed_conda {
        docker {
            enabled = false
        }
        process {
            withLabel:wf_human_sv {
                conda = "/path/to/my/conda/environment"
            }
            shell = ['/bin/bash', '-euo', 'pipefail']
        }
    }
}
```

and running nextflow by setting the profile to `fixed_conda`:

```
nextflow run epi2me-labs/wf-human-sv \
    -c my_config.cfg \
    -profile fixed_conda \
    ...
```

## Updating the workflow

Periodically when running the workflow, users may find that a message is displayed
indicating that an update to the workflow is available.

To update the workflow simply run:

    nextflow pull epi2me-labs/wf-human-sv

## Building the docker container from source

The docker image used for running the `wf-human-sv` workflow is available on
[dockerhub](https://hub.docker.com/repository/docker/ontresearch/wf-human-sv).
The image is built from the Dockerfile present in the git repository. Users
wishing to modify and build the image can do so with:

```
CONTAINER_TAG=ontresearch/wf-human-sv:latest

git clone https://github.com/epi2me-labs/wf-human-sv
cd wf-human-sv

docker build \
    -t ${CONTAINER_TAG} -f Dockerfile \
    --build-arg BASEIMAGE=ontresearch/base-workflow-image:v0.1.0 \
    .
```

In order to run the workflow with this new image it is required to give
`nextflow` the `--wfversion` parameter:

```
nextflow run epi2me-labs/wf-human-sv \
    --wfversion latest
```

## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
