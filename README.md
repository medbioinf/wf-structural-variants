# medbioinf/pangenomesv


[![GitHub Actions CI Status](https://github.com/medbioinf/wf-structural-variants/actions/workflows/nf-test.yml/badge.svg)](https://github.com/medbioinf/wf-structural-variants/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/medbioinf/wf-structural-variants/actions/workflows/linting.yml/badge.svg)](https://github.com/medbioinf/wf-structural-variants/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/medbioinf/wf-structural-variants)

## Introduction

**medbioinf/pangenomesv** is a bioinformatics pipeline that ...

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/community/brand/workflow-schematics#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

> [!NOTE]
> This repository contains a Nextflow-adapted version of [Swave](https://github.com/songbowang125/Swave). It will be migrated to a dedicated repository with a published container image. Currently, the pipeline only runs with Docker, and the image must be built locally.

### 1. Build the Docker Container Locally

Before running the pipeline, you must build the required Swave environment image locally. The tag **must** match the following name exactly so Nextflow can recognize it:

```bash
docker build -t quay.io/swave:latest .
```

### 2. Run with Included Test Data

The repository comes with small, pre-configured test data (found under `assets/testdata/`). You can perform a minimal test run using the test profile (configured under `conf/test.config`):

```bash
nextflow run main.nf -profile test,docker --outdir test_results
```

### 3. Run with Your Own Data

To execute the pipeline with custom data, you need to provide a reference genome and an input samplesheet (`.csv`).

#### 3.1. Prepare Input Files

By convention, it is recommended to organize your input files in a structured root directory:

- Save your assemblies (FASTA format) inside `data/assemblies/`.
- Place your reference genome FASTA anywhere accessible (e.g., directly under `data/` or in a `data/reference/` directory).

(The `data/` directory is gitignored and must be created locally.)

#### 3.2. Generate the Samplesheet

The pipeline requires a three-column samplesheet (`sample,haplotype,fasta`).

**Option 1: Automatic Generation**

You can run the provided automated Python script to scan your assembly directory and create the samplesheet:

```bash
python3 scripts/create_samplesheet.py --dir data/assemblies --out data/samplesheet.csv
```

(`--dir` and `--out` default to these values. Add `--allow-unphased` if your assemblies are unphased, to treat them as haplotype 0.)

**Option 2: Manual Creation**

Alternatively, create a `samplesheet.csv` manually with the following format:

```
sample,haplotype,fasta
assembly1,1,/path/to/assemblies/assembly1_hap1.fa
assembly1,2,/path/to/assemblies/assembly1_hap2.fa
assembly2,1,/path/to/assemblies/assembly2_hap1.fa
assembly2,2,/path/to/assemblies/assembly2_hap2.fa
...
```

#### 3.3. Run the Workflow

Run the pipeline by passing the paths to your generated samplesheet and reference genome:

```bash
nextflow run main.nf \
   -profile docker \
   --input <path_to_samplesheet.csv> \
   --fasta <path_to_reference_fasta> \
   --outdir <output_directory>
```

(`--outdir` defaults to `results` if not specified.)

For more detailed information on all available pipeline parameters, run:

```bash
nextflow run main.nf --help
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

medbioinf/pangenomesv was originally written by Jonah Kapski.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use medbioinf/pangenomesv for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
