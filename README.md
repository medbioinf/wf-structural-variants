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

**medbioinf/pangenomesv** is a bioinformatics pipeline for population-level structural variant (SV) detection using pangenome graphs. It takes long-read assemblies or raw long-read BAM files as input, builds a pangenome graph (Minigraph, PGGB, or Cactus), and applies [Swave](https://github.com/songbowang125/Swave), a sequence-to-image, deep-learning-based method, to classify simple and complex SVs directly from the graph.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/community/brand/workflow-schematics#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

> [!NOTE]
> This pipeline is built around [Swave](https://github.com/songbowang125/Swave), implemented as a modular set of tools in [structural-variants-swave](LINK_ZUM_SWAVE_REPO), containerized and currently available on Docker Hub as [`jonahkps/panswave`](https://hub.docker.com/r/jonahkps/panswave). Docker, Apptainer/Singularity, and Conda profiles are all supported.

### 1. Run with Included Test Data

The repository comes with small, pre-configured test data (found under `assets/testdata/`). You can perform a minimal test run using the test profile (configured under `conf/test.config`):

```bash
nextflow run main.nf -profile test,docker --outdir test_results
```

### 2. Run with Your Own Data

To execute the pipeline with custom data, you need to provide a reference genome and an input samplesheet (`.csv`).

#### 2.1. Prepare Input Files

It is recommended to organize your input files in a structured root directory:

- Save your assemblies (FASTA format) inside `data/assemblies/`, or raw long-read BAM directories inside `data/bams/`, if you want the pipeline to assemble them automatically with hifiasm.
- Place your reference genome FASTA anywhere accessible (e.g., directly under `data/` or in a `data/reference/` directory).

(The `data/` directory is gitignored and must be created locally.)

#### 2.2. Generate the Samplesheet

The pipeline requires a four-column samplesheet (`sample,haplotype,fasta,bam_dir`).

**Option 1: Automatic Generation**

You can run the provided automated Python script to scan your assembly and BAM directories and create the samplesheet:

```bash
python3 scripts/generate_samplesheet.py --assemblies_dir data/assemblies --bams_dir data/bams --out data/samplesheet.csv
```

(All paths default to the values shown above. Add `--allow_unphased` to treat files with unrecognized haplotypes as unphased (assigned to `0`) instead of skipping them. Use `--exclude_bams` or `--exclude_assemblies` to scan only one of the two directories.)

**Option 2: Manual Creation**

Alternatively, create a `samplesheet.csv` manually with the following format:
```bash
sample,haplotype,fasta,bam_dir
assembly1,1,/path/to/assemblies/assembly1_hap1.fa,
assembly1,2,/path/to/assemblies/assembly1_hap2.fa,
sample3,,,/path/to/bams/sample2_dir
```

#### 2.3. Run the Workflow

Run the pipeline by passing the paths to your generated samplesheet and reference genome:

```bash
nextflow run main.nf \
   -profile <docker/singularity/apptainer/conda> \
   --input <path_to_samplesheet.csv> \
   --fasta <path_to_reference_fasta> \
   --outdir <output_directory>
```

(`--outdir` defaults to `results` if not specified.)

### 3. Additional Notes

**Graph construction tools.** The pipeline supports three tools for pangenome graph construction, selected via `--graph_construction_tool`: `minigraph` (default, fast and incremental building possible), `pggb`, and `cactus` (both more fine-grained, also capturing SNPs and small indels, but slower and always require a full rebuild when adding new samples).

**Long-read assembly.** Since assembling raw long reads with hifiasm is resource-intensive, it might make sense to use `--assembly_only` to run this step on its own, then start the actual pipeline run pointing at the resulting assemblies.

**Incremental Minigraph graphs.** Minigraph does not track which assemblies are already represented in a given graph. When using `--minigraph_incremental` with `--gfa`, only the new assemblies to be added should be listed in the samplesheet. Re-including assemblies already present in the graph would result in minigraph re-aligning every provided assembly regardless of whether it's already incorporated, i.e. providing no benefit over rebuilding the graph from scratch Once the extended graph has been built, the actual pipeline run can then be started with the full samplesheet, passing the resulting GFA via `--gfa`.

**ANNOVAR annotation.** Gene/exon annotation via ANNOVAR is optional and requires your own, individually registered installation (`--annovar_dir`, pointing to a directory with `table_annovar.pl` and a populated `humandb/` folder, plus `--annovar_db` matching the downloaded database, e.g. `hg38` or `hs1`). ANNOVAR's license does not permit redistribution, so it cannot be bundled with this pipeline. See [ANNOVAR's registration page](https://annovar.openbioinformatics.org) to obtain your own copy. If `--annovar_dir` is not set, the annotation step is skipped entirely. Downloaded RefGene databases commonly reference chromosomes by RefSeq accession (e.g. `NC_060925.1`) rather than by name. Since this pipeline's VCFs use PanSN-formatted chromosome names (e.g. `CHM13#0#chr1`), the chromosome column of the downloaded database file may need to be remapped (e.g. `NC_060925.1` → `CHM13#0#chr1`) before annotation will produce meaningful results, otherwise ANNOVAR will report all variants as `intergenic`.

For more detailed information on all available pipeline parameters, run:

```bash
nextflow run main.nf --help
```

## Credits

medbioinf/pangenomesv was originally written by Jonah Kapski.

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
