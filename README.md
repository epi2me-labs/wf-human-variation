# Human variation workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow
for analysing variation in human genomic data. Specifically this workflow can
perform the following:

* basecalling of FAST5 (or POD5) sequencing data
* diploid variant calling
* structural variant calling
* aggregation of modified base counts
* copy number variant calling
* short tandem repeat (STR) expansion genotyping

The wf-human-variation workflow consolidates the small variant calling from the
previous wf-human-snp, structural variant calling from wf-human-sv, CNV calling from wf-cnv (all
of which are now deprecated), as well as performing STR expansion genotyping. This pipeline performs the steps of the four
pipelines simultaneously and the results are generated and output in the same
way as they would have been had the pipelines been run separately.




## Introduction

This workflow uses [Clair3](https://www.github.com/HKU-BAL/Clair3) for calling small
variants from long reads. Clair3 makes the best of two methods: pileup (for fast
calling of variant candidates in high confidence regions), and full-alignment
(to improve precision of calls of more complex candidates).

This workflow uses [sniffles2](https://github.com/fritzsedlazeck/Sniffles) for
calling structural variants.

This workflow uses [modbam2bed](https://github.com/epi2me-labs/modbam2bed) to
aggregate modified base counts into a [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/) file.

This workflow uses [Dorado](https://github.com/nanoporetech/dorado/tree/master/dorado)
for basecalling `pod5` or `fast5` signal data.

This workflow uses [QDNAseq](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html) for calling copy number variants.

This workflow uses a fork of [Straglr](https://github.com/philres/straglr) for genotyping short tandem repeat expansions.




## Quickstart

The workflow uses [Nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such Nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://sylabs.io/singularity/) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit our website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-human-variation --help
```

to see the options for the workflow.

**Download demonstration data**

A small test dataset is provided for the purposes of testing the workflow software,
it can be downloaded using:

```
wget -O demo_data.tar.gz \
    https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-human-variation/demo_data.tar.gz
tar -xzvf demo_data.tar.gz
```

The basecalling, SNP, SV, 5mC aggregation, and CNV workflows are all independent and can be
run in isolation or together using options to activate them. The STR workflow can also be run independently but will trigger the SNP workflow to run first, as a phased VCF is required to haplotag the input BAM file in order to successfully perform STR genotyping.

The SNP and SV workflows can be run with the demonstration data using:

```
OUTPUT=output
nextflow run epi2me-labs/wf-human-variation \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --snp --sv \
    --bam demo_data/demo.bam \
    --bed demo_data/demo.bed \
    --ref demo_data/demo.fasta \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0'  \
    --sample_name MY_SAMPLE \
    --out_dir ${OUTPUT}
```

Each subworkflow is enabled with a command line option:

* Basecalling: `--fast5_dir <input_dir>`
* SNP calling: `--snp`
* SV calling: `--sv`
* Methylation aggregation: `--methyl`
* CNV calling: `--cnv`
* STR genotyping: `--str`

Subworkflows where the relevant option is omitted will not be run.

Some subworkflows have additional required options:

* The SV workflow takes an optional `--tr_bed` option to specify tandem
repeats in the reference sequence --- see the [sniffles](https://github.com/fritzsedlazeck/Sniffles)
documentation for more information.

* The STR workflow takes a required `--sex` option which is `male` or `female`. If `--sex` is not specified, the workflow will default to `female`. Please be aware that incorrect sex assignment will result in the wrong number of calls for all repeats on chrX. 

To enable the 5mC aggregation step use the `--methyl` option. For this
step to produce meaningful output the input BAM file must have been produced
by a basecaller capable of emitting the 5mC calls.

This brings us to activating the basecalling workflow. To run all the above
including basecalling:

```
OUTPUT=output
nextflow run epi2me-labs/wf-human-variation \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --snp --sv --methyl \
    --fast5_dir path/to/fast5/dir \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0'  \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v4.1.0_5mCG_5hmCG@v2' \
    --bed path/to.bed \
    --ref path/to.fasta \
    --out_dir ${OUTPUT}
```

**Genome build compatibility**

The workflow carries out a check to determine the version of the human genome build used during alignment, as certain subworkflows
are only compatible with specific genome versions:

* By default, `--snp`, `--sv`, and `--phase_methyl` require either hg19/GRCh37 or hg38/GRCh38 to enable generation of annotations using SnpEff.
However, by disabling annotations with `--skip_annotation`, these subworkflows can be run with other human genome builds (and non-human genomes).
* `--str`: requires genome build hg38/GRCh38.
* `--cnv`: requires genome builds hg19/GRCh37 or hg38/GRCh38.


**Workflow outputs**

The primary outputs of the workflow include:

* a gzipped VCF file containing annotated SNPs found in the dataset (`--snp`)
* a gzipped VCF file containing annotated SVs called from the dataset (`--sv`)
* a gzipped [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/) file aggregating modified CpG base counts (`--methyl`)
* a VCF of CNV calls, QDNAseq-generated plots, and BED files of both raw read counts per bin and corrected, normalised, and smoothed read counts per bin (`--cnv`)
* a gzipped VCF file containing STRs found in the dataset, TSV Straglr output containing reads spanning STRs, and a haplotagged BAM (`--str`)
* an HTML report detailing the primary findings of the workflow, for SNP, SV, CNV calling, and STR genotyping
* if basecalling and alignment was conducted, the workflow will output two sorted, indexed CRAMs of basecalls aligned to the provided references, with reads separated by their quality score:
    * `<sample_name>.pass.cram` contains reads with `qscore >= threshold` (only pass reads are used to make downstream variant cals)
    * `<sample_name>.fail.cram` contains reads with `qscore < threshold`
* if unaligned reads were provided, the workflow will output a CRAM file containing the alignments used to make the downstream variant calls

The secondary outputs of the workflow include:
* `{sample_name}.mapula.csv` and `{sample_name}.mapula.json` provide basic alignment metrics (primary and secondary counts, read N50, median accuracy)
* `mosdepth` outputs include:
    * `{sample_name}.mosdepth.global.dist.txt`: a cumulative distribution indicating the proportion of total bases for each and all reference sequences [more info](https://github.com/brentp/mosdepth#distribution-output)
    * `{sample_name}.regions.bed.gz`: reports the mean coverage for each region in the provided BED file
    * `{sample_name}.thresholds.bed.gz`: reports the number of bases in each region that are covered at or above each threshold value (1, 10, 20, 30X) [more info](https://github.com/brentp/mosdepth#thresholds)
* `{sample_name}.readstats.tsv.gz`: a gzipped TSV summarising per-alignment statistics produced by [`bamstats`](https://github.com/epi2me-labs/fastcat#bamstats)

**Workflow tips**

- Users familiar with `wf-human-snp` and `wf-human-sv` are recommended to familiarise themselves with any parameter changes by using `--help`, in particular:
    - All arms of the variation calling workflow use `--ref` (not `--reference`) and `--bed` (not `--target_bedfile`)
- Annotations for `--snp` and `--sv` are generated using [SnpEff](https://pcingola.github.io/SnpEff/), with additional [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) annotations displayed in the report where available (please note, the report will not display any variants classified as 'Benign' or 'Likely benign', however these variants will be present in the
output VCF).
- Specifying a suitable [tandem repeat BED for your reference](https://raw.githubusercontent.com/fritzsedlazeck/Sniffles/master/annotations/) with `--tr_bed` can improve the accuracy of SV calling.
- Aggregation of methylation calls with `--methyl` requires data to be basecalled with a model that includes base modifications, providing the `MM` and `ML` BAM tags
- Refer to the [Dorado documentation](https://github.com/nanoporetech/dorado#available-basecalling-models) for a list of available basecalling models
- Take care to retain the input reference when basecalling or alignment has been performed as CRAM files cannot be read without the corresponding reference!
- Refer to our [blogpost](https://labs.epi2me.io/copy-number-calling/) and [CNV workflow documentation](https://github.com/epi2me-labs/wf-cnv) for more information on running the copy number calling subworkflow.
- The STR workflow performs genotyping of specific repeats, which can be found [here](https://github.com/epi2me-labs/wf-human-variation/blob/master/data/wf_str_repeats.bed).
- The workflow can perform physical phasing of SNP, Indels and SVs using with the `--joint_phasing` option.

### Support for basecalling on GPU

This section will be kept up to date with latest advice for running our workflows on the GPU.

#### Prerequisites

Basecalling with `Dorado` requires an NVIDIA GPU with [Volta architecture or newer](https://www.nvidia.com/en-gb/technologies/) and at least 8 GB of vRAM.

#### Windows

Windows should not be considered as a supported operating systems for wf-basecalling as we do not directly support configuration of accelerated computing through WSL2 and Docker.
Although we do not offer support, it is possible to set up Docker to use GPUs for most versions of Windows 11 and some versions of Windows 10 and we direct users to the [CUDA on WSL User Guide](https://docs.nvidia.com/cuda/wsl-user-guide/index.html).
Users should take note of the support constraints section to ensure their environment is suitable before following the guidance. **Do not install an NVIDIA driver into your WSL2 environment**.
Users are encouraged to download Dorado for Windows from the [Dorado GitHub repository](https://github.com/nanoporetech/dorado#installation).

#### MacOS

MacOS should not be considered as a supported operating systems for wf-basecalling as we do not support accelerated computing through Docker on MacOS.
On MacOS, GPU support through Docker remains in technical infancy. In addition, the containers we provide will not be able to leverage the M1 and M2 architecture and will not run as performantly as if Dorado had been run natively.
Users are encouraged to download Dorado for MacOS directly from the [Dorado GitHub repository](https://github.com/nanoporetech/dorado#installation).

#### Linux

When using Docker for accelerated computing on Linux, you will need the `nvidia-container-toolkit` installed.
If you observe the error "could not select device driver with capabilities gpu", you should follow the instructions to install `nvidia-container-toolkit` [here](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html#setting-up-nvidia-container-toolkit). You will need to follow the steps to:

- Setup the package repository and the GPG key (ignore the box about experimental releases)
- Update package listings
- Install nvidia-container-toolkit
- Configure the Docker daemon to recognize the NVIDIA Container Runtime
- Restart the Docker daemon to complete the installation after setting the default runtime

By default, workflows are configured to run GPU tasks in serial. That is, only one basecalling task will be run at a time. This is to prevent the GPU from running out of memory on local execution.
When running workflows on a cluster, or in a cloud where GPU resources are isolated from one another, users should specify `-profile discrete_gpus` as part of the command invocation. This will allow for parallel execution of GPU tasks.
You should ask your system administrator if you need to configure any additional options to leverage GPUs on your cluster. For example, you may need to provide a special string to the workflow's `--cuda_device` option to ensure tasks use the GPU assigned to them by the job scheduler.




## Useful links

* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/products/docker-desktop)
* [Singularity](https://sylabs.io/singularity/)

