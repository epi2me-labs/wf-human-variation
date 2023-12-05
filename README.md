# Human variation workflow

Basecalling; SNV, SV and CNV calling, modified base calling, and STR genotyping of human samples.



## Introduction

This repository contains a [nextflow](https://www.nextflow.io/) workflow
for analysing variation in human genomic data. Specifically this workflow can
perform the following:

* basecalling of FAST5 (or POD5) sequencing data
* diploid variant calling
* structural variant calling
* analysis of modified base calls
* copy number variant calling
* short tandem repeat (STR) expansion genotyping

The wf-human-variation workflow consolidates the small variant calling from the previous wf-human-snp, structural variant calling from wf-human-sv, CNV calling from wf-cnv (all of which are now deprecated), as well as performing STR expansion genotyping.
This pipeline performs the steps of the four pipelines simultaneously and the results are generated and output in the same
way as they would have been had the pipelines been run separately.




## Compute requirements

Recommended requirements:

+ CPUs = 32
+ Memory = 128GB

Minimum requirements:

+ CPUs = 12
+ Memory = 32GB

Approximate run time: Variable depending on whether it is targeted sequencing or whole genome sequencing, as well as coverage and the individual analyses requested. For instance, a 90X human sample run (options: `--snp --sv --mod --str --cnv --phased_vcf --joint_phasing --phase_mod --sex male`) takes less than 8h with recommended resources.

ARM processor support: False




## Install and run

These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).  

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore nextflow will need to be installed before attempting to run the workflow. 

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or [singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of the required software.
Both methods are automated out-of-the-box provided either docker or singularity is installed.
This is controlled by the [`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles) parameter as exemplified in the example below. 

It is not required to clone or download the git repository in order to run the workflow. 
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This  will pull the repository in to the assets folder of nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-human-variation –help 
```
A demo dataset is provided for testing of the workflow. It can be downloaded using: 
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-human-variation/wf-human-variation-demo.tar.gz 
tar -xzvf wf-human-variation-demo.tar.gz
```
The workflow can be run with the demo data using: 
```
nextflow run epi2me-labs/wf-human-variation \ 
            --bam 'wf-human-variation-demo/demo.bam' \
            --basecaller_cfg 'clair3:dna_r10.4.1_e8.2_400bps_hac_prom' \
            --mod \
            --ref 'wf-human-variation-demo/demo.fasta' \
            --sample_name 'DEMO' \
            --snp \
            --sv \
            -profile standard 
```
For further information about running a workflow on the command line see [our workflow quick start guide](https://labs.epi2me.io/wfquickstart/).



## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Inputs

### Workflow Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sv | boolean | Call for structural variants. | If this option is selected, structural variant calling will be carried out using Sniffles2. | False |
| snp | boolean | Call for small variants | If this option is selected, small variant calling will be carried out using Clair3. | False |
| cnv | boolean | Call for copy number variants. | If this option is selected, copy number variant calling will be carried out using QDNAseq. | False |
| str | boolean | Enable Straglr to genotype STR expansions. | If this option is selected, genotyping of STR expansions will be carried out using Straglr. This sub-workflow is only compatible with genome build 38. | False |
| mod | boolean | Enable output of modified calls to a bedMethyl file [requires input BAM with Ml and Mm tags] | This option is automatically selected and aggregation of modified calls with be carried out using modkit if Ml and Mm tags are found. Disable this option to prevent output of a bedMethyl file. | False |


### Main options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_name | string | Sample name to be displayed in workflow outputs. |  | SAMPLE |
| fast5_dir | string | Directory containing FAST5 signal for basecalling. | This directory will be searched recursively. All FAST5 or POD5 files (depending on which extension you select in the Basecalling Options) in this directory or any subdirectory (no matter how deep) will be basecalled. You may choose to provide a FAST5 directory, or a BAM/CRAM, but not both. |  |
| bam | string | Path to a BAM (or CRAM) containing aligned or unaligned reads. | You may choose to provide a BAM/CRAM, or a FAST5 directory for basecalling, but not both. |  |
| ref | string | Path to a reference FASTA file. | Reference against which to compare reads for variant calling. |  |
| old_ref | string | Reference FASTA file for CRAM input (only required if the CRAM requires realignment) | You do not need to provide this unless the workflow specifically asks you to. If your input CRAM headers do not match the metadata of the input reference, the workflow will assume you want to realign your reads to the new input reference. CRAM files are compressed using the reference, so the read sequences cannot be realigned without the old reference. |  |
| basecaller_cfg | string | Name of the model to use for converting signal and selecting a small variant calling model. | Required for basecalling and/or small variant calling. The basecaller configuration is used to automatically select the appropriate small variant calling model. The model list shows all models that are compatible for small variant calling with this workflow. Only a subset of these models can be used for basecalling with the workflow, these are shown at the top of the list. Models that begin with 'clair3:' cannot be used for basecalling and are included to allow SNP calling on existing datasets. You should select 'custom' to override the basecaller_cfg with basecaller_model_path. | dna_r10.4.1_e8.2_400bps_sup@v4.1.0 |
| bed | string | An optional BED file enumerating regions to process for variant calling. |  |  |
| annotation | boolean | SnpEff annotation. | If this option is unselected, VCFs will not be annotated with SnpEff. | True |
| out_dir | string | Directory for output of all workflow results. |  | output |


### Structural variant calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| tr_bed | string | Input BED file containing tandem repeat annotations for the reference genome. | Providing a tandem repeat BED can improve calling in repetitive regions. An appropriate tandem repeat BED can be downloaded for your reference genome [from the Sniffles2 repository](https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations). |  |


### Structural variant benchmarking options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sv_benchmark | boolean | Benchmark called structural variants. | If this option is selected, automated benchmarking of structural variant calls will be carried out using Truvari. | False |


### Small variant calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| phase_vcf | boolean | Output phasing information in VCF (SNP and SV workflows). | By default the final VCF outputs from SNP and SV workflows are not phased, although some intermediates are approximately phased. | False |
| include_all_ctgs | boolean | Call for variants on all sequences in the reference, otherwise small variants will only be called on chr{1..22,X,Y}. | Enabling this option will call for variants on all contigs of the input reference sequence. Typically this option is not required as standard human reference sequences contain decoy and unplaced contigs that are usually omitted for the purpose of variant calling. This option might be useful for non-standard reference sequence databases. | False |


### Modified base calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| phase_mod | boolean | Generate phased modification data. | Enabling this will output three [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/) files of modified sites for each haplotype (named SAMPLE.N, where N can be 1, 2 or untagged) by leveraging phased reads generated using the small variant calling subworkflow. This option will automatically enable the small variant and modified subworkflows. | False |
| force_strand | boolean | Require modkit to call strand-aware modifications. |  | False |


### Basecalling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| remora_cfg | string | Name of the model to use for calling modified bases. | Required for calling modified bases while basecalling. The model list only shows models that are compatible with this workflow. You should select 'custom' to override the remora_cfg with remora_model_path. |  |
| dorado_ext | string | File extension for Dorado input. | Set this to fast5 if you have not converted your fast5 to pod5. It is recommended to [convert existing fast5 files to pod5 for use with Dorado](https://github.com/nanoporetech/pod5-file-format/blob/master/python/README.md#pod5-convert-from-fast5). | fast5 |
| basecaller_basemod_threads | number | Number of threads to use for base modification calling. | You must set this to > 0 when using a modbase aware model. Modbase calling does not require much additional CPU and should be set carefully when using GPU servers with a small number of CPUs per GPU. | 2 |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| depth_intervals | boolean | Output a bedGraph file with entries for each genomic interval featuring homogeneous depth. | The output [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) file will have an entry for each genomic interval in which all positions have the same alignment depth. By default this workflow outputs summary depth information from your aligned reads. Per-base depth outputs are slower to generate but may be required for some downstream applications. | False |
| GVCF | boolean | Enable to output a gVCF file in addition to the VCF outputs (experimental). | By default the the workflow outputs a VCF file containing only records where a variant has been detected. Enabling this option will output additionally a gVCF with records spanning all reference positions regardless of whether a variant was detected in the sample. | False |
| joint_phasing | boolean | Jointly phase SNPs, small indels and, optionally, SVs. | By default, the phasing is performed for SNPs and SVs separately. This option allows the joint phasing of both variant classes at the end of the process, allowing consistency in the method and of the resulting phased sites. This step uses exclusively longphase. | False |
| downsample_coverage | boolean | Downsample the coverage to along the genome. | This options will trigger a downsampling of the read alignments to the target coverage specified by --downsample_coverage_target. Downsampling will make the workflow run faster but could lead to non-deterministic variant calls. | False |
| downsample_coverage_target | number | Average coverage or reads to use for the analyses. | This options will set the target coverage for the downsampling stage, if downsampling has been enabled. | 60 |


### Advanced basecalling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| qscore_filter | number | Mean qscore by which to filter reads. Inclusive such that reads with score >= qscore_filter are kept. | The mean qscore of reads is calculated by dorado and rounded to an integer by dorado and stored as a tag in dorado's SAM output. The pipeline separates reads into pass and fail categories based on this SAM tag. | 10 |
| basecaller_chunk_size | number | Number of input files to basecall in each basecalling process. |  | 25 |
| cuda_device | string | GPU device to use for basecalling [cuda:all] | For local execution this can be used to pin GPU tasks to one (or more) specific GPU devices. Use cuda:all to use all available GPU devices, or cuda:<idx>[,idx,...] where idx is an index number of the GPU device to use. | cuda:all |
| basecaller_args | string | Additional command line arguments to pass to the basecaller process. To provide custom arguments to `dorado` from command line proceed as follow: `--basecaller_args="--no-trim"` |  |  |


### Multiprocessing Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Set max number of threads to use for more intense processes (limited by config executor cpus) |  | 4 |
| ubam_map_threads | integer | Set max number of threads to use for aligning reads from uBAM (limited by config executor cpus) |  | 8 |
| ubam_sort_threads | integer | Set max number of threads to use for sorting and indexing aligned reads from uBAM (limited by config executor cpus) |  | 3 |
| ubam_bam2fq_threads | integer | Set max number of threads to use for uncompressing uBAM and generating FASTQ for alignment (limited by config executor cpus) |  | 1 |
| merge_threads | integer | Set max number of threads to use for merging alignment files (limited by config executor cpus) |  | 4 |
| modkit_threads | integer | Total number of threads to use in modkit modified base calling (limited by config executor cpus) |  | 4 |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |






## Outputs

Outputs files may be aggregated including information for all             samples or provided per sample. Per sample files             will be prefixed with respective aliases and represented             below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Report of the alignment statistics | {{ alias }}.wf-human-alignment-report.html | Report summarising the results of the alignment statistics for the sample. | per-sample |
| Report of the SNP workflow | {{ alias }}.wf-human-snp-report.html | Report summarising the results of the SNP subworkflow for the sample. | per-sample |
| Report of the SV workflow | {{ alias }}.wf-human-sv-report.html | Report summarising the results of the SV subworkflow for the sample. | per-sample |
| Report of the CNV workflow | {{ alias }}.wf-human-cnv-report.html | Report summarising the results of the CNV subworkflow for the sample. | per-sample |
| Report of the STR workflow | {{ alias }}.wf-human-str-report.html | Report summarising the results of the short tandem repeat subworkflow for the sample. | per-sample |
| Short variant VCF | {{ alias }}.wf_snp.vcf.gz | VCF file with the SNPs for the sample. | per-sample |
| Structural variant VCF | {{ alias }}.wf_sv.vcf.gz | VCF file with the SVs for the sample. | per-sample |
| SNP and SV phased VCF | {{ alias }}.wf_human_variation.phased.vcf.gz | VCF file with the jointly phased SNPs and SVs for the sample. | per-sample |
| Copy number variants VCF | {{ alias }}.wf_cnv.vcf.gz | VCF file with the CNV for the sample. | per-sample |
| Modified bases BEDMethyl | {{ alias }}.wf_mods.bedmethyl.gz | BED file with the aggregated modification counts for the sample. | per-sample |
| Short tandem repeat VCF | {{ alias }}.wf_str.vcf.gz | VCF file with the STR sites for the sample. | per-sample |
| Alignment file | {{ alias }}.{{ format }} | BAM or CRAM file with the aligned reads for the sample, generated when the input file is unaligned. | per-sample |
| Alignment file index | {{ alias }}.{{ format }}.{{ index_format }} | The index of the resulting BAM or CRAM file with the reads for the sample, generated when the input file is unaligned. | per-sample |
| Haplotagged alignment file | {{ alias }}.haplotagged.{{ format }} | BAM or CRAM file with the haplotagged reads for the sample. | per-sample |
| Haplotagged alignment file index | {{ alias }}.haplotagged.{{ format }}.{{ index_format }} | The index of the resulting BAM or CRAM file with the haplotagged reads for the sample. | per-sample |




## Pipeline overview

The workflow is composed of 6 distinct subworkflows, each enabled by a command line option:

* Basecalling: `--fast5_dir <input_dir>`
* SNP calling: `--snp`
* SV calling: `--sv`
* Analysis of modified bases: `--mod`
* CNV calling: `--cnv`
* STR genotyping: `--str`

Subworkflows where the relevant option is omitted will not be run.

### 1. Input and data preparation

The workflow relies on three primary input files:
1. A reference genome in [FASTA format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
2. Sequencing data for the sample in the form of:
    1. A single [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf) for the sample (either aligned or unaligned) or
    2. [POD5](https://github.com/nanoporetech/pod5-file-format)/[FAST5](https://github.com/nanoporetech/ont_fast5_api) files

The BAM file can be generated from:
1. [POD5](https://github.com/nanoporetech/pod5-file-format)/[FAST5](https://github.com/nanoporetech/ont_fast5_api) files using the [wf-basecalling](https://github.com/epi2me-labs/wf-basecalling) workflow, or
2. [fastq](https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#fastq) files using [wf-alignment](https://github.com/epi2me-labs/wf-alignment).
Both workflows will generate aligned BAM file that is ready to be used with `wf-human-variation`.

### 2. Basecalling
The basecalling is performed automatically when the user provides `pod5`/`fast5` inputs through the `--fast5_dir`.
Basecalling with [Dorado](https://github.com/nanoporetech/dorado) requires an NVIDIA GPU with [Volta architecture or newer](https://www.nvidia.com/en-gb/technologies/) and at least 8 GB of vRAM.

#### Windows
Windows should not be considered as a supported operating system for wf-basecalling as we do not directly support configuration of accelerated computing through WSL2 and Docker.
Although we do not offer support, it is possible to set up Docker to use GPUs for most versions of Windows 11 and some versions of Windows 10 and we direct users to the [CUDA on WSL User Guide](https://docs.nvidia.com/cuda/wsl-user-guide/index.html).
Users should take note of the support constraints section to ensure their environment is suitable before following the guidance. **Do not install an NVIDIA driver into your WSL2 environment**.
Users are encouraged to download Dorado for Windows from the [Dorado GitHub repository](https://github.com/nanoporetech/dorado#installation).

#### MacOS
MacOS should not be considered as a supported operating system for wf-basecalling as we do not support accelerated computing through Docker on MacOS.
On MacOS, GPU support through Docker remains in technical infancy. In addition, the containers we provide will not be able to leverage the M1 and M2 architecture and will not run as performantly as if Dorado had been run natively.
Users are encouraged to download Dorado for MacOS directly from the [Dorado GitHub repository](https://github.com/nanoporetech/dorado#installation).

#### Linux
When using Docker for accelerated computing on Linux, you will need the `nvidia-container-toolkit` installed.
If you observe the error "could not select device driver with capabilities gpu", you should follow the instructions to install `nvidia-container-toolkit` [here](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html#setting-up-nvidia-container-toolkit).
You will need to follow the steps to:

- Setup the package repository and the GPG key (ignore the box about experimental releases)
- Update package listings
- Install nvidia-container-toolkit
- Configure the Docker daemon to recognize the NVIDIA Container Runtime
- Restart the Docker daemon to complete the installation after setting the default runtime

By default, workflows are configured to run GPU tasks in serial. That is, only one basecalling task will be run at a time. This is to prevent the GPU from running out of memory on local execution.
When running workflows on a cluster, or in a cloud where GPU resources are isolated from one another, users should specify `-profile discrete_gpus` as part of the command invocation. This will allow for parallel execution of GPU tasks.

#### Options
Ensure to provide the appropriate model for your data with `--basecaller_cfg`. If you are interested in analyzing the modifications, you should also provide the appropriate remora model with `--remora_cfg`. Refer to the [Dorado documentation](https://github.com/nanoporetech/dorado#available-basecalling-models) for a list of available basecalling models.

If basecalling and alignment were conducted, the workflow will output two sorted, indexed CRAMs of basecalls aligned to the provided references, with reads separated by their quality score:
    * `<sample_name>.pass.cram` contains reads with `qscore >= threshold` (only pass reads are used to make downstream variant calls)
    * `<sample_name>.fail.cram` contains reads with `qscore < threshold`
Take care to retain the input reference when basecalling or alignment has been performed as CRAM files cannot be read without the corresponding reference!

### 3. Data QC and pre-processing
The workflow starts by performing multiple checks of the input BAM file, as well as computing:
1. depth of sequencing with [mosdepth](https://github.com/brentp/mosdepth);
2. read alignment statistics with [fastcat](https://github.com/epi2me-labs/fastcat).

After computing the coverage, the workflow will check that the input BAM file has a depth greater than `--bam_min_coverage`.
In case the user specify `--bam_min_coverage 0`, the check will be skipped and the workflow will proceed directly to the downstream analyses.
Some components work better withing certain ranges of coverage, and the user might achieve better results by providing a target coverage to downsample to. The user can set `--downsample_coverage true` to enable the downsampling of the reads, and `--downsample_coverage_target {{ X }}` to specify the target coverage (default: 60x).

### 4. Small variant calling with Clair3

The workflow implements a deconstructed version of [Clair3](https://github.com/HKU-BAL/Clair3) (v1.0.4) to call germline variants.
This workflow takes advantage of the parallel nature of Nextflow, providing optimal efficiency in high-performance, distributed systems. The workflow will automatically call small variants (SNPs and indels), collect statistics, annotate them with [SnpEff](https://pcingola.github.io/SnpEff/) (and additionally for SNPs, ClinVar details), and create a report summarising the findings.

If desired, the workflow can perform phasing of structural variants by using the `--phase_vcf` option. This will lead the workflow to use [longphase](https://github.com/twolinin/longphase) to perform phasing of the variants, with the option to use [whatshap](https://whatshap.readthedocs.io/) instead by setting `--use_longphase false --use_longphase_intermediate false`. The phasing will also generate a GFF file with the annotation of the phase blocks, facilitating the detection of these within genome visualizers.

### 5. Structural variant (SV) calling with Sniffles2

The workflow allows for calling of SVs using long-read sequencing data with [Sniffles2](https://github.com/fritzsedlazeck/Sniffles).
The workflow will perform SV calling, filtering and generation of a report.
Optionally, the workflow can also evaluate calls on HG002 against a truth set (provided the input data was aligned to HG19).
The SV workflow takes an optional `--tr_bed` option to specify tandem repeats in the reference sequence --- see the [sniffles](https://github.com/fritzsedlazeck/Sniffles) documentation for more information.
SVs can be phased using `--phase_vcf`. However, this will cause the workflow to run SNP analysis, as SV phasing relies on the haplotagged reads generated in this stage.

### 6. Modified base calling with modkit

Modified base calling can be performed by specifying `--mod`. The workflow will call modified bases using [modkit](https://github.com/nanoporetech/modkit). 
The workflow will automatically check whether the files contain the appropriate `MM`/`ML` tags, required for running [modkit pileup](https://nanoporetech.github.io/modkit/intro_bedmethyl.html). If the tags are not found, the workflow will not run the individual analysis, but will still run the other subworkflows requested by the user.
The default behaviour of the workflow is to run modkit with the `--cpg --combine-strands` options set. It is possible to report strand-aware modifications by providing `--force_strand`, which will trigger modkit to run in default mode. The resulting bedMethyl will include modifications for each site on each strand separately.
The modkit run can be fully customized by providing `--modkit_args`. This will override any preset, and allow full control over the run of modkit.
Haplotype-resolved aggregated counts of modified bases can be obtained with the `--phase_mod` option. This will generate three distinct BEDMethyl files with the naming pattern `{{ alias }}_{{ haplotype }}.wf_mods.bedmethyl.gz`, where `haplotype` can be `1`, `2` or `ungrouped`.

### 7. Copy number variants (CNV) calling with QDNASeq

CNV calling is performed using [QDNAseq](https://github.com/ccagc/QDNAseq). This workflow is compatible with genome builds hg19/GRCh37 or hg38/GRCh38.
In addition to the VCF of CNV calls, the workflow emits QDNAseq-generated plots and BED files of both raw read counts per bin and corrected, normalised, and smoothed read counts per bin.

### 8. Short tandem repeat (STR) genotyping with Straglr

STR genotyping is performed using a fork of [straglr](https://github.com/philres/straglr). This workflow is compatible with genome build hg38/GRCh38.
The STR workflow takes a required `--sex` option which is `male` or `female`. If `--sex` is not specified, the workflow will default to `female`. Please be aware that incorrect sex assignment will result in the wrong number of calls for all repeats on chrX.
In addition to a gzipped VCF file containing STRs found in the dataset, the workflow emits a TSV straglr output containing reads spanning STRs, and a haplotagged BAM. 

### 9. Joint phasing of SNP, small Indels and SV
The workflow can perform joint physical phasing with `longphase` of SNP, Indels and SVs by setting the `--joint_phasing` option.
This will trigger the generation of a single phased VCF file, that will include all the variants from both the `--snp` and `--sv`
workflows.

### 9. Variant annotation
Annotation will be performed automatically by the SNP and SV subworkflows, and can be disabled by the user with `--annotation false`. The workflow will annotate the variants using [SnpEff](https://pcingola.github.io/SnpEff/), and currently only support the human hg19 and hg38 genomes. Additionally, the workflow will add the [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) annotations for the SNP variants.

Running the workflow on non-human samples will require this option to be disabled.




## Troubleshooting

+ You should ask your system administrator if you need to configure any additional options to leverage GPUs on your cluster. For example, you may need to provide a special string to the workflow's `--cuda_device` option to ensure tasks use the GPU assigned to them by the job scheduler.
+ Annotations for `--snp` and `--sv` are generated using [SnpEff](https://pcingola.github.io/SnpEff/). For `--snp`, additional [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) annotations are displayed in the report where available (please note, the report will not display any variants classified as 'Benign' or 'Likely benign', however these variants will be present in the
output VCF).
+ Specifying a suitable [tandem repeat BED for your reference](https://raw.githubusercontent.com/fritzsedlazeck/Sniffles/master/annotations/) with `--tr_bed` can improve the accuracy of SV calling.
+ Aggregation of modified calls with `--mod` requires data to be basecalled with a model that includes base modifications, providing the `MM` and `ML` BAM tags
+ Refer to the [Dorado documentation](https://github.com/nanoporetech/dorado#available-basecalling-models) for a list of available basecalling models
+ CRAM files generated within the workflow cannot be read without the corresponding reference
+ The STR workflow performs genotyping of specific repeats, which can be found [here](https://github.com/epi2me-labs/wf-human-variation/blob/master/data/wf_str_repeats.bed).
+ While designed to work on human genomes, the workflow can be run on non-human species by setting `--cnv false --str false --annotation false`.
+ Ensure that the provided reference and BED files use the same chromosome coding (for example, that they both have the `chr` prefix, or they both to not have it).
+ If unaligned reads were provided, the workflow will output a CRAM file (or BAM if the user runs the `--cnv` option) containing the alignments used to make the downstream variant calls



## FAQ's

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-template/issues) page or start a discussion on the [community](https://community.nanoporetech.com/). 



## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)
+ [blogpost](https://labs.epi2me.io/copy-number-calling/) and [CNV workflow documentation](https://github.com/epi2me-labs/wf-cnv) for more information on running the copy number calling subworkflow.
+ [blogpost](https://labs.epi2me.io/human-targeted-analysis/) on how to run a targeted BRCA Gene Analysis with `wf-human-variation`
+ [blogpost](https://labs.epi2me.io/giab-2023.05/) on the generation of four Genome in a Bottle samples processed with `wf-human-variation`

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



