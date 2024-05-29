# Human variation workflow

SNV, SV and CNV calling, modified base calling, and STR genotyping of human samples.



## Introduction

This repository contains a [nextflow](https://www.nextflow.io/) workflow
for analysing variation in human genomic data. Specifically this workflow can
perform the following:

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

+ CPUs = 16
+ Memory = 32GB

Approximate run time: Variable depending on whether it is targeted sequencing or whole genome sequencing, as well as coverage and the individual analyses requested. For instance, a 90X human sample run (options: `--snp --sv --mod --str --cnv --phased --sex male`) takes less than 8h with recommended resources.

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
nextflow run epi2me-labs/wf-human-variation -–help 
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
            --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac_prom' \
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



## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
The `--bam` input parameter for this workflow accepts a path to a single BAM file, or a folder containing multiple BAM files for the sample. A sample name can be supplied with `--sample`.

```
(i)                     (ii)    
input_reads.bam     ─── input_directory
                        ├── reads0.bam
                        └── reads1.bam
```



## Input parameters

### Workflow Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sv | boolean | Call for structural variants. | If this option is selected, structural variant calling will be carried out using Sniffles2. | False |
| snp | boolean | Call for small variants | If this option is selected, small variant calling will be carried out using Clair3. | False |
| cnv | boolean | Call for copy number variants. | If this option is selected, copy number variant calling will be carried out with either Spectre (default) or QDNAseq. To use QDNAseq instead of Spectre, use the option `--use_qdnaseq`. Spectre is only compatible with genome build hg38, and if QDNAseq is used, it is only compatible with genome builds hg37 and hg38. | False |
| str | boolean | Enable Straglr to genotype STR expansions. | If this option is selected, genotyping of STR expansions will be carried out using Straglr. This sub-workflow is only compatible with genome build hg38. | False |
| mod | boolean | Enable output of modified calls to a bedMethyl file [requires input BAM with Ml and Mm tags] | This option is automatically selected and aggregation of modified calls with be carried out using modkit if Ml and Mm tags are found. Disable this option to prevent output of a bedMethyl file. | False |


### Main options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_name | string | Sample name to be displayed in workflow outputs. |  | SAMPLE |
| bam | string | BAM or unaligned BAM (uBAM) files for the sample to use in the analysis. | This accepts one of two cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files. A sample name can be supplied with `--sample`. |  |
| ref | string | Path to a reference FASTA file. | Reference against which to compare reads for variant calling. |  |
| basecaller_cfg | string | Name of the model to use for selecting a small variant calling model. | The workflow will attempt to find the basecaller model from the headers of your input data. If the model cannot be found in the header, it must be provided with this option as the basecaller model is required for small variant calling. The basecaller model is used to automatically select the appropriate small variant calling model. The model list shows all models that are compatible for small variant calling with this workflow. You should select 'custom' to override the basecaller_cfg with clair3_model_path. | dna_r10.4.1_e8.2_400bps_sup@v4.1.0 |
| bam_min_coverage | number | Minimum read coverage required to run analysis. |  | 20 |
| bed | string | An optional BED file enumerating regions to process for variant calling. |  |  |
| annotation | boolean | SnpEff annotation. | If this option is unselected, VCFs will not be annotated with SnpEff. | True |
| phased | boolean | Perform phasing. | This option enables phasing of SV, SNP and modifications, depending on which sub-workflow has been chosen; see [README](README.md#9-phasing-variants) for more details. | False |
| include_all_ctgs | boolean | Call for variants on all sequences in the reference, otherwise small and structural variants will only be called on chr{1..22,X,Y,MT}. | Enabling this option will call for variants on all contigs of the input reference sequence. Typically this option is not required as standard human reference sequences contain decoy and unplaced contigs that are usually omitted for the purpose of variant calling. This option might be useful for non-standard reference sequence databases. | False |
| output_gene_summary | boolean | If set to true, the workflow will generate gene-level coverage summaries. | If set to true, a 4-column BED file must be supplied, where column 4 is the gene label. The workflow will generate a list of all genes in the BED and their percentage coverage at a range of thresholds (1x, 10x, 15x, 20x, and 30x), as well as the average coverage of each gene. | False |
| out_dir | string | Directory for output of all workflow results. |  | output |


### Structural variant calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| tr_bed | string | Input BED file containing tandem repeat annotations for the reference genome. | Providing a tandem repeat BED can improve calling in repetitive regions. An appropriate tandem repeat BED can be downloaded for your reference genome [from the Sniffles2 repository](https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations). |  |


### Structural variant benchmarking options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sv_benchmark | boolean | Benchmark called structural variants. | If this option is selected, automated benchmarking of structural variant calls will be carried out using Truvari. | False |


### Copy number variant calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| use_qdnaseq | boolean | Use QDNAseq for CNV calling. | Set this to true to use QDNASeq for CNV calling instead of Spectre. QDNAseq is better suited to shorter reads such as those generated from adaptive sampling experiments. | False |
| qdnaseq_bin_size | integer | Bin size for QDNAseq in kbp. | Pre-computed bin annotations are available for a range of bin sizes. Larger sizes reduce noise, however this may result in reduced sensitivity. | 500 |


### Modified base calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| force_strand | boolean | Require modkit to call strand-aware modifications. | By default strand calls are collapsed (strand reported as '.'). Enabling this will force stranding to be considered when calling modifications, creating one output per modification per strand and the report will be tabulated by both modification and strand. | False |


### Short tandem repeat expansion genotyping options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sex | string | Sex (XX or XY) to be passed to Straglr-genotype. | The sex determines how many calls will be obtained for all repeats on chrX. If not specified, the workflow will naively attempt to infer whether the sample carries XX or XY based on relative coverage of the allosomes. |  |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| depth_intervals | boolean | Output a bedGraph file with entries for each genomic interval featuring homogeneous depth. | The output [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) file will have an entry for each genomic interval in which all positions have the same alignment depth. By default this workflow outputs summary depth information from your aligned reads. Per-base depth outputs are slower to generate but may be required for some downstream applications. | False |
| GVCF | boolean | Enable to output a gVCF file in addition to the VCF outputs (experimental). | By default the the workflow outputs a VCF file containing only records where a variant has been detected. Enabling this option will output additionally a gVCF with records spanning all reference positions regardless of whether a variant was detected in the sample. | False |
| downsample_coverage | boolean | Downsample the coverage to along the genome. | This options will trigger a downsampling of the read alignments to the target coverage specified by --downsample_coverage_target. Downsampling will make the workflow run faster but could lead to non-deterministic variant calls. | False |
| downsample_coverage_target | number | Average coverage or reads to use for the analyses. | This options will set the target coverage for the downsampling stage, if downsampling has been enabled. | 60 |
| output_xam_fmt | string | Desired file format of alignment files created by alignment and phasing. | This setting controls the file format of (1) alignment files created by aligning or re-aligning an input BAM and (2) alignment files with haplotag information created during phasing of an input BAM. If using QDNASeq for CNV calling, the setting will be ignored for alignment or realignment as QDNASeq requires BAM input. | cram |


### Multiprocessing Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Set max number of threads to use for more intense processes (limited by config executor cpus) |  | 4 |
| ubam_map_threads | integer | Set max number of threads to use for aligning reads from uBAM (limited by config executor cpus) |  | 8 |
| ubam_sort_threads | integer | Set max number of threads to use for sorting and indexing aligned reads from uBAM (limited by config executor cpus) |  | 3 |
| ubam_bam2fq_threads | integer | Set max number of threads to use for uncompressing uBAM and generating FASTQ for alignment (limited by config executor cpus) |  | 1 |
| modkit_threads | integer | Total number of threads to use in modkit modified base calling (limited by config executor cpus) |  | 4 |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Report of the alignment statistics | {{ alias }}.wf-human-alignment-report.html | Report summarising the results of the alignment statistics for the sample. | per-sample |
| JSON file of some base statistics | {{ alias }}.stats.json | This JSON file contains base statistics on the reads, mappings, SNPs and SVs for the sample. | per-sample |
| Report of the SNP workflow | {{ alias }}.wf-human-snp-report.html | Report summarising the results of the SNP subworkflow for the sample. | per-sample |
| Report of the SV workflow | {{ alias }}.wf-human-sv-report.html | Report summarising the results of the SV subworkflow for the sample. | per-sample |
| Report of the CNV workflow | {{ alias }}.wf-human-cnv-report.html | Report summarising the results of the CNV subworkflow for the sample. | per-sample |
| Report of the STR workflow | {{ alias }}.wf-human-str-report.html | Report summarising the results of the short tandem repeat subworkflow for the sample. | per-sample |
| Short variant VCF | {{ alias }}.wf_snp.vcf.gz | VCF file with the SNPs for the sample. | per-sample |
| Structural variant VCF | {{ alias }}.wf_sv.vcf.gz | VCF file with the SVs for the sample. | per-sample |
| Copy number variants VCF | {{ alias }}.wf_cnv.vcf.gz | VCF file with the CNV for the sample. | per-sample |
| Modified bases BEDMethyl | {{ alias }}.wf_mods.bedmethyl.gz | BED file with the aggregated modification counts for the sample. | per-sample |
| Modified bases BEDMethyl (haplotype 1) | {{ alias }}.wf_mods.1.bedmethyl.gz | BED file with the aggregated modification counts for haplotype 1 of the sample. | per-sample |
| Modified bases BEDMethyl (haplotype 2) | {{ alias }}.wf_mods.2.bedmethyl.gz | BED file with the aggregated modification counts for haplotype 2 of the sample. | per-sample |
| Modified bases BEDMethyl (ungrouped) | {{ alias }}.wf_mods.ungrouped.bedmethyl.gz | BED file with the aggregated modification counts of non-haplotagged reads for the sample. | per-sample |
| Short tandem repeat VCF | {{ alias }}.wf_str.vcf.gz | VCF file with the STR sites for the sample. | per-sample |
| Alignment file | {{ alias }}.cram | CRAM or BAM file with the aligned reads for the sample, generated when the input file is unaligned. | per-sample |
| Alignment file index | {{ alias }}.cram.crai | The index of the resulting CRAM or BAM file with the reads for the sample, generated when the input file is unaligned. | per-sample |
| Haplotagged alignment file | {{ alias }}.haplotagged.cram | CRAM or BAM file with the haplotagged reads for the sample. | per-sample |
| Haplotagged alignment file index | {{ alias }}.haplotagged.cram.crai | The index of the resulting CRAM or BAM file with the haplotagged reads for the sample. | per-sample |
| Mean coverage for each region | {{ alias }}.regions.bed.gz | The mean coverage in the individual regions of the genome in BED format. | per-sample |
| Coverage per region above the given thresholds | {{ alias }}.thresholds.bed.gz | The BED reporting the number of bases in each region that are covered at or above each threshold values (1x, 10x, 20x and 30x). | per-sample |
| Distribution of the proportion of total bases covered by a given coverage value | {{ alias }}.mosdepth.global.dist.txt | The cumulative distribution indicating the proportion of total bases covered by a given coverage value, both genome-wide and by sequence. | per-sample |
| Mean coverage per sequence and target region | {{ alias }}.mosdepth.summary.txt | The summary of mean depths per chromosome and within specified regions per chromosome. | per-sample |
| BEDgraph of the single-base coverage | {{ alias }}.per-base.bedgraph.gz | The single-base coverage of the genome in BED graph format. | per-sample |
| Gene level coverage summary | SAMPLE.gene_summary.tsv | A table where each gene of the input BED file has columns describing the percentage of positions along the gene region that are covered to a given threshold, and a column with the average coverage. | per-sample |
| Haplocheck contamination summary | {{ alias }}.haplocheck.tsv | A table generated by [haplocheck](https://mitoverse.readthedocs.io/haplocheck/haplocheck/), with estimate of contamination from the MT genome. | per-sample |




## Pipeline overview

The workflow is composed of 6 distinct subworkflows, each enabled by a command line option:

* [SNP calling](#3-small-variant-calling-with-clair3): `--snp`
* [SV calling](#4-structural-variant-sv-calling-with-sniffles2): `--sv`
* [Analysis of modified bases](#5-modified-base-calling-with-modkit): `--mod`
* [CNV calling](#6-copy-number-variants-cnv-calling-with-qdnaseq): `--cnv`
* [STR genotyping](#7-short-tandem-repeat-str-genotyping-with-straglr): `--str`

Subworkflows where the relevant option is omitted will not be run.

### 1. Input and data preparation

The workflow relies on two primary input files:
1. A reference genome in [FASTA format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
2. Sequencing data for the sample in the form of a single [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf) or a folder of BAM files, either aligned or unaligned.

When analysing human data, we recommend using [human_g1k_v37.fasta.gz (FTP link)](ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz) or [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz (FTP link)](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz). For more information see [this blog post](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) which outlines potential pitfalls with the various flavours of human references.

The input BAM file can be generated using the [wf-basecalling](https://github.com/epi2me-labs/wf-basecalling/) workflow, which is up to date with the current dorado releases and models.

### 2. Data QC and pre-processing
The workflow starts by performing multiple checks of the input BAM file, as well as computing:
1. depth of sequencing with [mosdepth](https://github.com/brentp/mosdepth);
2. read alignment statistics with [fastcat](https://github.com/epi2me-labs/fastcat).

After computing the coverage, the workflow will check that the input BAM file has a depth greater than `--bam_min_coverage`.
In case the user specify `--bam_min_coverage 0`, the check will be skipped and the workflow will proceed directly to the downstream analyses.
Some components work better withing certain ranges of coverage, and the user might achieve better results by providing a target coverage to downsample to. The user can set `--downsample_coverage true` to enable the downsampling of the reads, and `--downsample_coverage_target {{ X }}` to specify the target coverage (default: 60x).

An optional gene summary may also be generated by specifying `--output_gene_summary`. A 4-column BED file is required for this, where the fourth column contains the gene label.

### 3. Small variant calling with Clair3

The workflow implements a deconstructed version of [Clair3](https://github.com/HKU-BAL/Clair3) (v1.0.4) to call germline variants. The appropriate model can be provided with the `--basecaller_cfg` option. To decide on the appropriate model you can check out the Dorado documentation for a list of available basecalling models.
This workflow takes advantage of the parallel nature of Nextflow, providing optimal efficiency in high-performance, distributed systems. The workflow will automatically call small variants (SNPs and indels), collect statistics, annotate them with [SnpEff](https://pcingola.github.io/SnpEff/) (and additionally for SNPs, ClinVar details), and create a report summarising the findings.

If desired, the workflow can perform phasing of structural variants by using the `--phased` option. This will lead the workflow to use [whatshap](https://whatshap.readthedocs.io/) to perform phasing of the variants, with the option to use [longphase](https://github.com/twolinin/longphase) instead by setting `--use_longphase true`. The phasing will also generate a GFF file with the annotation of the phase blocks, enabling the visualisation of these blocks in genome browsers.

### 4. Structural variant (SV) calling with Sniffles2

The workflow allows for calling of SVs using long-read sequencing data with [Sniffles2](https://github.com/fritzsedlazeck/Sniffles).
The workflow will perform SV calling, filtering and generation of a report.
Optionally, the workflow can also evaluate calls on HG002 against a truth set (provided the input data was aligned to HG19).
The SV workflow takes an optional `--tr_bed` option to specify tandem repeats in the reference sequence --- see the [sniffles](https://github.com/fritzsedlazeck/Sniffles) documentation for more information.
SVs can be phased using `--phased`. However, this will cause the workflow to run SNP analysis, as SV phasing relies on the haplotagged reads generated in this stage.

### 5. Modified base calling with modkit

Modified base calling can be performed by specifying `--mod`. The workflow will call modified bases using [modkit](https://github.com/nanoporetech/modkit). 
The workflow will automatically check whether the files contain the appropriate `MM`/`ML` tags, required for running [modkit pileup](https://nanoporetech.github.io/modkit/intro_bedmethyl.html). If the tags are not found, the workflow will not run the individual analysis, but will still run the other subworkflows requested by the user.
The default behaviour of the workflow is to run modkit with the `--cpg --combine-strands` options set. It is possible to report strand-aware modifications by providing `--force_strand`, which will trigger modkit to run in default mode. The resulting bedMethyl will include modifications for each site on each strand separately.
The modkit run can be fully customized by providing `--modkit_args`. This will override any preset, and allow full control over the run of modkit.
Haplotype-resolved aggregated counts of modified bases can be obtained with the `--phased` option. This will generate three distinct BEDMethyl files with the naming pattern `{{ alias }}.wf_mods.{{ haplotype }}.bedmethyl.gz`, where `haplotype` can be `1`, `2` or `ungrouped`.

### 6a. Copy number variants (CNV) calling with Spectre

CNV calling is performed using a fork of [Spectre](https://github.com/fritzsedlazeck/Spectre/tree/ont-dev), using the `--cnv` flag. Spectre is the default CNV caller in the workflow, and is compatible with hg38/GRCh38. The output of this workflow is a VCF of CNV calls, annotated with SnpEff.

### 6b. Copy number variants (CNV) calling with QDNASeq

CNV calling may alternatively be performed using [QDNAseq](https://github.com/ccagc/QDNAseq), using `--cnv --use_qdnaseq`. This workflow is compatible with genome builds hg19/GRCh37 or hg38/GRCh38, and is recommended for shallow WGS or adaptive sampling data. In addition to the VCF of CNV calls, the workflow emits QDNAseq-generated plots and BED files of both raw read counts per bin and corrected, normalised, and smoothed read counts per bin. Please note that QDNAseq was the default CNV caller until version 1.11.0 of the workflow, and the additional `--use_qdnaseq` flag is now required to use it.

### 7. Short tandem repeat (STR) genotyping with Straglr

STR genotyping is performed using a fork of [straglr](https://github.com/philres/straglr). This workflow is compatible with genome build hg38/GRCh38.
The number of calls for repeats on chrX is dependent on the sample's genetic sex which should be provided if known with `--sex XX` or `--sex XY`.
If `--sex` is not specified, the workflow will attempt to infer the genetic sex from coverage of the allosomes, falling back to `XX` if a determination is unclear.
Please be aware that incorrect sex assignment will result in the wrong number of calls for all repeats on chrX.
In addition to a gzipped VCF file containing STRs found in the dataset, the workflow emits a TSV straglr output containing reads spanning STRs, and a haplotagged BAM.

### 8. Phasing variants
Variant phasing is switched on simply using the `--phased` option.
By default, the workflow uses [whatshap](https://whatshap.readthedocs.io/) to perform phasing of the variants, with the option to use [longphase](https://github.com/twolinin/longphase) instead by setting `--use_longphase true`.
The workflow will automatically turn on the necessary phasing processes based on the selected subworkflows.
The behaviour of the phasing is summarised in the below table:

|         |        |         |            | Phased SNP VCF | Phased SV VCF | Phased bedMethyl |
|---------|--------|---------|------------|----------------|---------------|------------------|
| `--snp` | `--sv` | `--mod` | `--phased` |     &check;    |     &check;   |       &check;    |
| `--snp` | `--sv` |         | `--phased` |     &check;    |     &check;   |                  |
| `--snp` |        |         | `--phased` |     &check;    |               |                  |
|         | `--sv` |         | `--phased` |                |     &check;   |                  |
|         |        | `--mod` | `--phased` |                |               |       &check;    |

Using `--GVCF` together with `--phased` will generate a phased GVCF, created by reflecting the phased genotype and the phase set annotation in the VCF file. This operation is performed using `bcftools annotate`, targeting the `GT` and `PS` fields.

Running the phasing is a compute intensive process. Running the workflow in phasing mode doubles the runtime, and significantly increases the storage requirements to the order of terabytes.

### 9. Variant annotation
Annotation will be performed automatically by the SNP and SV subworkflows, and can be disabled by the user with `--annotation false`. The workflow will annotate the variants using [SnpEff](https://pcingola.github.io/SnpEff/), and currently only support the human hg19 and hg38 genomes. Additionally, the workflow will add the [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) annotations for the SNP variants.

Running the workflow on non-human samples will require this option to be disabled. For more detail, see Section 10 below.

### 10. Genome compatibility and running the workflow on non-human genomes
Some of the sub-workflows in wf-human-variation are restricted to certain genome builds, which means they will not be executable on non-human genomes or human genome builds outside hg19/GRCh37 and hg38/GRCh38. The following table summarises which subworkflows and options are available (or required) for a desired input genome:

|    Genome    | `--snp`  | `--sv`  | `--mod` | `--cnv` | `--cnv --use_qdnaseq` | `--str` | `--annotation false` | `--include_all_ctgs` |
|--------------|----------|---------|---------|---------|-----------------------|---------|----------------------|----------------------|
| hg19/GRCh37  | &check;  | &check; | &check; |         |        &check;        |         |         \*           |                      |
| hg38/GRCh38  | &check;  | &check; | &check; | &check; |        &check;        | &check; |         \*           |                      |
| Other human  | &check;  | &check; | &check; |         |                       |         |       &check;        |                      |
| Non human    | &check;  | &check; | &check; |         |                       |         |       &check;        |       &check;        |

\* As noted above, annotation is performed by default but can be switched off for hg19/GRCh37 and hg38/GRCh38.




## Troubleshooting

+ Annotations for `--snp` and `--sv` are generated using [SnpEff](https://pcingola.github.io/SnpEff/). For `--snp`, additional [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) annotations are displayed in the report where available (please note, the report will not display any variants classified as 'Benign' or 'Likely benign', however these variants will be present in the
output VCF).
+ Specifying a suitable [tandem repeat BED for your reference](https://raw.githubusercontent.com/fritzsedlazeck/Sniffles/master/annotations/) with `--tr_bed` can improve the accuracy of SV calling.
+ Aggregation of modified calls with `--mod` requires data to be basecalled with a model that includes base modifications, providing the `MM` and `ML` BAM tags
+ CRAM files generated within the workflow cannot be read without the corresponding reference
+ The STR workflow performs genotyping of specific repeats, which can be found [here](https://github.com/epi2me-labs/wf-human-variation/blob/master/data/wf_str_repeats.bed).
+ By default, SNPs and SVs will be called on chromosomes 1-22, X, Y and MT; to call sites on other sequences enable the `--include_all_ctgs` option.
+ While designed to work on human genomes, the workflow can be run on non-human species by setting `--cnv false --str false --annotation false --include_all_ctgs true`.
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



