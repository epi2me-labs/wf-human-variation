# Human variation workflow

Human SNV, SV, CNV, genotyping STR and modified base calling.



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

The tools embedded in individual sub-workflows within wf-human-variation are specifically designed for use with whole-genome Oxford Nanopore Technologies sequencing data. While 20x average coverage is the absolute minimum requirement for the workflow to run, we recommend an average coverage above 30x to ensure optimal performance. Usage below the minimum coverage may cause the workflow to terminate with an error, or yield unexpected outcomes.




## Compute requirements

Recommended requirements:

+ CPUs = 32
+ Memory = 128GB

Minimum requirements:

+ CPUs = 16
+ Memory = 32GB

Approximate run time: Variable depending on whether it is targeted sequencing or whole genome sequencing, as well as coverage and the individual analyses requested. For instance, a 90X human sample run (options: `--snp --sv --mod --str --cnv --phased --sex XY`) takes less than 8h with recommended resources.

ARM processor support: False




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://docs.docker.com/get-started/)
or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-human-variation --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull epi2me-labs/wf-human-variation
```

A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-human-variation/hg002%2bmods.v1/wf-human-variation-demo.tar.gz
tar -xzvf wf-human-variation-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-human-variation \
	--bam 'wf-human-variation-demo/demo.bam' \
	--ref 'wf-human-variation-demo/demo.fasta' \
	--bed 'wf-human-variation-demo/demo.bed' \
	--sample_name 'DEMO' \
	--snp \
	--sv \
	--mod \
	--phased \
	-profile standard
```

For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
The `--bam` input parameter for this workflow accepts a path to a single BAM file, or a folder containing multiple BAM files for the sample. A sample name can be supplied with `--sample_name`.

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
| cnv | boolean | Call for copy number variants. | If this option is selected, copy number variant calling will be carried out with either Spectre (default) or QDNAseq. To use QDNAseq instead of Spectre, use the option --use_qdnaseq. | False |
| str | boolean | Enable Straglr to genotype STR expansions. | If this option is selected, genotyping of STR expansions will be carried out using Straglr. This option will also automatically enable the SNP subworkflow as the STR caller requires a haplotagged BAM. This sub-workflow is only compatible with genome build hg38. | False |
| mod | boolean | Enable output of modified calls to a bedMethyl file [requires input BAM with Ml and Mm tags] | This option is automatically selected and aggregation of modified calls with be carried out using modkit if Ml and Mm tags are found. Disable this option to prevent output of a bedMethyl file. | False |


### Main options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_name | string | Sample name to be displayed in workflow outputs. | When providing a MinKNOW experiment folder, the sample name must match one of the sample sub-folders in the input folder to specify which sample to analyse. |  |
| bam | string | BAM or unaligned BAM (uBAM) files for the sample to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a folder containing BAM files; (iii) the path to a MinKNOW experiment folder containing sub-folders for each sample in the experiment. For cases (i) and (ii) an optional sample name can be supplied with `--sample_name`. For case (iii), a sample name must be provided to select the sample for analysis from the experiment folder. |  |
| ref | string | Path to a reference FASTA file. | Reference against which to compare reads for variant calling. |  |
| bam_min_coverage | number | Minimum read coverage required to run analysis. |  | 20 |
| bed | string | An optional BED file enumerating regions to process for variant calling. If a fourth column is present, this will be used to generate coverage metrics for the provided regions. |  |  |
| coverage_bed | string | An optional BED file enumerating regions to process for coverage metrics only. Requires a fourth BED column specifying region name in order to be valid. |  |  |
| annotation | boolean | SnpEff annotation. | If this option is unselected, VCFs will not be annotated with SnpEff. | True |
| phased | boolean | Perform phasing. | This option enables phasing of SV, SNP and modifications, depending on which sub-workflow has been chosen; see [README](README.md#9-phasing-variants) for more details. | False |
| include_all_ctgs | boolean | Call for variants on all sequences in the reference, otherwise small and structural variants will only be called on chr{1..22,X,Y,MT}. | Enabling this option will call for variants on all contigs of the input reference sequence. Typically this option is not required as standard human reference sequences contain decoy and unplaced contigs that are usually omitted for the purpose of variant calling. This option might be useful for non-standard reference sequence databases. | False |
| igv | boolean | Visualize outputs in the EPI2ME IGV visualizer. | Enabling this option will visualize the alignments and VCF files in the EPI2ME desktop app IGV visualizer. | False |
| out_dir | string | Directory for output of all workflow results. |  | output |
| partner | string | Prepare outputs for tertiary analyses with partners (geneyx or fabric). |  |  |


### Copy number variant calling options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| use_qdnaseq | boolean | Use QDNAseq for CNV calling. | Set this to true to use QDNAseq for CNV calling instead of Spectre. QDNAseq is better suited to shorter reads such as those generated from adaptive sampling experiments. | False |
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
| output_xam_fmt | string | Desired file format of alignment files created by alignment and phasing. | This setting controls the file format of (1) alignment files created by aligning or re-aligning an input BAM and (2) alignment files with haplotag information created during phasing of an input BAM. If using QDNAseq for CNV calling, the setting will be ignored for alignment or realignment as QDNAseq requires BAM input. | cram |
| modkit_args | string | The additional options for modkit. | This is an advanced option to allow running modkit with custom settings. The arguments specified in this option will fully override all options set by the workflow. To provide custom arguments to [modkit pileup](https://nanoporetech.github.io/modkit/advanced_usage.html#pileup) from command line proceed as follows: `--modkit_args="--preset traditional"` |  |
| sniffles_args | string | Additional command line arguments to pass to the Sniffles2 process | The additional command line arguments will be passed directly to [Sniffles2](https://github.com/fritzsedlazeck/Sniffles/tree/v2.6.3); ensure to use the right commands for the version and from command line provide them as follows: `--sniffles_args="--mosaic"`. |  |
| spectre_args | string | Additional command line arguments to pass to the Spectre process | The additional command line arguments will be passed directly to [Spectre](https://github.com/epi2me-labs/ont-spectre); from the command line provide them as follows: `--spectre_args="--min-cnv-len 80000"`. |  |
| override_basecaller_cfg | string | Name of the model to use for selecting a small variant calling model. | The workflow will attempt to find the basecaller model from the headers of your input data, providing a value for this option will override the model found in the data. If the model cannot be found in the header, it must be provided with this option as the basecaller model is required for small variant calling. The basecaller model is used to automatically select the appropriate small variant calling model. The model list shows all models that are compatible for small variant calling with this workflow. You should select 'custom' to override the basecaller_cfg with clair3_model_path. |  |
| tr_bed | string | Input BED file containing tandem repeat annotations for the reference genome. | Providing a tandem repeat BED can improve calling in repetitive regions. The workflow will attempt to select an appropriate hg19 or hg38 TR BED depending on the genome detected, but this can be overridden with a custom TR BED with this parameter. |  |
| alignment_report_coverage_threshold | number | If a BED file is provided to the workflow, the "Regions below target coverage" section in the alignment report tabulates BED regions with mean coverage below this threshold. | This parameter is used to identify target regions where the mean coverage falls below a specified threshold | 20 |


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
| Structural variant SNF | {{ alias }}.wf_sv.snf | SNF file with the SVs for the sample, for onward multi-sample SV calling. | per-sample |
| Copy number variants VCF | {{ alias }}.wf_cnv.vcf.gz | VCF file with the CNV for the sample. | per-sample |
| ClinVar variant VCF | {{ alias }}.wf_snp_clinvar.vcf.gz | VCF file with ClinVar annotations. | per-sample |
| Modified bases BEDMethyl | {{ alias }}.wf_mods.bedmethyl.gz | BED file with the aggregated modification counts for the sample. | per-sample |
| Modified bases BEDMethyl (haplotype 1) | {{ alias }}.wf_mods.1.bedmethyl.gz | BED file with the aggregated modification counts for haplotype 1 of the sample. | per-sample |
| Modified bases BEDMethyl (haplotype 2) | {{ alias }}.wf_mods.2.bedmethyl.gz | BED file with the aggregated modification counts for haplotype 2 of the sample. | per-sample |
| Modified bases BEDMethyl (ungrouped) | {{ alias }}.wf_mods.ungrouped.bedmethyl.gz | BED file with the aggregated modification counts of non-haplotagged reads for the sample. | per-sample |
| Short tandem repeat VCF | {{ alias }}.wf_str.vcf.gz | VCF file with the STR sites for the sample. | per-sample |
| Alignment file | {{ alias }}.cram | CRAM or BAM file with the aligned reads for the sample, generated when the input file is unaligned. | per-sample |
| Alignment file index | {{ alias }}.cram.crai | The index of the resulting CRAM or BAM file with the reads for the sample, generated when the input file is unaligned. | per-sample |
| Haplotagged alignment file | {{ alias }}.haplotagged.cram | CRAM or BAM file of all input reads with haplotags added by phasing. | per-sample |
| Haplotagged alignment file index | {{ alias }}.haplotagged.cram.crai | The index of the resulting CRAM or BAM file produced when haplotags have been added by phasing. | per-sample |
| Mean coverage for each region | {{ alias }}.regions.bed.gz | The mean coverage in the individual regions of the genome in BED format. | per-sample |
| Coverage per region above the given thresholds | {{ alias }}.thresholds.bed.gz | The BED reporting the number of bases in each region that are covered at or above each threshold values (1x, 10x, 20x and 30x). | per-sample |
| Distribution of the proportion of total bases covered by a given coverage value | {{ alias }}.mosdepth.global.dist.txt | The cumulative distribution indicating the proportion of total bases covered by a given coverage value, both genome-wide and by sequence. | per-sample |
| Mean coverage per sequence and target region | {{ alias }}.mosdepth.summary.txt | The summary of mean depths per chromosome and within specified regions per chromosome. | per-sample |
| BEDgraph of the single-base coverage | {{ alias }}.per-base.bedgraph.gz | The single-base coverage of the genome in BED graph format. | per-sample |
| BED summary | {{ alias }}.bed_summary.tsv | Tab separated table with the average coverage and percentage of positions that are covered to a given threshold, for each region specified in the `--bed` BED file. Only generated if the BED file has at least 4 columns, where the fourth column is the region label. | per-sample |
| Coverage BED summary | {{ alias }}.coverage_bed_summary.tsv | Tab separated table with the average coverage and percentage of positions that are covered to a given threshold, for each region specified in the `--coverage_bed` BED file. Only generated if the BED file has at least 4 columns, where the fourth column is the region label. | per-sample |
| Haplocheck contamination summary | {{ alias }}.haplocheck.tsv | A table generated by [haplocheck](https://mitoverse.readthedocs.io/haplocheck/haplocheck/), with estimate of contamination from the MT genome. | per-sample |
| FAI index of the reference FASTA file | {{ ref }}.fai | FAI Index of the reference FASTA file. | aggregated |
| GZI index of the reference FASTA file | {{ ref }}.gzi | GZI Index of the reference FASTA file. | aggregated |




## Pipeline overview

The workflow is composed of 6 distinct subworkflows, each enabled by a command line option:

* [SNP calling](#3-small-variant-calling-with-clair3): `--snp`
* [SV calling](#4-structural-variant-sv-calling-with-sniffles2): `--sv`
* [Analysis of modified bases](#5-modified-base-calling-with-modkit): `--mod`
* [CNV calling (Spectre)](#6a-copy-number-variants-cnv-calling-with-spectre): `--cnv`
* [CNV calling (QDNAseq)](#6b-copy-number-variants-cnv-calling-with-qdnaseq): `--cnv --use_qdnaseq`
* [STR genotyping](#7-short-tandem-repeat-str-genotyping-with-straglr): `--str`

Subworkflows where the relevant option is omitted will not be run.

### 1. Input and data preparation

The workflow relies on two primary input files:
1. A reference genome in [FASTA format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
2. Sequencing data for the sample in the form of a single [BAM file](https://samtools.github.io/hts-specs/SAMv1.pdf) or a folder of BAM files, either aligned or unaligned.

When analysing human data, we recommend using [human_g1k_v37.fasta.gz](https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/ref/human_g1k_v37.fasta.gz) or [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz](https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz). For more information see [this blog post](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) which outlines potential pitfalls with the various flavours of human references.

The input BAM file can be generated using the [wf-basecalling](https://github.com/epi2me-labs/wf-basecalling/) workflow, which is up to date with the current dorado releases and models.

### 2. Data QC and pre-processing
The workflow starts by performing multiple checks of the input BAM file, as well as computing:
1. depth of sequencing with [mosdepth](https://github.com/brentp/mosdepth);
2. read alignment statistics with [fastcat](https://github.com/epi2me-labs/fastcat).

After computing the coverage, the workflow will check that the input BAM file has a depth greater than `--bam_min_coverage`.
In case the user specify `--bam_min_coverage 0`, the check will be skipped and the workflow will proceed directly to the downstream analyses.
Some components work better withing certain ranges of coverage, and the user might achieve better results by providing a target coverage to downsample to. The user can set `--downsample_coverage true` to enable the downsampling of the reads, and `--downsample_coverage_target {{ X }}` to specify the target coverage (default: 60x).

Two optional BED files can be provided to customise workflow behaviour. The `--bed` file may be used to restrict variant calling in the `--snp` and `--sv` subworkflows. If the file provided to `--bed` has 4 or more columns, a summary of coverage describing the percentage of each region covered at various read depths will automatically be made available.

Secondly, a `--coverage_bed` file may be used to specify additional genomic regions for which to similarly generate a summary of read coverage in the alignment report, without affecting analysis. The `--coverage_bed` file is used exclusively for generating coverage summaries. The coverage BED must contain at least four columns (as the name field is required) otherwise the workflow will terminate with an error.

### 3. Small variant calling with Clair3

The workflow implements a deconstructed version of [Clair3](https://github.com/HKU-BAL/Clair3) (v1.0.4) to call germline variants.
The workflow will select an appropriate Clair3 model by detecting the basecall model from the input data.
If the input data does not have the required information to determine the basecall model, the workflow will require the basecall model to be provided explicitly with the `--override_basecaller_cfg` option.
This workflow takes advantage of the parallel nature of Nextflow, providing optimal efficiency in high-performance, distributed systems. The workflow will automatically call small variants (SNPs and indels), collect statistics, annotate them with [SnpEff](https://pcingola.github.io/SnpEff/) (and additionally for SNPs, ClinVar details), and create a report summarising the findings.

If desired, the workflow can perform phasing of structural variants by using the `--phased` option. The workflow will use [whatshap](https://whatshap.readthedocs.io/) to perform phasing of the variants. Phasing will also generate a GFF file with the annotation of the phase blocks, enabling the visualisation of these blocks in genome browsers.

### 4. Structural variant (SV) calling with Sniffles2

The workflow allows for calling of SVs using long-read sequencing data with [Sniffles2](https://github.com/fritzsedlazeck/Sniffles).
The workflow will perform SV calling, filtering and generation of a report.
The SV workflow can accept a tandem repeat annotations BED file to improve calling in repetitive regions --- see the [sniffles](https://github.com/fritzsedlazeck/Sniffles) documentation for more information. The workflow will attempt to select an appropriate hg19 or hg38 TR BED but you can override the BED file used with the `--tr_bed` parameter.
SVs can be phased using `--phased`. However, this will cause the workflow to run SNP analysis, as SV phasing relies on the haplotagged reads generated in this stage.

### 5. Modified base calling with modkit

Modified base calling can be performed by specifying `--mod`. The workflow will call modified bases using [modkit](https://github.com/nanoporetech/modkit). 
The workflow will automatically check whether the files contain the appropriate `MM`/`ML` tags, required for running [modkit pileup](https://nanoporetech.github.io/modkit/intro_pileup.html). If the tags are not found, the workflow will not run the individual analysis, but will still run the other subworkflows requested by the user.
The default behaviour of the workflow is to run modkit with the `--cpg --combine-strands` options set. It is possible to report strand-aware modifications by providing `--force_strand`, which will trigger modkit to run in default mode. The resulting bedMethyl will include modifications for each site on each strand separately.
The modkit run can be fully customized by providing `--modkit_args`. This will override any preset, and allow full control over the run of modkit.
Haplotype-resolved aggregated counts of modified bases can be obtained with the `--phased` option. This will generate three distinct BEDMethyl files with the naming pattern `{{ alias }}.wf_mods.{{ haplotype }}.bedmethyl.gz`, where `haplotype` can be `1`, `2` or `ungrouped`.

### 6a. Copy number variants (CNV) calling with Spectre

CNV calling is performed using an ONT implementation of [Spectre](https://github.com/epi2me-labs/ont-spectre/) using the `--cnv` flag. Spectre is the default CNV caller in the workflow and is compatible with genome builds hg19/GRCh37 or hg38/GRCh38. The output of this workflow is a VCF of CNV calls annotated with SnpEff. Advanced users may wish to override default parameters using `--spectre_args`, e.g. `--spectre_args "--min-cnv-len 100000"`. We don't recommend setting `--min-cnv-len` to less than N50 * 5 to avoid false positive calls.

### 6b. Copy number variants (CNV) calling with QDNASeq

CNV calling may alternatively be performed using [QDNAseq](https://github.com/ccagc/QDNAseq), using `--cnv --use_qdnaseq`. This workflow is compatible with genome builds hg19/GRCh37 or hg38/GRCh38, and is recommended for shallow WGS or adaptive sampling data. In addition to the VCF of CNV calls, the workflow emits QDNAseq-generated plots and BED files of both raw read counts per bin and corrected, normalised, and smoothed read counts per bin. Please note that QDNAseq was the default CNV caller until version 1.11.0 of the workflow, and the additional `--use_qdnaseq` flag is now required to use it.

### 7. Short tandem repeat (STR) genotyping with Straglr

STR genotyping is performed using a fork of [straglr](https://github.com/philres/straglr). This workflow is compatible with genome build hg38/GRCh38.
The number of calls for repeats on chrX is dependent on the sample's genetic sex which should be provided if known with `--sex XX` or `--sex XY`.
If `--sex` is not specified, the workflow will attempt to infer the genetic sex from coverage of the allosomes, falling back to `XX` if a determination is unclear.
Please be aware that incorrect sex assignment will result in the wrong number of calls for all repeats on chrX.
The STR subworkflow requires a haplotagged BAM file to accurately determine repeat expansions per haplotype.
To generate this haplotagged BAM, the workflow automatically enables the `--snp` subworkflow to produce a phased VCF, which is then used to assign haplotype tags to reads.
It is therefore not possible to genotype STRs without running the SNP subworkflow.
In addition to a gzipped VCF file containing STRs found in the dataset, the workflow emits a TSV straglr output containing reads spanning STRs, and a haplotagged BAM.

### 8. Phasing variants
Variant phasing is switched on simply using the `--phased` option.
The workflow uses [whatshap](https://whatshap.readthedocs.io/) to perform phasing of the variants. The workflow will also add phase information to variants based on haplotagged reads using [whatshap haplotagphase](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotagphase).
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
| hg19/GRCh37  | &check;  | &check; | &check; | &check; |        &check;        |         |         \*           |                      |
| hg38/GRCh38  | &check;  | &check; | &check; | &check; |        &check;        | &check; |         \*           |                      |
| Other human  | &check;  | &check; | &check; |         |                       |         |       &check;        |                      |
| Non human    | &check;  | &check; | &check; |         |                       |         |       &check;        |       &check;        |

\* As noted above, annotation is performed by default but can be switched off for hg19/GRCh37 and hg38/GRCh38.

> Please note that while running the workflow is possible on non-human genomes by following the guidance above, this is not a supported use-case of wf-human-variation. Even when following the suggestions in this section, the workflow may terminate with an error or yield unexpected outcomes on non-human inputs.



## Troubleshooting

+ Annotations for `--snp` and `--sv` are generated using [SnpEff](https://pcingola.github.io/SnpEff/). For `--snp`, additional [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) annotations are displayed in the report where available (please note, the report will not display any variants classified as 'Benign' or 'Likely benign', however these variants will be present in the
output VCF).
+ Aggregation of modified calls with `--mod` requires data to be basecalled with a model that includes base modifications, providing the `MM` and `ML` BAM tags
+ CRAM files generated within the workflow cannot be read without the corresponding reference
+ The STR workflow performs genotyping of specific repeats, which can be found [here](https://github.com/epi2me-labs/wf-human-variation/blob/master/data/wf_str_repeats.bed).
+ By default, SNPs and SVs will be called on chromosomes 1-22, X, Y and MT; to call sites on other sequences enable the `--include_all_ctgs` option.
+ While designed to work on human genomes, the workflow can be run on non-human species by setting `--cnv false --str false --annotation false --include_all_ctgs true`.
+ Ensure that the provided reference and BED files use the same chromosome coding (for example, that they both have the `chr` prefix, or they both to not have it).
+ If unaligned reads were provided, the workflow will output a CRAM file (or BAM if the user runs the `--cnv` option) containing the alignments used to make the downstream variant calls.
+ Renaming, moving or deleting the input BAM, reference genome or the output directory from the location provided at runtime will cause IGV not to load.
+ The workflow expects either an uncompressed or [`bgzip`](https://www.htslib.org/doc/bgzip.html)-compressed reference. If the user provides a reference compressed not with `bgzip`, the workflow will run to completion, but won't be able to generate the necessary indexes to visualize the outputs in IGV.




## FAQ's

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-template/issues) page or start a discussion on the [community](https://community.nanoporetech.com/). 

+ *The number of SNVs and indels in the report do not sum up to the number of record, is that normal?* - Yes; this can be due to some multiallelic sites carrying a mixture of SNV and indel alleles.



## Related blog posts

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)
+ [blogpost](https://labs.epi2me.io/copy-number-calling/) and [CNV workflow documentation](https://github.com/epi2me-labs/wf-cnv) for more information on running the copy number calling subworkflow.
+ [blogpost](https://labs.epi2me.io/human-targeted-analysis/) on how to run a targeted BRCA Gene Analysis with `wf-human-variation`
+ [blogpost](https://labs.epi2me.io/giab-2023.05/) on the generation of four Genome in a Bottle samples processed with `wf-human-variation`

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



