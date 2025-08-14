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