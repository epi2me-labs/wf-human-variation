# Human variation workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow
for analysing variation in human genomic data. Specifically this workflow performs:

* basecalling of FAST5 (or POD5) sequencing data
* diploid variant calling
* structural variant calling
* aggregation of modified base counts

The wf-human-variation workflow consolidates the small variant calling
from wf-human-snp with the structual variant calling from wf-human-sv. This pipeline
performs the steps of the two pipelines simultaneously and the results are generated
and output in the same way as they would have been had the pipelines been run separately.

## Introduction

This workflow uses [Clair3](https://www.github.com/HKU-BAL/Clair3) for calling small
variants from long reads. Clair3 makes the best of two methods: pileup (for fast
calling of variant candidates in high confidence regions), and full-alignment
(to improve precision of calls of more complex candidates).

This workflow uses [sniffles2](https://github.com/fritzsedlazeck/Sniffles) for
calling structural variants.

This workflow uses [modbam2bed](https://github.com/epi2me-labs/modbam2bed) for
aggregate modified base counts into a [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/) file.

This workflow uses [Dorado](https://github.com/nanoporetech/dorado/tree/master/dorado)
for basecalling `pod5` or `fast5` signal data.
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
    --model demo_data/ont_r104_e81_sup_g5015 \
    --sample_name MY_SAMPLE \
    --out_dir ${OUTPUT}
```

**Additional examples**

***SNP only***
```
OUTPUT=output
nextflow run epi2me-labs/wf-human-variation \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --snp \
    --bam path/to.bam \
    --bed path/to.bed \
    --ref path/to.fasta \
    --model path/to/clair3/model/dir \
    --out_dir ${OUTPUT}
```

***SV only***
```
OUTPUT=output
nextflow run epi2me-labs/wf-human-variation \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --sv \
    --bam path/to.bam \
    --bed path/to.bed \
    --ref path/to.fasta \
    --tr_bed path/to/tandem_repeat.bed
    --out_dir ${OUTPUT}
```

***Methyl***
```
OUTPUT=output
nextflow run epi2me-labs/wf-human-variation \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --methyl \
    --bam path/to.bam \
    --ref path/to.fasta \
    --out_dir ${OUTPUT}
```

***Basecalling, SNP and SV***
```
OUTPUT=output
nextflow run epi2me-labs/wf-human-variation \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --snp --sv \
    --fast5_dir path/to/fast5/dir \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v3.5.2'  \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v3.5.2_5mCG@v2' \
    --basecaller_basemod_threads 2 \
    --bed path/to.bed \
    --ref path/to.fasta \
    --model path/to/clair3/model/dir \
    --out_dir ${OUTPUT}
```

**Workflow outputs**

The primary outputs of the workflow include:

* a gzipped VCF file containing SNPs found in the dataset (`--snp`)
* a gzipped VCF file containing the SVs called from the dataset (`--sv`)
* a gzipped [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/) file aggregating modified CpG base counts (`--methyl`)
* an HTML report detailing the primary findings of the workflow, for both the SNP and SV calling
* if unaligned reads were provided, a CRAM file containing the alignments used to make the calls

The secondary outputs of the workflow include:
* `{sample_name}.mapula.csv` and `{sample_name}.mapula.json` provide basic alignment metrics (primary and secondary counts, read N50, median accuracy)
* `mosdepth` outputs include:
    * `{sample_name}.mosdepth.global.dist.txt`: a cumulative distribution indicating the proportion of total bases for each and all reference sequences [more info](https://github.com/brentp/mosdepth#distribution-output)
    * `{sample_name}.regions.bed.gz`: reports the mean coverage for each region in the provided BED file
    * `{sample_name}.thresholds.bed.gz`: reports the number of bases in each region that are covered at or above each threshold value (1, 10, 20, 30X) [more info](https://github.com/brentp/mosdepth#thresholds)

**Workflow tips**

- Users familiar with `wf-human-snp` and `wf-human-sv` are recommended to familiarise themselves with any parameter changes by using `--help`, in particular:
    - All arms of the variation calling workflow use `--ref` (not `--reference`) and `--bed` (not `--target_bedfile`)
- Specifying a suitable [tandem repeat BED for your reference](https://raw.githubusercontent.com/fritzsedlazeck/Sniffles/master/annotations/) with `--tr_bed` can improve the accuracy of SV calling.
- Aggregation of methylation calls with `--methyl` requires data to be basecalled with a model that includes base modifications, providing the `MM` and `ML` BAM tags
- Refer to the [Dorado documentation](https://github.com/nanoporetech/dorado#available-basecalling-models) for a list of available basecalling models
## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
