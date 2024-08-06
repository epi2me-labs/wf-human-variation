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

Please note, the tools embedded in individual sub-workflows within wf-human-variation are intended for use with 20x whole-genome Oxford Nanopore Technologies sequencing data (with the exception of QDNAseq, please see [this section](#6b-copy-number-variants-cnv-calling-with-qdnaseq) for more information). Usage outside of this (e.g. with adaptive sampling data, or using lower coverage inputs) may cause the workflow to terminate with an error, or produce unexpected results.
