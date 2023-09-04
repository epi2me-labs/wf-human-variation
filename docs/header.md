# Human variation workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow
for analysing variation in human genomic data. Specifically this workflow can
perform the following:

* basecalling of FAST5 (or POD5) sequencing data
* diploid variant calling
* structural variant calling
* analysis of modified base calls
* copy number variant calling
* short tandem repeat (STR) expansion genotyping

The wf-human-variation workflow consolidates the small variant calling from the
previous wf-human-snp, structural variant calling from wf-human-sv, CNV calling from wf-cnv (all
of which are now deprecated), as well as performing STR expansion genotyping. This pipeline performs the steps of the four
pipelines simultaneously and the results are generated and output in the same
way as they would have been had the pipelines been run separately.
