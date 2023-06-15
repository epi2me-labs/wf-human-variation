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
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v3.5.2'  \
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
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_hac@v3.5.2'  \
    --remora_cfg 'dna_r10.4.1_e8.2_400bps_hac@v3.5.2_5mCG@v2' \
    --bed path/to.bed \
    --ref path/to.fasta \
    --out_dir ${OUTPUT}
```

**Workflow outputs**

The primary outputs of the workflow include:

* a gzipped VCF file containing SNPs found in the dataset (`--snp`)
* a gzipped VCF file containing the SVs called from the dataset (`--sv`)
* a gzipped [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/) file aggregating modified CpG base counts (`--methyl`)
* a VCF of CNV calls, QDNAseq-generated plots, and BED files of both raw read counts per bin and corrected, normalised, and smoothed read counts per bin (`--cnv`)
* a gzipped VCF file containing STRs found in the dataset, TSV Straglr output containing reads spanning STRs, and a haplotagged BAM (`--str`)
* an HTML report detailing the primary findings of the workflow, for SNP, SV, CNV calling, and STR genotyping.
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
- Specifying a suitable [tandem repeat BED for your reference](https://raw.githubusercontent.com/fritzsedlazeck/Sniffles/master/annotations/) with `--tr_bed` can improve the accuracy of SV calling.
- Aggregation of methylation calls with `--methyl` requires data to be basecalled with a model that includes base modifications, providing the `MM` and `ML` BAM tags
- Refer to the [Dorado documentation](https://github.com/nanoporetech/dorado#available-basecalling-models) for a list of available basecalling models
- Take care to retain the input reference when basecalling or alignment has been performed as CRAM files cannot be read without the corresponding reference!
- Refer to our [blogpost](https://labs.epi2me.io/copy-number-calling/) and [CNV workflow documentation](https://github.com/epi2me-labs/wf-cnv) for more information on running the copy number calling subworkflow.
- The STR workflow performs genotyping of specific repeats, which can be found [here](https://github.com/epi2me-labs/wf-human-variation/blob/master/data/wf_str_repeats.bed), and is only compatible with genome build 38.

### Support for basecalling on GPU

Currently, neither Windows nor macOS should be considered as supported operating systems for wf-basecalling as neither can efficiently perform basecalling through Docker. Users are encouraged to download Dorado for Windows and MacOS directly from the [Dorado GitHub repository](https://github.com/nanoporetech/dorado#installation).

We do not directly support configuration of accelerated computing through WSL2 and Docker on Windows, but it is possible to set up for most versions of Windows 11 and some versions of Windows 10 and we direct users to the [CUDA on WSL User Guide](https://docs.nvidia.com/cuda/wsl-user-guide/index.html). Take note of the support constraints section.

On MacOS, GPU support through Docker remains in technical infancy. In addition, the containers we provide will not be able to leverage the M1 and M2 architecture and will not run as performantly as if Dorado had been run natively.

Note that when using Docker for accelerated computing on Linux, you will need the `nvidia-container-toolkit` installed. If you observe the error "could not select device driver with capabilities gpu", you should follow the instructions to install `nvidia-container-toolkit` [here](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html#setting-up-nvidia-container-toolkit). You will need to follow the steps to:

- Setup the package repository and the GPG key (ignore the box about experimental releases)
- Update package listings
- Install nvidia-container-toolkit
- Configure the Docker daemon to recognize the NVIDIA Container Runtime
- Restart the Docker daemon to complete the installation after setting the default runtime

This section will be kept up to date with latest advice for running our workflows on the GPU.
