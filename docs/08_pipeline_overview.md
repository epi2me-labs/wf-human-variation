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
