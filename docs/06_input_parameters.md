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


### Multiprocessing Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Set max number of threads to use for more intense processes (limited by config executor cpus) |  | 4 |
| ubam_map_threads | integer | Set max number of threads to use for aligning reads from uBAM (limited by config executor cpus) |  | 8 |
| ubam_sort_threads | integer | Set max number of threads to use for sorting and indexing aligned reads from uBAM (limited by config executor cpus) |  | 3 |
| ubam_bam2fq_threads | integer | Set max number of threads to use for uncompressing uBAM and generating FASTQ for alignment (limited by config executor cpus) |  | 1 |
| merge_threads | integer | Set max number of threads to use for merging alignment files (limited by config executor cpus) |  | 4 |
| modkit_threads | integer | Total number of threads to use in modkit modified base calling (limited by config executor cpus) |  | 4 |


