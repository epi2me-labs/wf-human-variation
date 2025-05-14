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
| bam | string | BAM or unaligned BAM (uBAM) files for the sample to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a folder containing BAM files; (iii) the path to a MinKNOW experiment folder containing sub-folders for each sample in the experiment. For cases (i) and (ii) an optional sample name can be supplied with `--sample`. For case (iii), a sample name must be provided to select the sample for analysis from the experiment folder. |  |
| ref | string | Path to a reference FASTA file. | Reference against which to compare reads for variant calling. |  |
| bam_min_coverage | number | Minimum read coverage required to run analysis. |  | 20 |
| bed | string | An optional BED file enumerating regions to process for variant calling. If a fourth column is present, this will be used to generate coverage metrics for the provided regions. |  |  |
| coverage_bed | string | An optional BED file enumerating regions to process for coverage metrics only. Requires a fourth BED column specifying region name in order to be valid. |  |  |
| annotation | boolean | SnpEff annotation. | If this option is unselected, VCFs will not be annotated with SnpEff. | True |
| phased | boolean | Perform phasing. | This option enables phasing of SV, SNP and modifications, depending on which sub-workflow has been chosen; see [README](README.md#9-phasing-variants) for more details. | False |
| include_all_ctgs | boolean | Call for variants on all sequences in the reference, otherwise small and structural variants will only be called on chr{1..22,X,Y,MT}. | Enabling this option will call for variants on all contigs of the input reference sequence. Typically this option is not required as standard human reference sequences contain decoy and unplaced contigs that are usually omitted for the purpose of variant calling. This option might be useful for non-standard reference sequence databases. | False |
| igv | boolean | Visualize outputs in the EPI2ME IGV visualizer. | Enabling this option will visualize the alignments and VCF files in the EPI2ME desktop app IGV visualizer. | False |
| out_dir | string | Directory for output of all workflow results. |  | output |


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
| sniffles_args | string | Additional command line arguments to pass to the Sniffles2 process | The additional command line arguments will be passed directly to [Sniffles2](https://github.com/fritzsedlazeck/Sniffles/tree/v2.0.7); ensure to use the right commands for the version and from command line provide them as follows: `--sniffles_args="--non-germline"`. |  |
| spectre_args | string | Additional command line arguments to pass to the Spectre process | The additional command line arguments will be passed directly to [Spectre](https://github.com/epi2me-labs/ont-spectre); from the command line provide them as follows: `--spectre_args="--min-cnv-len 80000"`. |  |
| override_basecaller_cfg | string | Name of the model to use for selecting a small variant calling model. | The workflow will attempt to find the basecaller model from the headers of your input data, providing a value for this option will override the model found in the data. If the model cannot be found in the header, it must be provided with this option as the basecaller model is required for small variant calling. The basecaller model is used to automatically select the appropriate small variant calling model. The model list shows all models that are compatible for small variant calling with this workflow. You should select 'custom' to override the basecaller_cfg with clair3_model_path. |  |
| tr_bed | string | Input BED file containing tandem repeat annotations for the reference genome. | Providing a tandem repeat BED can improve calling in repetitive regions. The workflow will attempt to select an appropriate hg19 or hg38 TR BED depending on the genome detected, but this can be overridden with a custom TR BED with this parameter. |  |


### Multiprocessing Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Set max number of threads to use for more intense processes (limited by config executor cpus) |  | 4 |
| ubam_map_threads | integer | Set max number of threads to use for aligning reads from uBAM (limited by config executor cpus) |  | 8 |
| ubam_sort_threads | integer | Set max number of threads to use for sorting and indexing aligned reads from uBAM (limited by config executor cpus) |  | 3 |
| ubam_bam2fq_threads | integer | Set max number of threads to use for uncompressing uBAM and generating FASTQ for alignment (limited by config executor cpus) |  | 1 |
| modkit_threads | integer | Total number of threads to use in modkit modified base calling (limited by config executor cpus) |  | 4 |


