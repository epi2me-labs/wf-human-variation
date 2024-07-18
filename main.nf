#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { snp; report_snp } from './workflows/wf-human-snp'
include { lookup_clair3_model } from './modules/local/wf-human-snp'

include { bam as sv } from './workflows/wf-human-sv'
include { output_sv } from './modules/local/wf-human-sv'

include { str } from './workflows/wf-human-str'
include { output_str } from './modules/local/wf-human-str'

include { cnv as cnv_spectre } from './workflows/wf-human-cnv'

include { cnv as cnv_qdnaseq } from './workflows/wf-human-cnv-qdnaseq'

include {
    index_ref_gzi;
    index_ref_fai;
    cram_cache;
    decompress_ref;
    mosdepth as mosdepth_input;
    mosdepth as mosdepth_downsampled;
    readStats;
    getAllChromosomesBed;
    publish_artifact;
    get_region_coverage;
    failedQCReport; 
    makeAlignmentReport; 
    getParams; 
    getVersions;
    getGenome; 
    eval_downsampling;
    downsampling;
    annotate_vcf as annotate_snp_vcf;
    concat_vcfs as concat_snp_vcfs;
    concat_vcfs as concat_refined_snp;
    sift_clinvar_vcf as sift_clinvar_snp_vcf;
    bed_filter;
    sanitise_bed;
    combine_metrics_json;
    output_cnv;
    infer_sex;
    haplocheck;
} from './modules/local/common'

include {
    ingress;
    cram_to_bam;
} from './lib/_ingress.nf'

include {
    prepare_reference;
} from './lib/reference.nf'

include {
    refine_with_sv;
    vcfStats;
    output_snp;
} from "./modules/local/wf-human-snp.nf"

include { 
    mod;
    validate_modbam;
    sample_probs;
} from './workflows/methyl'



// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Map colors = NfcoreTemplate.logColours(params.monochrome_logs)

    can_start = true
    // Check for deprecated options
    if (params.containsKey('methyl')) {
        log.error (colors.red + "The workflow now uses modkit instead of the deprecated modbam2bed. Please use --mod instead of --methyl to enable modkit." + colors.reset)
        can_start = false
    }
    if (params.containsKey('phase_methyl') || params.containsKey('phase_mod') || params.containsKey('phase_vcf')) {
        log.error (colors.red + "phase_methyl, phase_mod and phase_vcf are deprecated. Please use --phased instead to enable phasing of modkit results." + colors.reset)
        can_start = false
    }

    if (!params.snp && !params.sv && !params.mod && !params.cnv && !params.str) {
        log.error (colors.red + "No work to be done! Choose one or more workflows to run from [--snp, --sv, --cnv, --str, --mod]" + colors.reset)
        can_start = false
    }

    if (params.containsKey("ubam_threads")) {
        log.error (colors.red + "--ubam_threads is deprecated. Use `nextflow run ${workflow.manifest.name} --help` to see the parameter list." + colors.reset)
        can_start = false
    }

    if (params.containsKey("ubam")) {
        log.error (colors.red + "--ubam is deprecated as this workflow can determine whether (re)alignment is required automatically, use --bam instead." + colors.reset)
        can_start = false
    }

    if (params.containsKey("fast5_dir")) {
        log.error (colors.red + "--fast5_dir is deprecated as this workflow does not run basecalling anymore. Use wf-basecalling to generate a valid BAM file instead." + colors.reset)
        can_start = false
    }

    // Check if it is in genotyping mode
    if (params.snp && params.vcf_fn) {
        if (params.bed){
            throw new Exception(colors.red + "Clair3 cannot run with both --vcf_fn and --bed." + colors.reset)
        }
        log.warn ("Running Clair3 in genotyping mode with --vcf_fn will override --snp_min_af and --indel_min_af to 0.0.")
    }

    // check SV calling will be done when benchmarking SV calls
    if(params.sv_benchmark && !params.sv) {
        throw new Exception(colors.red + "Cannot benchmark SV subworkflow without running SV subworkflow! Enable the SV subworkflow with --sv." + colors.reset)
    }

    // If downsampling is required, check that the requested coverage is above the min threshold
    if(params.downsample_coverage) {
        if (params.downsample_coverage_target < params.bam_min_coverage){
            log.error (colors.red + "Downsampling target ${params.downsample_coverage_target} is lower than the minimum BAM coverage requested of ${params.bam_min_coverage}" + colors.reset)
            can_start = false
        }
    }

    // If gene summaries are requested, check a BED is provided and warn if input BED doesn't have 4 columns
    // Set gene_summary_bed accordingly so we can avoid running mosdepth on incompatible BED files
    def gene_summary_bed = false
    if(params.output_gene_summary) {
        if (!params.bed) {
            log.warn ("A BED file has not been provided, and therefore a gene summary will not be generated.")
        }
        else {
            col_size = file(params.bed).splitCsv(sep: '\t').first().size
            if (col_size < 4){
                log.warn ("The input BED file has fewer than 4 columns, and therefore a gene summary will not be generated.")
            }
            else {
                gene_summary_bed = true
            }
        }
    }

    // Programmatically define chromosome codes.
    // note that we avoid interpolation (eg. "${chr}N") to ensure that values
    // are Strings and not GStringImpl, ensuring that .contains works.
    ArrayList chromosome_codes = []
    ArrayList chromosomes = [1..22] + ["X", "Y", "M", "MT"]
    for (N in chromosomes.flatten()){
        chromosome_codes += ["chr" + N, "" + N]
    }

    // Trigger haplotagging
    def run_haplotagging = params.str || params.phased

    // Trigger CRAM to BAM conversion (for qdnaseq)
    // This will:
    // - cause downsampling to always be emitted as BAM
    // - OR if not downsampling, cause (re)alignment to always be emitted as BAM
    // - OR if not downsampling or (re)aligning, explicitly convert input CRAM to BAM
    def convert_cram_to_bam = params.cnv && params.use_qdnaseq

    // User desired alignment extentions
    def desired_xam_ext = params.output_xam_fmt == "cram" ? ["cram", "crai"] : ["bam", "bai"]

    // Determine what extentions should be output by ingress
    // Note that ingress does not handle downsampling so we carve that case out here
    if (convert_cram_to_bam && !params.downsample_coverage) {
        // Force BAM if not downsampling and BAM is needed downstream
        ingress_ext = ['bam', 'bai']
    }
    else {
        // No need to force a BAM - do what the user wants
        ingress_ext = desired_xam_ext
    }

    // Set extensions for the final haplotagged XAM
    // CNV is run on the ingressed BAM channel,
    //   and STR is run on the intermediate phased BAM,
    //   so we are free to output CRAM here, if desired.
    def haplotagged_output_fmt = desired_xam_ext

    // Notify users that QDNAseq usage will override the format of output XAM
    if (convert_cram_to_bam && params.output_xam_fmt == "cram") {
        log.warn "CNV calling subworkflow using QDNAseq does not support CRAM, but you have selected CRAM for your output file format."
        log.warn "You do not need to do anything, but any alignment or realignment will ignore your CRAM selection and be written as BAM to maintain compatibility with QDNAseq."
    }

    // Trigger the SNP workflow based on a range of different conditions:
    def run_snp = params.snp || run_haplotagging || (params.cnv && !params.use_qdnaseq)

    // Trigger gene summary if: gene summary requested, BED provided, and BED compatible
    def create_gene_summary = params.output_gene_summary && params.bed && gene_summary_bed

    reference = prepare_reference([
        "input_ref": params.ref,
        "output_cache": true,
        "output_mmi": false
    ])
    ref = reference.ref
    ref_index = reference.ref_idx
    ref_cache = reference.ref_cache

    // canonical ref and BAM channels to pass around to all processes
    ref_channel = ref
    | concat(ref_index)
    | concat(ref_cache)
    | flatten
    | buffer(size: 4)

    // Otherwise handle (u)BAM/CRAM
    if (!params.bam) {
        throw new Exception(colors.red + "Missing required --bam input argument." + colors.reset)
    }

    // ************************************************************************
    // Bail from the workflow for a reason we should have already specified
    if (!can_start){
        throw new Exception("The workflow could not be started.")
    }
    // ************************************************************************

    // Dummy optional file
    // TODO should be a channel?
    OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")

    Pinguscript.ping_start(nextflow, workflow, params)

    // Determine if (re)alignment is required for input BAM
    bam_channel = ingress(
        ref,
        ref_index,
        params.bam,
        ingress_ext,
    )

    // enforce_genome_build determines if getGenome should be run
    //   and can be used later to determine if a genome build was enforced
    // NOTE Logic for whether humvar should make a decision as to continue
    //   based on the genome build should be activated only by this boolean
    def enforce_genome_build = \
        // always check genome build for CNV and STR subworkflows
        // getGenome will take care of checking which build is required for the CNV flavours
        (params.cnv || params.str) \
        // or if annotating, check genome build when using SNP, SV or phasing
        // as SnpEff annotations are only provided for hg19 and hg38
        || (params.annotation && (params.snp || params.sv || params.phased))

    // Check if the genome build in the BAM is suitable for any workflows that have restrictions
    // NOTE getGenome will cause the workflow to terminate if the build is neither hg19 or hg38
    //   so it shouldn't be called if annotation is skipped to allow other genomes (including non-human)
    if (enforce_genome_build) {
        genome_build = getGenome(bam_channel)
    }
    else {
        genome_build = null
    }

    // Check for contamination, if MT is present.
    if (params.haplocheck){
        // First, let's get the mitogenome code.
        if (params.mitogenome){
            // Ensure that the given chromosome code is in the reference genome
            mt_code = ref_index
            | splitCsv(sep:'\t', header: false)
            | map{ it[0] }
            | filter{it == params.mitogenome}
            | ifEmpty{
                throw new Exception(colors.red + "Mitochondrial genome ${params.mitogenome} not present in the reference." + colors.reset)
            }
        } else {
            default_mt_codes = Channel.of(['chrM', 'Mt', 'MT']) | flatten
            mt_code = ref_index
            | splitCsv(sep:'\t', header: false)
            | map{ it[0] }
            | cross(default_mt_codes)
            | map{it[0]}
        }
        // Do not run if there are multiple mitochondrial codes.
        n_mt_codes = mt_code
        | count
        | subscribe {
            if (it != 1){
                throw new Exception(colors.red + "Unexpected number of mitochondrial chromosome found: ${it}." + colors.reset) 
            }
        }
        hap_check = haplocheck(bam_channel, ref_channel.collect(), mt_code)
        | ifEmpty{
            log.warn "Haplocheck failed to run. The workflow will continue, but will not output a contamination determination."
            file("$projectDir/data/OPTIONAL_FILE")
        }
    } else {
        // If haplocheck is not needed, use the predefined NV file.
        hap_check = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")
    }

    // Set BED (and create the default all chrom BED if necessary)
    // Make a second bed channel that won't be filtered based on coverage,
    // to be used as a final ROI filter
    bed = null
    using_user_bed = false
    if(params.bed){
        using_user_bed = true
        // Sanitise the input BED file
        input_bed = Channel.fromPath(params.bed, checkIfExists: true)

        bed = sanitise_bed(input_bed, ref_channel)
        roi_filter_bed = bed
    }
    else {
        bed = getAllChromosomesBed(ref_channel).all_chromosomes_bed
    }

    // mosdepth for depth traces -- passed into wf-snp :/

    mosdepth_input(bam_channel, bed, ref_channel, params.depth_window_size, create_gene_summary)
    mosdepth_stats = mosdepth_input.out.mosdepth_tuple
    mosdepth_summary = mosdepth_input.out.summary
    if (params.depth_intervals){
        mosdepth_perbase = mosdepth_input.out.perbase
    } else {
        mosdepth_perbase = Channel.empty()
    }

    if (create_gene_summary){
        coverage_summary = mosdepth_input.out.gene_summary
    }
    else {
        coverage_summary = Channel.empty()
    }

    // Determine if the coverage threshold is met to perform analysis.
    // If too low, it creates an empty input channel, 
    // avoiding the subsequent processes to do anything
    software_versions = getVersions()
    workflow_params = getParams()
    if (params.bam_min_coverage > 0){
        if (params.bed){
            // Filter out the data based on the individual region's coverage
            coverage_check = get_region_coverage(bed, mosdepth_stats)
            bed = coverage_check.filt_bed
            mosdepth_stats = coverage_check.mosdepth_tuple
        }
    } 
    bam_channel.set{pass_bam_channel}
    discarded_bams = Channel.empty()

    // Check and perform downsampling if needed.
    if (params.downsample_coverage){
        // Define reduction rate
        eval_downsampling(
            mosdepth_input.out.summary,
            params.bed ? mosdepth_stats.map{it[0]} : OPTIONAL
        )
        eval_downsampling.out.downsampling_ratio
            .splitCsv()
            .branch{
                subset: it[0] == 'true'
                ready: it[0] == 'false'
            }
            .set{ratio}

        // Define extension based on whether we are asking for CNV. If so,
        // use BAM, otherwise use what the user wants.
        downsampling_ext = pass_bam_channel.map{
            xam, xai, meta -> 
            convert_cram_to_bam ? ['bam', 'bai'] : desired_xam_ext
        }
        downsampling(pass_bam_channel, ref_channel, ratio.subset, downsampling_ext)

        // prepare ready files
        ratio.ready
            .combine(pass_bam_channel)
            .map{ready, ratio, xam, xai, meta -> [xam, xai, meta]}
            .branch{
                xam, xai, meta ->
                cram: xam.name.endsWith('.cram')
                bam: xam.name.endsWith('.bam')
            }
            .set{branched_bam_channel}

        // Convert aligned CRAMs that could not be downsampled to BAM if needed and mix with other ingested BAMs
        // Avoid issues with BAM being passed to `cram_to_bam`.
        ready_bam_channel = cram_to_bam(
            branched_bam_channel.cram,
            ref_channel.map { ref, index, cache, path -> [ref, index] }
        )
        | map { xam, xai, meta -> [xam, xai, meta + [output: false, is_cram: false]] }
        | mix(branched_bam_channel.bam)


        // Join allowing a remainder, so that only one for each is retained.
        // we drop all null, and due to the structure the joined channel can only be:
        // - [meta, null, xam, xai], or
        // - [meta, xam, xai, null]
        // Using it - null removes the inputs from the wrong channel, retaining 
        // Before merging properly, we first check that the merged channel size is not malformed
        downsampling.out.xam
            .join(ready_bam_channel, by:2, remainder: true)
            .filter{it.size() > 4}
            .subscribe{
                throw new Exception(colors.red + "Unexpected channel size when merging." + colors.reset) 
            }
        // If this passes, then we can create the proper channel.
        downsampling.out.xam
            .join(ready_bam_channel, by:2, remainder: true)
            .map{it - null}
            .map{meta, xam, xai -> [xam, xai, meta]}
            .set{pass_bam_channel}

        // Prepare the output files for mosdepth.
        // First, we compute the depth for the downsampled files, if it
        // exists 
        mosdepth_downsampled(downsampling.out, bed, ref_channel, params.depth_window_size, false)
        // Then, choose which output will be used in the report. 
        // If it needs to be subset, then the combined output exists, whereas 
        // the original mosdepth file is merged with the empty ready channel, leaving 
        // the correct file to output. Otherwise, the reverse happens and it emits 
        // the original mosdepth files. 
        mosdepth_summary = 
            mosdepth_downsampled.out.summary
                .combine(ratio.subset)
                .map{it[0]}
                .join(
                    mosdepth_input.out.summary
                        .combine(ratio.ready)
                        .map{it[0]}
                    , remainder: true
                    )
        mosdepth_stats = 
            mosdepth_downsampled.out.mosdepth_tuple
                .combine(ratio.subset)
                .map{[it[0], it[1], it[2]]}
                .join(
                    mosdepth_input.out.mosdepth_tuple
                        .combine(ratio.ready)
                        .map{[it[0], it[1], it[2]]}
                    , remainder: true
                    )
                .map{it - null}
        if (params.depth_intervals){
            mosdepth_perbase = 
                mosdepth_downsampled.out.perbase
                    .combine(ratio.subset)
                    .map{it[0]}
                    .join(
                        mosdepth_input.out.perbase
                            .combine(ratio.ready)
                            .map{it[0]}
                        , remainder: true
                        )
                    .map{it - null}
        } else {
            mosdepth_perbase = Channel.empty()
        }
    }

    // TODO downsampling should be incorporated to ingress to avoid
    //      call to bootleg readStats here
    // Run readStats depending on the downsampling, if requested.
    if (params.downsample_coverage){
        readStats(pass_bam_channel, bed, ref_channel)
    // Otherwise, use input bam
    } else {
        readStats(bam_channel, bed, ref_channel)
    }
    bam_stats = readStats.out.read_stats
    bam_flag = readStats.out.flagstat
    bam_hists = readStats.out.histograms
    // populate output json with ingressed runids and models
    bam_runids = readStats.out.runids
    bam_basecallers = readStats.out.basecallers
    ArrayList ingressed_run_ids = []
    bam_runids.splitText().subscribe(
        onNext: {
            ingressed_run_ids += it.strip()
        },
        onComplete: {
            params.wf["ingress.run_ids"] = ingressed_run_ids
        }
    )

    // Define depth_pass channel
    if (params.bam_min_coverage > 0){
        // If bam_min_coverage is > 0, then check the coverage
        if (params.bed){
            // Count the number of lines in the file to ensure that
            // there are intervals with enough coverage for downstream
            // analyses.
            n_lines = mosdepth_stats
            | map{ it[0] }
            | countLines()

            // Ensure that the data have enough region coverage
            // and intervals in the output coverage BED file.
            // First, load and split the summary file, keeping only
            // the `total_region` value (`total_region` and `total`
            // are identical in absence of a BED file).
            depth_pass = mosdepth_summary
                | splitCsv(sep: "\t", header: true)
                | filter{it -> it.chrom == "total_region"}
                // Extract the mean coverage as floating value
                | map{
                    it -> 
                    float mean = it.mean as float
                    [mean]}
                // Add line number in the coverage BED file
                | combine(n_lines)
                // Check if the coverage is appropriate
                | map {
                    mean, n_lines_v -> 
                    int n_lines = n_lines_v as int
                    boolean pass = mean > params.bam_min_coverage && n_lines > 0
                    [pass, mean]
                }

        // Without a BED, use summary values for the region
        } else {
            depth_pass = mosdepth_summary
                | splitCsv(sep: "\t", header: true)
                | filter{it -> it.chrom == "total_region"}
                | map{
                    it -> 
                    float mean = it.mean as float
                    boolean pass = mean > params.bam_min_coverage
                    [pass, mean]}
        }
    } else {
        // Otherwise, set all BAM to pass.
        depth_pass = bam_channel
            | map{ it -> [true, null] }
    }


    // Implement the BAM stats barrier after the pre-processing.
    // This will use the reads after the downsampling when requested.
    // Currently, it works using only the BAM coverage, but in the
    // future will allow to easily implement additional thresholds.
    filter = depth_pass
        .combine(pass_bam_channel)
        .branch{
            dp_pass, dp_val_env, bam, bai, meta ->
            pass: dp_pass && meta.has_mapped_reads
            not_pass: true
            }
    // Create the pass_bam_channel  channel when they pass
    filter.pass
        .map{it ->
            it.size > 0 ? [it[-3], it[-2], it[-1]] : it
        }
        .set{pass_bam_channel}

    // If it doesn't pass the minimum depth required, 
    // emit a bam channel of discarded bam files.
    filter.not_pass
        .subscribe {
            dp_pass, dp, bam, bai, meta ->
            // check where it failed
            def fail_depth = !meta.has_mapped_reads ? "No mapped reads." : dp < params.bam_min_coverage ? "Depth: ${dp} < ${params.bam_min_coverage}" : "Unknown."
            // Log where it failed
            log.error "File ${bam.getName()} will not be processed by the workflow because:\n - ${fail_depth}\n"
        }
    filter.not_pass
        .map{it ->
            it.size > 0 ? [it[-3], it[-2], it[-1]] : it
        }
        .set{discarded_bams}

    // Set biological sex to the user-provided sex if given
    // Otherwise, attempt to infer if genome_build is set (as we're likely going to need it)
    // NOTE You may feel compelled to add sex to the bam channel meta but then
    //      this will block any downstream access to the bam channel on infer_sex!
    if (params.sex) {
        sex = Channel.of(params.sex)
    }
    else if (genome_build) {
        log.warn "Inferring genetic sex of sample as params.sex was not provided."
        sex = infer_sex(mosdepth_summary)
    }
    else {
        sex = Channel.of(null)
    }

    // Create reports for pass and fail channels
    if (params.output_report){
        // Create passing bam report
        report_pass = pass_bam_channel
            .combine(bam_stats)
            .combine(bam_flag)
            .combine(bam_hists)
            .combine(mosdepth_stats.map{it[0]})
            .combine(mosdepth_summary)
            .combine(ref_channel)
            .combine(software_versions.collect())
            .combine(workflow_params)
            .flatten()
            .collect() | makeAlignmentReport
        // Create failing bam report
        report_fail = discarded_bams
            .combine(bam_stats)
            .combine(bam_flag)
            .combine(bam_hists)
            .combine(mosdepth_stats.map{it[0]})
            .combine(mosdepth_summary)
            .combine(ref_channel)
            .combine(software_versions.collect())
            .combine(workflow_params)
            .flatten()
            .collect() | failedQCReport
    } else {
        report_pass = Channel.empty()
        report_fail = Channel.empty()
    }

    // Set up BED for wf-human-snp, wf-human-str or run_haplotagging
    // CW-2383: we first call the SNPs to generate an haplotagged bam file for downstream analyses
    if (run_snp) {
        if(using_user_bed) {
            snp_bed = bed
        }
        else {
            // wf-human-snp uses OPTIONAL_FILE for empty bed for legacy reasons
            snp_bed = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE", checkIfExists: true)
        }

        if(params.clair3_model_path) {
            log.warn "Overriding Clair3 model with ${params.clair3_model_path}."
            clair3_model = Channel.fromPath(params.clair3_model_path, type: "dir", checkIfExists: true)
        }
        else {
            // Add back basecaller models, if available.
            // Combine each BAM channel with the appropriate basecaller file
            // Fetch the unique basecaller models and, if these are more than the
            // ones in the metadata, add them in there.
            // We do it in the snv scope as it is the only workflow relying on the
            // model, and given it has to wait for the readStats process, we try
            // minimizing the waits
            pass_bam_channel = pass_bam_channel
            | combine(bam_basecallers)
            | map{
                xam, xai, meta, bc ->
                def models = bc.splitText().collect { it.strip() }
                [xam, xai, meta + [basecall_models: models]]
            }

            // map basecalling model to clair3 model
            lookup_table = Channel.fromPath("${projectDir}/data/clair3_models.tsv", checkIfExists: true)

            // attempt to pull out basecaller_cfg from metadata
            def observed_pass_bam = 0
            metamap_basecaller_cfg = pass_bam_channel
                | map { xam, bai, meta ->
                    observed_pass_bam++ // keep count of observed BAMs to guard against emitting basecaller_cfg logging later on failed BAM
                    meta["basecall_models"]
                }
                | flatten  // squash lists

            // check returned basecaller list cardinality
            metamap_basecaller_cfg
                | count
                | map { int n_models ->
                    // n_models of 0 may indicate an empty pass_bam_channel, so
                    // we keep a count of observed_pass_bam to determine whether
                    // we should handle basecaller_cfg errors
                    if (n_models == 0 && observed_pass_bam > 0){
                        if (params.override_basecaller_cfg) {
                            log.warn "Found zero basecall_model in the input alignment header, falling back to the model provided with --override_basecaller_cfg: ${params.override_basecaller_cfg}"
                        }
                        else {
                            String input_data_err_msg = '''\
                            ################################################################################
                            # INPUT DATA PROBLEM
                            Your input alignment does not indicate the basecall model in the header and you
                            did not provide an alternative with --override_basecaller_cfg.

                            wf-human-variation requires the basecall model in order to automatically select
                            an appropriate SNP calling model.

                            ## Next steps
                            You must re-run the workflow specifying the basecaller model with the
                            --override_basecaller_cfg option.
                            ################################################################################
                            '''.stripIndent()
                            error input_data_err_msg
                        }
                    }
                    else if (n_models > 1){
                        String input_data_err_msg = '''\
                        ################################################################################
                        # INPUT DATA PROBLEM
                        Your input data contains reads basecalled with more than one basecaller model.

                        Our workflows automatically select appropriate configuration and models for
                        downstream tools for a given basecaller model. This cannot be done reliably when
                        reads with different basecaller models are mixed in the same data set.

                        ## Next steps
                        To use this workflow you must separate your input files, making sure all reads
                        are have been basecalled with the same basecaller model.
                        ################################################################################
                        '''.stripIndent()
                        error input_data_err_msg
                    }
                }

            // use params.override_basecaller_cfg as basecaller_cfg if provided, regardless of what was detected
            // we'll have exploded by now if we have no idea what the config is
            if (params.override_basecaller_cfg) {
                metamap_basecaller_cfg.map {
                    log.info "Detected basecaller_model: ${it}"
                    log.warn "Overriding basecaller_model: ${params.override_basecaller_cfg}"
                }
                basecaller_cfg = Channel.of(params.override_basecaller_cfg)
            }
            else {
                basecaller_cfg = metamap_basecaller_cfg
                    | map { log.info "Detected basecaller_model: ${it}"; it }
                    | ifEmpty(params.override_basecaller_cfg)
                    | map { log.info "Using basecaller_model: ${it}"; it }
                    | first  // unpack from list
            }

            clair3_model = lookup_clair3_model(lookup_table, basecaller_cfg).map {
                log.info "Autoselected Clair3 model: ${it[0]}" // use model name for log message
                it[1] // then just return the path to match the interface above
            }
        }

        clair_vcf = snp(
            pass_bam_channel,
            snp_bed,
            ref_channel,
            clair3_model,
            genome_build,
            haplotagged_output_fmt,
            run_haplotagging,
            using_user_bed,
            chromosome_codes
        )
    }
    
    // wf-human-sv
    // CW-2383: we then call SVs using either the pass bam or haplotagged bam, depending on the settings
    if(params.sv) {
        // If haplotagged bam is available and phase_snv is required, then phase.
        // Otherwise, use pass_bam_file (passing a haplotagged bam and not requiring phase_snv would
        // cause the workflow to wait for the tagged reads, but not enable phasing of sv since --phase 
        // won't be set; hence skip it if not required).
        if (run_haplotagging){
            sv_bam = clair_vcf.haplotagged_xam
        } else {
            sv_bam = pass_bam_channel
        }
        results_sv = sv(
            sv_bam,
            ref_channel,
            bed,
            mosdepth_input.out.summary,
            OPTIONAL,
            genome_build,
            chromosome_codes
        )
        artifacts = results_sv.report.flatten()
        sniffles_vcf = results_sv.sniffles_vcf
        json_sv = results_sv.sv_stats_json
        output_sv(artifacts)
    } else {
        json_sv = Channel.empty()
        sniffles_vcf = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE", checkIfExists: true)
    }

    // Then, we finish working on the SNPs by refining with SVs and annotating them. This is needed to
    // maximise the interaction between Clair3 and Sniffles.
    if (run_snp){
        // Channel of results.
        // We drop the raw .vcf(.tbi) file from Clair3 in it to then add back the files in the 
        // final_vcf channel, allowing for the latest file to be emitted.
        // Channel structure is
        /*  [
        *   [CRAM, CRAI]
        *   [vcf, tbi]
        *   [gvcf, tbi] (optional)
        *   haploblocks (optional)
         ] */
        // If first element ends with .vcf.gz, then discard it
        clair_vcf.clair3_results
            .filter{
                !it[0].name.endsWith('.vcf.gz')
            }
            .collect()
            .set{clair3_results}

        // Define which bam to use for final refinement
        if (run_haplotagging){
            snp_refinement_xam = clair_vcf.haplotagged_xam
        } else {
            snp_refinement_xam = pass_bam_channel
        }

        // Refine the SNP phase using SVs from Sniffles
        if (params.refine_snp_with_sv && params.sv){
            // Run by chromosome to reduce memory usage
            // Use collect on the reference, the SNP VCF
            // and the SV VCFs to ensure running on each contig. 
            refined_snps = refine_with_sv(
                ref_channel.collect(),
                clair_vcf.vcf_files.combine(clair_vcf.contigs),
                snp_refinement_xam | first,
                sniffles_vcf.collect()
            )
            final_snp_vcf = concat_refined_snp(
                refined_snps.collect(),
                "${params.sample_name}.wf_snp"
            )
        } else {
            // If refine_with_sv not requested, passthrough
            final_snp_vcf = clair_vcf.vcf_files
        }

        // Filter by BED, if provided
        if (params.bed) {
            final_snp_vcf_filtered = bed_filter(final_snp_vcf, roi_filter_bed, "snp", "vcf").filtered
        }
        else {
            final_snp_vcf_filtered = final_snp_vcf
        }

        // Run annotation, when requested.
        if (!params.annotation) {
            final_vcf = final_snp_vcf_filtered
            // no ClinVar VCF, pass empty VCF to makeReport
            clinvar_vcf = Channel.fromPath("${projectDir}/data/empty_clinvar.vcf")
        }
        else {
            // do annotation and get a list of ClinVar variants for the report
            // snpeff is slow so we'll just pass the whole VCF but annotate per contig
            annotations = annotate_snp_vcf(
                final_snp_vcf_filtered.combine(clair_vcf.contigs), genome_build.first(), "snp"
            )
            final_vcf = concat_snp_vcfs(annotations.collect(), "${params.sample_name}.wf_snp").final_vcf
            clinvar_vcf = sift_clinvar_snp_vcf(final_vcf, genome_build, "snp").final_vcf_clinvar
        }

        // Run vcf statistics on the final VCF file
        vcf_stats = vcfStats(final_vcf)

        // Prepare the report
        snp_reporting = report_snp(vcf_stats[0], clinvar_vcf)
        json_snp = snp_reporting.snp_stats_json
        if (params.output_report){
            snp_report = snp_reporting.report
        } else {
            snp_report = Channel.empty()
        }

        // Output for SNP
        snp_report
            .concat(clair3_results)
            .concat(final_vcf)
            .concat(clinvar_vcf)
            .flatten() | output_snp
    } else {
        json_snp = Channel.empty()
    }

    // wf-human-mod
    // Validate modified bam
    if (params.mod){
        // Perform validation on the initial BAM, to allow running on the
        // fragmented BAMs when phasing is required
        validate_modbam(pass_bam_channel, ref_channel)

        // Warn of input without modified base tags
        validate_modbam.out.branch{
            stdbam: it[-1] == '65'
            modbam: it[-1] == '0'
            }.set{validated_modbam}
        // Log warn if it is not modbam
        validated_modbam.stdbam.subscribe{
            it -> log.warn "Input ${it[0]} does not contain modified base tags. Was a modified basecalling model selected when basecalling this data?"
        }
        modbam_ch = validated_modbam.modbam
                .map{cram, crai, meta, code -> [cram, crai, meta]}

        // Compute the probabilities on the valid modbam
        modkit_probs = sample_probs(modbam_ch, ref_channel)

        // Save the other as input, keeping only the necessary elements
        if (run_haplotagging){
            modkit_bam = clair_vcf.str_bams
        } else {
            modkit_bam = modbam_ch
        }

        // If the input is not modBAM, the workflow won't process anything because the
        // filtering probabilities are not calculated, preventing downstream processes.
        results = mod(
            modkit_bam,  // Input BAM for modkit
            bam_flag,  // Flagstats used to define chromosomes to analyse
            chromosome_codes,  // Accepted chromosome codes for the human genome
            modkit_probs,  // modkit probabilities for filtering
            ref_channel,
            run_haplotagging  // Define if the data are haplotagged.
        )
        mod_stats = results.modkit.flatten()
    } else {
        mod_stats = Channel.empty()
    }

    // wf-human-cnv
    if (params.cnv) {
        // cnv calling with qdnaseq
        if (params.use_qdnaseq) {
            results_cnv = cnv_qdnaseq(
                pass_bam_channel,
                bam_stats,
                genome_build
            )
        // cnv calling with spectre
        } else {
            results_cnv = cnv_spectre(
                pass_bam_channel,
                ref_channel,
                clair_vcf.vcf_files,
                bed
            )
        }
        output_cnv(results_cnv)
    }

    // wf-human-str
    if (params.str) {
        // use haplotagged bam from snp() as input to str()
        bam_channel_str = clair_vcf.str_bams

        results_str = str(
          bam_channel_str,
          ref_channel,
          bam_stats,
          sex
        )
        output_str(results_str)
    }

    // Combine into a final JSON of analyses stats
    analyses_jsons = Channel.empty()
        | mix(
            json_snp,
            json_sv
        )
        | collect
        | ifEmpty(OPTIONAL)

    final_json = combine_metrics_json(
        analyses_jsons,
        bam_flag,
        bam_hists,
        mosdepth_stats,
        mosdepth_summary,
        hap_check,
        sex,
    )

    publish_artifact(
        // emit bams with the "to_align" meta tag
        // but only if haplotagging is not on
        bam_channel
        | filter( { it[2].to_align && !run_haplotagging} )
        | mix(
            bam_stats.flatten(),
            bam_flag.flatten(),
            mosdepth_stats.flatten(),
            mosdepth_summary.flatten(),
            mosdepth_perbase.flatten(),
            mod_stats.flatten(),
            report_pass.flatten(),
            report_fail.flatten(),
            final_json.flatten(),
            coverage_summary.flatten(),
            hap_check.flatten()
        )
        | filter{it.name != 'OPTIONAL_FILE'}
    )

}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
