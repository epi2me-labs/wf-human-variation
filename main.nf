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

include { partners } from './workflows/partners'

include {
    mosdepth as mosdepth_input;
    mosdepth as mosdepth_downsampled;
    mosdepth as mosdepth_coverage;
    readStats;
    getAllChromosomesBed;
    publish_artifact;
    get_region_coverage;
    failedQCReport; 
    makeAlignmentReport; 
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
    sanitise_bed as sanitise_coverage_bed;
    combine_metrics_json;
    output_cnv;
    infer_sex;
    haplocheck;
} from './modules/local/common'

include {
    getParams;
} from './lib/common.nf'

include {
    detect_basecall_model
} from './lib/model.nf'

include {
    ingress;
    cram_to_bam;
} from './lib/_ingress.nf'

include {
    igv
} from './lib/igv.nf'

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

    // If coverage summaries are requested, check if BED files are provided and warn if they don't have 4 columns
    // Set create_bed_summary and create_coverage_bed_summary accordingly so we can avoid running mosdepth on incompatible BED files

    def create_bed_summary = false
    def create_coverage_bed_summary = false

    // check if a BED file has at least 4 columns
    def check_bed_has_name_col = { bed_file ->
        if (bed_file) {
            def col_size = file(bed_file).splitCsv(sep: '\t').first().size
            if (col_size < 4) {
                return false
            }
            return true
        }
        return false
    }

    // Check and set create_bed_summary if --bed provided
    if (params.bed) {
        create_bed_summary = check_bed_has_name_col(params.bed)
    }

    // check and set create_coverage_bed_summary if --coverage_bed provided
    // exit workflow if this BED file doesn't have at least 4 columns
    if (params.coverage_bed) {
        if (check_bed_has_name_col(params.coverage_bed)) {
            create_coverage_bed_summary = true
        }
        else {
            log.error (colors.red + "The provided BED file (${params.coverage_bed}) has fewer than 4 columns, and therefore a coverage summary can not be generated." + colors.reset)
            can_start = false
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
    
    // Combine data for partners
    def run_partners = params.partner?: false
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

    reference = prepare_reference([
        "input_ref": params.ref,
        "output_cache": true,
        "output_mmi": false
    ])
    ref = reference.ref
    ref_index = reference.ref_idx
    ref_cache = reference.ref_cache
    ref_gzindex = reference.ref_gzidx
    // canonical ref and BAM channels to pass around to all processes
    ref_channel = ref
    | concat(ref_index)
    | concat(ref_cache)
    | flatten
    | buffer(size: 4)

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
        // getGenome will take care of checking which build is required for STR
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

    // Set coverage_bed, used only to generate coverage metrics
    coverage_bed = null
    if (params.coverage_bed) {
        // Sanitise the coverage BED file
        input_coverage_bed = Channel.fromPath(params.coverage_bed, checkIfExists: true)

        coverage_bed = sanitise_coverage_bed(input_coverage_bed, ref_channel)

    }
    else {
        coverage_bed = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")
    }

    // mosdepth for depth traces -- passed into wf-snp :/

    mosdepth_input(bam_channel, bed, ref_channel, params.depth_window_size, create_bed_summary, "bed")
    mosdepth_stats = mosdepth_input.out.mosdepth_tuple
    mosdepth_summary = mosdepth_input.out.summary
    if (params.depth_intervals){
        mosdepth_perbase = mosdepth_input.out.perbase
    } else {
        mosdepth_perbase = Channel.empty()
    }
    

    if (create_bed_summary){
        bed_summary = mosdepth_input.out.bed_summary
    }
    else {
        bed_summary = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")
    }

    // if requested, run mosdepth again to generate coverage summary for `--coverage_bed`
    if (create_coverage_bed_summary){
        mosdepth_coverage(bam_channel, coverage_bed, ref_channel, params.depth_window_size, create_coverage_bed_summary, "coverage_bed")
        coverage_bed_summary = mosdepth_coverage.out.bed_summary
    }
    else {
        coverage_bed_summary = Channel.fromPath("$projectDir/data/OPTIONAL_FILE")
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
            params.bed ? mosdepth_stats.map{it[1]} : OPTIONAL
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
        mosdepth_downsampled(downsampling.out, bed, ref_channel, params.depth_window_size, false, "bed")
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
                .map{[it[0], it[1], it[2], it[3]]}
                .join(
                    mosdepth_input.out.mosdepth_tuple
                        .combine(ratio.ready)
                        .map{[it[0], it[1], it[2], it[3]]}
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
    // Also check if using_user_bed is true, in which case pass the sanitised 
    // BED to readStats, rather than the filtered BED
    if (params.downsample_coverage) {
        readStats(
            pass_bam_channel,
            using_user_bed ? roi_filter_bed : bed,
            ref_channel
        )
    } else {
        readStats(
            bam_channel,
            using_user_bed ? roi_filter_bed : bed,
            ref_channel
        )
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
            | map{ it[1] }
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
            .combine(mosdepth_stats.map{it[1]})
            .combine(mosdepth_summary)
            .combine(ref_channel)
            .combine(software_versions.collect())
            .combine(workflow_params)
            .combine(Channel.value(using_user_bed))
            .combine(bed_summary)
            .combine(coverage_bed_summary)
            .flatten()
            .collect() | makeAlignmentReport
        // Create failing bam report
        report_fail = discarded_bams
            .combine(bam_stats)
            .combine(bam_flag)
            .combine(bam_hists)
            .combine(mosdepth_stats.map{it[1]})
            .combine(mosdepth_summary)
            .combine(ref_channel)
            .combine(software_versions.collect())
            .combine(workflow_params)
            .combine(Channel.value(using_user_bed))
            .combine(bed_summary)
            .combine(coverage_bed_summary)
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
            detect_basecall_model(pass_bam_channel, bam_basecallers)
            basecaller_cfg = detect_basecall_model.out.basecaller_cfg
            pass_bam_channel = detect_basecall_model.out.bam_channel

            // Get Clair3 model
            clair3_model = lookup_clair3_model(
                Channel.fromPath("${projectDir}/data/clair3_models.tsv", checkIfExists: true),
                basecaller_cfg
            )
            | map {
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
            chromosome_codes,
            workflow_params
        )
        artifacts = results_sv.report.flatten()
        sniffles_vcf = results_sv.sniffles_vcf
        json_sv = results_sv.sv_stats_json
        sv_vcf = results_sv.for_phasing
        output_sv(artifacts)
    } else {
        json_sv = Channel.empty()
        sv_vcf = Channel.empty()
        sniffles_vcf = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE", checkIfExists: true)
    }

    // Then, we finish working on the SNPs by refining with SVs and annotating them. This is needed to
    // maximise the interaction between Clair3 and Sniffles.
    if (run_snp){
        // Channel of results.
        // We drop the raw .vcf(.tbi) file from Clair3 in it to then add back the files in the 
        // snp_vcf channel, allowing for the latest file to be emitted.
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
                sniffles_vcf.map{meta, vcf -> vcf}.collect()
            )
            final_snp_vcf = concat_refined_snp(
                refined_snps.map{ meta, vcf, tbi -> [meta, vcf]}.groupTuple(),
                "wf_snp"
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
            snp_vcf = final_snp_vcf_filtered
            // no ClinVar VCF, pass empty VCF to makeReport
            clinvar_vcf = Channel.fromPath("${projectDir}/data/empty_clinvar.vcf")
        }
        else {
            // do annotation and get a list of ClinVar variants for the report
            // snpeff is slow so we'll just pass the whole VCF but annotate per contig
            annotations = annotate_snp_vcf(
                final_snp_vcf_filtered.combine(clair_vcf.contigs), genome_build.first(), "snp"
            )
            snp_vcf = concat_snp_vcfs(annotations.map{ meta, vcf, tbi -> [meta,vcf]}.groupTuple(), "wf_snp").final_vcf
            
            sift_clinvar_snp_vcf(snp_vcf, genome_build, "snp")
            clinvar_vcf = sift_clinvar_snp_vcf.out.final_vcf_clinvar.map{ vcf, tbi -> vcf }
        }

        // Run vcf statistics on the final VCF file
        vcf_stats = vcfStats(snp_vcf)

        // Prepare the report
        snp_reporting = report_snp(vcf_stats, clinvar_vcf, workflow_params)
        json_snp = snp_reporting.snp_stats_json
        if (params.output_report){
            snp_report = snp_reporting.report
        } else {
            snp_report = Channel.empty()
        }

        // Output for SNP
        snp_report
            .concat(clair3_results)
            .concat(snp_vcf.map{meta, vcf, tbi -> [vcf, tbi]})
            .flatten() | output_snp
    } else {
        json_snp = Channel.empty()
        snp_vcf = Channel.empty()
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
                genome_build,
                workflow_params
            )
        // cnv calling with spectre
        } else {
            results_cnv = cnv_spectre(
                pass_bam_channel,
                ref_channel,
                clair_vcf.vcf_files,
                bed,
                genome_build,
                workflow_params
            )
        }
        cnv_vcf = results_cnv.cnv_vcf
        output_cnv(results_cnv.output)
    } else {
        cnv_vcf = Channel.empty()
    }

    // wf-human-str
    if (params.str) {
        // use haplotagged bam from snp() as input to str()
        bam_channel_str = clair_vcf.str_bams

        results_str = str(
          bam_channel_str,
          ref_channel,
          bam_stats,
          sex,
          workflow_params
        )
        str_vcf = results_str.str_vcf
        output_str(results_str.output)
    } else {
        str_vcf = Channel.empty()
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

    // If the workflow is set to run for the partner, then execute the merging
    if (run_partners){
        partners(
            snp_vcf,
            sv_vcf,
            cnv_vcf,
            str_vcf,
            pass_bam_channel,
            run_haplotagging
        )
    }

    // Prepare IGV viewer
    if (params.igv){
        // Indicate which output files should be displayed
        // Note igv() is not responsible for publishing these files
        igv_out = ref_channel
            // Add gzipped reference indexes
            | combine(ref_gzindex | ifEmpty([null, null, null]))
            | map {
                fasta, fai, cache, path_env, gzref, gzfai, gzi ->
                if (gzref){
                    [gzref, gzfai, gzi]
                } else {
                    [fasta, fai]
                }
            }
            | mix(
                snp_vcf | map { meta, vcf, tbi -> [vcf, tbi] },
                sv_vcf | map { meta, vcf, tbi -> [vcf, tbi] },
                str_vcf | map { meta, vcf, tbi -> [vcf, tbi] },
                cnv_vcf | map { meta, vcf, tbi -> [vcf, tbi] },
                // set correct BAM for IGV depending on whether haplotagging requested
                // or alignment carried out - if neither then fall back to the original
                // unchanged BAM
                (run_haplotagging 
                    ? clair_vcf.haplotagged_xam | map { xam, xai, meta -> [xam, xai] } 
                    : bam_channel | map { 
                        xam, xai, meta -> [
                            meta.to_align ? xam : file(meta.src_xam),
                            meta.to_align || !meta.src_xai ? xai : file(meta.src_xai)
                        ] 
                    }
                ),
                // haplotagged bedMethyl
                mod_stats.take(2)
            )
            | igv
    } else {
        igv_out = Channel.empty()
    }

    publish_artifact(
        // emit bams with the "to_align" meta tag
        // but only if haplotagging is not on
        // drop meta which otherwise creates a spurious input.1 file
        bam_channel
        | filter( { it[2].to_align && !run_haplotagging} )
        | map { xam, xai, meta -> [xam, xai] }
        // Emit fasta or fai if they were changed from the input
        // (i.e. decompressed for fasta, generated for the fai)
        // if they are required for use with IGV
        | mix(
            ref_channel
            | map {
                fasta, fai, cache, path_env -> [fasta, fai]
            }
            | flatten
            | filter{
                it.toString().startsWith("${workflow.workDir}") && params.igv
            }
        )
        | mix(
            bam_stats.flatten(),
            bam_flag.flatten(),
            mosdepth_stats.map{ meta, bed, dist, threshold -> [bed, dist, threshold]}.flatten(),
            mosdepth_summary.flatten(),
            mosdepth_perbase.flatten(),
            mod_stats.flatten(),
            report_pass.flatten(),
            report_fail.flatten(),
            final_json.flatten(),
            bed_summary.flatten(),
            coverage_bed_summary.flatten(),
            hap_check.flatten(),
            igv_out.flatten()
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
