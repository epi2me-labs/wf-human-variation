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
    configure_jbrowse;
    get_coverage; 
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
    output_cnv
    } from './modules/local/common'

include {
    ingress;
    cram_to_bam
} from './lib/_ingress.nf'

include {
    refine_with_sv;
    vcfStats;
    output_snp
} from "./modules/local/wf-human-snp.nf"

include { mod; validate_modbam} from './workflows/methyl'



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

    // check snp has basecaller config for clair3 model lookup
    if(params.snp || params.phased) {
        if(!params.basecaller_cfg && !params.clair3_model_path) {
            throw new Exception(colors.red + "You must provide a basecaller profile with --basecaller_cfg <profile> to ensure the right Clair3 model is chosen!" + colors.reset)
        }
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

    // switch workflow to BAM if calling CNV
    // CW-3324: Prevent ingress from running the CRAM to BAM conversion if
    //   downsampling is on, as the downsampling step can output BAM instead.
    if (params.cnv && params.use_qdnaseq) {
        log.warn "CNV calling subworkflow using QDNAseq does not support CRAM. You don't need to do anything, but we're just letting you know that:"
        log.warn "- If your input file is CRAM, it will be converted to a temporary BAM inside the workflow automatically."
        log.warn "- If your input requires alignment, the outputs will be saved to your output directory as BAM instead of CRAM."
        output_bam = params.downsample_coverage ? false : true
    }
    else {
        output_bam = false
    }

    // Define chromosome codes
        // Programmatically define chromosome codes.
    ArrayList chromosome_codes = []
    ArrayList chromosomes = [1..22] + ["X", "Y", "M", "MT"]
    for (N in chromosomes.flatten()){
        chromosome_codes += ["chr${N}", "${N}"]
    }

    // Trigger haplotagging
    def run_haplotagging = params.str || params.phased

    // Trigger CRAM to BAM conversion
    def convert_cram_to_bam = params.cnv && params.use_qdnaseq

    // Trigger the SNP workflow based on a range of different conditions:
    def run_snp = params.snp || run_haplotagging || (params.cnv && !params.use_qdnaseq)

    // Check ref and decompress if needed
    ref = null
    ref_index_fp = null
    if (params.ref.toLowerCase().endsWith("gz")) {
        // gzipped ref not supported by some downstream tools (e.g. cram_cache)
        // easier to just decompress and pass it around rather than confusing the user
        decompress_ref(file(params.ref))
        ref = decompress_ref.out.decompressed_ref
    }
    else {
        ref = Channel.fromPath(params.ref, checkIfExists: true)
        ref_index_fp = file(params.ref + ".fai")
    }

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

    // Create ref index if required
    if (!ref_index_fp || !ref_index_fp.exists()) {
        index_ref = index_ref_fai(ref)
        ref_index = index_ref.reference_index
    }
    else {
        ref_index = Channel.of(ref_index_fp)
    }

    Pinguscript.ping_start(nextflow, workflow, params)

    // Define extension based on presence/absence of downsampling.
    // This is only relevant when reads are remapped.
    def ingress_ext = convert_cram_to_bam && !params.downsample_coverage ? ['bam', 'bai'] : ['cram', 'crai']

    // Determine if (re)alignment is required for input BAM
    bam_channel = ingress(
        ref,
        ref_index,
        params.bam,
        ingress_ext,
    )
    
    // Check if the genome build in the BAM is suitable for any workflows that have restrictions
    // NOTE getGenome will exit non-zero if the build is neither hg19 or hg38, so it shouldn't be called
    // if annotation is skipped for snp, sv and phased, to allow other genomes (including non-human)
    // to be processed

    // always getGenome for CNV and STR
    if (params.cnv || params.str) {
        genome_build = getGenome(bam_channel)
    }
    // getGenome for SNP, SV and phased as long as annotation is not disabled
    else if ((params.snp || params.sv || params.phased) && params.annotation) {
        genome_build = getGenome(bam_channel)
    }
    else {
        genome_build = null
    }

    // Build ref cache for CRAM steps that do not take a reference
    cram_cache(ref)
    ref_cache = cram_cache.out.ref_cache
    ref_path = cram_cache.out.ref_path
    // canonical ref and BAM channels to pass around to all processes
    ref_channel = ref.concat(ref_index).concat(ref_cache).concat(ref_path).buffer(size: 4)

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

    mosdepth_input(bam_channel, bed, ref_channel, params.depth_window_size)
    mosdepth_stats = mosdepth_input.out.mosdepth_tuple
    mosdepth_summary = mosdepth_input.out.summary
    if (params.depth_intervals){
        mosdepth_perbase = mosdepth_input.out.perbase
    } else {
        mosdepth_perbase = Channel.from("$projectDir/data/OPTIONAL_FILE")
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
        } else {
            // Define if a dataset passes or not the filtering
            coverage_check = get_coverage(mosdepth_input.out.summary)
        }
        // Combine with the bam and branch by passing the depth filter
        coverage_check.pass
            .combine(bam_channel)
            .branch{ 
                pass: it[0] == "true" 
                not_pass: it[0] == "false" 
                }
            .set{bamdepth_filter}
        // Create the pass_bam_channel  channel when they pass
        bamdepth_filter.pass
            .map{it ->
                it.size > 0 ? [it[2], it[3], it[4]] : it
            }
            .set{pass_bam_channel}
        // If it doesn't pass the minimum depth required, 
        // emit a bam channel of discarded bam files.
        bamdepth_filter.not_pass
            .subscribe {
                log.error "ERROR: File ${it[2].getName()} will not be processed by the workflow as the detected coverage of ${it[1]}x is below the minimum coverage threshold of ${params.bam_min_coverage}x required for analysis."
            }
        bamdepth_filter.not_pass
            .map{it ->
                it.size > 0 ? [it[2], it[3], it[4]] : it
            }
            .set{discarded_bams}
    } else {
        // If the bam_min_depth is 0, then create alignment report for everything.
        bam_channel.set{pass_bam_channel}
        discarded_bams = Channel.empty()
    }
    // Set extensions for downstream analyses based on the input type
    // This will affect only the haplotagging.
    extensions = pass_bam_channel.map{
        xam, xai, meta -> 
        meta.is_cram ? ['cram', 'crai'] : ['bam', 'bai']
    }

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
        // use BAM, otherwise CRAM.
        downsampling_ext = pass_bam_channel.map{
            xam, xai, meta -> 
            convert_cram_to_bam ? ['bam', 'bai'] : ['cram', 'crai']
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
        mosdepth_downsampled(downsampling.out, bed, ref_channel, params.depth_window_size)
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
            mosdepth_perbase = Channel.from("$projectDir/data/OPTIONAL_FILE")
        }
    }

    // Run readStats depending on the downsampling, if requested.
    if (params.downsample_coverage){
        readStats(pass_bam_channel, bed, ref_channel)
    // Otherwise, use input bam
    } else {
        readStats(bam_channel, bed, ref_channel)
    }
    bam_stats = readStats.out.read_stats
    bam_flag = readStats.out.flagstat

    // Create reports for pass and fail channels
    // Create passing bam report
    report_pass = pass_bam_channel
                .combine(bam_stats)
                .combine(bam_flag)
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
                .combine(mosdepth_stats.map{it[0]})
                .combine(mosdepth_summary)
                .combine(ref_channel)
                .combine(software_versions.collect())
                .combine(workflow_params)
                .flatten()
                .collect() | failedQCReport

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
            // map basecalling model to clair3 model
            lookup_table = Channel.fromPath("${projectDir}/data/clair3_models.tsv", checkIfExists: true)
            // TODO basecaller_model_path
            clair3_model = lookup_clair3_model(lookup_table, params.basecaller_cfg).map {
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
            extensions,
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
        results = sv(
            sv_bam,
            ref_channel,
            bed,
            mosdepth_input.out.summary,
            OPTIONAL,
            genome_build,
            chromosome_codes
        )
        artifacts = results.report.flatten()
        sniffles_vcf = results.sniffles_vcf
        output_sv(artifacts)
    } else {
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

        // Define which bam to use
        if (run_haplotagging){
            snp_bam = clair_vcf.haplotagged_xam
        } else {
            snp_bam = pass_bam_channel
        }

        // Refine the SNP phase using SVs from Sniffles
        if (params.refine_snp_with_sv && params.sv){
            // Run by chromosome to reduce memory usage
            // Use collect on the reference, the SNP VCF
            // and the SV VCFs to ensure running on each contig. 
            refined_snps = refine_with_sv(
                ref_channel.collect(),
                clair_vcf.vcf_files.combine(clair_vcf.contigs),
                snp_bam.collect(),
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
        snp_report = report_snp(vcf_stats[0], clinvar_vcf)

        // Output for SNP
        snp_report
            .concat(clair3_results)
            .concat(final_vcf)
            .concat(clinvar_vcf)
            .flatten() | output_snp
    }

    // wf-human-mod
    // Validate modified bam
    if (params.mod){
        if (run_haplotagging){
            validate_modbam(clair_vcf.haplotagged_xam, ref_channel)
        } else {
            validate_modbam(pass_bam_channel, ref_channel)
        }

        // Warn of input without modified base tags
        validate_modbam.out.branch{
            stdbam: it[-1] == '65'
            modbam: it[-1] == '0'
            }.set{validation_results}
        // Log warn if it is not modbam
        validation_results.stdbam.subscribe{
            it -> log.warn "Input ${it[0]} does not contain modified base tags. Was a modified basecalling model selected when basecalling this data?"
        }

        // Save the other as input, keeping only the necessary elements
        validated_bam = validation_results.modbam.map{cram, crai, meta, code -> [cram, crai, meta]}

        results = mod(validated_bam, ref_channel)
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
                bed,
                chromosome_codes
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
          bam_stats
        )
        output_str(results_str)
    }

    jb_conf = configure_jbrowse(
        ref_channel,
        bam_channel,
    )

    publish_artifact(
        // CW-1033: remove environment variable from output
        ref_channel.map{it[0..2]}.flatten().mix(
            // emit bams with the "output" meta tag
            bam_channel.filter( { it[2].output } ),
            bam_stats.flatten(),
            bam_flag.flatten(),
            mosdepth_stats.flatten(),
            mosdepth_summary.flatten(),
            mosdepth_perbase.flatten(),
            mod_stats.flatten(),
            jb_conf.flatten(),
            report_pass.flatten(),
            report_fail.flatten()
        )
    )

}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
