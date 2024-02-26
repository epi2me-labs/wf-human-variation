
include {
    make_chunks;
    pileup_variants;
    aggregate_pileup_variants;
    select_het_snps;
    phase_contig;
    cat_haplotagged_contigs;
    get_qual_filter;
    create_candidates;
    evaluate_candidates;
    aggregate_full_align_variants;
    merge_pileup_and_full_vars;
    post_clair_phase_contig;
    post_clair_contig_haplotag;
    aggregate_all_variants;
    hap;
    getParams;
    getVersions;
    makeReport;
} from "../modules/local/wf-human-snp.nf"

include { 
    haploblocks as haploblocks_snp
} from '../modules/local/common.nf'

// workflow module
workflow snp {
    take:
        bam_channel
        bed
        ref
        model
        genome_build
        extensions
        run_haplotagging
        using_user_bed
        chromosome_codes
    main:

        // truncate bam channel to remove meta to keep compat with snp pipe
        bam = bam_channel.map{ it -> tuple(it[0], it[1]) }

        // Run preliminaries to find contigs and generate regions to process in
        // parallel.
        // > Step 0
        make_chunks(bam, ref, bed, model, chromosome_codes)
        chunks = make_chunks.out.chunks_file
            .splitText(){ 
                cols = (it =~ /(.+)\s(.+)\s(.+)/)[0]
                ["contig": cols[1], "chunk_id":cols[2], "total_chunks":cols[3]]}
        contigs = make_chunks.out.contigs_file.splitText() { it.trim() }
        cmd_file = make_chunks.out.cmd_file
        // use clair3 split beds if BED was provided
        if (using_user_bed) {
            split_beds = make_chunks.out.split_beds
        }
        else {
            split_beds = Channel.from("$projectDir/data/OPTIONAL_FILE").collect()
        }
        // Run the "pileup" caller on all chunks and collate results
        // > Step 1 
        pileup_variants(chunks, bam, ref, model, bed, cmd_file, split_beds)
        aggregate_pileup_variants(
            ref, pileup_variants.out.pileup_vcf_chunks.collect(),
            make_chunks.out.contigs_file, cmd_file)

        // Filter collated results to produce per-contig SNPs for phasing.
        // > Step 2
        select_het_snps(
            contigs,
            aggregate_pileup_variants.out.pileup_vcf,
            aggregate_pileup_variants.out.phase_qual)

        // Perform phasing for each contig.
        // `each` doesn't work with tuples, so we have to make the product ourselves
        phase_inputs = select_het_snps.out.het_snps_vcf
            .combine(bam).combine(ref)
        // > Step 3
        // > Step 4 (haplotagging is now done at the end of the workflow, rather than here)
        phase_contig(phase_inputs)
        phase_contig.out.phased_bam_and_vcf.set { phased_bam_and_vcf }

        // Find quality filter to select variants for "full alignment"
        // processing, then generate bed files containing the candidates.
        // > Step 5
        get_qual_filter(aggregate_pileup_variants.out.pileup_vcf)
        create_candidates(
            contigs, ref, 
            aggregate_pileup_variants.out.pileup_vcf,
            get_qual_filter.out.full_qual)

        // Run the "full alignment" network on candidates. Have to go through a
        // bit of a song and dance here to generate our input channels here
        // with various things duplicated (again because of limitations on 
        // `each` and tuples).
        // > Step 6
        candidate_beds = create_candidates.out.candidate_bed.flatMap {
            x ->
                // output globs can return a list or single item
                y = x[1]; if(! (y instanceof java.util.ArrayList)){y = [y]}
                // effectively duplicate chr for all beds - [chr, bed]
                y.collect { [x[0], it] } }
        // produce something emitting: [[chr, bam, bai, vcf], [chr20, bed], [ref, fai, cache], model]
        bams_beds_and_stuff = phased_bam_and_vcf
            .cross(candidate_beds)
            .combine(ref.map {it->[it]})
            .combine(model)
            .combine(cmd_file)
        // take the above and destructure it for easy reading
        bams_beds_and_stuff.multiMap {
            it ->
                bams: it[0]
                candidates: it[1]
                ref: it[2]
                model: it[3]
                cmd_file: it[4]
            }.set { mangled }
        // phew! Run all-the-things

        evaluate_candidates(
            mangled.bams, mangled.candidates, mangled.ref, mangled.model, mangled.cmd_file)

        // merge and sort all files for all chunks for all contigs
        // gvcf is optional, stuff an empty file in, so we have at least one
        // item to flatten/collect and tthis stage can run.
        gvcfs = pileup_variants.out.pileup_gvcf_chunks
            .flatten()
            .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
            .collect()
        pileup_variants.out.pileup_gvcf_chunks.flatten().collect()
        aggregate_full_align_variants(
            ref,
            evaluate_candidates.out.full_alignment.collect(),
            make_chunks.out.contigs_file,
            gvcfs,
            cmd_file)

        // merge "pileup" and "full alignment" variants, per contig
        // note: we never create per-contig VCFs, so this process
        //       take the whole genome VCFs and the list of contigs
        //       to produce per-contig VCFs which are then finally
        //       merge to yield the whole genome results.

        // First merge whole-genome results from pileup and full_alignment
        //   for each contig...
        // note: the candidate beds aren't actually used by the program for ONT
        // > Step 7
        non_var_gvcf = aggregate_full_align_variants.out.non_var_gvcf
            .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
        merge_pileup_and_full_vars(
            contigs, ref,
            aggregate_pileup_variants.out.pileup_vcf,
            aggregate_full_align_variants.out.full_aln_vcf,
            non_var_gvcf,
            candidate_beds.map {it->it[1] }.collect())

        // phased requires haplotagged bam to perform appropriate phasing
        // perform internal phasing only if snp+phase is requested, but not sv.
        // Otherwise use final joint phasing only.
        if (run_haplotagging) {
            post_clair_phase = merge_pileup_and_full_vars.out.merged_vcf
                .combine(bam)
                .combine(ref) |
                post_clair_phase_contig
            post_clair_phase.vcf
                .map { it -> [it[1]] }
                .set { final_vcfs }
            post_clair_phase.for_tagging | post_clair_contig_haplotag

            // intermediate ctg BAMs can flow to STR
            haplotagged_ctg_bams = post_clair_contig_haplotag.out.phased_bam
            // meanwhile cat all the intermediate ctg BAMs to desired XAM format for user output dir
            haplotagged_cat_xam = cat_haplotagged_contigs(
                haplotagged_ctg_bams.collect{it[1]}, // only need the BAMs themselves
                ref,
                extensions
            )

        } else {
            merge_pileup_and_full_vars.out.merged_vcf
                .map { it -> [it[1]] }
                .set { final_vcfs }
            // SNP only so we don't need these
            haplotagged_ctg_bams = Channel.empty()
            haplotagged_cat_xam = Channel.empty()
        }

        // ...then collate final per-contig VCFs for whole genome results
        gvcfs = merge_pileup_and_full_vars.out.merged_gvcf
            .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
        clair_final = aggregate_all_variants(
            ref,
            final_vcfs.collect(),
            gvcfs.collect(),
            params.phased,
            make_chunks.out.contigs_file,
            cmd_file)

        if (params.phased && (!params.sv || params.output_separate_phased)){
            hp_snp_blocks = haploblocks_snp(clair_final.final_vcf, 'snp')
        } else {
            hp_snp_blocks = Channel.empty()
        }

        // Define clair3 results, adding GVCF if needed
        if (params.GVCF){
            clair3_results = haplotagged_cat_xam.concat(clair_final.final_vcf).concat(clair_final.final_gvcf).concat(hp_snp_blocks)
        } else {
            clair3_results = haplotagged_cat_xam.concat(clair_final.final_vcf).concat(hp_snp_blocks)
        }

    emit:
        clair3_results = clair3_results
        str_bams = haplotagged_ctg_bams // intermediate haplotagged contigs used for STR
        vcf_files = clair_final.final_vcf
        haplotagged_xam = haplotagged_cat_xam.combine(bam_channel.map{it[2]}) // haplotagged XAM with meta appended
        contigs = contigs
}


// Reporting workflow
workflow report_snp {
    take:
        vcf_stats
        clinvar_vcf

    main:

        // reporting
        software_versions = getVersions()
        workflow_params = getParams()

        // Create report
        report = makeReport(
            vcf_stats, software_versions.collect(), workflow_params, clinvar_vcf)

    emit:
        report = report
}
