
include {
    make_chunks;
    pileup_variants;
    aggregate_pileup_variants;
    select_het_snps;
    phase_contig;
    get_qual_filter;
    create_candidates;
    evaluate_candidates;
    aggregate_full_align_variants;
    merge_pileup_and_full_vars;
    post_clair_phase_contig;
    aggregate_all_variants;
    hap;
    mosdepth;
    readStats;
    getParams;
    getVersions;
    vcfStats;
    makeReport;
} from "../modules/wf-human-snp.nf"

// workflow module
workflow clair3 {
    take:
        bam
        bed
        ref
        model
    main:

        // Run preliminaries to find contigs and generate regions to process in
        // parallel.
        // > Step 0
        make_chunks(bam, ref, bed)
        chunks = make_chunks.out.chunks_file
            .splitText(){ 
                cols = (it =~ /(.+)\s(.+)\s(.+)/)[0]
                ["contig": cols[1], "chunk_id":cols[2], "total_chunks":cols[3]]}
        contigs = make_chunks.out.contigs_file.splitText() { it.trim() }

        // Run the "pileup" caller on all chunks and collate results
        // > Step 1 
        pileup_variants(chunks, bam, ref, model)
        aggregate_pileup_variants(
            ref, pileup_variants.out.pileup_vcf_chunks.collect(),
            make_chunks.out.contigs_file)

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
        // > Step 3 (Step 4 haplotag step removed in clair3 v0.1.11)
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
        // produce something emitting: [[chr, bam, bai, vcf], [chr20, bed], [ref, fai], model]
        bams_beds_and_stuff = phased_bam_and_vcf
            .cross(candidate_beds)
            .combine(ref.map {it->[it]})
            .combine(model)
        // take the above and destructure it for easy reading
        bams_beds_and_stuff.multiMap {
            it ->
                bams: it[0]
                candidates: it[1]
                ref: it[2]
                model: it[3]
            }.set { mangled }
        // phew! Run all-the-things
        evaluate_candidates(
            mangled.bams, mangled.candidates, mangled.ref, mangled.model)

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
            gvcfs)

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

        if (params.phase_vcf) {
            data = merge_pileup_and_full_vars.out.merged_vcf
                .combine(bam).combine(ref).view()
            post_clair_phase_contig(data)
                .map { it -> [it[1]] }
                .set { final_vcfs }
        } else {
            merge_pileup_and_full_vars.out.merged_vcf
                .map { it -> [it[1]] }
                .set { final_vcfs }
        }

        // ...then collate final per-contig VCFs for whole genome results
        gvcfs = merge_pileup_and_full_vars.out.merged_gvcf
            .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
        clair_final = aggregate_all_variants(
            ref,
            final_vcfs.collect(),
            gvcfs.collect(),
            params.phase_vcf,
            make_chunks.out.contigs_file)

        // reporting
        read_stats = readStats(bam)
        software_versions = getVersions()
        workflow_params = getParams()
        vcf_stats = vcfStats(clair_final)
        read_depth = mosdepth(bam, bed)
        report = makeReport(
            read_stats.stats, mosdepth.out.mosdepth_dist, vcf_stats[0],
            software_versions.collect(), workflow_params)
        telemetry = workflow_params

    emit:
        clair_final.concat(report)
        telemetry
}
