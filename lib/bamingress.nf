include {
    minimap2_ubam;
} from '../modules/local/common'

import ArgumentParser

def create_metamap(Map arguments) {
    def parser = new ArgumentParser(
        args: [],
        kwargs:[
            "output": false,
            "is_cram": false,
        ],
        name:"create_metamap",
    )
    return parser.parse_args(arguments)
}


process check_for_alignment {

    input:
        tuple path(reference), path(ref_idx)
        tuple path(xam), path(xam_idx)
    output:
        tuple env(realign), path(xam), path(xam_idx)
    script:
        """
        realign=0
        workflow-glue check_sq_ref --xam ${xam} --ref ${reference} || realign=\$?

        # Allow EX_OK and EX_DATAERR, otherwise explode
        if [ \$realign -ne 0 ] && [ \$realign -ne 65 ]; then
            exit 1
        fi
        """
}

process samtools_index {
    input:
        path(xam)
    output:
        path("${xam}.*ai")
    script:
    """
    samtools index $xam
    """
}


workflow bam_ingress {

    take:
        ref_file
        ref_idx_file
        bam_file_fp
    main:
        is_cram = false

        // load bam as channel
        bam_file = Channel.fromPath(bam_file_fp, checkIfExists: true)

        if (bam_file_fp.toLowerCase().endsWith("cram")) {
            bam_idx_fp = params.bam + ".crai"
        }
        else {
            // Just assume BAM/BAI
            bam_idx_fp = params.bam + ".bai"
        }

        // load index as channel
        // index input if needed
        // TODO bit pointless to do it for ubam but we can come back to this
        if(file(bam_idx_fp).exists()) {
            bam_idx_file = Channel.fromPath(bam_idx_fp)
        }
        else {
            samtools_index(bam_file)
            bam_idx_file = samtools_index.out
        }

        check_bam_channel = bam_file.concat(bam_idx_file).buffer(size: 2)

        check_ref = ref_file.concat(ref_idx_file).buffer(size: 2) // don't wait for cram_cache to perform check_for_alignment
        checked_alignments = check_for_alignment(check_ref, check_bam_channel)

        // non-zero exit code from realignment step will be interpreted as alignment required
        // fork BAMs into realign (for (re)alignment) and noalign subchannels
        checked_alignments.branch {
            realign: it[0] == '65'
            noalign: it[0] == '0'
        }.set{alignment_fork}

        already_aligned_bams = alignment_fork.noalign.map{
            // map already aligned bam to (xam_path, xam_index, xam_meta) tuple
            // setting the meta.output to false because the bam has not been realigned
            it -> tuple(it[1], it[2], create_metamap([
                output: false,
                is_cram: is_cram,
            ]))
        }

        // Check old ref
        // TODO unaligned CRAM will trigger a request for old_ref which doesn't really make sense
        int to_realign = 0
        alignment_fork.realign.subscribe onNext: { to_realign++ }, onComplete: {
            if(to_realign > 0 && is_cram && !params.old_ref) {
                log.error(colors.red + "Input CRAM requires (re)alignment. You must provide a path to the original reference with '--old_ref <path>'" + colors.reset)
                throw new Exception("--old_ref not provided.")
            }
        }
        old_ref = Channel.fromPath("${projectDir}/data/OPTIONAL_FILE", checkIfExists: true)
        if(is_cram && to_realign > 0) {
            old_ref = Channel.fromPath(params.old_ref, checkIfExists: true)
        }

        // call minimap on bams that require (re)alignment
        // first map alignment_fork.realign to remove the exit status from the fork
        // then map the result of minimap2_ubam to the canonical (reads, index, meta) tuple
        new_mapped_bams = minimap2_ubam(ref_file, old_ref, alignment_fork.realign.map{ it -> tuple(it[1], it[2]) }).map{
            // map newly aligned bam to (xam_path, xam_index, xam_meta) tuple
            // setting the meta.output to true because the bam been (re)aligned
            it -> tuple(it[0], it[1], create_metamap([
                output: true,
                is_cram: is_cram,
            ]))
        }

        // mix realign and noalign forks back to canonical bam_channel with (reads, reads_idx, meta) format
        bam_channel = already_aligned_bams.mix(new_mapped_bams)

    emit:
        bam_channel
}
