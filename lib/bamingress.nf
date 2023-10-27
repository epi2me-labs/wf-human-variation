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

/**
 * Parse input arguments for `bam_ingress`.
 *
 * @param arguments: map with input arguments
 * @return: map of parsed arguments
 */
Map parse_arguments(Map arguments) {
    ArgumentParser parser = new ArgumentParser(
        args: [],
        kwargs: [
            "output_bam": false,
        ],
        name: "bam_ingress")
    return parser.parse_args(arguments)
}


// convert user CRAM to BAM for CNV
//  QDNASeq is stuck at an impasse as Rsamtools has not been formally updated to support CRAM
//  last i checked this is stuck on resolving an API call they shouldn't have been using
// as you can imagine i am not happy about this
process cram_to_bam {
    cpus 2
    input:
        tuple path(cram), path(crai)
        tuple path(ref), path(ref_idx)
    output:
        tuple path("${cram.simpleName}.bam"), path("${cram.simpleName}.bam.bai")
    script:
    log.info "Converting input CRAM to BAM..."
    """
    samtools view -@1 --reference ${ref} -b -o ${cram.simpleName}.bam##idx##${cram.simpleName}.bam.bai --write-index ${cram}
    """
}


process minimap2_alignment {
    memory '16 GB' // needs to be around 12 for humvar, bumped for breathing room (CW-2574) could be smaller if user isnt doing WGA
    cpus {params.ubam_map_threads + params.ubam_sort_threads + params.ubam_bam2fq_threads}
    input:
        path reference
        tuple path(reads), path(reads_idx), path(old_reference)
        tuple val(align_ext), val(index_ext) // either [bam, bai] or [cram, crai]
    output:
        tuple path("${params.sample_name}.${align_ext}"), path("${params.sample_name}.${align_ext}.${index_ext}"), emit: alignment
    script:
    def bam2fq_ref = old_reference.name != "OPTIONAL_FILE" ? "--reference ${old_reference}" : ''
    """
    samtools bam2fq -@ ${params.ubam_bam2fq_threads} -T 1 ${bam2fq_ref} ${reads} | minimap2 -y -t ${params.ubam_map_threads} -ax map-ont ${reference} - \
    | samtools sort -@ ${params.ubam_sort_threads} --write-index -o ${params.sample_name}.${align_ext}##idx##${params.sample_name}.${align_ext}.${index_ext} -O ${align_ext} --reference ${reference} -
    """
}


process check_for_alignment {

    input:
        tuple path(reference), path(ref_idx)
        tuple path(xam), path(xam_idx)
    output:
        tuple env(realign), env(xam_sq_len), path(xam), path(xam_idx)
    script:
        """
        realign=0
        workflow-glue check_sq_ref --xam ${xam} --ref ${reference} || realign=\$?

        # Count SQ lines to help determine if this is unaligned BAM/CRAM later
        # permit grep to exit 1 if the count is zero, otherwise blow up here
        xam_sq_len=\$(samtools view -H ${xam} | { grep -c '^@SQ' || [[ \$? == 1 ]]; })

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
        arguments // take care not to use the word args, friends!
    main:
        // check arguments
        Map margs = parse_arguments(arguments)

        // allow override of CRAM output in favour of BAM to permit garbage subworkflows
        def align_ext = "cram"
        def index_ext = "crai"
        if (margs.output_bam) {
            align_ext = "bam"
            index_ext = "bai"
        }
        alignment_exts = Channel.of([align_ext, index_ext]).collect()

        // load bam as channel
        xam_file = Channel.fromPath(bam_file_fp, checkIfExists: true)

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
            xam_idx_file = Channel.fromPath(bam_idx_fp)
        }
        else {
            samtools_index(xam_file)
            xam_idx_file = samtools_index.out
        }

        // check alignments against ref
        check_ref = ref_file.concat(ref_idx_file).buffer(size: 2) // don't wait for cram_cache to perform check_for_alignment
        check_xam_channel = xam_file.concat(xam_idx_file).buffer(size: 2)
        checked_alignments = check_for_alignment(check_ref, check_xam_channel)

        // non-zero exit code from realignment step will be interpreted as alignment required
        // fork BAMs into realign (for (re)alignment) and noalign subchannels
        // use `return` to drop the alignment status as it is not needed downstream
        checked_alignments.branch {
            realign: it[0] == '65' // (sq_len, xam, xidx)
                return tuple(it[1], it[2], it[3])
            noalign: it[0] == '0' // (xam, xidx)
                return tuple(it[2], it[3])
        }.set{alignment_fork}

        // fork already aligned CRAMs for conversion if necessary
        alignment_fork.noalign.branch {
            convert: margs.output_bam && it[0].name.toLowerCase().endsWith("cram")
            noconvert: true
        }.set{convert_fork}
        // and convert
        converted_bams = cram_to_bam(convert_fork.convert, check_ref)

        already_aligned_bams = convert_fork.noconvert.mix(converted_bams).map{
            // map already aligned bam to (xam_path, xam_index, xam_meta) tuple
            // setting the meta.output to false because the bam has not been realigned
            it -> tuple(it[0], it[1], create_metamap([
                output: false,
                is_cram: it[0].name.toLowerCase().endsWith("cram"),
            ]))
        }

        // Check old ref for CRAM requiring realignment
        alignment_fork.realign.map {
            old_ref = "${projectDir}/data/OPTIONAL_FILE"
            // it[0] will be the sq_len, ignore unaligned xams as the old_ref will not be needed
            // TODO -- check_for_alignment should probably emit whether this is cram or not but this will do for now
            if (it[0].toInteger() > 0 && it[1].name.toLowerCase().endsWith("cram")) {
                if (!params.old_ref) {
                    log.error("Input CRAM requires (re)alignment. You must provide a path to the original reference with '--old_ref <path>'")
                    throw new Exception("--old_ref not provided, but is required to realign a CRAM.")
                }
                old_ref = params.old_ref
            }
            // emit for alignment (xam, xidx, old_ref)
            tuple(it[1], it[2], file(old_ref))
        }.set {to_realign}

        // call minimap on bams that require (re)alignment
        // first map alignment_fork.realign to remove the sq_len from the fork
        // then map the result of minimap2_alignment to the canonical (reads, index, meta) tuple
        // files are returned in the format indicated by alignment_exts channel
        new_mapped_bams = minimap2_alignment(
            ref_file,
            to_realign,
            alignment_exts,
        ).map{
            // map newly aligned bam to (xam_path, xam_index, xam_meta) tuple
            // setting the meta.output to true because the bam has been (re)aligned
            it -> tuple(it[0], it[1], create_metamap([
                output: true,
                is_cram: !margs.output_bam,
            ]))
        }

        // mix realign and noalign forks back to canonical bam_channel with (reads, reads_idx, meta) format
        bam_channel = already_aligned_bams.mix(new_mapped_bams)

    emit:
        bam_channel
}
