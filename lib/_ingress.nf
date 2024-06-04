/*
 * This sub-workflow acts as a wrapper around `ingress.nf`,
 * extending its functionality to do additional input preprocessing:
 *  - check whether the aligned BAM files match the input reference
 *  - (re-)align the input BAM if unaligned or not matching the reference
 */
include { xam_ingress } from './ingress.nf'

// convert user CRAM to BAM for CNV
//  QDNASeq is stuck at an impasse as Rsamtools has not been formally updated to support CRAM
//  last i checked this is stuck on resolving an API call they shouldn't have been using
// as you can imagine i am not happy about this
process cram_to_bam {
    cpus 2
    memory 4.GB
    input:
        tuple path(cram), path(crai), val(meta)
        tuple path(ref), path(ref_idx)
    output:
        tuple path("${cram.simpleName}.bam"), path("${cram.simpleName}.bam.bai"), val(meta)
    script:
    log.info "Converting input CRAM to BAM..."
    """
    samtools view -@1 --reference ${ref} -b -o ${cram.simpleName}.bam##idx##${cram.simpleName}.bam.bai --write-index ${cram}
    """
}

// Minimap2 mapping
process minimap2_alignment {
    cpus {params.ubam_map_threads + params.ubam_sort_threads + params.ubam_bam2fq_threads}
    memory { (32.GB * task.attempt) - 1.GB }
    maxRetries 1
    errorStrategy = {task.exitStatus in [137,140] ? 'retry' : 'finish'}
    input:
        path reference
        tuple val(meta), path(reads), path(reads_idx)
        tuple val(align_ext), val(index_ext) // either [bam, bai] or [cram, crai]

    output:
        tuple val(meta), env(has_maps), path("${params.sample_name}.${align_ext}"), path("${params.sample_name}.${align_ext}.${index_ext}"), emit: alignment
    script:
    """
    samtools bam2fq -@ ${params.ubam_bam2fq_threads} -T 1 ${reads} \
        | minimap2 -y -t ${params.ubam_map_threads} -ax map-ont --cap-kalloc 100m --cap-sw-mem 50m \
            ${reference} - \
        | samtools sort -@ ${params.ubam_sort_threads} \
            --write-index -o ${params.sample_name}.${align_ext}##idx##${params.sample_name}.${align_ext}.${index_ext} \
            -O ${align_ext} --reference ${reference} -

    # Check that the first line is not unmapped
    REF_PATH=${reference} workflow-glue check_mapped_reads --xam ${params.sample_name}.${align_ext} > env.vars
    source env.vars
    """
}

// Check if the XAM header refers to the input reference.
process check_for_alignment {
    label "wf_common"
    cpus 2
    memory 4.GB
    input:
        tuple path(reference), path(ref_idx)
        tuple val(meta), path(xam), path(xam_idx)
    output:
        tuple env(to_align), env(has_maps), val(meta), path(xam), path(xam_idx)
    script:
        """
        to_align=0
        workflow-glue check_sq_ref --xam ${xam} --ref ${reference} || to_align=\$?

        # Count SQ lines to help determine if this is unaligned BAM/CRAM later
        # permit grep to exit 1 if the count is zero, otherwise blow up here
        xam_sq_len=\$(samtools view -H ${xam} | { grep -c '^@SQ' || [[ \$? == 1 ]]; })

        # Check that the first line is not unmapped
        workflow-glue check_mapped_reads --xam ${xam} > env.vars
        source env.vars

        # Allow EX_OK and EX_DATAERR, otherwise explode
        if [ \$to_align -ne 0 ] && [ \$to_align -ne 65 ]; then
            exit 1
        fi
        """
}


// Create ingress workflow, wrapping functions from `ingress.nf` and
// performing additional processings.
workflow ingress {
    take:
        ref_file
        ref_idx_file
        bam_file_fp
        alignment_exts
    main:
        // load bam as channel
        // We do not want to perform statistics here as we will do them downstream.
        // We also want to keep all unaligned reads.
        ingressed_bam = xam_ingress([
            "input":bam_file_fp,
            "sample":params.sample_name,
            "sample_sheet":null,
            "analyse_unclassified":true,
            "keep_unaligned": true,
            "stats": false,
            "watch_path": false
        ])
        // Check that we have a single BAM/folder with BAMs in it.
        // by counting how many entries are in the channel.
        // If there are more than 1, then throw an error.
        ingressed_bam
            .count()
            .subscribe { int n_samples -> 
                if (n_samples > 1){
                    error "Too many samples found: (${n_samples}) in ${bam_file_fp}.\nPlease, ensure you provide a single folder with all the BAM files for a single individual."
                }
            }

        // Prepare reference channel
        check_ref = ref_file.combine(ref_idx_file) // don't wait for cram_cache to perform check_for_alignment
        
        // Index the input BAM, and then check if the header matches the reference.
        // Add also if it is a CRAM (for downstream compatibility) and move if they need realignment
        // to meta.
        checked_bam = check_for_alignment(
                check_ref,
                ingressed_bam.map{ it - null }
            ) |
            map{ to_align, has_maps, meta, xam, xai ->
                [meta + [to_align: to_align != '0', has_mapped_reads: has_maps == '1'], xam, xai]
            }

        // fork BAMs into to_align and noalign subchannels
        checked_bam.branch {
            meta, xam, xai -> 
            to_align: meta.is_unaligned || meta.to_align
            noalign: true
        }.set{alignment_fork}

        // call minimap on bams that require (re)alignment
        // then map the result of minimap2_alignment to the canonical (reads, index, meta) tuple
        new_mapped_bams = minimap2_alignment(ref_file, alignment_fork.to_align, alignment_exts).map{
            meta, has_maps, xam, xai -> 
            [meta + [has_mapped_reads: has_maps == '1'], xam, xai]
        }

        // map to (xam_path, xam_index, xam_meta) tuple
        // check for `is_cram` as the file might have been (re)aligned
        bam_channel = alignment_fork.noalign
            .mix(new_mapped_bams)
            .map{
                meta, xam, xai ->
                [xam, xai, meta + [is_cram: xam.name.endsWith('.cram')]]
            }

    emit:
        bam_channel
}
