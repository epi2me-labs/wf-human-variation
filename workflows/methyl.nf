process modbam2bed {
    label "wf_human_methyl"
    cpus params.threads
    input:
        tuple path(alignment), path(alignment_index)
        tuple path(reference), path(reference_index), path(reference_cache)
    output:
        path "${params.sample_name}.cpg.bed.gz", emit: methyl_bed

    script:
    def ref_path = "${reference_cache}/%2s/%2s/%s:" + System.getenv("REF_PATH")
    """
    export REF_PATH=${ref_path}
    modbam2bed \
        -e -m 5mC --cpg -t ${task.cpus} \
        ${reference} \
        ${alignment} \
        | bgzip -c > ${params.sample_name}.cpg.bed.gz
    """
}

// NOTE uses base image
//process bigwig {
//    input:
//        path(methylbed)
//        tuple path(reference), path(reference_index), path(reference_cache)
//    output:
//        path "${params.sample_name}.cpg.bw", emit: methyl_bigwig
//    """
//    gzip -fd ${methylbed}
//    cut -f1,2 ${reference_index} > chrom.sizes
//    cut -f1,2,3,11 ${methylbed.baseName} | grep -v nan | sort -k1,1 -k2,2n > in.bedgraph
//    bedGraphToBigWig in.bedgraph chrom.sizes ${params.sample_name}.cpg.bw
//    """
//}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_methyl {
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}


workflow methyl {
    take:
        alignment
        reference
    main:
        out = modbam2bed(alignment, reference)
        //bw = bigwig(out.methyl_bed, reference)
    emit:
        bed = out.methyl_bed
        //bw = bw.methyl_bigwig
}
