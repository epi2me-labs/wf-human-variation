process modbam2bed {
    label "wf_human_methyl"
    cpus params.threads
    input:
        tuple path(alignment), path(alignment_index), val(alignment_meta), val(hap)
        tuple path(ref), path(ref_idx), path(ref_cache)
    output:
        path "${params.sample_name}*.methyl.*", emit: methyl_outputs

    script:
    def ref_path = "${ref_cache}/%2s/%2s/%s:" + System.getenv("REF_PATH")
    def modbam2bed_args = params.modbam2bed_args ?: ''
    // Define the haplotype, if required
    def is_phased = hap != 0 ? "--haplotype ${hap}" : ""
    def haptag = hap != 0 ? ".hap${hap}" : ""
    """
    export REF_PATH=${ref_path}
    modbam2bed \
        -e -m 5mC --cpg -t ${task.cpus} \
        ${ref} \
        ${alignment} \
        --aggregate \
        --prefix ${params.sample_name}${haptag}.methyl \
        ${modbam2bed_args} ${is_phased} \
        | bgzip -c > ${params.sample_name}${haptag}.methyl.cpg.bed.gz
    # also compress the aggregate
    bgzip ${params.sample_name}${haptag}.methyl.cpg.acc.bed
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


workflow methyl {
    take:
        alignment
        reference
    main:
        // CW-2457: modbam2bed with --haplotype ignores untagged reads, causing 
        // to not call regions that are not phased. Adding 0 as haplotype allows to
        // call sites on all regions. 
        haps = params.phase_methyl ? Channel.of(0, 1, 2) : Channel.of(0)
        out = modbam2bed(alignment.combine(haps), reference.collect())
        //bw = bigwig(out.methyl_bed, reference)
    emit:
        modbam2bed = out
        //bw = bw.methyl_bigwig
}
