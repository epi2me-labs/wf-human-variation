// Combine VCFs for Geneyx
process publish_geneyx {
    publishDir "${params.outdirNonPhi}/integrations/geneyx/", mode: 'copy', pattern: "*"
    cpus 2

    input:
        tuple val(meta), path('snv.vcf.gz'), path('snv.vcf.gz.tbi')
        tuple val(meta), path('sv.vcf.gz'), path('sv.vcf.gz.tbi')
        tuple val(meta), path('cnv.vcf.gz'), path('cnv.vcf.gz.tbi')
        tuple val(meta), path('str.vcf.gz'), path('str.vcf.gz.tbi')

    output:
        tuple path("${meta.alias}.wf_snp.geneyx.vcf.gz"), path("${meta.alias}.wf_snp.geneyx.vcf.gz.tbi"), optional: true
        tuple path("${meta.alias}.wf_sv.geneyx.vcf.gz"), path("${meta.alias}.wf_sv.geneyx.vcf.gz.tbi"), optional: true

    script:
    // Define switches to prepare output data
    prep_snp = params.snp ?: false
    prep_sv = params.sv || params.cnv || params.str ?: false
    // Inputs for creation of unified VCF of SV/STR/CNV
    def sv = params.sv ? "-s sv.vcf.gz" : ""
    def cnv = params.cnv ? "-c cnv.vcf.gz" : ""
    def str = params.str ? "-r str.vcf.gz" : ""
    """
    if ${prep_sv}; then
        workflow-glue unify_vcf ${sv} ${cnv} ${str} -o ${meta.alias}.wf_sv.geneyx.vcf && \
            bcftools sort -O z -o ${meta.alias}.wf_sv.geneyx.vcf.gz ${meta.alias}.wf_sv.geneyx.vcf && \
            bcftools index -t ${meta.alias}.wf_sv.geneyx.vcf.gz && rm ${meta.alias}.wf_sv.geneyx.vcf
    fi

    # Prepare SNPs if available
    if ${prep_snp}; then
        cp snv.vcf.gz ${meta.alias}.wf_snp.geneyx.vcf.gz && bcftools index -t ${meta.alias}.wf_snp.geneyx.vcf.gz
    fi
    """
}

process publish_fabric {
    publishDir "${params.outdirNonPhi}/integrations/fabric/", mode: 'copy', pattern: "*"
    cpus 2

    input:
        tuple val(alias), path(vcfs), path(tbis)

    output:
        tuple path("${alias}.fabric.vcf.gz"), path("${alias}.fabric.vcf.gz.tbi")

    script:
    """
    bcftools concat ${vcfs} --rm-dups exact -a -O u | \
        bcftools sort -O z > ${alias}.fabric.vcf.gz && \
        bcftools index -t ${alias}.fabric.vcf.gz
    """
}

workflow partners {
    take:
        snv
        sv
        cnv
        str
        bam
        run_haplotagging
    main:
        combined_vcf_ch = Channel.empty()
        // Placeholder channel enables to run the process regardless of the analysis been run.
        placeholder_ch = bam
        | map{
            xam, xai, meta -> 
            [meta, file("${projectDir}/data/OPTIONAL_FILE"), file("${projectDir}/data/OPTIONAL_FILE")]
        }

        if (params.partner == "geneyx"){
            combined_vcf_ch = publish_geneyx(
                params.snp || run_haplotagging ? snv : placeholder_ch,
                params.sv ? sv : placeholder_ch,
                params.cnv ? cnv : placeholder_ch,
                params.str ? str : placeholder_ch
            )
        }
        if (params.partner == "fabric") {
            combined_vcf_ch = snv
            | mix(sv, cnv, str)
            // Group using sample alias
            | map {
                meta, vcf, tbi -> 
                [meta.alias, vcf, tbi]
            }
            | groupTuple
            | publish_fabric
        }
}