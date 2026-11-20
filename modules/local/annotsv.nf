process ANNOTSV_SV {

    tag "${variant_vcf.name}"

    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/annotsv:3.5.9--hdfd78af_0' :
        'quay.io/biocontainers/annotsv:3.5.9--hdfd78af_0'}"

    publishDir "${params.out_dir}",
        mode: 'copy',
        saveAs: { filename ->
            filename.endsWith('.SVs.annotated.tsv')  ? "annot_SVs/${filename}" :
            filename.endsWith('.CNVs.annotated.tsv') ? "annot_CNVs/${filename}" :
            null
        }

    input:
    tuple val(meta), path(variant_vcf)

    output:
    tuple val(meta), path("*.annotated.tsv"), emit: annotated_SVs
    path "versions.yml", emit: versions

    script:
    def variant_type = meta.variant_source?.toString()

    if (!(variant_type in ['SVs', 'CNVs'])) {
        error "ANNOTSV_SV: meta.variant_source debe ser 'SVs' o 'CNVs'. Valor recibido: ${variant_type}"
    }

    def prefix = variant_vcf.name
        .replaceFirst(/\.vcf\.gz$/, '')
        .replaceFirst(/\.vcf$/, '')
        .replaceFirst(/\.(wf_)?(sv|cnv)$/, '')

    """
    AnnotSV \\
        -SVinputFile ${variant_vcf} \\
        -outputDir . \\
        -outputFile ${prefix}.${variant_type}.annotated.tsv \\
        -SVinputInfo 1 \\
        -genomeBuild GRCh38 \\
        -annotationsDir "${params.annotsv_annotations}" \\
        -txFile "${params.annotsv_txfile}"

    cat <<-END_VERSIONS > versions.yml
    "ANNOTSV_SV":
        annotsv: "3.5.9"
    END_VERSIONS
    """
}