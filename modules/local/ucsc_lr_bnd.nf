process GENERATE_UCSC_LR_BND_CANDIDATES {

    tag "${annotated_tsv.name}"

    label 'process_medium'

    container "${params.ucsc_lr_container}"

    publishDir "${params.out_dir}",
        mode: 'copy',
        saveAs: { filename ->
            filename.endsWith('.SVs.bnd_ucsc_match_summary.tsv') ||
            filename.endsWith('.SVs.bnd_ucsc_candidate_pairs.tsv') ?
                "annot_SVs/with_ucsc_lr/bnd/${filename}" :
            filename.endsWith('.CNVs.bnd_ucsc_match_summary.tsv') ||
            filename.endsWith('.CNVs.bnd_ucsc_candidate_pairs.tsv') ?
                "annot_CNVs/with_ucsc_lr/bnd/${filename}" :
            null
        }

    input:
    tuple val(meta), path(annotated_tsv)
    path reference_bed
    path generator_script

    output:
    tuple val(meta),
        path("*.annotated.UCSC_LR.tsv", includeInputs: true),
        path("*.bnd_ucsc_match_summary.tsv"),
        path("*.bnd_ucsc_candidate_pairs.tsv"),
        emit: bnd_candidates

    script:
    def variant_type = meta.variant_source?.toString()

    if (!(variant_type in ['SVs', 'CNVs'])) {
        throw new IllegalArgumentException(
            "GENERATE_UCSC_LR_BND_CANDIDATES: tipo no válido: ${variant_type}"
        )
    }

    def suffix = ".${variant_type}.annotated.UCSC_LR.tsv"

    if (!annotated_tsv.name.endsWith(suffix)) {
        throw new IllegalArgumentException(
            "El TSV no termina en ${suffix}: ${annotated_tsv.name}"
        )
    }

    def sample_prefix = annotated_tsv.name.substring(
        0,
        annotated_tsv.name.length() - suffix.length()
    )

    def summary_name =
        "${sample_prefix}.${variant_type}.bnd_ucsc_match_summary.tsv"

    def candidates_name =
        "${sample_prefix}.${variant_type}.bnd_ucsc_candidate_pairs.tsv"

    """
    set -euo pipefail

    python "${generator_script}" \
        --input "${annotated_tsv}" \
        --reference "${reference_bed}" \
        --summary "${summary_name}" \
        --candidates "${candidates_name}" \
        --max-distance 100

    test -s "${summary_name}"
    test -s "${candidates_name}"
    """
}

process FIX_UCSC_LR_BND {

    tag "${annotated_tsv.name}"

    label 'process_medium'

    container "${params.ucsc_lr_container}"

    publishDir "${params.out_dir}",
        mode: 'copy',
        saveAs: { filename ->
            filename.endsWith('.SVs.annotated.UCSC_LR_FINAL.tsv') ?
                "annot_SVs/with_ucsc_lr/${filename}" :
            filename.endsWith('.CNVs.annotated.UCSC_LR_FINAL.tsv') ?
                "annot_CNVs/with_ucsc_lr/${filename}" :
            filename.endsWith('.SVs.bnd_ucsc_candidates_final.tsv') ?
                "annot_SVs/with_ucsc_lr/bnd/${filename}" :
            filename.endsWith('.CNVs.bnd_ucsc_candidates_final.tsv') ?
                "annot_CNVs/with_ucsc_lr/bnd/${filename}" :
            null
        }

    input:
    tuple val(meta),
        path(annotated_tsv),
        path(summary_tsv),
        path(candidates_tsv)
    path reference_bed
    path fixer_script

    output:
    tuple val(meta),
        path("*.annotated.UCSC_LR_FINAL.tsv"),
        emit: annotated_ucsc_lr_final

    path("*.bnd_ucsc_candidates_final.tsv"),
        emit: candidate_report_final

    script:
    def variant_type = meta.variant_source?.toString()

    if (!(variant_type in ['SVs', 'CNVs'])) {
        throw new IllegalArgumentException(
            "FIX_UCSC_LR_BND: tipo de variante no válido: ${variant_type}"
        )
    }

    def suffix = ".${variant_type}.annotated.UCSC_LR.tsv"

    if (!annotated_tsv.name.endsWith(suffix)) {
        throw new IllegalArgumentException(
            "El TSV no termina en ${suffix}: ${annotated_tsv.name}"
        )
    }

    def sample_prefix = annotated_tsv.name.substring(
        0,
        annotated_tsv.name.length() - suffix.length()
    )

    def output_name =
        "${sample_prefix}.${variant_type}.annotated.UCSC_LR_FINAL.tsv"

    def candidate_output =
        "${sample_prefix}.${variant_type}.bnd_ucsc_candidates_final.tsv"

    """
    set -euo pipefail

    python "${fixer_script}" \
        --input "${annotated_tsv}" \
        --summary "${summary_tsv}" \
        --candidates "${candidates_tsv}" \
        --reference "${reference_bed}" \
        --output "${output_name}" \
        --candidate-output "${candidate_output}"

    test -s "${output_name}"
    test -s "${candidate_output}"
    """
}