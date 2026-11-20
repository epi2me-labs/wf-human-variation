process svanna_prioritize {

    tag { sample_prefix }

    container "${params.svanna_container}"

    cpus 2
    memory 64.GB
    time '4h'

    publishDir "${params.out_dir}",
        mode: 'copy',
        saveAs: { filename -> "svanna/${filename}" }

    input:
        tuple val(meta), path(vcf), val(sample_prefix), path(hpo_path)
        path jar_path
        val db_path
        path normalizer

    output:
        path("${sample_prefix}.html", emit: html)
        path("${sample_prefix}.csv.gz", emit: csv)
        path("${sample_prefix}.vcf.gz", emit: vcf)

    script:
    """
    set -euo pipefail

    echo "Normalizando contigs BND para SvAnna"
    python3 "${normalizer}" "${vcf}" svanna_input.vcf

    tr -d '\\r' < "${hpo_path}" > hpo_terms_clean.txt

    HPO_ARGS=()
    while IFS= read -r hp || [[ -n "\$hp" ]]; do
        [[ -z "\$hp" ]] && continue
        [[ "\${hp:0:1}" == "#" ]] && continue
        HPO_ARGS+=(--phenotype-term "\$hp")
    done < hpo_terms_clean.txt

    if [[ "\${#HPO_ARGS[@]}" -eq 0 ]]; then
        echo "ERROR: no hay términos HPO válidos" >&2
        exit 2
    fi

    java -version

    java -Xms4g -Xmx48g -jar "${jar_path}" prioritize -d "${db_path}" --vcf svanna_input.vcf "\${HPO_ARGS[@]}" --output-format html,csv,vcf --out-dir "\$PWD" --prefix "${sample_prefix}"

    test -s "${sample_prefix}.html"
    test -s "${sample_prefix}.csv.gz"
    test -s "${sample_prefix}.vcf.gz"
    """
}
