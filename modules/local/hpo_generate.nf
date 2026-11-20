
process GENERATE_HPO {

    tag "${meta.alias}"

    label 'process_low'

     container "${workflow.containerEngine == 'singularity' ?
        'docker://python:3.11.18-bookworm' :
        'python:3.11.18-bookworm'}"

    secret 'GPAP_AUTH'

    input:
    val(meta)
    path(hp_json)

    output:
    tuple val(meta), path('hpo_terms.txt'), emit: hpo_terms
    tuple val(meta), path('gpap_case.json'), emit: gpap_json
    tuple val(meta), path('hpo_terms_report.tsv'), emit: hpo_report

    script:
    def sampleId = meta.alias.toString()

    assert (sampleId ==~ /.+-\d+$/) :
        "No numeric suffix found in sample ID: ${sampleId}"

    def caseId = sampleId.replaceFirst(/-\d+$/, '')

    """
    set -euo pipefail

    echo "Consulting GPAP for case: ${caseId}"

    if [[ -z "\${GPAP_AUTH:-}" ]]; then
        echo "ERROR: GPAP_AUTH is not available in the task environment" >&2
        exit 2
    fi

    curl --fail --silent --show-error \
        -H "Authorization: Basic \$GPAP_AUTH" \
        "https://gpap.nasertic.es/impactgenomica/casos/${caseId}" \
        -o gpap_case.json

    python3 - "${hp_json}" <<'PY'
import json
import sys
import re
import csv

# Cargar respuesta de GPAP
with open("gpap_case.json") as f:
    data = json.load(f)

# Cargar ontologia HPO utilizada por SvAnna
with open(sys.argv[1]) as f:
    hp_data = json.load(f)

graph = hp_data["graphs"][0]
nodes = graph.get("nodes", [])

# Convertir identificadores OBO a formato HP:XXXXXXX
def to_curie(value):
    if value is None:
        return ""
    value = str(value).strip()
    local_id = value.rsplit("/", 1)[-1]

    if local_id.startswith("HP_"):
        return local_id.replace("HP_", "HP:", 1)

    return value

# Crear diccionario de terminos HPO
nodes_by_id = {}

for node in nodes:
    term_id = to_curie(node.get("id", ""))

    if re.fullmatch(r"HP:\\d{7}", term_id):
        nodes_by_id[term_id] = node

# Obtener fenotipos positivos de GPAP
positive_terms = set()

for feature in data.get("features", []):
    term_id = to_curie(feature.get("id"))
    presence = str(feature.get("presencia", "")).strip().lower()

    if term_id and presence in ("sí", "si", "yes"):
        positive_terms.add(str(term_id).strip())

# Resolver terminos obsoletos
REPLACED_BY = "http://purl.obolibrary.org/obo/IAO_0100001"

def resolve_hpo(term_id):

    current = term_id
    visited = set()

    while True:

        if current in visited:
            return None, "replacement_cycle"

        visited.add(current)

        node = nodes_by_id.get(current)

        if node is None:
            return None, "not_present_in_svanna_hpo"

        metadata = node.get("meta") or {}

        deprecated = (
            metadata.get("deprecated") is True
            or str(metadata.get("deprecated", "")).lower() == "true"
        )

        # El termino esta vigente
        if not deprecated:
            status = "kept" if current == term_id else "replaced"
            return current, status

        # Buscar reemplazo oficial
        replacements = set()

        for prop in metadata.get("basicPropertyValues", []):
            if prop.get("pred") == REPLACED_BY:
                replacements.add(to_curie(prop.get("val")))

        if len(replacements) == 0:
            return None, "obsolete_without_replacement"

        if len(replacements) > 1:
            return None, "multiple_replacements"

        # Comprobar tambien el estado del termino sustituto
        current = next(iter(replacements))


# Normalizar todos los terminos positivos
normalized_terms = set()
report = []

for term in sorted(positive_terms):

    if not re.fullmatch(r"HP:\\d{7}", term):
        report.append([term, "", "discarded", "invalid_hpo_id"])
        continue

    resolved, status = resolve_hpo(term)

    if resolved:
        normalized_terms.add(resolved)
        report.append([term, resolved, status, ""])
    else:
        report.append([term, "", "discarded", status])

# Generar fichero HPO definitivo
with open("hpo_terms.txt", "w") as f:
    for term in sorted(normalized_terms):
        f.write(term + "\\n")

# Guardar trazabilidad de sustituciones y descartes
with open("hpo_terms_report.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\\t")
    writer.writerow(["original_id", "final_id", "status", "reason"])
    writer.writerows(report)

n_replaced = sum(1 for r in report if r[2] == "replaced")
n_discarded = sum(1 for r in report if r[2] == "discarded")

print(f"HPO positivos originales: {len(positive_terms)}")
print(f"HPO finales para SvAnna: {len(normalized_terms)}")
print(f"HPO sustituidos: {n_replaced}")
print(f"HPO descartados: {n_discarded}")

print("HPO finales:")
for term in sorted(normalized_terms):
    print(term)

PY

    echo "HPO generados:"
    cat hpo_terms.txt

    echo "Informe de normalizacion:"
    cat hpo_terms_report.tsv
    """
}
