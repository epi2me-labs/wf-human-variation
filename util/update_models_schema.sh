#!/usr/bin/env bash
# Generate nextflow_schema with updated basecaller enumerations
#
# This script uses `nextflow config` to obtain the basecaller container,
# creates JSON arrays of the models using the container's list-models script
# and injects them with jq to create nextflow_schema.json.new.
set -euo pipefail
LC_ALL=C

TARGET=$1
ENGINE=$2

if ! command -v nextflow &> /dev/null
then
    # we should be in CI, nextflow is installed right here
    NEXTFLOW="./nextflow"
else
    NEXTFLOW=`which nextflow`
fi

# work out how to inspect the container contents
DORADO_CONTAINER=$(${NEXTFLOW} config -flat | grep "process.'withLabel:wf_basecalling'.container" | awk -F'= ' '{print $2}' | sed "s,',,g")
echo "# DORADO_CONTAINER=${DORADO_CONTAINER}"
if [ "$ENGINE" = "simg" ]; then
    CMD_PREFIX="singularity exec docker://${DORADO_CONTAINER}"
else
    CMD_PREFIX="docker run ${DORADO_CONTAINER}"
fi

# Get basecaller_cfg with a valid Clair3 model
awk 'NR>1 && $1 == "dorado" && $3 != "-" {print $2}' data/clair3_models.tsv | sort | uniq > valid_clair3_dorado_models.ls
awk 'NR>1 && $1 == "dorado" && $3 == "-" {print $2}' data/clair3_models.tsv | sort | uniq > invalid_clair3_dorado_models.ls
awk 'NR>1 && $1 != "dorado" && $3 != "-" {print $2}' data/clair3_models.tsv | sort | uniq > valid_clair3_nondorado_models.ls

# Get basecaller_cfg with a valid Dorado model in the container (should be sorted but lets be sure)
${CMD_PREFIX} list-models --simplex --only-names | grep -v '^rna' | sort > simplex_models.ls
comm -12 valid_clair3_dorado_models.ls simplex_models.ls > basecalling_models.ls # basecalling models WITH a clair3 model

cat valid_clair3_dorado_models.ls invalid_clair3_dorado_models.ls | sort > mentioned_dorado_models.ls
comm -13 mentioned_dorado_models.ls simplex_models.ls > unmentioned_dorado_models.ls
if [[ $(wc -l < unmentioned_dorado_models.ls) -ne 0 ]]; then
    echo "There are Dorado models in the container that are not mapped to a Clair3 model:"
    sed 's,^,* ,' unmentioned_dorado_models.ls
    exit 1
fi
comm -23 valid_clair3_dorado_models.ls simplex_models.ls | cat - valid_clair3_nondorado_models.ls | sed 's,^,clair3:,' | sort > clair3_only_models.ls

# Convert model lists to JSON arrays
SIMPLEX_MODELS=$(cat simplex_models.ls | sed '$a\custom' | cat - clair3_only_models.ls | jq -Rn '[inputs]')
MODBASE_MODELS=$(${CMD_PREFIX} list-models --modbase --only-names | sed '$a\custom' | jq -Rn '[inputs]')

# Inject JSON arrays to relevant schema enum
jq \
    -j \
    --indent 4 \
    --argjson simplex_models "${SIMPLEX_MODELS}" \
    --argjson modbase_models "${MODBASE_MODELS}" \
    '(.definitions.input.properties.basecaller_cfg.enum) = $simplex_models |
    (.definitions.basecalling_options.properties.remora_cfg.enum) = $modbase_models' \
    ${TARGET}/nextflow_schema.json > ${TARGET}/nextflow_schema.json.new

echo "# Updated schema generated, you should inspect it before adopting it!"
echo "diff ${TARGET}/nextflow_schema.json ${TARGET}/nextflow_schema.json.new"
echo "mv ${TARGET}/nextflow_schema.json.new ${TARGET}/nextflow_schema.json"
