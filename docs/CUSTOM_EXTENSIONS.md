# Custom extensions

This document describes the extensions added to the custom-6mA-svanna
branch of wf-human-variation, based on upstream version 2.8.1.

## Requirements

- Nextflow >= 24.04.0.
- Apptainer/Singularity for container execution.
- A compatible GRCh38 reference for the human-specific extensions.
- External software, databases and annotation resources must be obtained
  separately by the user, subject to their respective licences.

## 6mA modified-base analysis

The modified-base workflow is configured to call 6mA using:

    modkit pileup --modified-bases 6mA

No CpG filtering or strand combining is applied. Strand-specific records
are retained in the bedMethyl output.

The `--force_strand` parameter controls the use of
`--negative-strand-values` when converting bedMethyl files to bigWig.
It does not control strand retention in bedMethyl.

The workflow produces sample-level bedMethyl files and 6mA bigWig tracks.
When phasing is enabled, haplotype-resolved bedMethyl and bigWig files
are also generated.

## GPAP and HPO generation

When `--svanna` is enabled, the workflow queries GPAP for each sample
that passes the input QC filter and generates a set of positive HPO terms.

The sample alias must end in a numeric suffix separated by a hyphen.
The final `-<digits>` suffix is removed to construct the GPAP case ID.
For example, `sample-03` is queried using case ID `sample`.

The GPAP authentication value must be configured as a Nextflow secret
named `GPAP_AUTH`. Do not put credentials in parameter files, configuration
files, scripts or version control.

Configure the secret using the same Nextflow installation and `NXF_HOME`
that will be used to launch the workflow:

    read -r -s -p "GPAP_AUTH: " GPAP_AUTH
    printf '\n'
    nextflow secrets set GPAP_AUTH "$GPAP_AUTH"
    unset GPAP_AUTH

The secret is injected into the `GENERATE_HPO` task as an environment
variable. The compute environment must be able to access the GPAP service.

HPO identifiers are normalised to CURIE format. Obsolete terms are
resolved using replacement terms from the HPO ontology. Terms that are
invalid, absent from the ontology, obsolete without a replacement, or
ambiguous are discarded.

SvAnna requires at least one valid positive HPO term. A sample without
usable positive HPO terms will fail its SvAnna prioritisation task.

## SvAnna prioritisation

Enable with `--svanna true`. At least one of `--sv` or `--cnv` must
also be enabled.

The user must provide:

- `--svanna_container`: container image containing SvAnna and Java.
- `--svanna_jar`: path to the SvAnna CLI JAR.
- `--svanna_db`: path to the SvAnna database directory, including
  `phenotype/hp.json`.

SvAnna receives the sample VCF and the corresponding GPAP-derived HPO terms.
SV and CNV VCFs are prioritised independently.

Outputs are published under `out_dir/svanna/`:
- `<sample>.html`
- `<sample>.csv.gz`
- `<sample>.vcf.gz`

The SvAnna JAR and database are not included in this repository.
Obtain them directly from the project and comply with its non-commercial
licence terms. Commercial use requires separate permission.

## AnnotSV annotation

Enable with `--annotsv true`. At least one of `--sv` or `--cnv` must
also be enabled.

The workflow uses AnnotSV 3.5.9 to annotate the structural and copy
number variant VCFs.

The user must provide the GRCh38 AnnotSV resources:

- `--annotsv_annotations`: AnnotSV annotations directory.
- `--annotsv_txfile`: AnnotSV transcript annotation file.

SV and CNV annotations are published separately:

- `out_dir/annot_SVs/`
- `out_dir/annot_CNVs/`

The AnnotSV container is retrieved from the Galaxy/Biocontainers
distribution when the workflow runs. Annotation resources must be
obtained separately and must be compatible with GRCh38.

## UCSC long-read population frequencies and BND matching

UCSC-LR annotation is enabled together with `--annotsv true`.

The user must provide `--ucsc_lr_container`, an Apptainer/Singularity
image containing Python, pandas, NumPy, curl and the UCSC
`bigBedToBed` utility.

The image can be built from the supplied recipe and environment file:

    apptainer build --fakeroot ucsc-lr-frequency.sif containers/ucsc-lr/Apptainer.def

The workflow downloads the hg38 Noyvert long-read structural variant
reference from UCSC at runtime. Network access to UCSC is therefore
required by the execution environment. A BED conversion and SHA-256
checksums are also generated.

The frequency matching scripts use the following criteria:

- DEL, DUP and INV: reciprocal overlap threshold of 0.70.
- INS: breakpoint distance up to 100 bp and relative length difference
  up to 0.50.
- BND: candidate breakends are examined within 100 bp, including
  orientation compatibility. A frequency is assigned only when the
  matching is unambiguous.

Annotated frequency tables are published under:

- `out_dir/annot_SVs/with_ucsc_lr/`
- `out_dir/annot_CNVs/with_ucsc_lr/`

BND candidate summaries and candidate-pair reports are published in
the corresponding `bnd/` subdirectories. Final tables use the
`.annotated.UCSC_LR_FINAL.tsv` suffix.

The downloaded reference and its checksums are published under
`out_dir/annot_SVs/with_ucsc_lr/reference/`.

## Mahi-mahi STR analysis

Enable with `--mahi_mahi true`.

Mahi-mahi is run using the haplotagged alignment generated by the
workflow. Enabling this extension automatically triggers the SNP
and haplotagging subworkflow.

The implementation uses the built-in `hg38` repeat catalogue and
therefore requires a compatible GRCh38 reference.

Mahi-mahi requires Oxford Nanopore R10 data. The supported basecalling
models documented by the project include MinKNOW v5.2 HAC/SUP and
Dorado v2.0. Check the official Mahi-mahi documentation for the
requirements of the specific release being used.

The user must provide:

- `--mahi_mahi_container`: Apptainer/Singularity image containing
  Mahi-mahi.
- `--mahi_mahi_repeat_id`: optional repeat identifier to restrict
  genotyping to a specific locus.
- `--sex XX` or `--sex XY`: optional biological sex assignment.
  If omitted, the workflow attempts to infer sex from chromosome
  coverage. The upstream inference falls back to XX when it cannot
  determine sex, so this inferred value should not be treated as
  independent biological confirmation.

The workflow passes the haplotagged BAM or CRAM, reference FASTA
and index, sample alias and sex to Mahi-mahi.

Outputs are published under `out_dir/mahi-mahi/`:

- `<sample>.vcf`
- `<sample>.tsv`
- `<sample>.tagged.bam`
- `<sample>.config.yaml`
- `<sample>.log`
- `<sample>_plots/` (optional plots)

Mahi-mahi is distributed separately by Oxford Nanopore Technologies.
Its software licence is the Oxford Nanopore Technologies PLC Public
License v1.0 and its intended use is research only. The container
image is not included in this repository.

## Licensing and third-party resources

The original workflow and covered modifications remain subject to
the Oxford Nanopore Technologies PLC Public License v1.0 included
in this repository. In particular, the licence defines permitted
use as Research Purposes. Users must read and comply with the full
licence before using or redistributing the workflow.

Third-party components retain their own licences and terms:

- **SvAnna:** non-commercial use only under the terms provided by
  The Jackson Laboratory. Commercial use requires separate permission.
  The JAR, database and container image are not distributed here.
- **Mahi-mahi:** Oxford Nanopore Technologies PLC Public License v1.0;
  research use only. The container image is not distributed here.
- **AnnotSV:** GNU General Public License GPL-3.0-only. The workflow
  retrieves the AnnotSV 3.5.9 container from the Biocontainers
  distribution. Annotation resources must be obtained separately.
- **UCSC Noyvert long-read SV reference:** downloaded from the UCSC
  Genome Browser at runtime. Consult the UCSC track documentation
  and the associated publication for provenance and attribution.

No credentials, private parameter files, reference databases, BAM/CRAM
files, or pre-built Apptainer images should be committed to this repository.
