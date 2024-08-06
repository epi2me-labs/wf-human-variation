+ Annotations for `--snp` and `--sv` are generated using [SnpEff](https://pcingola.github.io/SnpEff/). For `--snp`, additional [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) annotations are displayed in the report where available (please note, the report will not display any variants classified as 'Benign' or 'Likely benign', however these variants will be present in the
output VCF).
+ Aggregation of modified calls with `--mod` requires data to be basecalled with a model that includes base modifications, providing the `MM` and `ML` BAM tags
+ CRAM files generated within the workflow cannot be read without the corresponding reference
+ The STR workflow performs genotyping of specific repeats, which can be found [here](https://github.com/epi2me-labs/wf-human-variation/blob/master/data/wf_str_repeats.bed).
+ By default, SNPs and SVs will be called on chromosomes 1-22, X, Y and MT; to call sites on other sequences enable the `--include_all_ctgs` option.
+ While designed to work on human genomes, the workflow can be run on non-human species by setting `--cnv false --str false --annotation false --include_all_ctgs true`.
+ Ensure that the provided reference and BED files use the same chromosome coding (for example, that they both have the `chr` prefix, or they both to not have it).
+ If unaligned reads were provided, the workflow will output a CRAM file (or BAM if the user runs the `--cnv` option) containing the alignments used to make the downstream variant calls