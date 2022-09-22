# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.2.2]
### Fixed
* "No such property" when using the `minimap2_ubam` alignment step
* Slow performance on `minimap2_ubam` step when providing CRAM as `--ubam`
* Slow performance on `snp:readStats` process
### Removed
* "Missing reference index" warning was unnecessary

## [v0.2.0]
### Added
* An experimental methylation subworkflow has been integrated, using [`modbam2bed`](https://github.com/epi2me-labs/modbam2bed) to aggregate modified base counts (input BAM should have `MM` and `ML` tags), enable this with `--methyl`
### Changed
* Workflow experimentally supports CRAM as input, uBAM input uses CRAM for intermediate files
* Reference FAI is now created if it does not exist, rather than raising an error
* `--sniffles_args` may be used to provide custom arguments to the `sniffles2` process
* Output files are uniformly prefixed with `--sample_name`
### Fixed
* Existence of Clair3 model directory is checked before starting workflow
* `--GVCF` and `--include_all_ctgs` are correctly typed as booleans
    * `--GVCF` now outputs GVCF files to the output directory as intended
    * `--include_all_ctgs` no longer needs to be set as `--include_all_ctgs y`

## [v0.1.1]
### Added
* `--ubam_bam2fq_threads` and `--ubam_sort_threads` allow finer control over the resourcing of the alignment step
    * The number of CPU required for `minimap2_ubam` is sum(`ubam_bam2fq_threads`, `ubam_map_threads`, `ubam_sort_threads`)
### Changed
* `--ubam_threads` is now `--ubam_map_threads` and only sets the threads for the mapping specifically
* `--ref` is now required by the schema, preventing obscure errors when it is not provided
* Print helpful warning if neither `--snp` or `--sv` have been specified
* Fastqingress metadata map
* Disable dag creation by default to suppress graphviz warning
### Fixed
* Tandem repeat BED is now correctly set by `--tr_bed`
* Update `fastcat` dependency to 0.4.11 to catch cases where input FASTQ cannot be read
* Sanitize fastq intermittent null object error
### Note
- Switched to "new flavour" CI, meaning that containers are released from workflows independently

## [v0.1.0]
# Added
- Ported wf-human-snp (v0.3.2) to modules/wf-human-snp
- Ported wf-human-sv (v0.1.0) to modules/wf-human-sv

## [v0.0.0]
- Initialised wf-human-variation from wf-template #195cab5
