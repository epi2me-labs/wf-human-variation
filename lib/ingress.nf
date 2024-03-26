import java.nio.file.NoSuchFileException

import ArgumentParser

enum InputType {
    SingleFile,
    TopLevelDir,
    DirWithSubDirs,
}

N_OPEN_FILES_LIMIT = 128


/**
* Check if a file ends with one of the target extensions.
*
* @param file: path to the file in question
* @param extensions: list of valid file extensions
* @return: boolean whether the file has one of the provided extensions
*/
def is_target_file(Path file, List extensions) {
    extensions.any { ext -> file.name.endsWith(ext) }
}


/**
 * Take a channel of the shape `[meta, reads, path-to-stats-dir | null]` (or
 * `[meta, [reads, index], path-to-stats-dir | null]` in the case of XAM) and extract the
 * run IDs from the `run_ids` file in the stats directory into the metamap. If the path
 * to the stats dir is `null`, add an empty list.
 *
 * @param ch: input channel of shape `[meta, reads, path-to-stats-dir | null]`
 * @return: channel with a list of run IDs added to the metamap
 */
def add_run_IDs_to_meta(ch) {
    // HashSet for all observed run_ids
    Set<String> ingressed_run_ids = new HashSet<String>()

    // extract run_ids from fastcat stats / bamstats results and add to metadata as well
    // as `ingressed_run_ids`
    ch = ch | map { meta, reads, stats ->
        ArrayList run_ids = []
        if (stats) {
            run_ids = stats.resolve("run_ids").splitText().collect { it.strip() }
            ingressed_run_ids += run_ids
        }
        // `meta + [...]` returns a new map which is handy to avoid any
        // modifying-maps-in-closures weirdness
        // See https://github.com/nextflow-io/nextflow/issues/2660
        [meta + [run_ids: run_ids], reads, stats]
    }
    // put run_ids somewhere global for trivial access later
    // bit grim but decouples ingress metadata from workflow main.nf
    // additionally no need to use CWUtil as we're not overriding any user params
    ch | subscribe(onComplete: {
        params.wf["ingress.run_ids"] = ingressed_run_ids
    })
    return ch
}


/**
 * Take a channel of the shape `[meta, reads, path-to-stats-dir | null]` and do the
 * following:
 *  - For `fastcat`, extract the number of reads from the `n_seqs` file.
 *  - For `bamstats`, extract the number of primary alignments and unmapped reads from
 *    the `bamstats.flagstat.tsv` file.
 * Then, add add these metrics to the meta map. If the path to the stats dir is `null`,
 * set the values to 0 when adding them.
 *
 * @param ch: input channel of shape `[meta, reads, path-to-stats-dir | null]`
 * @return: channel with a list of number of reads added to the metamap
 */
def add_number_of_reads_to_meta(ch, String input_type_format) {
    // extract reads from fastcat stats / bamstats results and add to metadata
    ch = ch | map { meta, reads, stats ->
        // Check that stats directory is present.
        if (stats) {
            if (input_type_format == "fastq") {
                // Stats from fastcat
                Integer n_seqs = stats.resolve("n_seqs").splitText()[0] as Integer
                // `meta + [...]` returns a new map which is handy to avoid any
                // modifying-maps-in-closures weirdness
                // See https://github.com/nextflow-io/nextflow/issues/2660
                [meta + [n_seqs: n_seqs], reads, stats]
            } else {
                // or bamstats
                ArrayList stats_csv = stats.resolve("bamstats.flagstat.tsv").splitCsv(header: true, sep:'\t')
                // get primary alignments and unmapped and sum them
                Integer n_primary = stats_csv["primary"].collect{it as Integer}.sum()
                Integer n_unmapped = stats_csv["unmapped"].collect{it as Integer}.sum()
                // `meta + [...]` returns a new map which is handy to avoid any
                // modifying-maps-in-closures weirdness
                // See https://github.com/nextflow-io/nextflow/issues/2660
                [meta + [n_primary: n_primary, n_unmapped: n_unmapped], reads, stats]
            }
         } else {
            // return defaults if stats is not there
            if (input_type_format == "fastq") {
                [meta + [n_seqs: null], reads, stats]
            } else {
                [meta + [n_primary: null, n_unmapped: null], reads, stats]
            }
        }
    }
    return ch
}

/**
 * Take a map of input arguments, find valid FASTQ inputs, and return a channel
 * with elements of `[metamap, seqs.fastq.gz | null, path-to-fastcat-stats | null]`.
 * The second item is `null` for sample sheet entries without a matching barcode
 * directory. The last item is `null` if `fastcat` was not run (it is only run on
 * directories containing more than one FASTQ file or when `stats: true`).
 *
 * @param arguments: map with arguments containing
 *  - "input": path to either: (i) input FASTQ file, (ii) top-level directory containing
 *     FASTQ files, (iii) directory containing sub-directories which contain FASTQ
 *     files
 *  - "sample": string to name single sample
 *  - "sample_sheet": path to CSV sample sheet
 *  - "analyse_unclassified": boolean whether to keep unclassified reads
 *  - "stats": boolean whether to write the `fastcat` stats
 *  - "fastcat_extra_args": string with extra arguments to pass to `fastcat`
 *  - "required_sample_types": list of required sample types in the sample sheet
 *  - "watch_path": boolean whether to use `watchPath` and run in streaming mode
 * @return: channel of `[Map(alias, barcode, type, ...), Path|null, Path|null]`.
 *  The first element is a map with metadata, the second is the path to the
 *  `.fastq.gz` file with the (potentially concatenated) sequences and the third is
 *  the path to the directory with the `fastcat` statistics. The second element is
 *  `null` for sample sheet entries for which no corresponding barcode directory was
 *  found. The third element is `null` if `fastcat` was not run.
 */
def fastq_ingress(Map arguments)
{
    // check arguments
    Map margs = parse_arguments("fastq_ingress", arguments, ["fastcat_extra_args": ""])

    ArrayList fq_extensions = [".fastq", ".fastq.gz", ".fq", ".fq.gz"]

    // `watch_path` will be handled within `get_valid_inputs()`
    def input = get_valid_inputs(margs, fq_extensions)

    def ch_result
    if (margs.stats) {
        // run fastcat regardless of input type
        ch_result = fastcat(input.files.mix(input.dirs), margs["fastcat_extra_args"])
    } else {
        // run `fastcat` only on directories and rename / compress single files
        ch_result = fastcat(input.dirs, margs["fastcat_extra_args"])
        | mix(
            input.files
            | move_or_compress_fq_file
            | map { meta, path -> [meta, path, null] }
        )
    }
    // add sample sheet entries without barcode dirs to the results channel and extract
    // the run IDs into the metamaps before returning
    ch_result = ch_result.mix(input.missing.map { [*it, null] })
    ch_result_run_IDs = add_run_IDs_to_meta(ch_result)
    // add number of reads after potential filtering under the field n_seqs
    return add_number_of_reads_to_meta(ch_result_run_IDs, "fastq")
}


/**
 * Take a map of input arguments, find valid (u)BAM inputs, and return a channel
 * with elements of `[metamap, reads.bam | null, path-to-bamstats-results | null]`.
 * The second item is `null` for sample sheet entries without a matching barcode
 * directory or samples containing only uBAM files when `keep_unaligned` is `false`.
 * The last item is `null` if `bamstats` was not run (it is only run when `stats:
 * true`).
 *
 * @param arguments: map with arguments containing
 *  - "input": path to either: (i) input (u)BAM file, (ii) top-level directory
 *    containing (u)BAM files, (iii) directory containing sub-directories which contain
 *    (u)BAM files
 *  - "sample": string to name single sample
 *  - "sample_sheet": path to CSV sample sheet
 *  - "analyse_unclassified": boolean whether to keep unclassified reads
 *  - "stats": boolean whether to run `bamstats`
 *  - "keep_unaligned": boolean whether to include uBAM files
 *  - "return_fastq": boolean whether to convert to FASTQ (this will always run
 *    `fastcat`)
 *  - "fastcat_extra_args": string with extra arguments to pass to `fastcat`
 *  - "required_sample_types": list of required sample types in the sample sheet
 *  - "watch_path": boolean whether to use `watchPath` and run in streaming mode
 * @return: channel of `[Map(alias, barcode, type, ...), Path|null, Path|null]`.
 *  The first element is a map with metadata, the second is the path to the
 *  `.bam` file with the (potentially merged) sequences and the third is
 *  the path to the directory with the `bamstats` statistics. The second element is
 *  `null` for sample sheet entries for which no corresponding barcode directory was
 *  found and for samples with only uBAM files when `keep_unaligned: false`. The third
 *  element is `null` if `bamstats` was not run.
 */
def xam_ingress(Map arguments)
{
    // check arguments
    Map margs = parse_arguments(
        "xam_ingress",
        arguments,
        ["keep_unaligned": false, "return_fastq": false, "fastcat_extra_args": ""]
    )

    // we only accept BAM or uBAM for now (i.e. no SAM or CRAM)
    ArrayList xam_extensions = [".bam", ".ubam"]

    def input = get_valid_inputs(margs, xam_extensions)

    // check BAM headers to see if any samples are uBAM
    ch_result = input.dirs
    | map { meta, path -> [meta, get_target_files_in_dir(path, xam_extensions)] }
    | mix(input.files)
    | map{
        // If there is more than one BAM in each folder we ignore
        // the indices. For single BAM we add it as a string to the
        // metadata for later use. If then the BAM returns as position
        // sorted, the index will be used.
        meta, paths -> 
        boolean is_array = paths instanceof ArrayList
        String xai_fn
        // Using `.uri` or `.Uri()` leads to S3 paths to be prefixed with `s3:///`
        // instead of `s3://`, causing the workflow to not find the index file.
        // `.toUriString()` returns the correct path.
        if (!is_array){
            def xai = file(paths.toUriString() + ".bai")
            if (xai.exists()){
                xai_fn = xai.toUriString()
            }
        }
        [meta + [xai_fn: xai_fn], paths]
    }
    | checkBamHeaders
    | map { meta, paths, is_unaligned_env, mixed_headers_env, is_sorted_env ->
        // convert the env. variables from strings ('0' or '1') into bools
        boolean is_unaligned = is_unaligned_env as int as boolean
        boolean mixed_headers = mixed_headers_env as int as boolean
        boolean is_sorted = is_sorted_env as int as boolean
        // throw an error if there was a sample with mixed headers
        if (mixed_headers) {
            error "Found mixed headers in (u)BAM files of sample '${meta.alias}'."
        }
        // add `is_unaligned` to the metamap (note the use of `+` to create a copy of
        // `meta` to avoid modifying every item in the channel;
        // https://github.com/nextflow-io/nextflow/issues/2660)
        [meta + [is_unaligned: is_unaligned, is_sorted: is_sorted], paths]
    }
    | branch { meta, paths ->
        // set `paths` to `null` for uBAM samples if unallowed (they will be added to
        // the results channel in shape of `[meta, null]` at the end of the function
        // (alongside the sample sheet entries without matching barcode dirs)
        if (!margs["keep_unaligned"] && meta["is_unaligned"]){
            paths = null
        }
        // get the number of files (`paths` can be a list, a single path, or `null`)
        int n_files = paths instanceof List ? paths.size() : (paths ? 1 : 0)
        // Preparations finished; we can do the branching now. There will be 5 branches
        // depending on the number of files per sample and whether the reads are already
        // aligned:
        // * no files: no need to do anything
        // * indexed: a single sorted and indexed BAM file. Index will be validated.
        // * to_index: a single sorted, but not indexed, BAM file
        // * to_catsort: `samtools cat` into `samtools sort`
        //  - a single aligned file
        //  - more than one unaligned file
        //  - too many aligned files to safely and quickly merge (`samtools merge` opens
        //    all files at the same time and some machines might have low limits for
        //    open file descriptors)
        // * to_merge: flatMap > sort > group > merge
        //  - between 1 and `N_OPEN_FILES_LIMIT` aligned files
        no_files: n_files == 0
        indexed: \
            n_files == 1 && (meta["is_unaligned"] || meta["is_sorted"]) && meta["xai_fn"]
        to_index: 
            n_files == 1 && (meta["is_unaligned"] || meta["is_sorted"]) && !meta["xai_fn"]
        to_catsort: \
            (n_files == 1) || (n_files > N_OPEN_FILES_LIMIT) || meta["is_unaligned"]
        to_merge: true
    }

    if (margs["return_fastq"]) {
        // only run samtools fastq on samples with at least one file
        ch_to_fastq = ch_result.indexed.mix(
            ch_result.to_index,
            ch_result.to_merge,
            ch_result.to_catsort
        )
    
        // input.missing: sample sheet entries without barcode dirs
        ch_result = input.missing
        | mix(ch_result.no_files)
        | map { [*it, null] }
        | mix(bamToFastq(ch_to_fastq, margs["fastcat_extra_args"]))
        | map{
            meta, path, stats ->
            [meta.findAll { it.key !in ['xai_fn', 'is_sorted'] }, path, stats]
        }
        return add_number_of_reads_to_meta(add_run_IDs_to_meta(ch_result), "fastq")
    }

    // deal with samples with few-enough files for `samtools merge` first
    ch_merged = ch_result.to_merge
    | flatMap { meta, paths -> paths.collect { [meta, it] } }
    | sortBam
    | groupTuple
    | mergeBams

    // now handle samples with too many files for `samtools merge`
    ch_catsorted = ch_result.to_catsort
    | catSortBams

    // Validate the index of the input BAM.
    // If the input BAM index is invalid, regenerate it.
    // First separate the BAM from the null input channels.
    ch_to_validate = ch_result.indexed
    | map{
        meta, paths ->
        bai = paths && meta.xai_fn ? file(meta.xai_fn) : null
        [meta, paths, bai]
    }
    | branch {
        meta, paths, bai ->
        to_validate: paths && bai
        no_op_needed: true
    }

    // Validate non-null files with index
    ch_validated = ch_to_validate.to_validate
    | validateIndex
    | branch {
        meta, bam, bai, has_valid_index_env -> 
        boolean has_valid_index = has_valid_index_env as int as boolean
        // Split if it is a valid index
        valid_idx: has_valid_index
            return [meta, bam, bai]
        invalid_idx: true
            return [meta, bam]
    }

    // Create channel for no_op needed (null channels and valid indexes)
    ch_no_op = ch_validated.valid_idx
    | mix(ch_to_validate.no_op_needed)

    // Re-index sorted-not-indexed BAM file
    ch_indexed = ch_result.to_index
    | mix( ch_validated.invalid_idx )
    | samtools_index

    // Add extra null for the missing index to input.missing
    // as well as the missing metadata.
    // input.missing: sample sheet entries without barcode dirs
    ch_missing = input.missing
    | mix(
        ch_result.no_files,
    )
    | map{
        meta, paths ->
        [meta + [xai_fn: null, is_sorted: false], paths, null]
    }

    // Combine all possible inputs
    ch_result = ch_missing | mix(
        ch_no_op,
        ch_indexed,
        ch_merged,
        ch_catsorted,
    )

    // run `bamstats` if requested
    if (margs["stats"]) {
        // branch and run `bamstats` only on the non-`null` paths
        ch_result = ch_result.branch { meta, path, index ->
            has_reads: path
            is_null: true
        }
        ch_bamstats = bamstats(ch_result.has_reads)

        // the channel comes from xam_ingress also have the BAM index in it.
        // Handle this by placing them in a nested array, maintaining the structure 
        // from fastq_ingress. We do not use variable name as assigning variable
        // name with a tuple not matching (e.g. meta, bam, bai, stats <- [meta, bam, stats] )
        // causes the workflow to crash.
        ch_result = ch_bamstats
        | map{
            it[3] ? [it[0], [it[1], it[2]], it[3]] : it
        }
        | add_run_IDs_to_meta
        | map{
            it.flatten()
        }
        | mix(
            ch_result.is_null.map{it + [null]}
        )
    } else {
        // add `null` instead of path to `bamstats` results dir
        ch_result = ch_result | map { meta, bam, bai -> [meta, bam, bai, null] }
    }

    // Remove metadata that are unnecessary downstream:
    // meta.xai_fn: not needed, as it will be part of the channel as a file
    // meta.is_sorted: if data are aligned, they will also be sorted/indexed
    //
    // The output meta can contain the following flags:
    // [
    //     barcode: always present
    //     type: always present
    //     run_id: always present, but can be empty (i.e. `[]`)
    //     alias: always present
    //     n_primary: always present, but can be `null`
    //     n_unmapped: always present, but can be `null`
    //     is_unaligned: present if there is a (u)BAM file
    // ]
    ch_result = add_number_of_reads_to_meta(
        ch_result
            | map{
                meta, bam, bai, stats ->
                [meta.findAll { it.key !in ['xai_fn', 'is_sorted'] }, [bam, bai], stats]
            }, 
        "xam"
    )
    | map{
        it.flatten()
    }

    return ch_result
}

process bamToFastq {
    label "ingress"
    label "wf_common"
    cpus 4
    memory "2 GB"
    input:
        tuple val(meta), path(bams, stageAs: "input_dir/reads*.bam")
        val extra_args
    output: tuple val(meta), path("seqs.fastq.gz"), path("fastcat_stats")
    script:
    """
    mkdir fastcat_stats

    # Save file as compressed fastq
    fastcat \
        -s ${meta["alias"]} \
        -r >(bgzip -c > fastcat_stats/per-read-stats.tsv.gz) \
        -f fastcat_stats/per-file-stats.tsv \
        --histograms histograms \
        $extra_args \
        <( 
            samtools cat -b <(find input_dir -name 'reads*.bam') | \
            samtools fastq - -n -T '*' -o - -0 - 
        ) \
    | bgzip -c > seqs.fastq.gz

    mv histograms/* fastcat_stats

    # extract the run IDs and number of sequences (n_seqs) from the per-read stats
    csvtk freq -tf runid fastcat_stats/per-read-stats.tsv.gz \
    | csvtk del-header \
    | tee >(cut -f 1 | sort > "fastcat_stats/run_ids") \
    | awk 'BEGIN{n=0}; {n+=\$2}; END{print n}' > "fastcat_stats/n_seqs"
    """
}

process checkBamHeaders {
    label "ingress"
    label "wf_common"
    cpus 1
    memory "2 GB"
    input: tuple val(meta), path("input_dir/reads*.bam")
    output:
        // set the two env variables by `eval`-ing the output of the python script
        // checking the XAM headers
        tuple(
            val(meta),
            path("input_dir/reads*.bam", includeInputs: true),
            env(IS_UNALIGNED),
            env(MIXED_HEADERS),
            env(IS_SORTED),
        )
    script:
    """
    workflow-glue check_bam_headers_in_dir input_dir > env.vars
    source env.vars
    """
}


process validateIndex {
    label "ingress"
    label "wf_common"
    cpus 1
    memory "2 GB"
    input: tuple val(meta), path("reads.bam"), path("reads.bam.bai")
    output:
        // set the two env variables by `eval`-ing the output of the python script
        // checking the XAM headers
        tuple(
            val(meta),
            path("reads.bam", includeInputs: true),
            path("reads.bam.bai", includeInputs: true),
            env(HAS_VALID_INDEX)
        )
    script:
    """
    workflow-glue check_xam_index reads.bam > env.vars
    source env.vars
    """
}


process mergeBams {
    label "ingress"
    label "wf_common"
    cpus 3
    memory "4 GB"
    input: tuple val(meta), path("input_bams/reads*.bam"), path("input_bams/reads*.bam.bai")
    output: tuple val(meta), path("reads.bam"), path("reads.bam.bai")
    script:
    def merge_threads = Math.max(1, task.cpus - 1)
    """
    samtools merge -@ ${merge_threads} \
        -b <(find input_bams -name 'reads*.bam') --write-index -o reads.bam##idx##reads.bam.bai
    """
}


process catSortBams {
    label "ingress"
    label "wf_common"
    cpus 4
    memory "4 GB"
    input: tuple val(meta), path("input_bams/reads*.bam")
    output: tuple val(meta), path("reads.bam"), path("reads.bam.bai")
    script:
    def sort_threads = Math.max(1, task.cpus - 2)
    """
    samtools cat -b <(find input_bams -name 'reads*.bam') \
    | samtools sort - -@ ${sort_threads} --write-index -o reads.bam##idx##reads.bam.bai
    """
}


process sortBam {
    label "ingress"
    label "wf_common"
    cpus 3
    memory "4 GB"
    input: tuple val(meta), path("reads.bam")
    output: tuple val(meta), path("reads.sorted.bam"), path("reads.sorted.bam.bai")
    script:
    def sort_threads = Math.max(1, task.cpus - 1)
    """
    samtools sort --write-index -@ ${sort_threads} reads.bam -o reads.sorted.bam##idx##reads.sorted.bam.bai
    """
}


process bamstats {
    label "ingress"
    label "wf_common"
    cpus 3
    memory "4 GB"
    input:
        tuple val(meta), path("reads.bam"), path("reads.bam.bai")
    output:
        tuple val(meta),
              path("reads.bam"),
              path("reads.bam.bai"),
              path("bamstats_results")
    script:
        def bamstats_threads = Math.max(1, task.cpus - 1)
    """
    mkdir bamstats_results
    bamstats reads.bam -s $meta.alias -u \
        -f bamstats_results/bamstats.flagstat.tsv -t $bamstats_threads \
        --histograms histograms \
    | bgzip > bamstats_results/bamstats.readstats.tsv.gz
    mv histograms/* bamstats_results/

    # extract the run IDs from the per-read stats
    csvtk cut -tf runid bamstats_results/bamstats.readstats.tsv.gz \
    | csvtk del-header | sort | uniq > bamstats_results/run_ids
    """
}
/**
 * Run `watchPath` on the input directory and return a channel of shape [metamap,
 * path-to-target-file]. The meta data is taken from the sample sheet in case one was
 * provided. Otherwise it only contains the `alias` (either `margs["sample"]` or the
 * name of the parent directory of the file).
 *
 * @param input: path to a directory to watch
 * @param margs: Map with parsed input arguments
 * @param extensions: list of valid extensions for the target file type
 * @return: Channel of [metamap, path-to-target-file]
 */
def watch_path(Path input, Map margs, ArrayList extensions) {
    // we have two cases to consider: (i) files being generated in the top-level
    // directory and (ii) files being generated in sub-directories. If we find files of
    // both kinds, throw an error.
    if (input.isFile()) {
        error "Input ($input) must be a directory when using `watch_path`."
    }
    // get existing target files first (look for relevant files in the top-level dir and
    // all sub-dirs)
    def ch_existing_input = Channel.fromPath(input)
    | concat(Channel.fromPath("$input/*", type: 'dir'))
    | map { get_target_files_in_dir(it, extensions) }
    | flatten
    // now get channel with files found by `watchPath`
    def ch_watched = Channel.watchPath("$input/**").until { it.name.startsWith('STOP') }
    // only keep target files
    | filter { is_target_file(it, extensions) }
    // merge the channels
    ch_watched = ch_existing_input | concat(ch_watched)
    // check if input is as expected; start by throwing an error when finding files in
    // top-level dir and sub-directories
    String prev_input_type
    ch_watched
    | map {
        String input_type = (it.parent == input) ? "top-level" : "sub-dir"
        if (prev_input_type && (input_type != prev_input_type)) {
            error "`watchPath` found input files in the top-level directory " +
                "as well as in sub-directories."
        }
        // if file is in a sub-dir, make sure it's not a sub-sub-dir
        if ((input_type == "sub-dir") && (it.parent.parent != input)) {
            error "`watchPath` found an input file more than one level of " +
                "sub-directories deep ('$it')."
        }
        // we also don't want files in the top-level dir when we got a sample sheet
        if ((input_type == "top-level") && margs["sample_sheet"]) {
            error "`watchPath` found input files in top-level directory even though " +
                "a sample sheet was provided ('${margs["sample_sheet"]}')."
        }
        prev_input_type = input_type
    }
    if (margs.sample_sheet) {
        // add metadata from sample sheet (we can't use join here since it does not work
        // with repeated keys; we therefore need to transform the sample sheet data into
        // a map with the barcodes as keys)
        def ch_sample_sheet = get_sample_sheet(file(margs.sample_sheet), margs.required_sample_types)
        | collect
        | map { it.collectEntries { [(it["barcode"]): it] } }
        // now we can use this channel to annotate all files with the corresponding info
        // from the sample sheet
        ch_watched = ch_watched
        | combine(ch_sample_sheet)
        | map { file_path, sample_sheet_map ->
            String barcode = file_path.parent.name
            Map sample_sheet_entry = sample_sheet_map[barcode]
            // throw error if the barcode was not in the sample sheet
            if (!sample_sheet_entry) {
                error "Sub-directory $barcode was not found in the sample sheet."
            }
            [create_metamap(sample_sheet_entry), file_path]
        }
    } else {
        ch_watched = ch_watched
        | map {
            // This file could be in the top-level dir or a sub-dir. In the first case
            // check if a sample name was provided. In the second case, the alias is
            // always the name of the sub-dir.
            String alias
            if (it.parent == input) {
                // top-level dir
                alias = margs["sample"] ?: it.parent.name
            } else {
                // sub-dir
                alias = it.parent.name
            }
            [create_metamap([alias: alias]), it]
        }
    }
    return ch_watched
}


process move_or_compress_fq_file {
    label "ingress"
    label "wf_common"
    cpus 1
    memory "2 GB"
    input:
        // don't stage `input` with a literal because we check the file extension
        tuple val(meta), path(input)
    output:
        tuple val(meta), path("seqs.fastq.gz")
    script:
        String out = "seqs.fastq.gz"
        if (input.name.endsWith('.gz')) {
            // we need to take into account that the file could already be named
            // "seqs.fastq.gz" in which case `mv` would fail
            """
            [ "$input" == "$out" ] || mv "$input" $out
            """
        } else {
            """
            cat "$input" | bgzip -@ $task.cpus > $out
            """
        }
}


process fastcat {
    label "ingress"
    label "wf_common"
    cpus 3
    memory "2 GB"
    input:
        tuple val(meta), path("input")
        val extra_args
    output:
        tuple val(meta),
              path("seqs.fastq.gz"),
              path("fastcat_stats")
    script:
        String out = "seqs.fastq.gz"
        String fastcat_stats_outdir = "fastcat_stats"
        """
        mkdir $fastcat_stats_outdir
        fastcat \
            -s ${meta["alias"]} \
            -r >(bgzip -c > $fastcat_stats_outdir/per-read-stats.tsv.gz) \
            -f $fastcat_stats_outdir/per-file-stats.tsv \
            --histograms histograms \
            $extra_args \
            input \
            | bgzip > $out

        mv histograms/* $fastcat_stats_outdir
        # extract the run IDs and number of sequences (n_seqs) from the per-read stats
        csvtk freq -tf runid $fastcat_stats_outdir/per-read-stats.tsv.gz \
        | csvtk del-header \
        | tee >(cut -f 1 | sort > "$fastcat_stats_outdir/run_ids") \
        | awk 'BEGIN{n=0}; {n+=\$2}; END{print n}' > "$fastcat_stats_outdir/n_seqs"
        """
}


/**
 * Parse input arguments for `fastq_ingress` or `xam_ingress`.
 *
 * @param arguments: map with input arguments (see the corresponding ingress function
 *  for details)
 * @param extra_kwargs: map of extra keyword arguments and their defaults (this allows
 *  the argument-parsing to be tailored to a particular ingress function)
 * @return: map of parsed arguments
 */
Map parse_arguments(String func_name, Map arguments, Map extra_kwargs=[:]) {
    ArrayList required_args = ["input"]
    Map default_kwargs = [
        "sample": null,
        "sample_sheet": null,
        "analyse_unclassified": false,
        "stats": true,
        "required_sample_types": [],
        "watch_path": false
    ]
    ArgumentParser parser = new ArgumentParser(
        args: required_args,
        kwargs: default_kwargs + extra_kwargs,
        name: func_name)
    return parser.parse_args(arguments)
}


/**
 * Find valid inputs based on the target extensions and return a branched channel with
 * branches `missing`, `files` and `dir`, which are of the shape `[metamap, input_path |
 * null]` (with `input_path` pointing to a target file or a directory containing target
 * files, respectively). `missing` contains sample sheet entries for which no
 * corresponding barcodes were found.
 * Unless `watchPath` was requested, the function checks whether the input is a single
 * target file, a top-level directory with target files, or a directory containing
 * sub-directories (usually barcodes) with target files.
 *
 * @param margs: parsed arguments (see `fastq_ingress` and `xam_ingress` for details)
 * @param extensions: list of valid extensions for the target file type
 * @return: branched channel with branches `missing`, `dir`, and `files`
*/
def get_valid_inputs(Map margs, ArrayList extensions){
    log.info "Searching input for $extensions files."
    Path input
    try {
        input = file(margs.input, checkIfExists: true)
    } catch (NoSuchFileException e) {
        error "Input path $margs.input does not exist."
    }
    // declare resulting input channel
    def ch_input
    // run `watchPath` if requested
    if (margs["watch_path"]) {
        ch_input = watch_path(input, margs, extensions)
    } else {
        // check which of the allowed input types (single file, top-lvl dir, dir with
        // sub-dirs) we got
        InputType input_type = determine_input_type(
            input, extensions, margs.analyse_unclassified
        )
        // handle case of `input` being a single file
        if (input_type == InputType.SingleFile) {
            ch_input = Channel.of(
                [create_metamap([alias: margs["sample"] ?: input.simpleName]), input])
        } else if (input_type == InputType.TopLevelDir) {
            // input is a directory containing target files
            ch_input = Channel.of(
                [create_metamap([alias: margs["sample"] ?: input.baseName]), input])
        } else {
            // input is a directory with sub-directories (e.g. barcodes) containing
            // target files --> find these sub-directories
            ArrayList sub_dirs_with_target_files = file(
                input.resolve('*'), type: "dir"
            ).findAll { get_target_files_in_dir(it, extensions) }
            // remove directories called 'unclassified' unless otherwise specified
            if (!margs.analyse_unclassified) {
                sub_dirs_with_target_files = sub_dirs_with_target_files.findAll {
                    it.baseName != "unclassified"
                }
            }
            // filter based on sample sheet in case one was provided
            if (margs.sample_sheet) {
                // get channel of entries in the sample sheet
                def ch_sample_sheet = get_sample_sheet(
                    file(margs.sample_sheet), margs.required_sample_types
                )
                // get the union of both channels (missing values will be replaced with
                // `null`)
                def ch_union = Channel.fromPath(sub_dirs_with_target_files).map {
                    [it.baseName, it]
                }.join(ch_sample_sheet.map{[it.barcode, it]}, remainder: true)
                // after joining the channels, there are three possible cases:
                // (i) valid input path and sample sheet entry are both present
                // (ii) there is a sample sheet entry but no corresponding input dir
                //      --> we'll emit `[metamap-from-sample-sheet-entry, null]`
                // (iii) there is a valid path, but the sample sheet entry is missing
                //      --> drop this entry and print a warning to the log
                ch_input = ch_union.map {barcode, path, sample_sheet_entry ->
                    if (sample_sheet_entry) {
                        [create_metamap(sample_sheet_entry), path]
                    } else {
                        log.warn "Input directory '$barcode' was found, but sample " +
                            "sheet '$margs.sample_sheet' has no such entry."
                    }
                }
            } else {
                // no sample sheet --> simply emit the sub-dirs with the target files
                ch_input = Channel.fromPath(sub_dirs_with_target_files).map {
                    [create_metamap([alias: it.baseName, barcode: it.baseName]), it]
                }
            }
        }
    }
    // finally, we "unwrap" directories containing only a single file and then split the
    // results channel into the three different output types (sample sheet entries
    // without corresponding barcodes -- i.e. with `path == null`, single files, and
    // dirs with multiple files)
    def ch_branched_results = ch_input.map { meta, path ->
        if (path && path.isDirectory()) {
            List fq_files = get_target_files_in_dir(path, extensions)
            if (fq_files.size() == 1) {
                path = fq_files[0]
            }
        }
        [meta, path]
    } .branch { meta, path ->
        missing: !path
        files: path.isFile()
        dirs: path.isDirectory()
    }
    return ch_branched_results
}

/**
 * Determine which of the allowed categories (single file, top-level directory, or
 * directory with sub-directory) an input path belongs to.
 *
 * @param margs: parsed arguments (see `fastq_ingress()` or `xam_ingress()` for details)
 * @param extensions: list of valid extensions for the target file type
 * @return: input type represented as an instance of the `InputType` enum
 */
InputType determine_input_type(
    Path input, ArrayList extensions, boolean analyse_unclassified
) {
    if (input.isFile()) {
        if (!is_target_file(input, extensions)) {
            error "Input file is not of required file type."
        }
        return InputType.SingleFile
    } else if (!input.isDirectory()){
        error "Input $input appears to be neither a file nor a directory."
    }
    // `input` is a directory --> we accept two cases: (i) a top-level directory with
    // target files and no sub-directories or (ii) a directory with one layer of
    // sub-directories containing target files. First, check if the directory contains
    // target files and find potential sub-directories (and sub-dirs with target files;
    // note that these lists can be empty)
    boolean dir_has_target_files = get_target_files_in_dir(input, extensions)
    ArrayList sub_dirs = file(input.resolve('*'), type: "dir")
    ArrayList sub_dirs_with_target_files = sub_dirs.findAll {
        get_target_files_in_dir(it, extensions)
    }.findAll { it.baseName != "unclassified" || analyse_unclassified }

    // define string to re-use in error messages below
    String target_files_str = \
        "target files (ending in ${extensions.collect{'\'' + it + '\''}.join(' / ')})"

    // check for target files in the top-level dir; if there are any, make sure there
    // are no sub-directories containing target files
    if (dir_has_target_files) {
        if (sub_dirs_with_target_files) {
            error "Input directory '$input' cannot contain $target_files_str " +
                "and also sub-directories with such files."
        }
        return InputType.TopLevelDir
    }

    // no target files in the top-level dir --> make sure there were sub-dirs with
    // target files
    if (!sub_dirs_with_target_files) {
        error "Input directory '$input' must contain either $target_files_str " +
            "or sub-directories containing such files (no more than one layer deep)."
    }
    // we don't allow sub-sub-directories with target files
    if (sub_dirs.any {
        ArrayList subsubdirs = file(it.resolve('*'), type: "dir")
        subsubdirs.any { get_target_files_in_dir(it, extensions) }
    }) {
        error "Input directory '$input' cannot contain more " +
            "than one level of sub-directories with $target_files_str."
    }
    return InputType.DirWithSubDirs
}


/**
 * Create a map that contains at least these keys: `[alias, barcode, type]`.
 * `alias` is required, `barcode` and `type` are filled with default values if
 * missing. Additional entries are allowed.
 *
 * @param arguments: map with input parameters; must contain `alias`
 * @return: map(alias, barcode, type, ...)
 */
Map create_metamap(Map arguments) {
    ArgumentParser parser = new ArgumentParser(
        args: ["alias"],
        kwargs: [
            "barcode": null,
            "type": "test_sample",
            "run_ids": [],
        ],
        name: "create_metamap",
    )
    def metamap = parser.parse_known_args(arguments)
    metamap['alias'] = metamap['alias'].replaceAll(" ","_")
    return metamap
}


/**
 * Get the target files in the directory (non-recursive).
 *
 * @param dir: path to the target directory
 * @param extensions: list of valid extensions for the target file type
 * @return: list of found target files
 */
ArrayList get_target_files_in_dir(Path dir, ArrayList extensions) {
    file(dir.resolve("*")).findAll { is_target_file(it, extensions) }
}


/**
 * Check the sample sheet and return a channel with its rows if it is valid.
 *
 * @param sample_sheet: path to the sample sheet CSV
 * @return: channel of maps (with values in sample sheet header as keys)
 */
def get_sample_sheet(Path sample_sheet, ArrayList required_sample_types) {
    // If `validate_sample_sheet` does not return an error message, we can assume that
    // the sample sheet is valid and parse it. However, because of Nextflow's
    // asynchronous magic, we might emit values from `.splitCSV()` before the
    // error-checking closure finishes. This is no big deal, but undesired nonetheless
    // as the error message might be overwritten by the traces of new nextflow processes
    // in STDOUT. Thus, we use the somewhat clunky construct with `concat` and `last`
    // below. This lets the CSV channel only start to emit once the error checking is
    // done.
    ch_err = validate_sample_sheet(sample_sheet, required_sample_types).map { stdoutput, sample_sheet_file ->
        // check if there was an error message
        if (stdoutput) error "Invalid sample sheet: ${stdoutput}."
        stdoutput
    }
    // concat the channel holding the path to the sample sheet to `ch_err` and call
    // `.last()` to make sure that the error-checking closure above executes before
    // emitting values from the CSV
    return ch_err.concat(Channel.fromPath(sample_sheet)).last().splitCsv(
        header: true, quote: '"'
    )
}


/**
 * Python script for validating a sample sheet. The script will write messages
 * to STDOUT if the sample sheet is invalid. In case there are no issues, no
 * message is emitted. The sample sheet will be published to the output dir.
 *
 * @param: path to sample sheet CSV
 * @param: list of required sample types (optional)
 * @return: string (optional)
 */
process validate_sample_sheet {
    publishDir params.out_dir, mode: 'copy', overwrite: true
    cpus 1
    label "ingress"
    label "wf_common"
    memory "2 GB"
    input:
        path "sample_sheet.csv"
        val required_sample_types
    output: 
        tuple stdout, path("sample_sheet.csv")
    script:
    String req_types_arg = required_sample_types ? "--required_sample_types "+required_sample_types.join(" ") : ""
    """
    workflow-glue check_sample_sheet sample_sheet.csv $req_types_arg
    """
}

// Generate an index for an input XAM file
process samtools_index {
    cpus 4
    label "ingress"
    label "wf_common"
    memory 4.GB
    input:
        tuple val(meta), path("reads.bam")
    output:
        tuple val(meta), path("reads.bam"), path("reads.bam.bai")
    script:
    """
    samtools index -@ $task.cpus reads.bam
    """
}
