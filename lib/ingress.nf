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
 * Take a channel of the shape `[meta, reads, path-to-stats-dir | null]` and extract the
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
    Map margs = parse_arguments(arguments, ["fastcat_extra_args": ""])

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
    return add_run_IDs_to_meta(ch_result)
}


/**
 * Take a map of input arguments, find valid (u)BAM inputs, and return a channel
 * with elements of `[metamap, reads.bam | null, path-to-bamstats-results | null]`.
 * The second item is `null` for sample sheet entries without a matching barcode
 * directory or samples containing only uBAM files when `keep_unaligned` is `false`.
 * The last item is `null` if `bamstats` was not run (it is only run when `stats: true`).
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
    Map margs = parse_arguments(arguments, ["keep_unaligned": false])

    // we only accept BAM or uBAM for now (i.e. no SAM or CRAM)
    ArrayList xam_extensions = [".bam", ".ubam"]

    def input = get_valid_inputs(margs, xam_extensions)

    ch_result = input.dirs
    | map { meta, path -> [meta, get_target_files_in_dir(path, xam_extensions)] }
    | mix(input.files)

    ch_is_unaligned = ch_result
    | checkBamHeaders
    | map { meta, is_unaligned_env, mixed_headers_env ->
        // convert the env. variables from strings ('0' or '1') into bools
        boolean is_unaligned = is_unaligned_env as int as boolean
        boolean mixed_headers = mixed_headers_env as int as boolean
        // throw an error if there was a sample with mixed headers
        if (mixed_headers) {
            error "Found mixed headers in (u)BAM files of sample '${meta.alias}'."
        }
        [meta, is_unaligned]
    }

    ch_result = ch_result | join(ch_is_unaligned)
    // add `is_unaligned` to the metamap (note the use of `+` to create a copy of `meta`
    // to avoid modifying every item in the channel;
    // https://github.com/nextflow-io/nextflow/issues/2660)
    | map { meta, paths, is_unaligned -> [meta + [is_unaligned: is_unaligned], paths] }
    | branch { meta, paths ->
        // set `paths` to `null` for uBAM samples if unallowed (they will be added to
        // the results channel in shape of `[meta, null]` at the end of the function
        // (alongside the sample sheet entries without matching barcode dirs)
        if (!margs["keep_unaligned"] && meta["is_unaligned"]){
            paths = null
        }
        // get the number of files (`paths` can be a list, a single path, or `null`)
        int n_files = paths instanceof List ? paths.size() : (paths ? 1 : 0)
        // Preparations finished; we can do the branching now. There will be 3 branches
        // depending on the number of files per sample and whether the reads are already
        // aligned:
        // * no_op_needed: no need to do anything; just add to the final results channel
        //   downstream
        //   - no files
        //   - a single unaligned file
        // * to_catsort: `samtools cat` into `samtools sort`
        //  - a single aligned file
        //  - more than one unaligned file
        //  - too many aligned files to safely and quickly merge (`samtools merge` opens
        //    all files at the same time and some machines might have low limits for
        //    open file descriptors)
        // * to_merge: flatMap > sort > group > merge
        //  - between 1 and `N_OPEN_FILES_LIMIT` aligned files
        no_op_needed: (n_files == 0) || (n_files == 1 && meta["is_unaligned"])
        to_catsort: \
            (n_files == 1) || (n_files > N_OPEN_FILES_LIMIT) || meta["is_unaligned"]
        to_merge: true
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

    ch_result = input.missing
    | mix(
        ch_result.no_op_needed,
        ch_merged,
        ch_catsorted,
    )

    // run `bamstats` if requested
    if (margs["stats"]) {
        // branch and run `bamstats` only on the non-`null` paths
        ch_result = ch_result.branch { meta, path ->
            has_reads: path
            is_null: true
        }
        ch_bamstats = bamstats(ch_result.has_reads)
        ch_result = add_run_IDs_to_meta(ch_bamstats) | mix(ch_result.is_null)
    } else {
        // add `null` instead of path to `bamstats` results dir
        ch_result = ch_result | map { meta, bam -> [meta, bam, null] }
    }
    return ch_result
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
        tuple val(meta), env(IS_UNALIGNED), env(MIXED_HEADERS)
    script:
    """
    workflow-glue check_bam_headers_in_dir input_dir > env.vars
    source env.vars
    """
}


process mergeBams {
    label "ingress"
    label "wf_common"
    cpus 3
    memory "4 GB"
    input: tuple val(meta), path("input_bams/reads*.bam")
    output: tuple val(meta), path("reads.bam")
    shell:
    """
    samtools merge -@ ${task.cpus - 1} \
        -b <(find input_bams -name 'reads*.bam') -o reads.bam
    """
}


process catSortBams {
    label "ingress"
    label "wf_common"
    cpus 4
    memory "4 GB"
    input: tuple val(meta), path("input_bams/reads*.bam")
    output: tuple val(meta), path("reads.bam")
    script:
    """
    samtools cat -b <(find input_bams -name 'reads*.bam') \
    | samtools sort - -@ ${task.cpus - 2} -o reads.bam
    """
}


process sortBam {
    label "ingress"
    label "wf_common"
    cpus 3
    memory "4 GB"
    input: tuple val(meta), path("reads.bam")
    output: tuple val(meta), path("reads.sorted.bam")
    script:
    """
    samtools sort -@ ${task.cpus - 1} reads.bam -o reads.sorted.bam
    """
}


process bamstats {
    label "ingress"
    label "wf_common"
    cpus 3
    memory "4 GB"
    input:
        tuple val(meta), path("reads.bam")
    output:
        tuple val(meta), path("reads.bam"), path("bamstats_results")
    script:
        def bamstats_threads = Math.max(1, task.cpus - 1)
    """
    mkdir bamstats_results
    bamstats reads.bam -s $meta.alias -u \
        -f bamstats_results/bamstats.flagstat.tsv -t $bamstats_threads \
    | bgzip > bamstats_results/bamstats.readstats.tsv.gz

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
        tuple val(meta), path("seqs.fastq.gz"), path("fastcat_stats")
    script:
        String out = "seqs.fastq.gz"
        String fastcat_stats_outdir = "fastcat_stats"
        """
        mkdir $fastcat_stats_outdir
        fastcat \
            -s ${meta["alias"]} \
            -r >(bgzip -c > $fastcat_stats_outdir/per-read-stats.tsv.gz) \
            -f $fastcat_stats_outdir/per-file-stats.tsv \
            $extra_args \
            input \
            | bgzip > $out

        # extract the run IDs from the per-read stats
        csvtk cut -tf runid $fastcat_stats_outdir/per-read-stats.tsv.gz \
        | csvtk del-header | sort | uniq > $fastcat_stats_outdir/run_ids
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
Map parse_arguments(Map arguments, Map extra_kwargs=[:]) {
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
        name: "fastq_ingress")
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
    ch_err = validate_sample_sheet(sample_sheet, required_sample_types).map {
        // check if there was an error message
        if (it) error "Invalid sample sheet: ${it}."
        it
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
 * message is emitted.
 *
 * @param: path to sample sheet CSV
 * @param: list of required sample types (optional)
 * @return: string (optional)
 */
process validate_sample_sheet {
    cpus 1
    label "ingress"
    label "wf_common"
    memory "2 GB"
    input:
        path "sample_sheet.csv"
        val required_sample_types
    output: stdout
    script:
    String req_types_arg = required_sample_types ? "--required_sample_types "+required_sample_types.join(" ") : ""
    """
    workflow-glue check_sample_sheet sample_sheet.csv $req_types_arg
    """
}
