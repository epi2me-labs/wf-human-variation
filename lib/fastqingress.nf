import java.nio.file.NoSuchFileException

import ArgumentParser

EXTENSIONS = ["fastq", "fastq.gz", "fq", "fq.gz"]

/**
 * Take a map of input arguments, find valid inputs, and return a channel
 * with elements of `[metamap, seqs.fastq.gz, path-to-fastcat-stats]`.
 * The last item is `null` if `fastcat` was not run. It is only run on directories
 * containing more than one FASTQ file or when `fastcat_stats: true`.
 *
 * @param arguments: map with arguments containing
 *  - "input": path to either: (i) input FASTQ file, (ii) top-level directory containing
 *     FASTQ files, (iii) directory containing sub-directories which contain FASTQ
 *     files
 *  - "sample": string to name single sample
 *  - "sample_sheet": path to CSV sample sheet
 *  - "analyse_unclassified": boolean whether to keep unclassified reads
 *  - "fastcat_stats": boolean whether to write the `fastcat` stats
 * @return Channel of `[Map(alias, barcode, type, ...), Path, Path|null]`.
 *  The first element is a map with metadata, the second is the path to the
 *  `.fastq.gz` file with the (potentially concatenated) sequences and the third is
 *  the path to the directory with the fastcat statistics (or `null` if `fastcat`
 *  wasn't run).
 */
def fastq_ingress(Map arguments)
{
    // check arguments
    Map margs = parse_arguments(arguments)
    // define the channel for holding the inputs [metamap, input_path]. It will be
    // either filled by `watchPath` (only emitting files) or by the data of the three
    // input types (single file or dir with fastq or subdirs with fastq).
    def ch_input
    // handle `watchPath` case
    if (margs["watch_path"]) {
        ch_input = watch_path(margs)
    } else {
        // create a channel with the inputs (single file / dir with fastq / subdirs
        // with fastq)
        ch_input = get_valid_inputs(margs)
    }
    // `ch_input` might contain elements of `[metamap, null]` if there were entries in
    // the sample sheet for which no FASTQ files were found. We put these into an extra
    // channel and combine with the result channel before returning.
    ch_input = ch_input.branch { meta, path ->
        reads_found: path as boolean
        no_reads_found: true
    }
    def ch_result
    if (margs.fastcat_stats) {
        // run fastcat regardless of input type
        ch_result = fastcat(ch_input.reads_found, margs["fastcat_extra_args"]).map {
            meta, reads, stats ->
                // extract run_ids parsed by fastcat into metadata
                ArrayList run_ids = stats.resolve("run_ids").splitText().collect {
                    it.strip()
                }
                // `meta + [...]` returns a new map which is handy to avoid any
                // modifying-maps-in-closures weirdness
                // See https://github.com/nextflow-io/nextflow/issues/2660
                [meta + [run_ids: run_ids], reads, stats]
        }
    } else {
        // the fastcat stats were not requested --> run fastcat only on directories with
        // more than one FASTQ file (and not on single files or directories with a
        // single file)
        def ch_branched = ch_input.reads_found.map {meta, path ->
            // find directories with only a single FASTQ file and "unwrap" the file
            if (path.isDirectory()) {
                List fq_files = get_fq_files_in_dir(path)
                if (fq_files.size() == 1) {
                    path = fq_files[0]
                }
            }
            [meta, path]
        } .branch { meta, path ->
            // now there can only be two cases:
            // (i) single FASTQ file (pass to `move_or_compress` later)
            // (ii) dir with multiple fastq files (pass to `fastcat` later)
            single_file: path.isFile()
            dir_with_fastq_files: true
        }
        // call the respective processes on both branches and return
        ch_result = fastcat(
            ch_branched.dir_with_fastq_files, margs["fastcat_extra_args"]
        ).concat(
            ch_branched.single_file | move_or_compress | map {
                meta, path -> [meta, path, null]
            }
        )
    }
    return ch_result.concat(ch_input.no_reads_found.map { [*it, null] })
}


/**
 * Run `watchPath` on the input directory and return a channel [metamap, path-to-fastq].
 * The meta data is taken from the sample sheet in case one was provided. Otherwise it
 * only contains the `alias` (either `margs["sample"]` or the name of the parent
 * directory of the file).
 *
 * @param margs: map with parsed input arguments
 * @return: Channel of [metamap, path-to-fastq]
 */
def watch_path(Map margs) {
    // we have two cases to consider: (i) files being generated in the top-level
    // directory and (ii) files being generated in sub-directories. If we find files of
    // both kinds, throw an error.
    Path input
    try {
        input = file(margs.input, checkIfExists: true)
    } catch (NoSuchFileException e) {
        error "Input path $margs.input does not exist."
    }
    if (input.isFile()) {
        error "Input ($input) must be a directory when using `watch_path`."
    }
    // get existing FASTQ files first (look for relevant files in the top-level dir and
    // all sub-dirs)
    def ch_existing_input = Channel.fromPath(input)
    | concat(Channel.fromPath("$input/*", type: 'dir'))
    | map { get_fq_files_in_dir(it) }
    | flatten
    // now get channel with files found by `watchPath`
    def ch_watched = Channel.watchPath("$input/**").until { it.name.startsWith('STOP') }
    // only keep FASTQ files
    | filter {
        for (ext in EXTENSIONS) {
            if (it.name.endsWith(ext)) return true
        }
        return false
    }
    // merge the channels
    ch_watched = ch_existing_input | concat(ch_watched)
    // check if input is as expected; start by throwing an error when finding files in
    // top-level dir and sub-directories
    String prev_input_type
    ch_watched
    | map {
        String input_type = (it.parent == input) ? "top-level" : "sub-dir"
        if (prev_input_type && (input_type != prev_input_type)) {
            error "`watchPath` found FASTQ files in the top-level directory " +
                "as well as in sub-directories."
        }
        // if file is in a sub-dir, make sure it's not a sub-sub-dir
        if ((input_type == "sub-dir") && (it.parent.parent != input)) {
            error "`watchPath` found a FASTQ file more than one level of " +
                "sub-directories deep ('$it')."
        }
        // we also don't want files in the top-level dir when we got a sample sheet
        if ((input_type == "top-level") && margs["sample_sheet"]) {
            error "`watchPath` found files in top-level directory even though a " +
                "sample sheet was provided ('${margs["sample_sheet"]}')."
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
            Map meta = sample_sheet_map[barcode]
            // throw error if the barcode was not in the sample sheet
            if (!meta) {
                error "Sub-directory $barcode was not found in the sample sheet."
            }
            [meta, file_path]
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


process move_or_compress {
    label "fastq_ingress"
    label "wf_common"
    cpus 1
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
    label "fastq_ingress"
    label "wf_common"
    cpus 3
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
            -r $fastcat_stats_outdir/per-read-stats.tsv \
            -f $fastcat_stats_outdir/per-file-stats.tsv \
            $extra_args \
            input \
            | bgzip -@ $task.cpus > $out
        csvtk cut -tf runid $fastcat_stats_outdir/per-read-stats.tsv | csvtk del-header | sort | uniq > $fastcat_stats_outdir/run_ids
        """
}


/**
 * Parse input arguments for `fastq_ingress`.
 *
 * @param arguments: map with input arguments (see `fastq_ingress` for details)
 * @return: map of parsed arguments
 */
Map parse_arguments(Map arguments) {
    ArgumentParser parser = new ArgumentParser(
        args:["input"],
        kwargs:["sample": null,
                "sample_sheet": null,
                "analyse_unclassified": false,
                "fastcat_stats": false,
                "fastcat_extra_args": "",
                "required_sample_types": [],
                "watch_path": false],
        name: "fastq_ingress")
    return parser.parse_args(arguments)
}


/**
 * Find valid inputs based on the input type.
 *
 * @param margs: parsed arguments (see `fastq_ingress` for details)
 * @return: channel of `[metamap, input-path]`; `input-path` can be the path to
 *  a single FASTQ file or to a directory containing FASTQ files
 */
def get_valid_inputs(Map margs){
    log.info "Checking fastq input."
    Path input
    try {
        input = file(margs.input, checkIfExists: true)
    } catch (NoSuchFileException e) {
        error "Input path $margs.input does not exist."
    }
    // declare resulting input channel and other variables needed in the outer scope
    def ch_input
    ArrayList sub_dirs_with_fastq_files
    // handle case of `input` being a single file
    if (input.isFile()) {
        // the `fastcat` process can deal with directories or single file inputs
        ch_input = Channel.of(
            [create_metamap([alias: margs["sample"] ?: input.simpleName]), input])
    } else if (input.isDirectory()) {
        // input is a directory --> we accept two cases: (i) a top-level directory with
        // fastq files and no sub-directories or (ii) a directory with one layer of
        // sub-directories containing fastq files
        boolean dir_has_fastq_files = get_fq_files_in_dir(input)
        // find potential sub-directories (and sub-dirs with FASTQ files; note that
        // these lists can be empty)
        ArrayList sub_dirs = file(input.resolve('*'), type: "dir")
        sub_dirs_with_fastq_files = sub_dirs.findAll { get_fq_files_in_dir(it) }
        // deal with first case (top-lvl dir with FASTQ files and no sub-directories
        // containing FASTQ files)
        if (dir_has_fastq_files) {
            if (sub_dirs_with_fastq_files) {
                error "Input directory '$input' cannot contain FASTQ " +
                    "files and sub-directories with FASTQ files."
            }
            ch_input = Channel.of(
                [create_metamap([alias: margs["sample"] ?: input.baseName]), input])
        } else {
            // deal with the second case (sub-directories with fastq data) --> first
            // check whether we actually found sub-directories
            if (!sub_dirs_with_fastq_files) {
                error "Input directory '$input' must contain either FASTQ files " +
                    "or sub-directories containing FASTQ files."
            }
            // make sure that there are no sub-sub-directories with FASTQ files and that
            // the sub-directories actually contain fastq files)
            if (sub_dirs.any {
                ArrayList subsubdirs = file(it.resolve('*'), type: "dir")
                subsubdirs.any { get_fq_files_in_dir(it) }
            }) {
                error "Input directory '$input' cannot contain more " +
                    "than one level of sub-directories with FASTQ files."
            }
            // remove directories called 'unclassified' unless otherwise specified
            if (!margs.analyse_unclassified) {
                sub_dirs_with_fastq_files = sub_dirs_with_fastq_files.findAll {
                    it.baseName != "unclassified"
                }
            }
            // filter based on sample sheet in case one was provided
            if (margs.sample_sheet) {
                // get channel of entries in the sample sheet
                def ch_sample_sheet = get_sample_sheet(file(margs.sample_sheet), margs.required_sample_types)
                // get the union of both channels (missing values will be replaced with
                // `null`)
                def ch_union = Channel.fromPath(sub_dirs_with_fastq_files).map {
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
                ch_input = Channel.fromPath(sub_dirs_with_fastq_files).map {
                    [create_metamap([alias: it.baseName, barcode: it.baseName]), it]
                }
            }
        }
    } else {
        error "Input $input appears to be neither a file nor a directory."
    }
    // a sample sheet only makes sense in the case of a directory with
    // sub-directories
    if (margs.sample_sheet && !sub_dirs_with_fastq_files) {
        error "Sample sheet was provided, but input does not contain " +
            "sub-directories with FASTQ files."
    }
    return ch_input
}


/**
 * Create a map that contains at least these keys: `[alias, barcode, type]`.
 * `alias` is required, `barcode` and `type` are filled with default values if
 * missing. Additional entries are allowed.
 *
 * @param kwargs: map with input parameters; must contain `alias`
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
 * Get the fastq files in the directory (non-recursive).
 *
 * @param dir: path to the target directory
 * @return: list of found fastq files
 */
ArrayList get_fq_files_in_dir(Path dir) {
    return EXTENSIONS.collect { file(dir.resolve("*.$it"), type: "file") } .flatten()
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
    label "fastq_ingress"
    label "wf_common"
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

