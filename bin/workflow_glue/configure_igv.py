"""Create an IGV config file."""

import json
from pathlib import Path
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


# Common variables
REF_EXTENSIONS = [".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz"]
DATA_TYPES_LISTS = {
    "bam": ["bam"],
    "bam_idx": ["bam.bai"],
    "cram": ["cram"],
    "cram_idx": ["cram.crai"],
    "vcf": ["vcf", "vcf.gz"],
    "vcf_idx": ["vcf.gz.tbi", "vcf.gz.csi"],
    "bcf": ["bcf"],
    "bcf_idx": ["bcf.csi"],
    "gtf": ["gtf", "gtf.gz"],
    "gtf_idx": ["gtf.gz.tbi"],
    "gff": ["gff", "gff.gz", "gff3", "gff3.gz"],
    "gff_idx": ["gff.gz.tbi", "gff3.gz.tbi"],
    "bed": ["bed", "bed.gz"],
    "bed_idx": ["bed.gz.tbi"],
    "bedmethyl": ["bedmethyl", "bedmethyl.gz"],
    "bedmethyl_idx": ["bedmethyl.gz.tbi"],
    "ref": REF_EXTENSIONS,
}
DATA_TYPES = {
    ext: ftype for ftype, extlist in DATA_TYPES_LISTS.items() for ext in extlist
}

# Data by idx
DATA_INDEXES_FMT = {
    fmt: f"{fmt}_idx" for fmt, dtype in DATA_TYPES.items() if "_idx" not in dtype
}

# Assign each format to its index
INDEX_PAIRS = {
    "bam": ("bai",),
    "cram": ("crai",),
    "vcf": ("tbi", "csi"),
    "bcf": ("csi",),
    "bed": ("tbi",),
    "bedmethyl": ("tbi",),
    "gff": ("tbi",),
    "gtf": ("tbi",),
}


class TrackBuilder:
    """Class that builds an IGV track."""

    def __init__(self):
        """Initialize properties for interval track."""
        # Reference properties
        self.ref = None
        self.fai = None
        self.gzi = None
        # Samples info
        self.samples = {}
        # Track properties
        self.igv_json = {"reference": {}, "tracks": []}
        self.track_type = {
            "bam": "alignment",
            "cram": "alignment",
            "bcf": "variant",
            "vcf": "variant",
            "bedmethyl": "annotation",
            "bed": "annotation",
            "gtf": "annotation",
            "gff": "annotation",
        }
        # Here we save aliases of file formats that IGV.js
        # wants and that do not match the input file extension.
        self.igv_fmt_alias = {"gff": "gff3"}
        # lookup of extra options for each data type
        self.extra_opts_lookups = {
            "bam": {},
            "cram": {},
            "bcf": {},
            "vcf": {},
            "bed": {},
            "bedmethyl": {},
            "gtf": {},
            "gff": {},
        }

    def add_ref(self, ref=None):
        """Add reference file, unless already defined."""
        if self.ref:
            raise Exception(
                f"Reference genome has already been set to {self.ref}.\n"
                "Only one reference FASTA file is expected."
            )
        else:
            self.ref = ref

    def add_ref_index(self, ref_index=None):
        """Add reference index if valid."""
        basename = Path(self.ref).name
        idx_basename = Path(ref_index).name
        if idx_basename == f"{basename}.fai":
            self.fai = ref_index
        if idx_basename == f"{basename}.gzi" and basename.endswith(".gz"):
            self.gzi = ref_index

    def parse_fnames(self, fofn):
        """Parse list with filenames and return them grouped.

        :param fofn: File with list of file names (one per line)
        """
        tmp_samples = {}
        with open(fofn, "r") as f:
            for line in f:
                # If the line contains the sample name, prepare the data structure
                if "," in line:
                    sample, fname = line.strip().split(",")
                    if sample not in tmp_samples:
                        tmp_samples[sample] = SampleBundle(sample=sample)
                    tmp_samples[sample].append(fname)
                else:
                    # Otherwise, assign everything to NO_SAMPLE
                    # Files will still be displayed, but in no specific order.
                    fname = line.strip()
                    if any(fname.endswith(ext) for ext in REF_EXTENSIONS):
                        self.add_ref(ref=fname)
                    elif fname.endswith(".fai") or fname.endswith(".gzi"):
                        self.add_ref_index(ref_index=fname)
                    else:
                        if "NO_SAMPLE" not in tmp_samples.keys():
                            tmp_samples["NO_SAMPLE"] = SampleBundle(sample="NO_SAMPLE")
                        tmp_samples["NO_SAMPLE"].append(fname)
        # Re-order samples in dict and add them to the list, leaving
        # NO_SAMPLE as last
        sorted_samples = (
            sorted([sample for sample in tmp_samples.keys() if sample != 'NO_SAMPLE'])
        )
        if 'NO_SAMPLE' in tmp_samples.keys():
            sorted_samples += ['NO_SAMPLE']
        for sample in sorted_samples:
            self.samples[sample] = tmp_samples[sample]

    def build_igv_json(self):
        """Ensure there is a reference genome."""
        if not self.ref:
            raise ValueError(
                "No reference file (i.e. file ending in one of "
                f"{REF_EXTENSIONS} was found)."
            )
        # Evaluate that a bgzipped reference has the appropriate index.
        if self.ref.endswith(".gz") and not self.gzi:
            raise ValueError(f"GZI reference index for {self.ref} not found.")

        # Create the base track if there is a reference genome.
        self.igv_json["reference"] = {
            "id": "ref",
            "name": "ref",
            "wholeGenomeView": False,
            "fastaURL": self.ref,
        }
        if self.fai:
            self.igv_json["reference"]["indexURL"] = self.fai
        if self.gzi:
            self.igv_json["reference"]["compressedIndexURL"] = self.gzi

        # Add samples data now
        for sample, bundle in self.samples.items():
            bundle.process_data()
            # Add the bundled data to the tracks
            for fname, index, file_fmt in bundle.data_bundles:
                self.add_track(
                    fname,
                    file_fmt,
                    sample_name=sample if sample != "NO_SAMPLE" else None,
                    index=index,
                    extra_opts=self.extra_opts_lookups[file_fmt],
                )

    def add_track(self, infile, file_fmt, sample_name=None, index=None, extra_opts={}):
        """Add a track to an IGV json.

        This function takes an input file, an optional index file, its
        file format and additional extra options for the track.

        :param infile: input file to create a track for
        :param file_fmt: input file track type
        :param sample_name: Name of the sample to display in the track name
        :param index: index for the input file
        :param extra_opts: dict of extra options for the track
        :return: dict with track options
        """
        # Define track name depending on whether the sample ID is provided
        track_name = Path(infile).name
        if sample_name:
            track_name = f"{sample_name}: {Path(infile).name}"
        track_dict = {
            "name": track_name,
            "type": self.track_type[file_fmt],
            "format": self.igv_fmt_alias.get(file_fmt, file_fmt),
            "url": infile,
        }
        # add the index, if present
        if index:
            track_dict["indexURL"] = index
        track_dict.update(extra_opts)
        self.igv_json["tracks"] += [track_dict]

    def add_locus(self, locus):
        """Add target locus to the json."""
        self.igv_json["locus"] = locus

    def add_extra_opts(
        self,
        extra_alignment_opts=None,
        extra_variant_opts=None,
        extra_interval_opts=None,
    ):
        """Import extra options from json files."""
        if extra_alignment_opts is not None:
            with open(extra_alignment_opts, "r") as f:
                extra_alignment_opts_json = json.load(f)
                for ftype in ["bam", "cram"]:
                    self.extra_opts_lookups[ftype] = extra_alignment_opts_json
        if extra_variant_opts is not None:
            with open(extra_variant_opts, "r") as f:
                extra_variant_opts_json = json.load(f)
                for ftype in ["vcf", "bcf"]:
                    self.extra_opts_lookups[ftype] = extra_variant_opts_json
        if extra_interval_opts is not None:
            with open(extra_interval_opts, "r") as f:
                extra_interval_opts_json = json.load(f)
                for ftype in ["bed", "bedmethyl", "gff", "gtf"]:
                    self.extra_opts_lookups[ftype] = extra_interval_opts_json


class SampleBundle:
    """Sample data class.

    This class stores the data for multiple tracks for a
    single sample, then is used to generate a collection of
    IGV.js tracks.
    """

    def __init__(self, sample):
        """Initialize properties for a sample."""
        self.sample = sample
        self.infiles = []
        self.data_bundles = []

    def append(self, fname):
        """Add a new raw file to the bundle."""
        self.infiles.append(fname)

    def process_data(self):
        """Process input files."""
        fbasenames = [Path(fname).name for fname in self.infiles]
        ftypes = [self.classify_files(bname) for bname in fbasenames]
        self.data_bundles = self.pair_file_with_index(self.infiles, fbasenames, ftypes)

    @staticmethod
    def classify_files(fname):
        """Classify inputs."""
        for extension, ftype in DATA_TYPES.items():
            if fname.endswith(f".{extension}"):
                return ftype

    @staticmethod
    def pair_file_with_index(infiles, fbasenames, ftypes):
        """Clump files with their indexes."""
        # Collect data by group type
        groups = {ftype: {"basenames": [], "paths": []} for ftype in set(ftypes)}
        # Group each file by its type and base name
        for ftype, fbasename, fname in zip(ftypes, fbasenames, infiles):
            groups[ftype]["basenames"] += [fbasename]
            groups[ftype]["paths"] += [fname]

        # Output bundles
        outputs = []
        # Start matching the variant files
        for ftype, itype in DATA_INDEXES_FMT.items():
            # Ignore file formats that are not present in the bundle.
            if ftype not in groups:
                continue
            # Make pairs of files.
            for fbasename, fpath in zip(
                groups[ftype]["basenames"], groups[ftype]["paths"]
            ):
                #  Construct potential index file names based on basename of input files
                idx_basenames = set(
                    [f"{fbasename}.{idx}" for idx in INDEX_PAIRS[ftype]]
                )
                # Find which indexes are available
                if itype in groups.keys():
                    idx_basenames = list(
                        idx_basenames.intersection(set(groups[itype]["basenames"]))
                    )
                    # Get the first index (if there are more than one,
                    # it doesn't matter)
                    bname = idx_basenames[0]
                    idx_fn = groups[itype]["paths"][
                        groups[itype]["basenames"].index(bname)
                    ]
                    outputs.append([fpath, idx_fn, ftype])
                # Otherwise, return only the simple file.
                else:
                    outputs.append([fpath, None, ftype])
        return outputs


def main(args):
    """Run the entry point."""
    logger = get_named_logger("configIGV")

    # parse the FOFN
    igv_builder = TrackBuilder()

    # Add the additional track configurations
    igv_builder.add_extra_opts(
        extra_alignment_opts=args.extra_alignment_opts,
        extra_variant_opts=args.extra_variant_opts,
        extra_interval_opts=args.extra_interval_opts
    )

    # Import files
    igv_builder.parse_fnames(args.fofn)

    # initialise the IGV options dict with the reference options
    igv_builder.build_igv_json()

    # Add locus information
    if args.locus is not None:
        igv_builder.add_locus(args.locus)

    json.dump(igv_builder.igv_json, sys.stdout, indent=4)

    logger.info("Printed IGV config JSON to STDOUT.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("configure_igv")
    parser.add_argument(
        "--fofn",
        required=True,
        help=(
            "File with list of names of reference / XAM / VCF files and indices "
            "(one filename per line)"
        ),
    )
    parser.add_argument(
        "--locus",
        help="Locus string to set initial genomic coordinates to display in IGV",
    )
    parser.add_argument(
        "--extra-alignment-opts",
        help="JSON file with extra options for alignment tracks",
    )
    parser.add_argument(
        "--extra-variant-opts",
        help="JSON file with extra options for variant tracks",
    )
    parser.add_argument(
        "--extra_interval_opts",
        help="JSON file with extra options for interval tracks",
    )
    return parser
