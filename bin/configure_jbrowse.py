#!/usr/bin/env python
"""configure_jbrowse."""
import json
import os
import sys

import pysam


class JbConfig:
    """JBrowse2 basic configuration container."""

    @staticmethod
    def construct_defaultsession_track(
        track_name,
        track_type,
        display_type,
        track_id=None,
    ):
        """Create a track to append to the defaultSession track list."""
        if not track_id:
            track_id = track_name
        return {
            "id": "%s-track" % track_name,
            "type": track_type,
            "configuration": track_id,
            "displays": [
                {
                    "id": "%s-track-display" % track_id,
                    "type": display_type,
                    "configuration": "%s-%s" % (track_id, display_type),
                }
            ]
        }

    def __init__(self):
        """Construct an empty LinearGenomeView configuration container."""
        self.assembly_names = []
        self.assembly = None  # single assembly support for now
        self.alignment = None
        self.tracks = {}  # structs indexed by track name

        self.default_sequence_name = None
        self.default_assembly_name = None
        self.default_view = {
            "id": "view",
            "type": "LinearGenomeView",
            "displayedRegions": [],
            "tracks": [],
            "location": None,  # outside of JB2 spec but used by LGV API
        }
        self.default_session = {
            "name": "epi2me-labs",
            "views": [self.default_view],
            "minimized": True,  # hide the track config drawer
        }

    def add_assembly(
        self, sequence_fp, index_fp, compressed_index_fp=None, name="ref",
        display_width=1000, actual_reference_fp=None
    ):
        """Add an assembly configuration."""
        if not actual_reference_fp:
            actual_reference_fp = sequence_fp
        the_ref = pysam.FastaFile(actual_reference_fp)
        self.default_sequence_name = the_ref.references[0]
        self.default_assembly_name = name

        assembly_type = "IndexedFastaAdapter"
        if compressed_index_fp:
            assembly_type = "BgzipFastaAdapter"

        assembly = {
            "name": name,
            "aliases": [],
            "sequence": {
                "type": "ReferenceSequenceTrack",
                "trackId": name,
                "adapter": {
                    "type": assembly_type,
                    "fastaLocation": {
                        "uri": sequence_fp,
                        "locationType": "UriLocation",
                    },
                    "faiLocation": {
                        "uri": index_fp,
                        "locationType": "UriLocation",
                    },
                },
                "rendering": {
                    "type": "DivSequenceRenderer"
                },
            }
        }
        if compressed_index_fp:
            assembly["sequence"]["adapter"]["gziLocation"] = {
                "uri": compressed_index_fp,
                "locationType": "UriLocation",
            }

        self.assembly = assembly
        self.assembly_names = [name]

        # configure default view nonsense
        # first add the default contig to the default display region
        self.default_view["displayedRegions"] = [
            {
                "assemblyName": self.default_assembly_name,
                "refName": self.default_sequence_name,
                "start": 0,
                "end": the_ref.lengths[0],
                "reversed": False,
            }
        ]
        # now add the reference as a track
        self.default_view["tracks"].append(
            JbConfig.construct_defaultsession_track(
                name,
                "ReferenceSequenceTrack",
                "LinearReferenceSequenceDisplay",
            )
        )
        # finally give the LGV plugin a heads up
        self.default_view["location"] = "%s:1..%d" % (
            self.default_sequence_name,
            display_width
        )

    def add_alignment(self, track_name, assembly_name, alignment_fp, index_fp):
        """Add an alignment configuration."""
        if assembly_name not in self.assembly_names:
            raise Exception(
                "Cannot add variant track for unknown assembly '%s'"
                % assembly_name
            )

        cram = False
        if "cram" in alignment_fp.lower():
            cram = True

        if cram:
            adapter = {
                "type": "CramAdapter",
                "cramLocation": {
                    "uri": alignment_fp,
                    "locationType": "UriLocation"
                },
                "craiLocation": {
                    "uri": index_fp,
                    "locationType": "UriLocation"
                }
            }
        else:
            adapter = {
                "type": "BamAdapter",
                "bamLocation": {
                    "uri": alignment_fp,
                    "locationType": "UriLocation"
                },
                "index": {
                    "indexType": "BAI",
                    "location": {
                        "uri": index_fp,
                        "locationType": "UriLocation"
                    }
                }
            }

        self.tracks[track_name] = {
            "trackId": track_name,
            "name": track_name,
            "assemblyNames": [assembly_name],
            "type": "AlignmentsTrack",
            "adapter": adapter
        }
        defaultsession_track = JbConfig.construct_defaultsession_track(
            track_name,
            "AlignmentsTrack",
            "LinearAlignmentsDisplay",
        )
        # statements by the utterly deranged
        defaultsession_track["displays"][0]["PileupDisplay"] = {
            "id": "%s-track-display-pileup" % track_name,
            "type": "LinearPileupDisplay",
            "configuration": {
                "type": "LinearPileupDisplay",
                "displayId": "%s-LinearAlignmentsDisplay_pileup_xyz" % track_name  # noqa: E501
            }
        }
        self.default_view["tracks"].append(defaultsession_track)

    def add_variant(self, name, assembly_name, variant_fp, index_fp):
        """Append a variant track to the configuration."""
        if assembly_name not in self.assembly_names:
            raise Exception(
                "Cannot add variant track for unknown assembly '%s'"
                % assembly_name
            )

        track = {
            "type": "VariantTrack",
            "trackId": name,
            "name": name,
            "assemblyNames": [assembly_name],
            "adapter": {
                "type": "VcfTabixAdapter",
                "vcfGzLocation": {
                    "uri": variant_fp,
                    "locationType": "UriLocation"
                },
                "index": {
                    "indexType": "TBI",
                    "location": {
                        "uri": index_fp,
                        "locationType": "UriLocation"
                    }
                }
            }
        }
        self.tracks[name] = track
        self.default_view["tracks"].append(
            JbConfig.construct_defaultsession_track(
                name,
                "VariantTrack",
                "LinearVariantDisplay",
            )
        )

    def as_dict(self):
        """Return configuration as dict."""
        if not self.assembly:
            raise Exception(
                "Cannot structify config without a configured assembly"
            )
        return {
            "assemblies": [self.assembly],
            "defaultSession": self.default_session,
            "tracks": [
                track_dat for track_name, track_dat in self.tracks.items()
            ],
            "configuration": {
                "disableAnalytics": True,
            }
        }

    def as_json(self):
        """Return configuration as JSON string."""
        return json.dumps(self.as_dict(), indent=4)


if __name__ == "__main__":
    import argparse
    from collections import namedtuple

    Reference = namedtuple(
        "Reference",
        "actual_ref_fp sequence_fp index_fp compressed_index_fp"
    )
    Alignment = namedtuple(
        "Alignment",
        "alignment_fp index_fp"
    )
    Variant = namedtuple(
        "Variant",
        "name vcf_fp index_fp"
    )

    # NOTE the alignment and variant file paths are passed through
    # to the output JSON and do not need to exist to run this script...
    # HOWEVER the default session must make reference to the first
    # reference sequence for display, and yet the reference path may
    # not be the one we wish to appear in the JSON (eg. from a workdir;
    # plus we cannot rely on reaching into the out_dir), so we get
    # around this for now by invoking --reference as:
    #   --reference <path_to_real_ref> \
    #               <ref_path_in_config> \
    #               <fai_path_in_config> \
    #               [gzi_path_in_config]
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", nargs='+', required=True)
    parser.add_argument("--alignment", nargs=2, action="append", default=[])
    parser.add_argument("--variant", action="append", nargs=3, default=[])
    args = parser.parse_args()

    if len(args.reference) == 3:
        ref = Reference(*args.reference, None)
    elif len(args.reference) == 4:
        ref = Reference(*args.reference)
    else:
        raise Exception("Incorrect number of args for --reference!")

    config = JbConfig()
    config.add_assembly(
        ref.sequence_fp,
        ref.index_fp,
        ref.compressed_index_fp,
        actual_reference_fp=ref.actual_ref_fp,
        name="ref"
    )
    for alignment_i, alignment_dat in enumerate(args.alignment):
        a = Alignment(*alignment_dat)
        config.add_alignment(
            "alignment-%s" % os.path.basename(a.alignment_fp).split('.')[0],
            "ref",
            a.alignment_fp,
            a.index_fp
        )
    for t in args.variant:
        track = Variant(*t)
        config.add_variant(
            track.name,
            "ref",
            track.vcf_fp,
            track.index_fp
        )

    sys.stdout.write(config.as_json())
