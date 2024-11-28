"""Reheader a SAM in a stream.

When using the bam2fq -> minimap2 pattern for (re)aligning BAM data, we
lose any existing RG and PG headers. This is particularly egregious when
handling basecalled data as lines related to dorado basecalling settings
as well as dorado RG headers are lost; orphaning RG tags in the reads.
This is problematic for downstream anaylses that would like to read the
XAM header to intelligently determine how to handle the reads based on
the basecaller model and basecaller configuration.

This script handles:
  - Inserting RG, PG and CO lines from an existing XAM header into the
    header of the SAM emitted from minimap2's alignment stream
  - Inserting a PG header to indicate that a call to bam2fq was made
  - Updating the first streamed PG.PP parent tag with the last PG.ID
    of the existing XAM header to maintain a chain of custody
  - Updating any streamed PG.ID (and PG.PP) tags to avoid collisions
    with inserted PG.ID

Handling collisions may seem like overkill but it is anticipated that
this script will be called immediately after minimap2, any previous
attempt to use minimap2 will lead to ambiguity. This would be the
expected case where users have used wf-basecalling or wf-alignment to
align a set of reads, only to realign them to another reference (eg.
via wf-human-variation). Arguably, we should remove older references to
minimap2 as they will have been invalidated by the call to bam2fq but
removing PG records and sticking the PG chain back together seems more
fraught with annoying future bugs than simply resolving conflicts.

This script will explode on a stream that contains:
  - PG lines in the original header where the last PG in the chain is
    ambiguous, or where the parent PP IDs are not injective
  - PG lines in the stream that do not appear in the order of their
    chain (that is if a PG.PP refers to a PG.ID that has not been
    encountered yet)

SQ lines are retained after an HD line. That is to say, the most recent
set of SQ lines observed after an HD will appear in the final output.
SQ, RG, PG and CO lines are emitted as a group together, with elements
written out in the order observed.

PG lines are naively appended to the last PG element in the chain. No
attempt is made to keep multiple program chains intact as this can lead
to bloated headers. Broken PG metadata is a known problem (see
samtools/hts-specs#275) but one that is preferable to headers that
become unwieldly large to process: there IS an upper limit to a SAM
header's size after all.

This script takes advantage of minimap2's SAM output to immediately
reheader the stream before any downstream calls to other programs pollute
the PG header. This script is a little overkill but attempts to be robust
with handling PG collisions and more obviously encapsulates reheadering
behaviour, and leaves some room to do more clever things as necessary.
"""
import sys

from .util import wf_parser  # noqa: ABS101


class SamHeader:
    """An overkill container to manage merging PG lines in SAM headers.

    Collision handling is simple. If a PG.ID is duplicated by the stream
    then we add a suffix to its name and keep an eye out for the
    corresponding PG.PP later. We assume that headers emitted by the
    stream are chronological because this script should not be called as
    part of any complicated pipework other than immediately following
    minimap2.
    """

    def __init__(self):
        """Initialise a collision aware PG container."""
        self.remapped_pgids = {}
        self.collision_suffix = 0

        # Default HD, in case the new stream does not provide one
        self.hd = "@HD\tVN:1.6\tSO:unknown"

        # We'll merge RG, CO and PG
        self.rg_records = []
        self.co_records = []
        self.pg_records = []

        # We keep the most recently observed block of SQ records by
        # resetting SQ on the first SQ seen after non-SQ. We cannot
        # rely on HD being emitted (as minimap2 does not do this!)
        self.sq_records = []
        self.reset_sq = False

        self.observed_rgids = set()
        self.observed_pgids = set()
        self.last_pgid = None

    @staticmethod
    def str_to_record(line):
        """Return an appropriate struct for a given string record."""
        try:
            record_type, record_data = line.strip().split('\t', 1)
        except ValueError:
            raise Exception(f"Record type could not be determined: {line}")

        if len(record_type) > 3:
            raise Exception(f"Record type malformed: {record_type}")

        record = {}
        if record_type in ["@HD", "@CO", "@SQ"]:
            return record_type, record_data
        elif record_type in ["@RG", "@PG"]:
            allowed_keys = {
                "@RG": ["ID", "BC", "CN", "DS", "DT", "FO", "KS", "LB", "PG", "PI", "PL", "PM", "PU", "SM"],  # noqa:E501
                "@PG": ["ID", "PN", "CL", "PP", "DS", "VN"]
            }
            for field in record_data.strip().split('\t'):
                k, v = field.split(':', 1)
                if k not in allowed_keys[record_type]:
                    raise Exception(f"{record_type} with bad key '{k}': {record_data}")
                record[k] = v
            if "ID" not in record:
                raise Exception(f"{record_type} with no ID: {record_data}")
            return record_type, record
        else:
            raise Exception(f"Unknown record type: {line}")

    @staticmethod
    def record_to_str(record_type, record_data):
        """Form a string from a header record."""
        if record_type in ["@PG", "@RG"]:
            tags = [f"{k}:{v}" for k, v in record_data.items()]
            return f"{record_type}\t" + '\t'.join(tags)
        elif record_type in ["@SQ", "@CO"]:
            return f"{record_type}\t{record_data}"

    @staticmethod
    def resolve_pg_chain(pg_dicts):
        """Check links between PG.ID and PP.ID, exploding if inconsistent."""
        links = {}
        # Document links between all ID and their PP parent
        pgids_without_ppid = 0
        for pgd in pg_dicts:
            pgid = pgd["ID"]
            pgpp = pgd.get("PP")
            links[pgid] = pgpp
            if pgpp is None:
                pgids_without_ppid += 1
        if len(links) > 0:
            # If there are links, exactly one should have a None parent
            # to indicate the first PG in the chain. Explode if we see
            # no head or multiple heads.
            if pgids_without_ppid == 0:
                raise Exception("PG chain does not have a head.")
            elif pgids_without_ppid > 1:
                raise Exception("PG chain has multiple heads.")
        for source in links:
            head = source
            path = [head]
            while True:
                head = links[head]
                if head is None:
                    break
                if head in path:
                    path.append(head)
                    raise Exception(f"PG chain appears to contain cycle: {path}")
                path.append(head)
        # This function is only really called to catch any explosions
        # but we'll return the links here as it is useful for testing
        return links

    def _bump_pg_collider(self):
        """Alter the collision suffix after determining a collision."""
        self.collision_suffix += 1

    def _uncollide_pgid(self, pgid):
        """Return an uncollided string for a given PG ID."""
        new_pgid = f"{pgid}-{self.collision_suffix}"
        self.remapped_pgids[pgid] = new_pgid
        self._bump_pg_collider()
        return new_pgid

    def add_line(self, line):
        """Add a header line to the header."""
        record_type, record = self.str_to_record(line)

        if record_type == "@HD":
            self.hd = f"@HD\t{record}"
        elif record_type == "@CO":
            self.co_records.append(record)
        elif record_type == "@SQ":
            if self.reset_sq:
                self.sq_records = []
                self.reset_sq = False
            self.sq_records.append(record)
        elif record_type == "@RG":
            rgid = record["ID"]
            if rgid not in self.observed_rgids:
                self.observed_rgids.add(rgid)
                self.rg_records.append(record)
            elif record not in self.rg_records:
                # if rgid has been seen before, abort if this record is different
                raise Exception(
                    f"Duplicate RG with ID '{rgid}' conflicts with previously seen RG with same ID."  # noqa:E501
                )
        elif record_type == "@PG":
            pgid = record["ID"]
            if pgid in self.observed_pgids:
                # collision, rewrite the pgid
                pgid = self._uncollide_pgid(pgid)
                record["ID"] = pgid
            else:
                self.observed_pgids.add(pgid)

            # maintain chain
            ppid = record.get("PP")
            if not ppid:
                # record has no parent, this is either
                # - the first record (last_pgid is None) so is the tail
                # - an inserted record that needs its parent to be the current tail
                if not self.last_pgid:
                    self.last_pgid = pgid
                else:
                    record["PP"] = self.last_pgid
                    self.last_pgid = pgid
            else:
                if ppid not in self.observed_pgids:
                    raise Exception(
                        f"Encountered PG.PP '{ppid}' before observing corresponding PG.ID"  # noqa:E501
                    )
                # remap parent id (if needed)
                record["PP"] = self.remapped_pgids.get(ppid, ppid)
                # set tail to this record
                self.last_pgid = pgid

            self.pg_records.append(record)

        if len(self.sq_records) > 0 and record_type != '@SQ':
            self.reset_sq = True

        return record

    def write_header(self, fh):
        """Write this header to a file handle."""
        self.resolve_pg_chain(self.pg_records)  # check PG header
        fh.write(f"{self.hd}\n")
        for sq in self.sq_records:
            fh.write(self.record_to_str("@SQ", sq) + '\n')
        for rg in self.rg_records:
            fh.write(self.record_to_str("@RG", rg) + '\n')
        for pg in self.pg_records:
            fh.write(self.record_to_str("@PG", pg) + '\n')
        for co in self.co_records:
            fh.write(self.record_to_str("@CO", co) + '\n')


def reheader_samstream(header_in, stream_in, stream_out, args):
    """Run reheader_samstream."""
    # read original header into container
    sh = SamHeader()
    for line in header_in:
        sh.add_line(line)

    # append user provided lines to container
    for line in args.insert:
        sh.add_line(line)

    # read the header portion of the minimap2 stream
    wrote_header = False
    for line in stream_in:
        if line[0] != '@':
            # write out header on first alignment
            sh.write_header(stream_out)
            wrote_header = True
            # and actually write the first alignment
            stream_out.write(line)
            break
        sh.add_line(line)

    # Pass through the rest of the alignments
    for line in stream_in:
        stream_out.write(line)
    # If there were no alignments, we won't have hit the != @ case in the first stdin,
    # and we won't have written the header out. Write a header if we haven't already.
    if not wrote_header:
        sh.write_header(stream_out)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("reheader_samstream")
    parser.add_argument("header_in")
    parser.add_argument("--insert", action="append", default=[])
    return parser


def main(args):
    """reheader_samstream default entry point."""
    with open(args.header_in) as header_in:
        reheader_samstream(header_in, sys.stdin, sys.stdout, args)
