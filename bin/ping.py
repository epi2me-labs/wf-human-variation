#!/usr/bin/env python
"""Send workflow ping."""

import argparse
import json
import uuid

from epi2melabs import ping


def get_uuid(val):
    """Construct UUID from string."""
    return uuid.UUID(str(val))


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--hostname", required=True, default=None,
        help="ping some meta")
    parser.add_argument(
        "--opsys", required=True, default=None,
        help="ping some meta")
    parser.add_argument(
        "--session", default=None,
        help="ping some meta")
    parser.add_argument(
        "--message", required=True, default=None,
        help="message to include in the ping")
    parser.add_argument(
        "--meta", default=None,
        help="JSON file of metadata to be included in the ping")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--disable", action='store_true',
        help="Run the script but don't send the ping")
    args = parser.parse_args()

    meta = None
    if args.meta:
        with open(args.meta, "r") as json_file:
            meta = json.load(json_file)

    if not args.disable:
        ping.Pingu(
            get_uuid(args.session),
            hostname=args.hostname,
            opsys=args.opsys
        ).send_workflow_ping(
            workflow='wf-template',
            message=args.message,
            revision=args.revision,
            commit=args.commit,
            meta=meta
        )


if __name__ == "__main__":
    main()
