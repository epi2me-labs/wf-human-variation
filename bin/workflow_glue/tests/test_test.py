"""A dummy test."""

import argparse

from workflow_glue import check_sample_sheet


def test():
    """Just showing that we can import using the workflow-glue."""
    assert isinstance(check_sample_sheet.argparser(), argparse.ArgumentParser)
