"""A dummy test."""

import argparse

from workflow_glue import report_al as report


def test():
    """Just showing that we can import using the workflow-glue."""
    assert isinstance(report.argparser(), argparse.ArgumentParser)
