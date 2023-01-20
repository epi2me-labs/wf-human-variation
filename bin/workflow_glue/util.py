"""The odd helper function."""

import argparse
import logging

_log_name = None


def get_main_logger(name):
    """Create the top-level logger."""
    global _log_name
    _log_name = name
    logging.basicConfig(
        format='[%(asctime)s - %(name)s] %(message)s',
        datefmt='%H:%M:%S', level=logging.INFO)
    return logging.getLogger(name)


def get_named_logger(name):
    """Create a logger with a name.

    :param name: name of logger.
    """
    name = name.ljust(10)[:10]  # so logging is aligned
    logger = logging.getLogger('{}.{}'.format(_log_name, name))
    return logger


def wf_parser(name):
    """Make an argument parser for a workflow command."""
    return argparse.ArgumentParser(
        name,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)


def _log_level():
    """Parser to set logging level and acquire software version/commit."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    modify_log_level = parser.add_mutually_exclusive_group()
    modify_log_level.add_argument(
        '--debug', action='store_const',
        dest='log_level', const=logging.DEBUG, default=logging.INFO,
        help='Verbose logging of debug information.')
    modify_log_level.add_argument(
        '--quiet', action='store_const',
        dest='log_level', const=logging.WARNING, default=logging.INFO,
        help='Minimal logging; warnings only.')

    return parser
