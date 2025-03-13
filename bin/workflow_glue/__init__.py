"""Workflow Python code."""
import argparse
import glob
import importlib
import itertools
import os
import sys

from .util import _log_level, get_main_logger  # noqa: ABS101


__version__ = "0.0.1"
_package_name = "workflow_glue"

HELPERS = "wfg_helpers"


def get_components(allowed_components=None):
    """Find a list of workflow command scripts."""
    logger = get_main_logger(_package_name)

    # gather all python files in the current directory and the wfg_helpers
    home_path = os.path.dirname(os.path.abspath(__file__))
    standard_lib = os.path.join(home_path, HELPERS)
    globs = itertools.chain.from_iterable((
        glob.glob(os.path.join(path, "*.py"))
        for path in (home_path, standard_lib)))

    components = dict()
    for fname in globs:
        name = os.path.splitext(os.path.basename(fname))[0]
        if name in ("__init__", "util"):
            continue
        if allowed_components is not None and name not in allowed_components:
            continue

        # leniently attempt to import module
        try:
            if HELPERS in fname:
                mod = importlib.import_module(f"{_package_name}.{HELPERS}.{name}")
            else:
                mod = importlib.import_module(f"{_package_name}.{name}")
        except ModuleNotFoundError as e:
            # if imports cannot be satisifed, refuse to add the component
            # rather than exploding
            logger.warn(f"Could not load {name} due to missing module {e.name}")
            continue

        # if theres a main() and and argparser() thats good enough for us.
        try:
            req = "main", "argparser"
            if all(callable(getattr(mod, x)) for x in req):
                components[name] = mod
        except Exception:
            pass
    return components


def cli():
    """Run workflow entry points."""
    logger = get_main_logger(_package_name)
    logger.info("Bootstrapping CLI.")
    parser = argparse.ArgumentParser(
        'wf-glue',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-v', '--version', action='version',
        version='%(prog)s {}'.format(__version__))

    subparsers = parser.add_subparsers(
        title='subcommands', description='valid commands',
        help='additional help', dest='command')
    subparsers.required = True

    # importing everything can take time, try to shortcut
    if len(sys.argv) > 1:
        components = get_components(allowed_components=[sys.argv[1]])
        if not sys.argv[1] in components:
            logger.warn("Importing all modules, this may take some time.")
            components = get_components()
    else:
        components = get_components()

    # add all module parsers to main CLI
    for name, module in components.items():
        p = subparsers.add_parser(
            name.split(".")[-1], parents=[module.argparser()])
        p.set_defaults(func=module.main)

    args = parser.parse_args()

    logger.info("Starting entrypoint.")
    args.func(args)
