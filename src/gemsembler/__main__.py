import argparse
import json
import sys

import yaml

from .__init__ import __version__
from .gathering import load_sbml_model
from .pathsfinding import pathsfinding
from .quickdiff import quickdiff


def pathsfinding_cli(args):
    desc = (
        "gemsembler will use metquest package to find a path "
        "between given metabolites in the metabolic network"
    )
    parser = argparse.ArgumentParser("pathsfinding", description=desc)
    parser.add_argument("path_to_model", help="Path to the model.")
    parser.add_argument(
        "path_to_config", help="Path to yaml configuration file for metquest."
    )
    parser.add_argument("output_dir", help="Output directory.")

    args = parser.parse_args(args)

    with open(args.path_to_config) as fh:
        config = yaml.safe_load(fh)

    pathsfinding(args.path_to_model, args.output_dir, **config)


def quickdiff_cli(args):
    desc = "Quick comparison between two cobrapy models."
    parser = argparse.ArgumentParser("quickdiff", description=desc)

    parser.add_argument("model1", help="Path to first model")
    parser.add_argument("model2", help="Path to second model")
    parser.add_argument("-o", "--output", help="Path to output json file.")

    args = parser.parse_args(args)

    model1 = load_sbml_model(args.model1)
    model2 = load_sbml_model(args.model2)

    diffs = quickdiff(model1, model2)

    if args.output is not None:
        if diffs is True:
            diffs = {}

        with open(args.output, "w") as fh:
            json.dump(diffs, fh)
    else:
        print(diffs)


def main():
    desc = (
        "CLI interface for gemsembler."
        "To get more information on each subcommand, type `gemsembler <subcommand> -h`"
        "1. gemsembler can run metquest to find paths within the metabolic network"
        "2. gemsembler can compare two cobrapy models with quickdiff subcommand"
    )

    parser = argparse.ArgumentParser(prog="gemsembler", description=desc)

    parser.add_argument("command", choices=["pathsfinding", "quickdiff"])
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__,
        help="Show version number",
    )
    # parse_args defaults to [1:] for args, but you need to
    # exclude the rest of the args too, or validation will fail
    args = parser.parse_args(sys.argv[1:2])

    if args.command == "pathsfinding":
        pathsfinding_cli(sys.argv[2:])
    elif args.command == "quickdiff":
        quickdiff_cli(sys.argv[2:])
