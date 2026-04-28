#!/usr/bin/env python3
"""Milestone 1 placeholder for future fresh run cleanup."""

import argparse


def build_parser():
    parser = argparse.ArgumentParser(
        description="Placeholder cleanup entry point; no files are deleted in Task 1."
    )
    parser.add_argument("--run-id", help="Future run identifier under fresh/runs/.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Accepted for future compatibility. Task 1 is always non-destructive.",
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    run_label = args.run_id or "<not provided>"
    print("Task 1 cleanup placeholder; no files deleted for run_id={0}.".format(run_label))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
