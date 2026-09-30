#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
import argparse
import sys
from ringtail import RingtailCore

def main():
    """Will generate a file config.json that contains all Ringtail options with their default values."""
    parser = argparse.ArgumentParser(
        prog="rt_generate_config_file",
        description="Write a config file template with the default values of the rt_process_vs options.",
    )
    parser.add_argument(
        "-o", "--output", default="config.json", help="file to write (default: config.json), must not exist"
    )
    args = parser.parse_args()
    try:
        print(f"Wrote {RingtailCore.generate_config_file_template(args.output)}")
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
