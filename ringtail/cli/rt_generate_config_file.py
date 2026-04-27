#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
import sys
from ringtail import RingtailCore

def main():
    """Will generate a file config.json that contains all Ringtail options with their default vlaues."""
    try:
        RingtailCore.generate_config_file_template()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
