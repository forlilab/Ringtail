#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
import sys
from ringtail import RingtailCore

def main():
    """Will generate a file config.json that contains all Ringtail options with their default vlaues."""
    RingtailCore.generate_config_file_template()
    return

if __name__ == "__main__":
    sys.exit(main())
