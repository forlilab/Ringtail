#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

from ringtail import RingtailCore

if __name__ == "__main__":
    """Will generate a file config.json that contains all Ringtail options with their default vlaues."""
    RingtailCore.generate_config_file_template()
