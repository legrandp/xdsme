#!/usr/bin/env python
# -*- coding: utf-8 -*-

"Small utility to extract XDS *.LP logfile information."

__version__ = "0.0.1"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "23-11-2008"
__copyright__ = "Copyright (c) 2009 Pierre Legrand"
__license__ = "New BSD"

import sys
import os

from XDS import XDSLogParser

for lpfile in sys.argv[1:]:
    if os.path.isfile(lpfile) and (".LP" in lpfile):
        res = XDSLogParser(lpfile, verbose=True).results
