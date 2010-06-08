#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from pprint import pprint
from xupy import XParam, xdsInp2Param, opWriteCl, \
                 saveLastVersion, LP_names, xdsinp_base, \
                 SPGlib, Lattice, resum_scaling, write_xscale_resum, \
		 get_BravaisToSpgs


__version__ = "0.3.8"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "14-10-2008"
__copyright__ = "Copyright (c) 2006-2008 Pierre Legrand"
__license__ = "LGPL"

sys.path.append("/data/bioxsoft/progs/PXPY/XDS")
    
def parse_spaceGroup(spg):
    try:
        _spg = int(spg)
    except ValueError:
        _spg = 0
        _spg_name = spg.upper()
        for spgn in SPGlib:
            if _spg_name in SPGlib[spgn]:
                _spg = spgn
                break
    if _spg:
        _spg_info = SPGlib[_spg]
        _spg_str = "  Imposed Space group:  %s,  number %d" % \
               (_spg_info[1], _spg)
    else:
        raise Exception, "\n ERROR: Unrecognised space group: %s\n" % spg
    return _spg, _spg_info, _spg_str

def select_strategy(idxref_results, xdsPar):
    selSpgN = xdsPar["SPACE_GROUP_NUMBER"]
    selAno =  xdsPar["FRIEDEL'S_LAW"] 
    R = idxref_results
    validInp = False
    Bravais_to_SPGs = get_BravaisToSpgs()
    # Select LATTICE
    while not validInp:
        defSel = 1
        if selSpgN != 0:
            # choose the lattice solution according to the selected spg.
            i = 0
            for LAT in R["lattices_table"]:
                if LAT.fit <= LATTICE_GEOMETRIC_FIT_CUTOFF:
                    i += 1
                    if selSpgN in Bravais_to_SPGs[LAT.Bravais_type]:
                        defSel = i
        selection = raw_input("\n Select a solution number [%d]: " % defSel)
        # If the selection is not compatible with the spg, set not valid
        _sel = selection.split()
        selnum = 1
        try:
            if len(_sel) == 1:
                selnum = int(_sel[0])
                validInp = True
            elif len(_sel) == 0:
                selnum = defSel
                validInp = True
            else:
                raise Exception, "Invalid selection input."
        except Exception, err:
            print "\n ERROR. ", err
    selLat = R["lattices_table"][selnum-1]
    if selSpgN == 0:
        selSpgN = selLat.symmetry_num
    validInp = False
    # Select SPACEGROUP
    print " Possible spacegroup for this lattice are:\n"
    for _sSpg in Bravais_to_SPGs[selLat.Bravais_type]:
        print "  %15s, number: %3d" % (SPGlib[_sSpg][1], _sSpg)
    while not validInp:
        selection = raw_input("\n Select the spacegroup [%s, %d]: "
                             % (SPGlib[selSpgN][1], selSpgN))
        _sel = selection.split()
        try:
            if len(_sel) == 1:
                selSpgN, _spg_info, _spg_str = parse_spaceGroup(_sel[0])
                selSpgS = _spg_info[1]
                validInp = True
            elif len(_sel) == 0:
                validInp = True
            else:
                raise Exception, "Invalid selection input."
            if selSpgN not in Bravais_to_SPGs[selLat.Bravais_type]:
                validInp = False
                msg = "Inconsistant combinaison of Bravais lattice"
                msg += " and spacegroup.\n For this Bravais Lattice"
                msg += " (%s), spacegroup should be one of these:\n\n" % \
                        (selLat.Bravais_type)
                for _sSpg in Bravais_to_SPGs[selLat.Bravais_type]:
                    msg += "  %15s, number: %3d\n" % (SPGlib[_sSpg][1], _sSpg)
                raise Exception, msg
        except Exception, err:
            print "\n ERROR. ", err
    validInp = False
    # Select ANOMALOUS
    while not validInp:
        if selAno:
            txt3 = "N/y"
        else:
            txt3 = "Y/n"
        selection = raw_input(" Anomalous [%s]: " % txt3)
        try:
            _ans =  selection.strip()
            if _ans == "":
                validInp = True    
            elif _ans[0] in "Yy":
                xdsPar["FRIEDEL'S_LAW"] = False
                validInp = True
            elif _ans[0] in "Nn":
                xdsPar["FRIEDEL'S_LAW"] = True
                validInp = True
            else:
                raise Exception, "Invalid answer [Y/N]."
        except Exception, err:
            print "\n ERROR. ", err
    print "\n Selected  cell paramters:  ", selLat
    if selSpgN > 2:
        selLat.idealize()
        print " Idealized cell parameters: ", selLat.prt()
        xdsPar["UNIT_CELL_CONSTANTS"] = selLat.prt()
    xdsPar["SPACE_GROUP_NUMBER"] = selSpgN
    # Select just the internal circle of the detector.
    xdsPar["TRUSTED_REGION"] = 0.0, 1.0
    return xdsPar

from XDS import XDSLogParser
LATTICE_GEOMETRIC_FIT_CUTOFF = 100
fmt_lat = "%3d  %5s %7.2f  %4d  %s"

for lpfile in sys.argv[1:]:
     if os.path.isfile(lpfile) and (".LP" in lpfile):
         res = XDSLogParser(lpfile, verbose=True).results
         i = 0
         print "    TABLE OF POSSIBLE LATTICES:\n"
         print " num  Symm  quality  mult     a       b       c",
         print "     alpha    beta   gamma"
         print " "+"-"*73
         if "lattices_table" in res:
            for LAT in res["lattices_table"]:
                if LAT.fit <= LATTICE_GEOMETRIC_FIT_CUTOFF:
                    i += 1
                    print fmt_lat % (i, LAT.symmetry_str1,
                                             LAT.fit, LAT.multiplicity, LAT)
         #print get_BravaisToSpgs()
	 print select_strategy(res, {"SPACE_GROUP_NUMBER":96, "FRIEDEL'S_LAW":True}, )
         print select_strategy(res, {"SPACE_GROUP_NUMBER":0, "FRIEDEL'S_LAW":True}, )
