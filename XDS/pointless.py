#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from xml.dom import minidom


def pointless(dir_name, hklinp="XDS_ASCII.HKL"):
    cmline = "pointless XDSIN %s XMLOUT XDS_pointless.xml"
    cmline += " HKLOUT XDS_pointless.mtz > XDS_pointless.log"
    os.chdir(dir_name)
    os.system(cmline % hklinp)

    xml_inp = "XDS_pointless.xml"
    likely_spacegroups = []
    prob_max = 0
    dom = minidom.parse(xml_inp)
    get_elem = lambda n, m, f: f(
                   n.getElementsByTagName(m)[0].childNodes[0].data.strip())
    cell = dom.getElementsByTagName('cell')[0]
    zone_list = dom.getElementsByTagName('ZoneScoreList')[0]
    spg_list = dom.getElementsByTagName('SpacegroupList')[0]

    cell_par = dict([(x, get_elem(cell, x, float)) for x in ('a', 'b', 'c')])
    # Looking at systematique extinctions
    print "\n      Systematique extinctions from pointless:"
    print "      Zone Type            axe len.   #obs     Condition    Prob."
    print "      "+59*"-"
    for node in zone_list.getElementsByTagName('Zone'):
        ztype = get_elem(node, 'ZoneType', str)
        nobs = get_elem(node, 'Nobs', int)
        prob = get_elem(node, 'Prob', float)
        condition = get_elem(node, 'Condition', str)
        axe = cell_par[ztype[ztype.index("[")+1:ztype.index("]")]]
        all_dat = (ztype, axe, nobs, condition, prob)
        print "%25s %8.1f %7d  %12s  %7.3f" % all_dat
    print "\n      Possible spacegroup from pointless:"
    print "      Symbol      num   TotalProb   SysAbsProb"
    print "      "+40*"-"
    # looking for most probable spacegroup
    for node in spg_list.getElementsByTagName('Spacegroup'):
        total_prob = get_elem(node, 'TotalProb', float)
        sys_abs_prob = get_elem(node, 'SysAbsProb', float)
        spg_name = get_elem(node, 'SpacegroupName', str)
        spg_num = get_elem(node, 'SGnumber', int)

        all_dat = (spg_name, spg_num, total_prob, sys_abs_prob)
        prob_max = max(total_prob, prob_max)
        print "%15s   #%d  %9.3f    %9.3f" % all_dat
        if total_prob == prob_max:
            likely_spacegroups.append(all_dat)
    os.chdir("..")
    print 
    return likely_spacegroups

if __name__ == '__main__':
    import os
    from xml.dom import minidom
    pointless(".")