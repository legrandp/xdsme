#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.8.3"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "26-06-2012"
__copyright__ = "Copyright (c) 2006-2012 Pierre Legrand"
__license__ = "New BSD http://www.opensource.org/licenses/bsd-license.php"

# Environemantal variable XDSHOME, if set, defines the place where the xds
# executables will be searched. The parallelized execs (xds_par, xscale_par)
# will be used be defaults.

import sys
import os

atom_names = ['Ag', 'Al', 'Ar', 'As', 'Au', 'B', 'Ba', 'Be', 'Bi', 'Br',
'C', 'Ca', 'Cd', 'Ce', 'Cl', 'Co', 'Cr', 'Cs', 'Cu', 'Dy', 'Er', 'Eu', 'F',
'Fe', 'Ga', 'Gd', 'Ge', 'H', 'He', 'Hf', 'Hg', 'Ho', 'I', 'In', 'Ir', 'K',
'Kr', 'La', 'Li', 'Lu', 'Mg', 'Mn', 'Mo', 'N', 'Na', 'Nb', 'Nd', 'Ne',
'Ni', 'O', 'Os', 'P', 'Pb', 'Pd', 'Pm', 'Pr', 'Pt', 'Pu', 'Rb', 'Re',
'Rh', 'Rn', 'Ru', 'S', 'Sb', 'Sc', 'Se', 'Si', 'Sm', 'Sn', 'Sr', 'Ta',
'Tb', 'Tc', 'Te', 'Th', 'Ti', 'Tl', 'Tm', 'U', 'V', 'W', 'Xe', 'Y',
'Yb', 'Zn', 'Zr']

options = ["CCP4","CCP4F","CNS","SHELX","SOLVE","EPMR","CRANK",
           "AMORE","SHARP","PHASER","REPLACE"]

usage   = """
    USAGE : %s FILE OPTIONS [free_hkl_to_inherit] [nSites] [atomType] EXPORT_MODE \n
        EXPORT_MODE can be one of these:\n
            %s\n
        FILE is a XDS or XSCALE reflection file (usually XDS_ASCII.HKL).\n
    OPTIONS:
        -a       force anomal output (Friedel's law = False)
        -n       force normal output (Friedel's law = True)
        -m       force merged output
        -u       force unmerged output
        -f       force generation of free reflection flag
        -nf      force no generation of free reflection flag
        Default is keeping the XDS input file settings.

        -l label, or -l=label
                 In case of CCP4 export, give a new label columns.
                 The defaults labels: FP, SIGFP, DANO, SIGDANO, ISYM with
                 -l pk or -l=pk becomes:
                 FP_pk, SIGFP_pk, DANO_pk ... in CCP4 mode
                 FP_pk, SIGFP_pk, F(+)_pk, SIGF(+)_pk, ... in PHASER mode

        free_hkl_to_inherit: is a reflection file containing a previously
                 selected set of free reflection to be kept in the newly
                 exported reflection file for calculation of unbiased Rfree.
                 The accepted file format are: SHELX, CNS and MTZ (with the
                 following labels:  FP=FP SIGFP=SIGFP FREE=FreeR_flag).\n
        nSites: integer describing the number of heavy atom sites expected.\n
        atomType: Symbol of the heavy atom type expected. Only one or two
                 letters symbols are recognised (like I, Se, S, Hg).\n
      NOTE:
          i)   Keywords free_hkl_to_inherit, nSites, atomeType and EXPORT_MODE
               can be given in any order.
          ii)  Cell parameters, space group number and wavelength are taken
               from the XDS reflection file header.
          iii) If Friedel's law == False and heavy atom type is not given, a
               guess is made base on the wavelength closest atome type edge.
          iv)  All the exported files are created in a new local directory named
               after the requested mode (./ccp4, ./phaser, ./solve...).
          v)   In most modes, custom bash scripts are created to allow a rapid
               interactive exploration of data processing.

      EXAMPLES:
          xdsconv.py XDS_ASCII.HKL shelx 12 Se
          xdsconv.py solve XDS_ASCII.HKL Se 12
          xdsconv.py 12 Se phaser XDS_ASCII.HKL
          xdsconv.py XDS_ASCII.HKL ccp4 -n FreeR_reference.mtz
          xdsconv.py XDS_ASCII.HKL ccp4 Se 12 -l=peak FreeR_reference.mtz
"""

fmt_inp = """
<== Input file:           %(file_name_in)s
<<< Space group number:      %(spgn_in)s
<<< Cell parameters:         %(cell_str)s
<<< Friedel's law:           %(friedel_in)s
<<< Merge:                   %(merge_in)s
<<< Dirname:                 %(dirname)s
<<< Wavelength               %(wavelength).4f
"""

fmt_outp = """==> Output file:          %(file_name_out)s
==> Output mode:             %(mode)s
>>> Friedel's law:           %(friedel_out)s
>>> Merge:                   %(merge_out)s
>>> Resolution limit:        %(res_low).2f - %(res_high).2f Angstroem
""" 

fmt_outp_ha = """
>>> From Wavelength          %(wavelength).4f
     ??? Guessed Atome Edge       %(ha_name)s
     ???         F'               %(fp).2f
     ???         F"               %(fpp).2f
"""

xdsconv_script = """UNIT_CELL_CONSTANTS=  %(cell_str)s
SPACE_GROUP_NUMBER=   %(spgn_in)s
INPUT_FILE=%(file_name_in)s %(file_type)s %(res_low).2f %(res_high).2f 
OUTPUT_FILE=%(dir_mode)s%(file_name_out)s %(mode_out)s FRIEDEL'S_LAW=%(friedel_out)s MERGE=%(merge_out)s
"""

shelxc_script = """SAD  %(file_name_out)s
CELL %(cell_str)s 
SPAG  %(spg_name)s
NTRY 1000
FIND %(num_sites)d
SFAC %(ha_name)s
MIND -3.5
MAXM 2
"""

shelxall_script = """#!/bin/sh

search_sites_shelx () {

ID=$1      # Identifier
RES=$2     # High resolution limit
NSITES=$3  # Number of sites
PATS=$4    # Patterson seeding option

shelxc ${ID} << EOFC > ${ID}_shelxc.log
SAD  %(file_name_out)s
CELL %(cell_str)s 
SPAG  %(spg_name)s
SHEL 999 ${RES}
NTRY 1000
FIND  ${NSITES}
SFAC %(ha_name)s
ESEL 1.4
MIND -3.5 1
MAXM 4
EOFC

if [ $PATS != "yes" ] ; then

#grep -v PATS ${ID}_fa.ins > ${ID}_nopats_fa.ins
sed -e 's/^PATS/WEED 0.3\\nSKIP 0.5/'  ${ID}_fa.ins > ${ID}_nopats_fa.ins
cp ${ID}.hkl ${ID}_nopats.hkl
cp ${ID}_fa.hkl ${ID}_nopats_fa.hkl

shelxd ${ID}_nopats_fa > ${ID}_nopats_shelxd.log

else

shelxd ${ID}_fa > ${ID}_shelxd.log

fi

}

id="aaa"
nsites=%(num_sites)d
pats="yes"

# loop to start batch process with different resoution limit
for res in 3.2 3.5 4.0 4.4 ; do

search_sites_shelx ${id}${res} ${res} ${nsites} ${pats} &

done
# end of the loop.
"""
sitcom_script ="""#!/bin/sh

do_sitcom () {
ID=$1      # Identifier
NSITES=$2  # Number of sites
PATS=$3    # Patterson seeding option

if [ $PATS = "yes" ] ; then

cat  << EOFS > ${ID}_sitcom.inp
unit_cell    %(cell_str)s 
space_group  %(spgn_in)s
str_name     ${ID}
deriv_atyp   %(ha_name)s
#
#        TAG     WEIGHT   FILE            N(SOL) N(SITES)
#
read_set ${ID}  1.0      ${ID}_fa.lst   10    ${NSITES}
#read_set SHELXD2  1.0      a4_fa.lst     10    ${nsites}
#read_sol HYSS    1.0      a3_fa.pdb
#read_sol SHELXD5  1.0      a2_fa.pdb
EOFS

sitcom < ${ID}_sitcom.inp > ${ID}_sitcom.log

cat  << EOFS >> tmp.all_sitcom.inp
read_sol ${ID}  1.0      ${ID}_fa.pdb   ${NSITES}
EOFS

else

cat  << EOFS > ${ID}_nopats_sitcom.inp
unit_cell    %(cell_str)s 
space_group  %(spgn_in)s
str_name     ${ID}_nopats
deriv_atyp   %(ha_name)s
#
#        TAG     WEIGHT   FILE            N(SOL) N(SITES)
#
read_set ${ID}np  1.0      ${ID}_nopats_fa.lst   10    ${NSITES}
#read_set SHELXD2  1.0      a4_fa.lst     10    ${nsites}
#read_sol HYSS    1.0      a3_fa.pdb
#read_sol SHELXD5  1.0      a2_fa.pdb
EOFS

sitcom < ${ID}_nopats_sitcom.inp > ${ID}_nopats_sitcom.log

cat  << EOFS >> tmp.all_sitcom.inp
read_sol ${ID}np  1.0      ${ID}_nopats_fa.pdb   ${NSITES}
EOFS

fi
}

id="aaa"
nsites=%(num_sites)d
pats="yes"

cat  << EOFS > tmp.all_sitcom.inp
unit_cell    %(cell_str)s 
space_group  %(spgn_in)s
str_name     ${id}
deriv_atyp   %(ha_name)s
#
#        TAG     WEIGHT   FILE            N(SOL) N(SITES)
#
EOFS

# loop to start batch process with different resoution limit
for res in 3.2 3.5 4.0 4.4 ; do

do_sitcom ${id}${res} ${nsites} ${pats}

done
# end of the loop

# compare all solutions
if [ $pats = "yes" ] ; then
mv tmp.all_sitcom.inp ${id}_all_sitcom.inp
sitcom < ${id}_all_sitcom.inp > ${id}_all_sitcom.log
else
mv tmp.all_sitcom.inp ${id}_all_np_sitcom.inp
sitcom < ${id}_all_np_sitcom.inp > ${id}_all_np_sitcom.log
fi
"""

xprep_script = """../%(file_name_in)s
X\nY\n\nS\nI
%(cns_sg)s
Y
A\nA\nA
100\n
%(ident)s\nY\n
%(ha_name)s\n4\n%(wavelength)s\n3.0\n
E\nP\nR\n15 3
P\n%(ident)s.pattX.ps\nX\n\n\n
P\n%(ident)s.pattY.ps\nY\n\n\n
P\n%(ident)s.pattZ.ps\nZ\n\n\n
E\nQ
"""

f2mtz_script = """
TITLE data from XDS 
FILE %(file_name_out)s
SYMMETRY    %(spgn_in)s
CELL    %(cell_str)s
LABOUT  H K L FP%(lbl)s SIGFP%(lbl)s %(cinp_ano)s %(free_lbl)s
CTYPOUT H H H  F   Q   %(cinp_ano2)s   %(free_code)s
NAME PROJECT %(ID)s CRYSTAL %(ID)s DATASET d%(ID)s
END
"""

f2mtz_phaser_script = """
TITLE data from XDS 
FILE %(file_name_out)s
SYMMETRY    %(spgn_in)s
CELL    %(cell_str)s
LABOUT  H K L FP%(lbl)s SIGFP%(lbl)s F(+)%(lbl)s SIGF(+)%(lbl)s F(-)%(lbl)s SIGF(-)%(lbl)s %(free_lbl)s
CTYPOUT H H H  F   Q    G     L      G     L %(free_code)s
END
"""

scala_script = """#!/bin/bash 

function run_scala() {
pointless XDSIN XDS_ASCII.HKL \\
          HKLOUT XDS_pointless_correct.mtz > XDS_pointless_correct.log
scala hklin XDS_pointless_correct.mtz hklout XDS_scala_correct.mtz \\
      scales     SCALA.scales \\
      rogues     SCALA.rogues \\
      rogueplot  SCALA.rogueplot \\
      correlplot SCALA.correlplot \\
      normplot   SCALA.norm \\
      anomplot   SCALA.anom \\
> XDS_scala.log << eof-1
cycles 0
bin 12
scales constant    # batch scaling is generally poorer than smoothed 
anomalous on 
eof-1
}

function run_aimless() {
pointless XDSIN XDS_ASCII.HKL \\
          HKLOUT XDS_pointless_correct.mtz > XDS_pointless_correct.log
aimless hklin XDS_pointless_correct.mtz hklout XDS_aimless_correct.mtz \\
      scales     AIMLESS.scales \\
      rogues     AIMLESS.rogues \\
      rogueplot  AIMLESS.rogueplot \\
      correlplot AIMLESS.correlplot \\
      normplot   AIMLESS.norm \\
      anomplot   AIMLESS.anom \\
> XDS_aimless.log << eof-2
cycles 0
bin 12
scales constant    # batch scaling is generally poorer than smoothed 
eof-2
}

#run_scala
run_aimless
"""

cad_crank_script = """
 cad HKLIN1 temp.mtz HKLOUT output_file_name.mtz<<EOF
 LABIN FILE 1 ALL
 XNAME FILE 1 ALL=XDS
 DWAVELENGTH FILE 1 XDS %(wavelength)s
 END
"""

f2mtz_sharp_script = """
TITLE data from XDS 
FILE %(file_name_out)s
SYMMETRY    %(spgn_in)s
CELL    %(cell_str)s
LABOUT  H K L FMID SMID DANO SANO ISYM
CTYPOUT H H H  F   Q    D     Q   Y
END
"""

solve_script = """#!/bin/sh
set noclobber

# run solve

solve_giant << eof-solve > solve.out
symfile %(spg_name)s.sym
cell  %(cell_str)s
resolution %(res_low).3f %(res_high).3f
logfile ./solve.logfile

readformatted
unmerged
read_intensities

checksolve
mad_atom %(ha_name)s

lambda 1 
label sad wavelength = %(wavelength)s
rawmadfile %(file_name_out)s
fixscattfactors
fprprv_mad %(fpp)s
nsolsite_deriv %(num_sites)d

# Add your known sites here
#xyz    0.3690   0.1194   0.0178
#xyz    0.9238   0.0377   0.0699
acceptance 0.30

!nres 200

#addsolve
SAD
eof-solve

# run resolve

resolve << eof-resolve > resolve.out
!solvent_content 0.30
!seq_file protein.seq
eof-resolve

fft HKLIN resolve.mtz MAPOUT resolve.ccp4map << eof-fft > resolve_fft.out 
TITLE    from resolve
LABIN   F1=FP PHI=PHIM W=FOMM
END
eof-fft
"""

phaser_script = """#!/bin/bash
# script written by xdsconv.py (pierre.legrand at synchrotron-soleil.fr)

label=""
scat="%(ha_name)s"
PARTIAL_MODEL_OPTION=""
solvent_content=0.5
parrot_cycles=5
parrot_resolution=1.0

function usage () {
  echo
  echo " $0 [-s x] ha_sites.pdb  [partial_model.pdb]"
  echo "    -l label, --label            mtz columns label"
  echo "    -a, --anom                   anomalous scatterer"
  echo "    -s solvc,   --solvent solvc  set the solvent content"
  echo "    -n ncycles, --parrot_cycles  number of parrot cycles"
  echo "    -r resol,   --resolution     parrot resolution cutoff"
  echo "    -h, --help                   print this help"
  echo
}

while [ $# -gt 0 ]; do
  case "$1" in
    -h | --help | -help )
      usage
      exit 0 ;;
    -l | --label )
      label="_$2"
      shift ; shift  ;;
    -a | --anom )
      scat=$2
      echo "INFO:  Using anomalous scatterer: $scat"
      shift; shift ;;
    -s | --solvent )
      solvent_content=$2
      echo "INFO:  Using a solvent fraction of: $solvent_content"
      shift; shift ;;
    -r | --resolution )
      parrot_resolution=$2
      echo "INFO:  Cuting high resolution for parrot to: $2"
      shift; shift ;;
    -n | --parrot_cycles )
      parrot_cycles=$2
      echo "INFO:  Number of parrot parrot cycles set to: $2"
      shift; shift ;;
    * )
      hatom_pdb=$1
      ID=`basename $hatom_pdb .pdb`
      break ;;
  esac
done

if [[ $# -eq 2 ]];then
   echo "Using partial model: $2"
   PARTIAL_MODEL_OPTION="PART PDB $2 ID 100"
   echo $PARTIAL_MODEL_OPTION
fi

function run_phaser() {
phaser << eof > phaser_auto_${ID}.log
MODE EP_AUTO
TITLe SAD phasing of ${ID} with %(num_sites)d ${scat}
HKLIn %(last_name)s
#COMPosition PROTein SEQ PROT.seq NUM 2
COMPosition BY SOLVent  
COMPosition PERCentage $solvent_content
CRYStal ${ID} DATAset sad LABIn F+=F(+)${label} SIG+=SIGF(+)${label} F-=F(-)${label} SIG-=SIGF(-)${label}
WAVElength %(wavelength)f
LLGComplete COMPLETE ON # CLASH 3.8
ATOM CRYStal ${ID} PDB $hatom_pdb SCATtering ${scat}
ATOM CHANge SCATterer ${scat}
ROOT ${ID}_auto
${PARTIAL_MODEL_OPTION}
eof
}

function run_parrot() {
hand=$1
if [ $hand = "ori" ] ; then 
parrID=${ID}_auto
elif [ $hand = "inv" ] ; then 
parrID=${ID}_auto.hand
fi
cparrot \\
-mtzin-wrk      ${parrID}.mtz \\
-pdbin-wrk-ha   ${parrID}.pdb \\
-colin-wrk-fo   "/*/*/[FP${label},SIGFP${label}]" \\
-colin-wrk-hl   "/*/*/[HLA,HLB,HLC,HLD]" \\
-colin-wrk-fc   "/*/*/[FWT,PHWT]" \\
-colin-wrk-free "/*/*/[FreeR_flag]" \\
-mtzout         ${parrID}_parrot_${solvent_content}_${parrot_cycles}.mtz \\
-colout         parrot \\
-solvent-flatten \\
-histogram-match \\
-cycles         ${parrot_cycles} \\
-resolution     ${parrot_resolution} \\
-solvent-content ${solvent_content} \\
-ncs-average \\
> cparrot_${parrID}_${solvent_content}_${parrot_cycles}.log
# -ncs-mask-filter-radius 22 \\
}

run_phaser
run_parrot ori
run_parrot inv

"""

replace_script = """#!/bin/bash

cat << eof > RF_SELF_3.inp
title ordinary self RF by slow RF from xdsme
!
print RF_SELF_3.prt
!
polar xzk
euler zyz
orthog axabz
!
acell %(cell_str)s
asymmetry %(spg_name)s
aobsfile data.hkl
acutoff 1.0 1.0 0.0
aformat 3i6, 2e10.3
! reading F power=2
! reading I power=1
apower 2
origin true
!
cutoff 0.25
!
resolution 40.0 3.5
boxsize 3 3 3
radius 20.0
geval 2
!
self true
cross false
fast true
!
sangle polar
oangle polar xyk
!
! This sets the search limits automatically
!
slimit 2 270 90 3
!
mapfile RF_SELF_3.map
peak 3.0 50
!
cntf RF_SELF_3.ps
cntl 3 9 0.5 1
cntl 1.5 2.5 0.5 4
!
stop
eof

glrf < RF_SELF_3.inp > RF_SELF.log

"""


crossec_script = """
ATOM %(ha_name)s
NWAV 1 %(wavelength)s
END
"""

cns_par = """a, b, c, alpha, beta, gamma = %(cns_cell)s
sg = %(cns_sg)s
low_res, high_res = %(res_low).3f %(res_high).3f
reflection_infile_1 = %(file_name_out)s
obs_f, obs_sigf, = "fobs", "sigma"
test_set, test_flag = "test", 1
"""

cad_script = """LABIN FILE 1 ALL
DWAVE FILE 1 %(ID)s d%(ID)s %(wavelength).5f\nEND"""
cad2_script = """LABIN FILE 1 ALL
LABIN FILE 2 %(cad_ano)s 
DWAVE FILE 2 %(ID)s d%(ID)s %(wavelength).5f\nEND"""

mtzutils_script = "END\n"

mtz2various_script = """#!/bin/bash
rm -f free_refl_shelx_F3.tmp
mtz2various hklin $1 hklout free_refl_shelx_F3.tmp << eof > mtz2shelx.log
LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag
OUTPUT SHELX
eof
grep -v "  0   0   0" free_refl_shelx_F3.tmp > free_refl_shelx_F3.hkl
echo "   0   0   0      0.    0.00   0" >> free_refl_shelx_F3.hkl
rm -f free_refl_shelx_F3.tmp
"""

HAd = {0.72227:("U ",92,"L3"),0.72766:("Y ",39,"K "),0.76973:("Sr",38,"K "),
0.81554:("Rb",37,"K "),0.86552:("Kr",36,"K "),0.92040:("Br",35,"K "),
0.92340:("Bi",83,"L3"),0.95073:("Pb",82,"L3"),0.97974:("Se",34,"K "),
1.00910:("Hg",80,"L3"),1.04000:("Au",79,"L3"),1.04500:("As",33,"K "),
1.07230:("Pt",78,"L3"),1.10580:("Ir",77,"L3"),1.11658:("Ge",32,"K "),
1.14080:("Os",76,"L3"),1.17730:("Re",75,"L3"),1.19580:("Ga",31,"K "),
1.21550:("W ",74,"L3"),1.25530:("Ta",73,"L3"),1.28340:("Zn",30,"K "),
1.29720:("Hf",72,"L3"),1.34050:("Lu",71,"L3"),1.38059:("Cu",29,"K "),
1.38620:("Yb",70,"L3"),1.43340:("Tm",69,"L3"),1.48350:("Er",68,"L3"),
1.48807:("Ni",28,"K "),1.53680:("Ho",67,"L3"),1.59160:("Dy",66,"L3"),
1.60815:("Co",27,"K "),1.64970:("Tb",65,"L3"),1.71170:("Gd",64,"L3"),
1.74346:("Fe",26,"K "),1.77610:("Eu",63,"L3"),1.84570:("Sm",62,"L3")}

# A dictionary containing all sg_number:sg_name from the spacegroup.lib
cns_sg_lib = {1:"P1",2:"P-1",3:"P2",4:"P2(1)",5:"C2",6:"PM",7:"PC",8:"CM",
9:"CC",10:"P2/M",11:"P2(1)/M",12:"C2/M",13:"P2/C",14:"P2(1)/C",15:"C2/C",
16:"P222",17:"P222(1)",18:"P2(1)2(1)2",19:"P2(1)2(1)2(1)",20:"C222(1)",
21:"C222",22:"F222",23:"I222",24:"I2(1)2(1)2(1)",25:"PMM2",26:"PMC2(1)",
27:"PCC2",28:"PMA2",29:"PCA2(1)",30:"PNC2",31:"PMN2(1)",32:"PBA2",
33:"PNA2(1)",34:"PNN2",35:"CMM2",36:"CMC2(1)",37:"CCC2",38:"AMM2",
39:"ABM2",40:"AMA2",41:"ABA2",42:"FMM2",43:"FDD2",44:"IMM2",45:"IBA2",
46:"IMA2",47:"PMMM",48:"PNNN",49:"PCCM",50:"PBAN",51:"PMMA",52:"PNNA",
53:"PMNA",54:"PCCA",55:"PBAM",56:"PCCN",57:"PBCM",58:"PNNM",59:"PMMN",
60:"PBCN",61:"PBCA",62:"PNMA",63:"CMCM",64:"CMCA",65:"CMMM",13065:"AMMM",
66:"CCCM",67:"CMMA",68:"CCCA",69:"FMMM",70:"FDDD",71:"IMMM",72:"IBAM",
73:"IBCA",74:"IMMA",75:"P4",76:"P4(1)",77:"P4(2)",78:"P4(3)",79:"I4",
80:"I4(1)",81:"P-4",82:"I-4",83:"P4/M",84:"P4(2)/M",85:"P4/N",86:"P4(2)/N",
87:"I4/M",88:"I4(1)/A",89:"P422",90:"P42(1)2",91:"P4(1)22",92:"P4(1)2(1)2",
93:"P4(2)22",94:"P4(2)2(1)2",95:"P4(3)22",96:"P4(3)2(1)2",97:"I422",
98:"I4(1)22",99:"P4MM",100:"P4BM",101:"P4(2)CM",102:"P4(2)NM",103:"P4CC",
104:"P4NC",105:"P4(2)MC",106:"P4(2)BC",107:"I4MM",108:"I4CM",109:"I4(1)MD",
110:"I4(1)CD",111:"P-42M",112:"P-42C",113:"P-42(1)M",114:"P-42(1)C",
115:"P-4M2",116:"P-4C2",117:"P-4B2",118:"P-4N2",119:"I-4M2",120:"I-4C2",
121:"I-42M",122:"I-42D",123:"P4/MMM",124:"P4/MCC",125:"P4/NBM",126:"P4/NNC",
127:"P4/MBM",128:"P4/MNC",129:"P4/NMM",130:"P4/NCC",131:"P4(2)/MMC",
132:"P4(2)/MCM",133:"P4(2)/NBC",134:"P4(2)/NNM",135:"P4(2)/MBC",
136:"P4(2)/MNM",137:"P4(2)/NMC",138:"P4(2)/NCM",139:"I4/MMM",140:"I4/MCM",
141:"I4(1)/AMD",142:"I4(1)/ACD",143:"P3",144:"P3(1)",145:"P3(2)",146:"R3",
20146:"R3R",147:"P-3",148:"R-3",20148:"R-3R",149:"P312",150:"P321",
151:"P3(1)12",152:"P3(1)21",153:"P3(2)12",154:"P3(2)21",155:"R32",
20155:"R32R",156:"P3M1",157:"P31M",158:"P3C1",159:"P31C",160:"R3M",
20160:"R3MR",161:"R3C",20161:"R3CR",162:"P-31M",163:"P-31C",164:"P-3M1",
165:"P-3C1",166:"R-3M",20166:"R-3MR",167:"R-3C",20167:"R-3CR",168:"P6",
169:"P6(1)",170:"P6(5)",171:"P6(2)",172:"P6(4)",173:"P6(3)",174:"P-6",
175:"P6/M",176:"P6(3)/M",177:"P622",178:"P6(1)22",179:"P6(5)22",
180:"P6(2)22",181:"P6(4)22",182:"P6(3)22",183:"P6MM",184:"P6CC",
185:"P6(3)CM",186:"P6(3)MC",187:"P-6M2",188:"P-6C2",189:"P-62M",
190:"P-62C",191:"P6/MMM",192:"P6/MCC",193:"P6(3)/MCM",194:"P6(3)/MMC",
195:"P23",196:"F23",197:"I23",198:"P2(1)3",199:"I2(1)3",200:"PM-3",
201:"PN-3",202:"FM-3",203:"FD-3",204:"IM-3",205:"PA-3",206:"IA-3",
207:"P432",208:"P4(2)32",209:"F432",210:"F4(1)32",211:"I432",
212:"P4(3)32",213:"P4(1)32",214:"I4(1)32",215:"P-43M",216:"F-43M",
217:"I-43M",218:"P-43N",219:"F-43C",220:"I-43D",221:"PM-3M",
222:"PN-3N",223:"PM-3N",224:"PN-3M",225:"FM-3M",226:"FM-3C",227:"FD-3M",
228:"FD-3C",229:"IM-3M",230:"IA-3D"}

amore_symops = {
1: ((1, 0, 'P1', 'PG1', 'TRICLINIC'),'x,y,z * end'),
2: ((2, 0, 'P-1', 'PG1BAR', 'TRICLINIC'),'x,y,z * -x,-y,-z * end'),
3: ((2, 0, 'P2', 'PG2', 'MONOCLINIC'),'x,y,z * -x,y,-z * end'),
4: ((2, 0, 'P21', 'PG2', 'MONOCLINIC'),'x,y,z * -x,1/2+y,-z * end'),
5: ((2, 1, 'C2', 'PG2', 'MONOCLINIC'),'x,y,z * -x,y,-z * 1/2,1/2,0 * end'),
6: ((2, 0, 'PM', 'PGM', 'MONOCLINIC'),'x,y,z * x,-y,z * end'),
7: ((2, 0, 'PC', 'PGM', 'MONOCLINIC'),'x,y,z * x,-y,1/2+z * end'),
8: ((2, 1, 'CM', 'PGM', 'MONOCLINIC'),'x,y,z * x,-y,z * 1/2,1/2,0 * end'),
9: ((2, 1, 'CC', 'PGM', 'MONOCLINIC'),'x,y,z * x,-y,1/2+z * 1/2,1/2,0 * end'),
10: ((4, 0, 'P2/M', 'PG2/M', 'MONOCLINIC'),'x,y,z * -x,y,-z * -x,-y,-z * x,-y,z * end'),
11: ((4, 0, 'P21/M', 'PG2/M', 'MONOCLINIC'),'x,y,z * -x,1/2+y,-z * -x,-y,-z * x,1/2-y,z * end'),
12: ((4, 1, 'C2/M', 'PG2/M', 'MONOCLINIC'),'x,y,z * -x,y,-z * -x,-y,-z * x,-y,z * 1/2,1/2,0 * end'),
13: ((4, 0, 'P2/C', 'PG2/M', 'MONOCLINIC'),'x,y,z * -x,y,1/2-z * -x,-y,-z * x,-y,1/2+z * end'),
14: ((4, 0, 'P21/C', 'PG2/M', 'MONOCLINIC'),'x,y,z * -x,-y,-z * -x,1/2+y,1/2-z * x,1/2-y,1/2+z * end'),
15: ((4, 1, 'C2/C', 'PG2/M', 'MONOCLINIC'),'x,y,z * -x,y,1/2-z * -x,-y,-z * x,-y,1/2+z * 1/2,1/2,0 * end'),
16: ((4, 0, 'P222', 'PG222', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * end'),
17: ((4, 0, 'P2221', 'PG222', 'ORTHORHOMBIC'),'x,y,z * -x,-y,1/2+z * -x,y,1/2-z * x,-y,-z * end'),
18: ((4, 0, 'P21212', 'PG222', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2-x,1/2+y,-z * 1/2+x,1/2-y,-z * end'),
19: ((4, 0, 'P212121', 'PG222', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * end'),
20: ((4, 1, 'C2221', 'PG222', 'ORTHORHOMBIC'),'x,y,z * -x,-y,1/2+z * -x,y,1/2-z * x,-y,-z * 1/2,1/2,0 * end'),
21: ((4, 1, 'C222', 'PG222', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * 1/2,1/2,0 * end'),
22: ((4, 3, 'F222', 'PG222', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
23: ((4, 1, 'I222', 'PG222', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * 1/2,1/2,1/2 * end'),
24: ((4, 1, 'I212121', 'PG222', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * 1/2,1/2,1/2 * end'),
25: ((4, 0, 'PMM2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * x,-y,z * -x,y,z * end'),
26: ((4, 0, 'PMC21', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,1/2+z * x,-y,1/2+z * -x,y,z * end'),
27: ((4, 0, 'PCC2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * x,-y,1/2+z * -x,y,1/2+z * end'),
28: ((4, 0, 'PMA2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2+x,-y,z * 1/2-x,y,z * end'),
29: ((4, 0, 'PCA21', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,1/2+z * 1/2+x,-y,z * 1/2-x,y,1/2+z * end'),
30: ((4, 0, 'PNC2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * x,1/2-y,1/2+z * -x,1/2+y,1/2+z * end'),
31: ((4, 0, 'PMN21', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,-y,1/2+z * 1/2+x,-y,1/2+z * -x,y,z * end'),
32: ((4, 0, 'PBA2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * end'),
33: ((4, 0, 'PNA21', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,1/2+z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,1/2+z * end'),
34: ((4, 0, 'PNN2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * end'),
35: ((4, 1, 'CMM2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * x,-y,z * -x,y,z * 1/2,1/2,0 * end'),
36: ((4, 1, 'CMC21', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,1/2+z * x,-y,1/2+z * -x,y,z * 1/2,1/2,0 * end'),
37: ((4, 1, 'CCC2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * x,-y,1/2+z * -x,y,1/2+z * 1/2,1/2,0 * end'),
38: ((4, 1, 'AMM2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * x,-y,z * -x,y,z * 0,1/2,1/2 * end'),
39: ((4, 1, 'ABM2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * x,1/2-y,z * -x,1/2+y,z * 0,1/2,1/2 * end'),
40: ((4, 1, 'AMA2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2+x,-y,z * 1/2-x,y,z * 0,1/2,1/2 * end'),
41: ((4, 1, 'ABA2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * 0,1/2,1/2 * end'),
42: ((4, 3, 'FMM2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * x,-y,z * -x,y,z * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
43: ((4, 3, 'FDD2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/4+x,1/4-y,1/4+z * 1/4-x,1/4+y,1/4+z * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
44: ((4, 1, 'IMM2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * x,-y,z * -x,y,z * 1/2,1/2,1/2 * end'),
45: ((4, 1, 'IBA2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * 1/2,1/2,1/2 * end'),
46: ((4, 1, 'IMA2', 'PGMM2', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2+x,-y,z * 1/2-x,y,z * 1/2,1/2,1/2 * end'),
47: ((8, 0, 'PMMM', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * end'),
48: ((8, 0, 'PNNN', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * 1/2-x,1/2-y,1/2-z * 1/2+x,1/2+y,1/2-z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * end'),
49: ((8, 0, 'PCCM', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,1/2-z * x,-y,1/2-z * -x,-y,-z * x,y,-z * x,-y,1/2+z * -x,y,1/2+z * end'),
50: ((8, 0, 'PBAN', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * 1/2-x,1/2-y,-z * 1/2+x,1/2+y,-z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * end'),
51: ((8, 0, 'PMMA', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,-y,z * -x,y,-z * 1/2+x,-y,-z * -x,-y,-z * 1/2+x,y,-z * x,-y,z * 1/2-x,y,z * end'),
52: ((8, 0, 'PNNA', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,-y,z * 1/2-x,1/2+y,1/2-z * x,1/2-y,1/2-z * -x,-y,-z * 1/2+x,y,-z * 1/2+x,1/2-y,1/2+z * -x,1/2+y,1/2+z * end'),
53: ((8, 0, 'PMNA', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,-y,1/2+z * 1/2-x,y,1/2-z * x,-y,-z * -x,-y,-z * 1/2+x,y,1/2-z * 1/2+x,-y,1/2+z * -x,y,z * end'),
54: ((8, 0, 'PCCA', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,-y,z * -x,y,1/2-z * 1/2+x,-y,1/2-z * -x,-y,-z * 1/2+x,y,-z * x,-y,1/2+z * 1/2-x,y,1/2+z * end'),
55: ((8, 0, 'PBAM', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2-x,1/2+y,-z * 1/2+x,1/2-y,-z * -x,-y,-z * x,y,-z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * end'),
56: ((8, 0, 'PCCN', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,1/2-y,z * -x,1/2+y,1/2-z * 1/2+x,-y,1/2-z * -x,-y,-z * 1/2+x,1/2+y,-z * x,1/2-y,1/2+z * 1/2-x,y,1/2+z * end'),
57: ((8, 0, 'PBCM', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,1/2+z * -x,1/2+y,1/2-z * x,1/2-y,-z * -x,-y,-z * x,y,1/2-z * x,1/2-y,1/2+z * -x,1/2+y,z * end'),
58: ((8, 0, 'PNNM', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2-x,1/2+y,1/2-z * 1/2+x,1/2-y,1/2-z * -x,-y,-z * x,y,-z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * end'),
59: ((8, 0, 'PMMN', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2-x,1/2+y,-z * 1/2+x,1/2-y,-z * 1/2-x,1/2-y,-z * 1/2+x,1/2+y,-z * x,-y,z * -x,y,z * end'),
60: ((8, 0, 'PBCN', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,1/2-y,1/2+z * -x,y,1/2-z * 1/2+x,1/2-y,-z * -x,-y,-z * 1/2+x,1/2+y,1/2-z * x,-y,1/2+z * 1/2-x,1/2+y,z * end'),
61: ((8, 0, 'PBCA', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * -x,-y,-z * 1/2+x,y,1/2-z * x,1/2-y,1/2+z * 1/2-x,1/2+y,z * end'),
62: ((8, 0, 'PNMA', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,-z * 1/2+x,1/2-y,1/2-z * -x,-y,-z * 1/2+x,y,1/2-z * x,1/2-y,z * 1/2-x,1/2+y,1/2+z * end'),
63: ((8, 1, 'CMCM', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,1/2+z * -x,y,1/2-z * x,-y,-z * -x,-y,-z * x,y,1/2-z * x,-y,1/2+z * -x,y,z * 1/2,1/2,0 * end'),
64: ((8, 1, 'CMCA', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,1/2-y,1/2+z * -x,1/2+y,1/2-z * x,-y,-z * -x,-y,-z * x,1/2+y,1/2-z * x,1/2-y,1/2+z * -x,y,z * 1/2,1/2,0 * end'),
65: ((8, 1, 'CMMM', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * 1/2,1/2,0 * end'),
66: ((8, 1, 'CCCM', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,1/2-z * x,-y,1/2-z * -x,-y,-z * x,y,-z * x,-y,1/2+z * -x,y,1/2+z * 1/2,1/2,0 * end'),
67: ((8, 1, 'CMMA', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,1/2-y,z * -x,1/2+y,-z * x,-y,-z * -x,-y,-z * x,1/2+y,-z * x,1/2-y,z * -x,y,z * 1/2,1/2,0 * end'),
68: ((8, 1, 'CCCA', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,1/2-y,z * -x,y,-z * 1/2+x,1/2-y,-z * -x,1/2-y,1/2-z * 1/2+x,y,1/2-z * x,1/2-y,1/2+z * 1/2-x,y,1/2+z * 1/2,1/2,0 * end'),
69: ((8, 3, 'FMMM', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
70: ((8, 3, 'FDDD', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * 1/4-x,1/4-y,1/4-z * 1/4+x,1/4+y,1/4-z * 1/4+x,1/4-y,1/4+z * 1/4-x,1/4+y,1/4+z * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
71: ((8, 1, 'IMMM', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * 1/2,1/2,1/2 * end'),
72: ((8, 1, 'IBAM', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,-y,z * 1/2-x,1/2+y,-z * 1/2+x,1/2-y,-z * -x,-y,-z * x,y,-z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * 1/2,1/2,1/2 * end'),
73: ((8, 1, 'IBCA', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * -x,-y,-z * 1/2+x,y,1/2-z * x,1/2-y,1/2+z * 1/2-x,1/2+y,z * 1/2,1/2,1/2 * end'),
74: ((8, 1, 'IMMA', 'PGMMM', 'ORTHORHOMBIC'),'x,y,z * -x,1/2-y,z * -x,1/2+y,-z * x,-y,-z * -x,-y,-z * x,1/2+y,-z * x,1/2-y,z * -x,y,z * 1/2,1/2,1/2 * end'),
75: ((4, 0, 'P4', 'PG4', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * end'),
76: ((4, 0, 'P41', 'PG4', 'TETRAGONAL'),'x,y,z * -x,-y,1/2+z * -y,x,1/4+z * y,-x,3/4+z * end'),
77: ((4, 0, 'P42', 'PG4', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,1/2+z * y,-x,1/2+z * end'),
78: ((4, 0, 'P43', 'PG4', 'TETRAGONAL'),'x,y,z * -x,-y,1/2+z * -y,x,3/4+z * y,-x,1/4+z * end'),
79: ((4, 1, 'I4', 'PG4', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * 1/2,1/2,1/2 * end'),
80: ((4, 1, 'I41', 'PG4', 'TETRAGONAL'),'x,y,z * 1/2-x,1/2-y,1/2+z * -y,1/2+x,1/4+z * 1/2+y,-x,3/4+z * 1/2,1/2,1/2 * end'),
81: ((4, 0, 'P-4', 'PG4BAR', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * end'),
82: ((4, 1, 'I-4', 'PG4BAR', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * 1/2,1/2,1/2 * end'),
83: ((8, 0, 'P4/M', 'PG4/M', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * -x,-y,-z * x,y,-z * y,-x,-z * -y,x,-z * end'),
84: ((8, 0, 'P42/M', 'PG4/M', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,1/2+z * y,-x,1/2+z * -x,-y,-z * x,y,-z * y,-x,1/2-z * -y,x,1/2-z * end'),
85: ((8, 0, 'P4/N', 'PG4/M', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,z * 1/2+y,1/2-x,z * 1/2-x,1/2-y,-z * 1/2+x,1/2+y,-z * y,-x,-z * -y,x,-z * end'),
86: ((8, 0, 'P42/N', 'PG4/M', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,1/2+z * 1/2+y,1/2-x,1/2+z * 1/2-x,1/2-y,1/2-z * 1/2+x,1/2+y,1/2-z * y,-x,-z * -y,x,-z * end'),
87: ((8, 1, 'I4/M', 'PG4/M', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * -x,-y,-z * x,y,-z * y,-x,-z * -y,x,-z * 1/2,1/2,1/2 * end'),
88: ((8, 1, 'I41/A', 'PG4/M', 'TETRAGONAL'),'x,y,z * 1/2-x,1/2-y,1/2+z * -y,1/2+x,1/4+z * 1/2+y,-x,3/4+z * -x,1/2-y,1/4-z * 1/2+x,y,3/4-z * y,-x,-z * 1/2-y,1/2+x,1/2-z * 1/2,1/2,1/2 * end'),
89: ((8, 0, 'P422', 'PG422', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * -x,y,-z * x,-y,-z * y,x,-z * -y,-x,-z * end'),
90: ((8, 0, 'P4212', 'PG422', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,z * 1/2+y,1/2-x,z * 1/2-x,1/2+y,-z * 1/2+x,1/2-y,-z * y,x,-z * -y,-x,-z * end'),
91: ((8, 0, 'P4122', 'PG422', 'TETRAGONAL'),'x,y,z * -x,-y,1/2+z * -y,x,1/4+z * y,-x,3/4+z * -x,y,-z * x,-y,1/2-z * y,x,3/4-z * -y,-x,1/4-z * end'),
92: ((8, 0, 'P41212', 'PG422', 'TETRAGONAL'),'x,y,z * -x,-y,1/2+z * 1/2-y,1/2+x,1/4+z * 1/2+y,1/2-x,3/4+z * 1/2-x,1/2+y,1/4-z * 1/2+x,1/2-y,3/4-z * y,x,-z * -y,-x,1/2-z * end'),
93: ((8, 0, 'P4222', 'PG422', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,1/2+z * y,-x,1/2+z * -x,y,-z * x,-y,-z * y,x,1/2-z * -y,-x,1/2-z * end'),
94: ((8, 0, 'P42212', 'PG422', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,1/2+z * 1/2+y,1/2-x,1/2+z * 1/2-x,1/2+y,1/2-z * 1/2+x,1/2-y,1/2-z * y,x,-z * -y,-x,-z * end'),
95: ((8, 0, 'P4322', 'PG422', 'TETRAGONAL'),'x,y,z * -x,-y,1/2+z * -y,x,3/4+z * y,-x,1/4+z * -x,y,-z * x,-y,1/2-z * y,x,1/4-z * -y,-x,3/4-z * end'),
96: ((8, 0, 'P43212', 'PG422', 'TETRAGONAL'),'x,y,z * -x,-y,1/2+z * 1/2-y,1/2+x,3/4+z * 1/2+y,1/2-x,1/4+z * 1/2-x,1/2+y,3/4-z * 1/2+x,1/2-y,1/4-z * y,x,-z * -y,-x,1/2-z * end'),
97: ((8, 1, 'I422', 'PG422', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * -x,y,-z * x,-y,-z * y,x,-z * -y,-x,-z * 1/2,1/2,1/2 * end'),
98: ((8, 1, 'I4122', 'PG422', 'TETRAGONAL'),'x,y,z * 1/2-x,1/2-y,1/2+z * -y,1/2+x,1/4+z * 1/2+y,-x,3/4+z * 1/2-x,y,3/4-z * x,1/2-y,1/4-z * 1/2+y,1/2+x,1/2-z * -y,-x,-z * 1/2,1/2,1/2 * end'),
99: ((8, 0, 'P4MM', 'PG4MM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * x,-y,z * -x,y,z * -y,-x,z * y,x,z * end'),
100: ((8, 0, 'P4BM', 'PG4MM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * 1/2-y,1/2-x,z * 1/2+y,1/2+x,z * end'),
101: ((8, 0, 'P42CM', 'PG4MM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,1/2+z * y,-x,1/2+z * x,-y,1/2+z * -x,y,1/2+z * -y,-x,z * y,x,z * end'),
102: ((8, 0, 'P42NM', 'PG4MM', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,1/2+z * 1/2+y,1/2-x,1/2+z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * -y,-x,z * y,x,z * end'),
103: ((8, 0, 'P4CC', 'PG4MM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * x,-y,1/2+z * -x,y,1/2+z * -y,-x,1/2+z * y,x,1/2+z * end'),
104: ((8, 0, 'P4NC', 'PG4MM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2+x,1/2+z * end'),
105: ((8, 0, 'P42MC', 'PG4MM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,1/2+z * y,-x,1/2+z * x,-y,z * -x,y,z * -y,-x,1/2+z * y,x,1/2+z * end'),
106: ((8, 0, 'P42BC', 'PG4MM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,1/2+z * y,-x,1/2+z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2+x,1/2+z * end'),
107: ((8, 1, 'I4MM', 'PG4MM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * x,-y,z * -x,y,z * -y,-x,z * y,x,z * 1/2,1/2,1/2 * end'),
108: ((8, 1, 'I4CM', 'PG4MM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * x,-y,1/2+z * -x,y,1/2+z * -y,-x,1/2+z * y,x,1/2+z * 1/2,1/2,1/2 * end'),
109: ((8, 1, 'I41MD', 'PG4MM', 'TETRAGONAL'),'x,y,z * 1/2-x,1/2-y,1/2+z * -y,1/2+x,1/4+z * 1/2+y,-x,3/4+z * x,-y,z * 1/2-x,1/2+y,1/2+z * -y,1/2-x,1/4+z * 1/2+y,x,3/4+z * 1/2,1/2,1/2 * end'),
110: ((8, 1, 'I41CD', 'PG4MM', 'TETRAGONAL'),'x,y,z * 1/2-x,1/2-y,1/2+z * -y,1/2+x,1/4+z * 1/2+y,-x,3/4+z * x,-y,1/2+z * 1/2-x,1/2+y,z * -y,1/2-x,3/4+z * 1/2+y,x,1/4+z * 1/2,1/2,1/2 * end'),
111: ((8, 0, 'P-42M', 'PG4BAR2M', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * -x,y,-z * x,-y,-z * -y,-x,z * y,x,z * end'),
112: ((8, 0, 'P-42C', 'PG4BAR2M', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * -x,y,1/2-z * x,-y,1/2-z * -y,-x,1/2+z * y,x,1/2+z * end'),
113: ((8, 0, 'P-421M', 'PG4BAR2M', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * 1/2-x,1/2+y,-z * 1/2+x,1/2-y,-z * 1/2-y,1/2-x,z * 1/2+y,1/2+x,z * end'),
114: ((8, 0, 'P-421C', 'PG4BAR2M', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * 1/2-x,1/2+y,1/2-z * 1/2+x,1/2-y,1/2-z * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2+x,1/2+z * end'),
115: ((8, 0, 'P-4M2', 'PG4BARM2', 'TETRAGONAL'),'x,y,z * -x,-y,z * y,-x,-z * -y,x,-z * x,-y,z * -x,y,z * y,x,-z * -y,-x,-z * end'),
116: ((8, 0, 'P-4C2', 'PG4BARM2', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * x,-y,1/2+z * -x,y,1/2+z * y,x,1/2-z * -y,-x,1/2-z * end'),
117: ((8, 0, 'P-4B2', 'PG4BARM2', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * 1/2+y,1/2+x,-z * 1/2-y,1/2-x,-z * end'),
118: ((8, 0, 'P-4N2', 'PG4BARM2', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * 1/2+y,1/2+x,1/2-z * 1/2-y,1/2-x,1/2-z * end'),
119: ((8, 1, 'I-4M2', 'PG4BARM2', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * x,-y,z * -x,y,z * y,x,-z * -y,-x,-z * 1/2,1/2,1/2 * end'),
120: ((8, 1, 'I-4C2', 'PG4BARM2', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * x,-y,1/2+z * -x,y,1/2+z * y,x,1/2-z * -y,-x,1/2-z * 1/2,1/2,1/2 * end'),
121: ((8, 1, 'I-42M', 'PG4BAR2M', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * -x,y,-z * x,-y,-z * -y,-x,z * y,x,z * 1/2,1/2,1/2 * end'),
122: ((8, 1, 'I-42D', 'PG4BAR2M', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,-z * y,-x,-z * 1/2-x,y,3/4-z * 1/2+x,-y,3/4-z * 1/2-y,-x,3/4+z * 1/2+y,x,3/4+z * 1/2,1/2,1/2 * end'),
123: ((16, 0, 'P4/MMM', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * -x,y,-z * x,-y,-z * y,x,-z * -y,-x,-z * -x,-y,-z * x,y,-z * y,-x,-z * -y,x,-z * x,-y,z * -x,y,z * -y,-x,z * y,x,z * end'),
124: ((16, 0, 'P4/MCC', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * -x,y,1/2-z * x,-y,1/2-z * y,x,1/2-z * -y,-x,1/2-z * -x,-y,-z * x,y,-z * y,-x,-z * -y,x,-z * x,-y,1/2+z * -x,y,1/2+z * -y,-x,1/2+z * y,x,1/2+z * end'),
125: ((16, 0, 'P4/NBM', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * -x,y,-z * x,-y,-z * y,x,-z * -y,-x,-z * 1/2-x,1/2-y,-z * 1/2+x,1/2+y,-z * 1/2+y,1/2-x,-z * 1/2-y,1/2+x,-z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * 1/2-y,1/2-x,z * 1/2+y,1/2+x,z * end'),
126: ((16, 0, 'P4/NNC', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * -x,y,-z * x,-y,-z * y,x,-z * -y,-x,-z * 1/2-x,1/2-y,1/2-z * 1/2+x,1/2+y,1/2-z * 1/2+y,1/2-x,1/2-z * 1/2-y,1/2+x,1/2-z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2+x,1/2+z * end'),
127: ((16, 0, 'P4/MBM', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * 1/2-x,1/2+y,-z * 1/2+x,1/2-y,-z * 1/2+y,1/2+x,-z * 1/2-y,1/2-x,-z * -x,-y,-z * x,y,-z * y,-x,-z * -y,x,-z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * 1/2-y,1/2-x,z * 1/2+y,1/2+x,z * end'),
128: ((16, 0, 'P4/MNC', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * 1/2-x,1/2+y,1/2-z * 1/2+x,1/2-y,1/2-z * 1/2+y,1/2+x,1/2-z * 1/2-y,1/2-x,1/2-z * -x,-y,-z * x,y,-z * y,-x,-z * -y,x,-z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2+x,1/2+z * end'),
129: ((16, 0, 'P4/NMM', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,z * 1/2+y,1/2-x,z * 1/2-x,1/2+y,-z * 1/2+x,1/2-y,-z * y,x,-z * -y,-x,-z * 1/2-x,1/2-y,-z * 1/2+x,1/2+y,-z * y,-x,-z * -y,x,-z * x,-y,z * -x,y,z * 1/2-y,1/2-x,z * 1/2+y,1/2+x,z * end'),
130: ((16, 0, 'P4/NCC', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,z * 1/2+y,1/2-x,z * 1/2-x,1/2+y,1/2-z * 1/2+x,1/2-y,1/2-z * y,x,1/2-z * -y,-x,1/2-z * 1/2-x,1/2-y,-z * 1/2+x,1/2+y,-z * y,-x,-z * -y,x,-z * x,-y,1/2+z * -x,y,1/2+z * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2+x,1/2+z * end'),
131: ((16, 0, 'P42/MMC', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,1/2+z * y,-x,1/2+z * -x,y,-z * x,-y,-z * y,x,1/2-z * -y,-x,1/2-z * -x,-y,-z * x,y,-z * y,-x,1/2-z * -y,x,1/2-z * x,-y,z * -x,y,z * -y,-x,1/2+z * y,x,1/2+z * end'),
132: ((16, 0, 'P42/MCM', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,1/2+z * y,-x,1/2+z * -x,y,1/2-z * x,-y,1/2-z * y,x,-z * -y,-x,-z * -x,-y,-z * x,y,-z * y,-x,1/2-z * -y,x,1/2-z * x,-y,1/2+z * -x,y,1/2+z * -y,-x,z * y,x,z * end'),
133: ((16, 0, 'P42/NBC', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,1/2+z * 1/2+y,1/2-x,1/2+z * -x,y,1/2-z * x,-y,1/2-z * 1/2+y,1/2+x,-z * 1/2-y,1/2-x,-z * 1/2-x,1/2-y,1/2-z * 1/2+x,1/2+y,1/2-z * y,-x,-z * -y,x,-z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * -y,-x,1/2+z * y,x,1/2+z * end'),
134: ((16, 0, 'P42/NNM', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,1/2+z * 1/2+y,1/2-x,1/2+z * -x,y,-z * x,-y,-z * 1/2+y,1/2+x,1/2-z * 1/2-y,1/2-x,1/2-z * 1/2-x,1/2-y,1/2-z * 1/2+x,1/2+y,1/2-z * y,-x,-z * -y,x,-z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * -y,-x,z * y,x,z * end'),
135: ((16, 0, 'P42/MBC', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,1/2+z * y,-x,1/2+z * 1/2-x,1/2+y,-z * 1/2+x,1/2-y,-z * 1/2+y,1/2+x,1/2-z * 1/2-y,1/2-x,1/2-z * -x,-y,-z * x,y,-z * y,-x,1/2-z * -y,x,1/2-z * 1/2+x,1/2-y,z * 1/2-x,1/2+y,z * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2+x,1/2+z * end'),
136: ((16, 0, 'P42/MNM', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,1/2+z * 1/2+y,1/2-x,1/2+z * 1/2-x,1/2+y,1/2-z * 1/2+x,1/2-y,1/2-z * y,x,-z * -y,-x,-z * -x,-y,-z * x,y,-z * 1/2+y,1/2-x,1/2-z * 1/2-y,1/2+x,1/2-z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * -y,-x,z * y,x,z * end'),
137: ((16, 0, 'P42/NMC', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,1/2+z * 1/2+y,1/2-x,1/2+z * 1/2-x,1/2+y,1/2-z * 1/2+x,1/2-y,1/2-z * y,x,-z * -y,-x,-z * 1/2-x,1/2-y,1/2-z * 1/2+x,1/2+y,1/2-z * y,-x,-z * -y,x,-z * x,-y,z * -x,y,z * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2+x,1/2+z * end'),
138: ((16, 0, 'P42/NCM', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * 1/2-y,1/2+x,1/2+z * 1/2+y,1/2-x,1/2+z * 1/2-x,1/2+y,-z * 1/2+x,1/2-y,-z * y,x,1/2-z * -y,-x,1/2-z * 1/2-x,1/2-y,1/2-z * 1/2+x,1/2+y,1/2-z * y,-x,-z * -y,x,-z * x,-y,1/2+z * -x,y,1/2+z * 1/2-y,1/2-x,z * 1/2+y,1/2+x,z * end'),
139: ((16, 1, 'I4/MMM', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * -x,y,-z * x,-y,-z * y,x,-z * -y,-x,-z * -x,-y,-z * x,y,-z * y,-x,-z * -y,x,-z * x,-y,z * -x,y,z * -y,-x,z * y,x,z * 1/2,1/2,1/2 * end'),
140: ((16, 1, 'I4/MCM', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * -x,-y,z * -y,x,z * y,-x,z * -x,y,1/2-z * x,-y,1/2-z * y,x,1/2-z * -y,-x,1/2-z * -x,-y,-z * x,y,-z * y,-x,-z * -y,x,-z * x,-y,1/2+z * -x,y,1/2+z * -y,-x,1/2+z * y,x,1/2+z * 1/2,1/2,1/2 * end'),
141: ((16, 1, 'I41/AMD', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * 1/2-x,1/2-y,1/2+z * -y,1/2+x,1/4+z * 1/2+y,-x,3/4+z * 1/2-x,y,3/4-z * x,1/2-y,1/4-z * 1/2+y,1/2+x,1/2-z * -y,-x,-z * -x,1/2-y,1/4-z * 1/2+x,y,3/4-z * y,-x,-z * 1/2-y,1/2+x,1/2-z * 1/2+x,1/2-y,1/2+z * -x,y,z * 1/2-y,-x,3/4+z * y,1/2+x,1/4+z * 1/2,1/2,1/2 * end'),
142: ((16, 1, 'I41/ACD', 'PG4/MMM', 'TETRAGONAL'),'x,y,z * 1/2-x,1/2-y,1/2+z * -y,1/2+x,1/4+z * 1/2+y,-x,3/4+z * 1/2-x,y,1/4-z * x,1/2-y,3/4-z * 1/2+y,1/2+x,-z * -y,-x,1/2-z * -x,1/2-y,1/4-z * 1/2+x,y,3/4-z * y,-x,-z * 1/2-y,1/2+x,1/2-z * 1/2+x,1/2-y,z * -x,y,1/2+z * 1/2-y,-x,1/4+z * y,1/2+x,3/4+z * 1/2,1/2,1/2 * end'),
143: ((3, 0, 'P3', 'PG3', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * end'),
144: ((3, 0, 'P31', 'PG3', 'TRIGONAL'),'x,y,z * -y,x-y,1/3+z * y-x,-x,2/3+z * end'),
145: ((3, 0, 'P32', 'PG3', 'TRIGONAL'),'x,y,z * -y,x-y,2/3+z * y-x,-x,1/3+z * end'),
146: ((3, 2, 'R3', 'PG3', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * 2/3,1/3,1/3 * 1/3,2/3,2/3 * end'),
146.1: ((3, 0, 'R3_R', 'PG3', 'RHOMBOHEDRAL'),'x,y,z * z,x,y * y,z,x * end'),
147: ((6, 0, 'P-3', 'PG3BAR', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,-z * y,y-x,-z * x-y,x,-z * end'),
148: ((6, 2, 'R-3', 'PG3BAR', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,-z * y,y-x,-z * x-y,x,-z * 2/3,1/3,1/3 * 1/3,2/3,2/3 * end'),
148.1: ((6, 0, 'R-3_R', 'PG3BAR', 'RHOMBOHEDRAL'),'x,y,z * z,x,y * y,z,x * -x,-y,-z * -z,-x,-y * -y,-z,-x * end'),
149: ((6, 0, 'P312', 'PG312', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -y,-x,-z * y-x,y,-z * x,x-y,-z * end'),
150: ((6, 0, 'P321', 'PG321', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * y,x,-z * x-y,-y,-z * -x,y-x,-z * end'),
151: ((6, 0, 'P3112', 'PG312', 'TRIGONAL'),'x,y,z * -y,x-y,1/3+z * y-x,-x,2/3+z * -y,-x,2/3-z * y-x,y,1/3-z * x,x-y,-z * end'),
152: ((6, 0, 'P3121', 'PG321', 'TRIGONAL'),'x,y,z * -y,x-y,1/3+z * y-x,-x,2/3+z * y,x,-z * x-y,-y,2/3-z * -x,y-x,1/3-z * end'),
153: ((6, 0, 'P3212', 'PG312', 'TRIGONAL'),'x,y,z * -y,x-y,2/3+z * y-x,-x,1/3+z * -y,-x,1/3-z * y-x,y,2/3-z * x,x-y,-z * end'),
154: ((6, 0, 'P3221', 'PG321', 'TRIGONAL'),'x,y,z * -y,x-y,2/3+z * y-x,-x,1/3+z * y,x,-z * x-y,-y,1/3-z * -x,y-x,2/3-z * end'),
155: ((6, 2, 'R32', 'PG32', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * y,x,-z * x-y,-y,-z * -x,y-x,-z * 2/3,1/3,1/3 * 1/3,2/3,2/3 * end'),
155.1: ((6, 0, 'R32_R', 'PG32', 'RHOMBOHEDRAL'),'x,y,z * z,x,y * y,z,x * -y,-x,-z * -x,-z,-y * -z,-y,-x * end'),
156: ((6, 0, 'P3M1', 'PG3M1', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -y,-x,z * y-x,y,z * x,x-y,z * end'),
157: ((6, 0, 'P31M', 'PG31M', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * y,x,z * x-y,-y,z * -x,y-x,z * end'),
158: ((6, 0, 'P3C1', 'PG3M1', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -y,-x,1/2+z * y-x,y,1/2+z * x,x-y,1/2+z * end'),
159: ((6, 0, 'P31C', 'PG31M', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * y,x,1/2+z * x-y,-y,1/2+z * -x,y-x,1/2+z * end'),
160: ((6, 2, 'R3M', 'PG3M', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -y,-x,z * y-x,y,z * x,x-y,z * 2/3,1/3,1/3 * 1/3,2/3,2/3 * end'),
160.1: ((6, 0, 'R3M_R', 'PG3M', 'RHOMBOHEDRAL'),'x,y,z * z,x,y * y,z,x * y,x,z * x,z,y * z,y,x * end'),
161: ((6, 2, 'R3C', 'PG3M', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -y,-x,1/2+z * y-x,y,1/2+z * x,x-y,1/2+z * 2/3,1/3,1/3 * 1/3,2/3,2/3 * end'),
161.1: ((6, 0, 'R3C_R', 'PG3M', 'RHOMBOHEDRAL'),'x,y,z * z,x,y * y,z,x * 1/2+y,1/2+x,1/2+z * 1/2+x,1/2+z,1/2+y * 1/2+z,1/2+y,1/2+x * end'),
162: ((12, 0, 'P-31M', 'PG3BAR1M', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -y,-x,-z * y-x,y,-z * x,x-y,-z * -x,-y,-z * y,y-x,-z * x-y,x,-z * y,x,z * x-y,-y,z * -x,y-x,z * end'),
163: ((12, 0, 'P-31C', 'PG3BAR1M', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -y,-x,1/2-z * y-x,y,1/2-z * x,x-y,1/2-z * -x,-y,-z * y,y-x,-z * x-y,x,-z * y,x,1/2+z * x-y,-y,1/2+z * -x,y-x,1/2+z * end'),
164: ((12, 0, 'P-3M1', 'PG3BARM1', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * y,x,-z * x-y,-y,-z * -x,y-x,-z * -x,-y,-z * y,y-x,-z * x-y,x,-z * -y,-x,z * y-x,y,z * x,x-y,z * end'),
165: ((12, 0, 'P-3C1', 'PG3BARM1', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * y,x,1/2-z * x-y,-y,1/2-z * -x,y-x,1/2-z * -x,-y,-z * y,y-x,-z * x-y,x,-z * -y,-x,1/2+z * y-x,y,1/2+z * x,x-y,1/2+z * end'),
166: ((12, 2, 'R-3M', 'PG3BARM', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * y,x,-z * x-y,-y,-z * -x,y-x,-z * -x,-y,-z * y,y-x,-z * x-y,x,-z * -y,-x,z * y-x,y,z * x,x-y,z * 2/3,1/3,1/3 * 1/3,2/3,2/3 * end'),
166.1: ((12, 0, 'R-3M_R', 'PG3BARM', 'RHOMBOHEDRAL'),'x,y,z * z,x,y * y,z,x * -y,-x,-z * -x,-z,-y * -z,-y,-x * -x,-y,-z * -z,-x,-y * -y,-z,-x * y,x,z * x,z,y * z,y,x * end'),
167: ((12, 2, 'R-3C', 'PG3BARM', 'TRIGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * y,x,1/2-z * x-y,-y,1/2-z * -x,y-x,1/2-z * -x,-y,-z * y,y-x,-z * x-y,x,-z * -y,-x,1/2+z * y-x,y,1/2+z * x,x-y,1/2+z * 2/3,1/3,1/3 * 1/3,2/3,2/3 * end'),
167.1: ((12, 0, 'R-3C_R', 'PG3BARM', 'RHOMBOHEDRAL'),'x,y,z * z,x,y * y,z,x * 1/2-y,1/2-x,1/2-z * 1/2-x,1/2-z,1/2-y * 1/2-z,1/2-y,1/2-x * -x,-y,-z * -z,-x,-y * -y,-z,-x * 1/2+y,1/2+x,1/2+z * 1/2+x,1/2+z,1/2+y * 1/2+z,1/2+y,1/2+x * end'),
168: ((6, 0, 'P6', 'PG6', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,z * y,y-x,z * x-y,x,z * end'),
169: ((6, 0, 'P61', 'PG6', 'HEXAGONAL'),'x,y,z * -y,x-y,1/3+z * y-x,-x,2/3+z * -x,-y,1/2+z * y,y-x,5/6+z * x-y,x,1/6+z * end'),
170: ((6, 0, 'P65', 'PG6', 'HEXAGONAL'),'x,y,z * -y,x-y,2/3+z * y-x,-x,1/3+z * -x,-y,1/2+z * y,y-x,1/6+z * x-y,x,5/6+z * end'),
171: ((6, 0, 'P62', 'PG6', 'HEXAGONAL'),'x,y,z * -y,x-y,2/3+z * y-x,-x,1/3+z * -x,-y,z * y,y-x,2/3+z * x-y,x,1/3+z * end'),
172: ((6, 0, 'P64', 'PG6', 'HEXAGONAL'),'x,y,z * -y,x-y,1/3+z * y-x,-x,2/3+z * -x,-y,z * y,y-x,1/3+z * x-y,x,2/3+z * end'),
173: ((6, 0, 'P63', 'PG6', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,1/2+z * y,y-x,1/2+z * x-y,x,1/2+z * end'),
174: ((6, 0, 'P-6', 'PG6BAR', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * x,y,-z * -y,x-y,-z * y-x,-x,-z * end'),
175: ((12, 0, 'P6/M', 'PG6/M', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,z * y,y-x,z * x-y,x,z * -x,-y,-z * y,y-x,-z * x-y,x,-z * x,y,-z * -y,x-y,-z * y-x,-x,-z * end'),
176: ((12, 0, 'P63/M', 'PG6/M', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,1/2+z * y,y-x,1/2+z * x-y,x,1/2+z * -x,-y,-z * y,y-x,-z * x-y,x,-z * x,y,1/2-z * -y,x-y,1/2-z * y-x,-x,1/2-z * end'),
177: ((12, 0, 'P622', 'PG622', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,z * y,y-x,z * x-y,x,z * y,x,-z * x-y,-y,-z * -x,y-x,-z * -y,-x,-z * y-x,y,-z * x,x-y,-z * end'),
178: ((12, 0, 'P6122', 'PG622', 'HEXAGONAL'),'x,y,z * -y,x-y,1/3+z * y-x,-x,2/3+z * -x,-y,1/2+z * y,y-x,5/6+z * x-y,x,1/6+z * y,x,1/3-z * x-y,-y,-z * -x,y-x,2/3-z * -y,-x,5/6-z * y-x,y,1/2-z * x,x-y,1/6-z * end'),
179: ((12, 0, 'P6522', 'PG622', 'HEXAGONAL'),'x,y,z * -y,x-y,2/3+z * y-x,-x,1/3+z * -x,-y,1/2+z * y,y-x,1/6+z * x-y,x,5/6+z * y,x,2/3-z * x-y,-y,-z * -x,y-x,1/3-z * -y,-x,1/6-z * y-x,y,1/2-z * x,x-y,5/6-z * end'),
180: ((12, 0, 'P6222', 'PG622', 'HEXAGONAL'),'x,y,z * -y,x-y,2/3+z * y-x,-x,1/3+z * -x,-y,z * y,y-x,2/3+z * x-y,x,1/3+z * y,x,2/3-z * x-y,-y,-z * -x,y-x,1/3-z * -y,-x,2/3-z * y-x,y,-z * x,x-y,1/3-z * end'),
181: ((12, 0, 'P6422', 'PG622', 'HEXAGONAL'),'x,y,z * -y,x-y,1/3+z * y-x,-x,2/3+z * -x,-y,z * y,y-x,1/3+z * x-y,x,2/3+z * y,x,1/3-z * x-y,-y,-z * -x,y-x,2/3-z * -y,-x,1/3-z * y-x,y,-z * x,x-y,2/3-z * end'),
182: ((12, 0, 'P6322', 'PG622', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,1/2+z * y,y-x,1/2+z * x-y,x,1/2+z * y,x,-z * x-y,-y,-z * -x,y-x,-z * -y,-x,1/2-z * y-x,y,1/2-z * x,x-y,1/2-z * end'),
183: ((12, 0, 'P6MM', 'PG6MM', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,z * y,y-x,z * x-y,x,z * -y,-x,z * y-x,y,z * x,x-y,z * y,x,z * x-y,-y,z * -x,y-x,z * end'),
184: ((12, 0, 'P6CC', 'PG6MM', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,z * y,y-x,z * x-y,x,z * -y,-x,1/2+z * y-x,y,1/2+z * x,x-y,1/2+z * y,x,1/2+z * x-y,-y,1/2+z * -x,y-x,1/2+z * end'),
185: ((12, 0, 'P63CM', 'PG6MM', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,1/2+z * y,y-x,1/2+z * x-y,x,1/2+z * -y,-x,1/2+z * y-x,y,1/2+z * x,x-y,1/2+z * y,x,z * x-y,-y,z * -x,y-x,z * end'),
186: ((12, 0, 'P63MC', 'PG6MM', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,1/2+z * y,y-x,1/2+z * x-y,x,1/2+z * -y,-x,z * y-x,y,z * x,x-y,z * y,x,1/2+z * x-y,-y,1/2+z * -x,y-x,1/2+z * end'),
187: ((12, 0, 'P-6M2', 'PG6BARM2', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * x,y,-z * -y,x-y,-z * y-x,-x,-z * -y,-x,z * y-x,y,z * x,x-y,z * -y,-x,-z * y-x,y,-z * x,x-y,-z * end'),
188: ((12, 0, 'P-6C2', 'PG6BARM2', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * x,y,1/2-z * -y,x-y,1/2-z * y-x,-x,1/2-z * -y,-x,1/2+z * y-x,y,1/2+z * x,x-y,1/2+z * -y,-x,-z * y-x,y,-z * x,x-y,-z * end'),
189: ((12, 0, 'P-62M', 'PG6BAR2M', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * x,y,-z * -y,x-y,-z * y-x,-x,-z * y,x,-z * x-y,-y,-z * -x,y-x,-z * y,x,z * x-y,-y,z * -x,y-x,z * end'),
190: ((12, 0, 'P-62C', 'PG6BAR2M', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * x,y,1/2-z * -y,x-y,1/2-z * y-x,-x,1/2-z * y,x,-z * x-y,-y,-z * -x,y-x,-z * y,x,1/2+z * x-y,-y,1/2+z * -x,y-x,1/2+z * end'),
191: ((24, 0, 'P6/MMM', 'PG6/MMM', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,z * y,y-x,z * x-y,x,z * y,x,-z * x-y,-y,-z * -x,y-x,-z * -y,-x,-z * y-x,y,-z * x,x-y,-z * -x,-y,-z * y,y-x,-z * x-y,x,-z * x,y,-z * y-x,-x,-z * -y,x-y,-z * -y,-x,z * y-x,y,z * x,x-y,z * y,x,z * x-y,-y,z * -x,y-x,z * end'),
192: ((24, 0, 'P6/MCC', 'PG6/MMM', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,z * y,y-x,z * x-y,x,z * y,x,1/2-z * x-y,-y,1/2-z * -x,y-x,1/2-z * -y,-x,1/2-z * y-x,y,1/2-z * x,x-y,1/2-z * -x,-y,-z * y,y-x,-z * x-y,x,-z * x,y,-z * y-x,-x,-z * -y,x-y,-z * -y,-x,1/2+z * y-x,y,1/2+z * x,x-y,1/2+z * y,x,1/2+z * x-y,-y,1/2+z * -x,y-x,1/2+z * end'),
193: ((24, 0, 'P63/MCM', 'PG6/MMM', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,1/2+z * y,y-x,1/2+z * x-y,x,1/2+z * y,x,1/2-z * x-y,-y,1/2-z * -x,y-x,1/2-z * -y,-x,-z * y-x,y,-z * x,x-y,-z * -x,-y,-z * y,y-x,-z * x-y,x,-z * x,y,1/2-z * y-x,-x,1/2-z * -y,x-y,1/2-z * -y,-x,1/2+z * y-x,y,1/2+z * x,x-y,1/2+z * y,x,z * x-y,-y,z * -x,y-x,z * end'),
194: ((24, 0, 'P63/MMC', 'PG6/MMM', 'HEXAGONAL'),'x,y,z * -y,x-y,z * y-x,-x,z * -x,-y,1/2+z * y,y-x,1/2+z * x-y,x,1/2+z * y,x,-z * x-y,-y,-z * -x,y-x,-z * -y,-x,1/2-z * y-x,y,1/2-z * x,x-y,1/2-z * -x,-y,-z * y,y-x,-z * x-y,x,-z * x,y,1/2-z * y-x,-x,1/2-z * -y,x-y,1/2-z * -y,-x,z * y-x,y,z * x,x-y,z * y,x,1/2+z * x-y,-y,1/2+z * -x,y-x,1/2+z * end'),
195: ((12, 0, 'P23', 'PG23', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * end'),
196: ((12, 3, 'F23', 'PG23', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
197: ((12, 1, 'I23', 'PG23', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * 1/2,1/2,1/2 * end'),
198: ((12, 0, 'P213', 'PG23', 'CUBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * z,x,y * 1/2+z,1/2-x,-y * 1/2-z,-x,1/2+y * -z,1/2+x,1/2-y * y,z,x * -y,1/2+z,1/2-x * 1/2+y,1/2-z,-x * 1/2-y,-z,1/2+x * end'),
199: ((12, 1, 'I213', 'PG23', 'CUBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * z,x,y * 1/2+z,1/2-x,-y * 1/2-z,-x,1/2+y * -z,1/2+x,1/2-y * y,z,x * -y,1/2+z,1/2-x * 1/2+y,1/2-z,-x * 1/2-y,-z,1/2+x * 1/2,1/2,1/2 * end'),
200: ((24, 0, 'PM-3', 'PGM3BAR', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * -z,-x,-y * -z,x,y * z,x,-y * z,-x,y * -y,-z,-x * y,-z,x * -y,z,x * y,z,-x * end'),
201: ((24, 0, 'PN-3', 'PGM3BAR', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * 1/2-x,1/2-y,1/2-z * 1/2+x,1/2+y,1/2-z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * 1/2-z,1/2-x,1/2-y * 1/2-z,1/2+x,1/2+y * 1/2+z,1/2+x,1/2-y * 1/2+z,1/2-x,1/2+y * 1/2-y,1/2-z,1/2-x * 1/2+y,1/2-z,1/2+x * 1/2-y,1/2+z,1/2+x * 1/2+y,1/2+z,1/2-x * end'),
202: ((24, 3, 'FM-3', 'PGM3BAR', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * -z,-x,-y * -z,x,y * z,x,-y * z,-x,y * -y,-z,-x * y,-z,x * -y,z,x * y,z,-x * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
203: ((24, 3, 'FD-3', 'PGM3BAR', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * 1/4-x,1/4-y,1/4-z * 1/4+x,1/4+y,1/4-z * 1/4+x,1/4-y,1/4+z * 1/4-x,1/4+y,1/4+z * 1/4-z,1/4-x,1/4-y * 1/4-z,1/4+x,1/4+y * 1/4+z,1/4+x,1/4-y * 1/4+z,1/4-x,1/4+y * 1/4-y,1/4-z,1/4-x * 1/4+y,1/4-z,1/4+x * 1/4-y,1/4+z,1/4+x * 1/4+y,1/4+z,1/4-x * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
204: ((24, 1, 'IM-3', 'PGM3BAR', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * -z,-x,-y * -z,x,y * z,x,-y * z,-x,y * -y,-z,-x * y,-z,x * -y,z,x * y,z,-x * 1/2,1/2,1/2 * end'),
205: ((24, 0, 'PA-3', 'PGM3BAR', 'CUBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * z,x,y * 1/2+z,1/2-x,-y * 1/2-z,-x,1/2+y * -z,1/2+x,1/2-y * y,z,x * -y,1/2+z,1/2-x * 1/2+y,1/2-z,-x * 1/2-y,-z,1/2+x * -x,-y,-z * 1/2+x,y,1/2-z * x,1/2-y,1/2+z * 1/2-x,1/2+y,z * -z,-x,-y * 1/2-z,1/2+x,y * 1/2+z,x,1/2-y * z,1/2-x,1/2+y * -y,-z,-x * y,1/2-z,1/2+x * 1/2-y,1/2+z,x * 1/2+y,z,1/2-x * end'),
206: ((24, 1, 'IA-3', 'PGM3BAR', 'CUBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * z,x,y * 1/2+z,1/2-x,-y * 1/2-z,-x,1/2+y * -z,1/2+x,1/2-y * y,z,x * -y,1/2+z,1/2-x * 1/2+y,1/2-z,-x * 1/2-y,-z,1/2+x * -x,-y,-z * 1/2+x,y,1/2-z * x,1/2-y,1/2+z * 1/2-x,1/2+y,z * -z,-x,-y * 1/2-z,1/2+x,y * 1/2+z,x,1/2-y * z,1/2-x,1/2+y * -y,-z,-x * y,1/2-z,1/2+x * 1/2-y,1/2+z,x * 1/2+y,z,1/2-x * 1/2,1/2,1/2 * end'),
207: ((24, 0, 'P432', 'PG432', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * y,x,-z * -y,-x,-z * y,-x,z * -y,x,z * x,z,-y * -x,z,y * -x,-z,-y * x,-z,y * z,y,-x * z,-y,x * -z,y,x * -z,-y,-x * end'),
208: ((24, 0, 'P4232', 'PG432', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * 1/2+y,1/2+x,1/2-z * 1/2-y,1/2-x,1/2-z * 1/2+y,1/2-x,1/2+z * 1/2-y,1/2+x,1/2+z * 1/2+x,1/2+z,1/2-y * 1/2-x,1/2+z,1/2+y * 1/2-x,1/2-z,1/2-y * 1/2+x,1/2-z,1/2+y * 1/2+z,1/2+y,1/2-x * 1/2+z,1/2-y,1/2+x * 1/2-z,1/2+y,1/2+x * 1/2-z,1/2-y,1/2-x * end'),
209: ((24, 3, 'F432', 'PG432', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * y,x,-z * -y,-x,-z * y,-x,z * -y,x,z * x,z,-y * -x,z,y * -x,-z,-y * x,-z,y * z,y,-x * z,-y,x * -z,y,x * -z,-y,-x * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
210: ((24, 3, 'F4132', 'PG432', 'CUBIC'),'x,y,z * -x,1/2-y,1/2+z * 1/2-x,1/2+y,-z * 1/2+x,-y,1/2-z * z,x,y * 1/2+z,-x,1/2-y * -z,1/2-x,1/2+y * 1/2-z,1/2+x,-y * y,z,x * 1/2-y,1/2+z,-x * 1/2+y,-z,1/2-x * -y,1/2-z,1/2+x * 3/4+y,1/4+x,3/4-z * 1/4-y,1/4-x,1/4-z * 1/4+y,3/4-x,3/4+z * 3/4-y,3/4+x,1/4+z * 3/4+x,1/4+z,3/4-y * 3/4-x,3/4+z,1/4+y * 1/4-x,1/4-z,1/4-y * 1/4+x,3/4-z,3/4+y * 3/4+z,1/4+y,3/4-x * 1/4+z,3/4-y,3/4+x * 3/4-z,3/4+y,1/4+x * 1/4-z,1/4-y,1/4-x * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
211: ((24, 1, 'I432', 'PG432', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * y,x,-z * -y,-x,-z * y,-x,z * -y,x,z * x,z,-y * -x,z,y * -x,-z,-y * x,-z,y * z,y,-x * z,-y,x * -z,y,x * -z,-y,-x * 1/2,1/2,1/2 * end'),
212: ((24, 0, 'P4332', 'PG432', 'CUBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * z,x,y * 1/2+z,1/2-x,-y * 1/2-z,-x,1/2+y * -z,1/2+x,1/2-y * y,z,x * -y,1/2+z,1/2-x * 1/2+y,1/2-z,-x * 1/2-y,-z,1/2+x * 1/4+y,3/4+x,3/4-z * 1/4-y,1/4-x,1/4-z * 3/4+y,3/4-x,1/4+z * 3/4-y,1/4+x,3/4+z * 1/4+x,3/4+z,3/4-y * 3/4-x,1/4+z,3/4+y * 1/4-x,1/4-z,1/4-y * 3/4+x,3/4-z,1/4+y * 1/4+z,3/4+y,3/4-x * 3/4+z,3/4-y,1/4+x * 3/4-z,1/4+y,3/4+x * 1/4-z,1/4-y,1/4-x * end'),
213: ((24, 0, 'P4132', 'PG432', 'CUBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * z,x,y * 1/2+z,1/2-x,-y * 1/2-z,-x,1/2+y * -z,1/2+x,1/2-y * y,z,x * -y,1/2+z,1/2-x * 1/2+y,1/2-z,-x * 1/2-y,-z,1/2+x * 3/4+y,1/4+x,1/4-z * 3/4-y,3/4-x,3/4-z * 1/4+y,1/4-x,3/4+z * 1/4-y,3/4+x,1/4+z * 3/4+x,1/4+z,1/4-y * 1/4-x,3/4+z,1/4+y * 3/4-x,3/4-z,3/4-y * 1/4+x,1/4-z,3/4+y * 3/4+z,1/4+y,1/4-x * 1/4+z,1/4-y,3/4+x * 1/4-z,3/4+y,1/4+x * 3/4-z,3/4-y,3/4-x * end'),
214: ((24, 1, 'I4132', 'PG432', 'CUBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * z,x,y * 1/2+z,1/2-x,-y * 1/2-z,-x,1/2+y * -z,1/2+x,1/2-y * y,z,x * -y,1/2+z,1/2-x * 1/2+y,1/2-z,-x * 1/2-y,-z,1/2+x * 3/4+y,1/4+x,1/4-z * 3/4-y,3/4-x,3/4-z * 1/4+y,1/4-x,3/4+z * 1/4-y,3/4+x,1/4+z * 3/4+x,1/4+z,1/4-y * 1/4-x,3/4+z,1/4+y * 3/4-x,3/4-z,3/4-y * 1/4+x,1/4-z,3/4+y * 3/4+z,1/4+y,1/4-x * 1/4+z,1/4-y,3/4+x * 1/4-z,3/4+y,1/4+x * 3/4-z,3/4-y,3/4-x * 1/2,1/2,1/2 * end'),
215: ((24, 0, 'P-43M', 'PG4BAR3M', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * y,x,z * -y,-x,z * y,-x,-z * -y,x,-z * x,z,y * -x,z,-y * -x,-z,y * x,-z,-y * z,y,x * z,-y,-x * -z,y,-x * -z,-y,x * end'),
216: ((24, 3, 'F-43M', 'PG4BAR3M', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * y,x,z * -y,-x,z * y,-x,-z * -y,x,-z * x,z,y * -x,z,-y * -x,-z,y * x,-z,-y * z,y,x * z,-y,-x * -z,y,-x * -z,-y,x * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
217: ((24, 1, 'I-43M', 'PG4BAR3M', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * y,x,z * -y,-x,z * y,-x,-z * -y,x,-z * x,z,y * -x,z,-y * -x,-z,y * x,-z,-y * z,y,x * z,-y,-x * -z,y,-x * -z,-y,x * 1/2,1/2,1/2 * end'),
218: ((24, 0, 'P-43N', 'PG4BAR3M', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * 1/2+y,1/2+x,1/2+z * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2-x,1/2-z * 1/2-y,1/2+x,1/2-z * 1/2+x,1/2+z,1/2+y * 1/2-x,1/2+z,1/2-y * 1/2-x,1/2-z,1/2+y * 1/2+x,1/2-z,1/2-y * 1/2+z,1/2+y,1/2+x * 1/2+z,1/2-y,1/2-x * 1/2-z,1/2+y,1/2-x * 1/2-z,1/2-y,1/2+x * end'),
219: ((24, 3, 'F-43C', 'PG4BAR3M', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * 1/2+y,1/2+x,1/2+z * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2-x,1/2-z * 1/2-y,1/2+x,1/2-z * 1/2+x,1/2+z,1/2+y * 1/2-x,1/2+z,1/2-y * 1/2-x,1/2-z,1/2+y * 1/2+x,1/2-z,1/2-y * 1/2+z,1/2+y,1/2+x * 1/2+z,1/2-y,1/2-x * 1/2-z,1/2+y,1/2-x * 1/2-z,1/2-y,1/2+x * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
220: ((24, 1, 'I-43D', 'PG4BAR3M', 'CUBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * z,x,y * 1/2+z,1/2-x,-y * 1/2-z,-x,1/2+y * -z,1/2+x,1/2-y * y,z,x * -y,1/2+z,1/2-x * 1/2+y,1/2-z,-x * 1/2-y,-z,1/2+x * 1/4+y,1/4+x,1/4+z * 1/4-y,3/4-x,3/4+z * 3/4+y,1/4-x,3/4-z * 3/4-y,3/4+x,1/4-z * 1/4+x,1/4+z,1/4+y * 3/4-x,3/4+z,1/4-y * 1/4-x,3/4-z,3/4+y * 3/4+x,1/4-z,3/4-y * 1/4+z,1/4+y,1/4+x * 3/4+z,1/4-y,3/4-x * 3/4-z,3/4+y,1/4-x * 1/4-z,3/4-y,3/4+x * 1/2,1/2,1/2 * end'),
221: ((48, 0, 'PM-3M', 'PGM3BARM', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * y,x,-z * -y,-x,-z * y,-x,z * -y,x,z * x,z,-y * -x,z,y * -x,-z,-y * x,-z,y * z,y,-x * z,-y,x * -z,y,x * -z,-y,-x * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * -z,-x,-y * -z,x,y * z,x,-y * z,-x,y * -y,-z,-x * y,-z,x * -y,z,x * y,z,-x * -y,-x,z * y,x,z * -y,x,-z * y,-x,-z * -x,-z,y * x,-z,-y * x,z,y * -x,z,-y * -z,-y,x * -z,y,-x * z,-y,-x * z,y,x * end'),
222: ((48, 0, 'PN-3N', 'PGM3BARM', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * y,x,-z * -y,-x,-z * y,-x,z * -y,x,z * x,z,-y * -x,z,y * -x,-z,-y * x,-z,y * z,y,-x * z,-y,x * -z,y,x * -z,-y,-x * 1/2-x,1/2-y,1/2-z * 1/2+x,1/2+y,1/2-z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * 1/2-z,1/2-x,1/2-y * 1/2-z,1/2+x,1/2+y * 1/2+z,1/2+x,1/2-y * 1/2+z,1/2-x,1/2+y * 1/2-y,1/2-z,1/2-x * 1/2+y,1/2-z,1/2+x * 1/2-y,1/2+z,1/2+x * 1/2+y,1/2+z,1/2-x * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2+x,1/2+z * 1/2-y,1/2+x,1/2-z * 1/2+y,1/2-x,1/2-z * 1/2-x,1/2-z,1/2+y * 1/2+x,1/2-z,1/2-y * 1/2+x,1/2+z,1/2+y * 1/2-x,1/2+z,1/2-y * 1/2-z,1/2-y,1/2+x * 1/2-z,1/2+y,1/2-x * 1/2+z,1/2-y,1/2-x * 1/2+z,1/2+y,1/2+x * end'),
223: ((48, 0, 'PM-3N', 'PGM3BARM', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * 1/2+y,1/2+x,1/2-z * 1/2-y,1/2-x,1/2-z * 1/2+y,1/2-x,1/2+z * 1/2-y,1/2+x,1/2+z * 1/2+x,1/2+z,1/2-y * 1/2-x,1/2+z,1/2+y * 1/2-x,1/2-z,1/2-y * 1/2+x,1/2-z,1/2+y * 1/2+z,1/2+y,1/2-x * 1/2+z,1/2-y,1/2+x * 1/2-z,1/2+y,1/2+x * 1/2-z,1/2-y,1/2-x * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * -z,-x,-y * -z,x,y * z,x,-y * z,-x,y * -y,-z,-x * y,-z,x * -y,z,x * y,z,-x * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2+x,1/2+z * 1/2-y,1/2+x,1/2-z * 1/2+y,1/2-x,1/2-z * 1/2-x,1/2-z,1/2+y * 1/2+x,1/2-z,1/2-y * 1/2+x,1/2+z,1/2+y * 1/2-x,1/2+z,1/2-y * 1/2-z,1/2-y,1/2+x * 1/2-z,1/2+y,1/2-x * 1/2+z,1/2-y,1/2-x * 1/2+z,1/2+y,1/2+x * end'),
224: ((48, 0, 'PN-3M', 'PGM3BARM', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * 1/2+y,1/2+x,1/2-z * 1/2-y,1/2-x,1/2-z * 1/2+y,1/2-x,1/2+z * 1/2-y,1/2+x,1/2+z * 1/2+x,1/2+z,1/2-y * 1/2-x,1/2+z,1/2+y * 1/2-x,1/2-z,1/2-y * 1/2+x,1/2-z,1/2+y * 1/2+z,1/2+y,1/2-x * 1/2+z,1/2-y,1/2+x * 1/2-z,1/2+y,1/2+x * 1/2-z,1/2-y,1/2-x * 1/2-x,1/2-y,1/2-z * 1/2+x,1/2+y,1/2-z * 1/2+x,1/2-y,1/2+z * 1/2-x,1/2+y,1/2+z * 1/2-z,1/2-x,1/2-y * 1/2-z,1/2+x,1/2+y * 1/2+z,1/2+x,1/2-y * 1/2+z,1/2-x,1/2+y * 1/2-y,1/2-z,1/2-x * 1/2+y,1/2-z,1/2+x * 1/2-y,1/2+z,1/2+x * 1/2+y,1/2+z,1/2-x * -y,-x,z * y,x,z * -y,x,-z * y,-x,-z * -x,-z,y * x,-z,-y * x,z,y * -x,z,-y * -z,-y,x * -z,y,-x * z,-y,-x * z,y,x * end'),
225: ((48, 3, 'FM-3M', 'PGM3BARM', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * y,x,-z * -y,-x,-z * y,-x,z * -y,x,z * x,z,-y * -x,z,y * -x,-z,-y * x,-z,y * z,y,-x * z,-y,x * -z,y,x * -z,-y,-x * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * -z,-x,-y * -z,x,y * z,x,-y * z,-x,y * -y,-z,-x * y,-z,x * -y,z,x * y,z,-x * -y,-x,z * y,x,z * -y,x,-z * y,-x,-z * -x,-z,y * x,-z,-y * x,z,y * -x,z,-y * -z,-y,x * -z,y,-x * z,-y,-x * z,y,x * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
226: ((48, 3, 'FM-3C', 'PGM3BARM', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * 1/2+y,1/2+x,1/2-z * 1/2-y,1/2-x,1/2-z * 1/2+y,1/2-x,1/2+z * 1/2-y,1/2+x,1/2+z * 1/2+x,1/2+z,1/2-y * 1/2-x,1/2+z,1/2+y * 1/2-x,1/2-z,1/2-y * 1/2+x,1/2-z,1/2+y * 1/2+z,1/2+y,1/2-x * 1/2+z,1/2-y,1/2+x * 1/2-z,1/2+y,1/2+x * 1/2-z,1/2-y,1/2-x * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * -z,-x,-y * -z,x,y * z,x,-y * z,-x,y * -y,-z,-x * y,-z,x * -y,z,x * y,z,-x * 1/2-y,1/2-x,1/2+z * 1/2+y,1/2+x,1/2+z * 1/2-y,1/2+x,1/2-z * 1/2+y,1/2-x,1/2-z * 1/2-x,1/2-z,1/2+y * 1/2+x,1/2-z,1/2-y * 1/2+x,1/2+z,1/2+y * 1/2-x,1/2+z,1/2-y * 1/2-z,1/2-y,1/2+x * 1/2-z,1/2+y,1/2-x * 1/2+z,1/2-y,1/2-x * 1/2+z,1/2+y,1/2+x * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
227: ((48, 3, 'FD-3M', 'PGM3BARM', 'CUBIC'),'x,y,z * -x,1/2-y,1/2+z * 1/2-x,1/2+y,-z * 1/2+x,-y,1/2-z * z,x,y * 1/2+z,-x,1/2-y * -z,1/2-x,1/2+y * 1/2-z,1/2+x,-y * y,z,x * 1/2-y,1/2+z,-x * 1/2+y,-z,1/2-x * -y,1/2-z,1/2+x * 3/4+y,1/4+x,3/4-z * 1/4-y,1/4-x,1/4-z * 1/4+y,3/4-x,3/4+z * 3/4-y,3/4+x,1/4+z * 3/4+x,1/4+z,3/4-y * 3/4-x,3/4+z,1/4+y * 1/4-x,1/4-z,1/4-y * 1/4+x,3/4-z,3/4+y * 3/4+z,1/4+y,3/4-x * 1/4+z,3/4-y,3/4+x * 3/4-z,3/4+y,1/4+x * 1/4-z,1/4-y,1/4-x * 1/4-x,1/4-y,1/4-z * 1/4+x,3/4+y,3/4-z * 3/4+x,3/4-y,1/4+z * 3/4-x,1/4+y,3/4+z * 1/4-z,1/4-x,1/4-y * 3/4-z,1/4+x,3/4+y * 1/4+z,3/4+x,3/4-y * 3/4+z,3/4-x,1/4+y * 1/4-y,1/4-z,1/4-x * 3/4+y,3/4-z,1/4+x * 3/4-y,1/4+z,3/4+x * 1/4+y,3/4+z,3/4-x * 1/2-y,-x,1/2+z * y,x,z * -y,1/2+x,1/2-z * 1/2+y,1/2-x,-z * 1/2-x,-z,1/2+y * 1/2+x,1/2-z,-y * x,z,y * -x,1/2+z,1/2-y * 1/2-z,-y,1/2+x * -z,1/2+y,1/2-x * 1/2+z,1/2-y,-x * z,y,x * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
228: ((48, 3, 'FD-3C', 'PGM3BARM', 'CUBIC'),'x,y,z * -x,1/2-y,1/2+z * 1/2-x,1/2+y,-z * 1/2+x,-y,1/2-z * z,x,y * 1/2+z,-x,1/2-y * -z,1/2-x,1/2+y * 1/2-z,1/2+x,-y * y,z,x * 1/2-y,1/2+z,-x * 1/2+y,-z,1/2-x * -y,1/2-z,1/2+x * 3/4+y,1/4+x,3/4-z * 1/4-y,1/4-x,1/4-z * 1/4+y,3/4-x,3/4+z * 3/4-y,3/4+x,1/4+z * 3/4+x,1/4+z,3/4-y * 3/4-x,3/4+z,1/4+y * 1/4-x,1/4-z,1/4-y * 1/4+x,3/4-z,3/4+y * 3/4+z,1/4+y,3/4-x * 1/4+z,3/4-y,3/4+x * 3/4-z,3/4+y,1/4+x * 1/4-z,1/4-y,1/4-x * 3/4-x,3/4-y,3/4-z * 3/4+x,1/4+y,1/4-z * 1/4+x,1/4-y,3/4+z * 1/4-x,3/4+y,1/4+z * 3/4-z,3/4-x,3/4-y * 1/4-z,3/4+x,1/4+y * 3/4+z,1/4+x,1/4-y * 1/4+z,1/4-x,3/4+y * 3/4-y,3/4-z,3/4-x * 1/4+y,1/4-z,3/4+x * 1/4-y,3/4+z,1/4+x * 3/4+y,1/4+z,1/4-x * -y,1/2-x,z * 1/2+y,1/2+x,1/2+z * 1/2-y,x,-z * y,-x,1/2-z * -x,1/2-z,y * x,-z,1/2-y * 1/2+x,1/2+z,1/2+y * 1/2-x,z,-y * -z,1/2-y,x * 1/2-z,y,-x * z,-y,1/2-x * 1/2+z,1/2+y,1/2+x * 0,1/2,1/2 * 1/2,0,1/2 * 1/2,1/2,0 * end'),
229: ((48, 1, 'IM-3M', 'PGM3BARM', 'CUBIC'),'x,y,z * -x,-y,z * -x,y,-z * x,-y,-z * z,x,y * z,-x,-y * -z,-x,y * -z,x,-y * y,z,x * -y,z,-x * y,-z,-x * -y,-z,x * y,x,-z * -y,-x,-z * y,-x,z * -y,x,z * x,z,-y * -x,z,y * -x,-z,-y * x,-z,y * z,y,-x * z,-y,x * -z,y,x * -z,-y,-x * -x,-y,-z * x,y,-z * x,-y,z * -x,y,z * -z,-x,-y * -z,x,y * z,x,-y * z,-x,y * -y,-z,-x * y,-z,x * -y,z,x * y,z,-x * -y,-x,z * y,x,z * -y,x,-z * y,-x,-z * -x,-z,y * x,-z,-y * x,z,y * -x,z,-y * -z,-y,x * -z,y,-x * z,-y,-x * z,y,x * 1/2,1/2,1/2 * end'),
230: ((48, 1, 'IA-3D', 'PGM3BARM', 'CUBIC'),'x,y,z * 1/2-x,-y,1/2+z * -x,1/2+y,1/2-z * 1/2+x,1/2-y,-z * z,x,y * 1/2+z,1/2-x,-y * 1/2-z,-x,1/2+y * -z,1/2+x,1/2-y * y,z,x * -y,1/2+z,1/2-x * 1/2+y,1/2-z,-x * 1/2-y,-z,1/2+x * 3/4+y,1/4+x,1/4-z * 3/4-y,3/4-x,3/4-z * 1/4+y,1/4-x,3/4+z * 1/4-y,3/4+x,1/4+z * 3/4+x,1/4+z,1/4-y * 1/4-x,3/4+z,1/4+y * 3/4-x,3/4-z,3/4-y * 1/4+x,1/4-z,3/4+y * 3/4+z,1/4+y,1/4-x * 1/4+z,1/4-y,3/4+x * 1/4-z,3/4+y,1/4+x * 3/4-z,3/4-y,3/4-x * -x,-y,-z * 1/2+x,y,1/2-z * x,1/2-y,1/2+z * 1/2-x,1/2+y,z * -z,-x,-y * 1/2-z,1/2+x,y * 1/2+z,x,1/2-y * z,1/2-x,1/2+y * -y,-z,-x * y,1/2-z,1/2+x * 1/2-y,1/2+z,x * 1/2+y,z,1/2-x * 1/4-y,3/4-x,3/4+z * 1/4+y,1/4+x,1/4+z * 3/4-y,3/4+x,1/4-z * 3/4+y,1/4-x,3/4-z * 1/4-x,3/4-z,3/4+y * 3/4+x,1/4-z,3/4-y * 1/4+x,1/4+z,1/4+y * 3/4-x,3/4+z,1/4-y * 1/4-z,3/4-y,3/4+x * 3/4-z,3/4+y,1/4-x * 3/4+z,1/4-y,3/4-x * 1/4+z,1/4+y,1/4+x * 1/2,1/2,1/2 * end')}

def opReadCl(filename):
    f = open(filename)
    r = f.read()
    f.close()
    return r

def opWriteCl(filename, _str):
    f = open(filename,"w")
    r = f.write(_str)
    f.close()

def guessHA(wavelength):
    "Trying to find the closest HA edge in HAd from wavelength"
    diffs = map(lambda x: abs(x - wavelength),HAd.keys())
    return HAd[HAd.keys()[diffs.index(min(diffs))]]

def get_crossec(P):
    idx = -1
    opWriteCl("crossec.inp", crossec_script % vars(P))
    try:
        os.system("crossec<crossec.inp>crossec.out")
        txt = opReadCl("crossec.out")
        idx = txt.index(" $$\n<B>")
    except: pass
    if idx != -1:
        try:
            fp, fpp = map(float, txt[idx-23:idx-2].split())
            return fp, fpp
        except: pass
    else: return 0., 0.

class Dumy:
    pass

class DoMode:
    def __init__(self, P):
        self.mode = P.mode
        self.dir_mode = self.mode.lower()+"/"
        P.dir_mode = self.dir_mode
        P.spg_name = (xupy.SPGlib[int(P.spgn_in)][0]).lower()
        try:
            # It will fail herre in case Numeric is not installed
            # Or for old version of XDS
            from xdsHKLinfos import get_info
            infos = get_info(P.file_name_in)
            P.res_high =  infos["res_high"]
            P.res_low = infos["res_low"]
        except: 
            P.res_high = 1.1
            P.res_low = 30

        if P.wavelength:
            self.wavelength = P.wavelength
        else:
            self.wavelength = 1.5418

        if P.ha_name:
            self.HA = P.ha_name
            P.ha_name = self.HA
            P.fp, P.fpp = get_crossec(P)
        elif P.wavelength:
            self.HA = guessHA(self.wavelength)
            P.ha_name = self.HA[0]
            P.fp, P.fpp = get_crossec(P)
        else:
            P.ha_name = "Unknown"
            P.fp, P.fpp = 0., 0.
            self.HA = "Unknown"

        if self.mode == "CNS":
            self.name_ext = ".cv"
            P.mode_out = self.mode

        elif self.mode == "SHELX":
            self.name_ext = "_F4.hkl"
            #P.free_out = "GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.05\n"
            P.mode_out = self.mode

        elif self.mode == "SOLVE":
            self.name_ext = ".hkli"
            P.mode_out = "NONE"
            P.merge_out == "FALSE"
            P.friedel_out = "TRUE"

        elif self.mode == "EPMR":
            self.name_ext = ".hkl"
            P.mode_out = "FALL"
            P.merge_out == "TRUE"
            P.friedel_out = "TRUE"

        elif self.mode == "AMORE":
            P.free_out = ""
            P.free_lbl = ""
            P.free_code = ""
            self.name_ext = ".HKL"
            P.mode_out = "CCP4"
            P.merge_out == "TRUE"
            P.friedel_out = "TRUE"

        elif self.mode == "REPLACE":
            P.free_out = ""
            P.free_lbl = ""
            P.free_code = ""
            self.name_ext = ".hkl"
            P.mode_out = "CCP4"
            P.merge_out == "TRUE"
            P.friedel_out = "TRUE"

        elif self.mode == "CRANK":
            self.name_ext = ".mtz"
            P.mode_out = "CCP4_I"
            P.merge_out == "TRUE"
            P.friedel_out = "FALSE"

        elif self.mode == "CCP4":
            self.name_ext = ".mtz"
            P.mode_out = self.mode
            if P.friedel_out == "FALSE":
                P.cinp_ano = "DANO%(lbl)s SIGDANO%(lbl)s ISYM%(lbl)s" \
                              % vars(P)
                P.cinp_ano2 = "D   Q   Y"
            else: P.cinp_ano = P.cinp_ano2 = ""
            P.merge_out = "TRUE"

        elif self.mode == "CCP4F":
            self.name_ext = ".mtz"
            P.mode_out = "CCP4"
            if P.friedel_out == "FALSE":
                P.cinp_ano = "DANO%(lbl)s SIGDANO%(lbl)s ISYM%(lbl)s" \
                    % vars(P)
                P.cad_ano = "E1=DANO%(lbl)s E2=SIGDANO%(lbl)s E3=ISYM%(lbl)s" \
                    % vars(P)
                P.cinp_ano2 = "D   Q   Y"
            else: P.cinp_ano = P.cad_ano = P.cinp_ano2 = ""
            P.merge_out = "TRUE"
            P.file_name_out = "F2MTZ.HKL.TMP2"
            opWriteCl("XDSCONV2.INP", xdsconv_script % vars(P) + P.free_out)
            opWriteCl("f2mtz.inp2", f2mtz_script % vars(P))
            opWriteCl("run_simple_scale.sh", scala_script)
            P.mode_out = "CCP4_F"
            if P.friedel_out == "FALSE":
                P.cinp_ano = "F(+)%(lbl)s SIGF(+)%(lbl)s F(-)%(lbl)s SIGF(-)%(lbl)s" \
                              % vars(P)
                P.cinp_ano2 = "G  L  G  L"
            else: P.cinp_ano = P.cinp_ano2 = ""
            P.merge_out = "TRUE"

# LABOUT  H K L FP SIGFP F(+) SIGF(+) F(-) SIGF(-) FreeRflag
# CTYPOUT H H H  F   Q    G     L      G     L         X
        elif self.mode == "PHASER":
            self.name_ext = ".mtz"
            P.mode_out = "CCP4_F"
            if P.friedel_out == "FALSE":
                P.cinp_ano = "F(+)%(lbl)s SIGF(+)%(lbl)s F(-)%(lbl)s SIGF(-)%(lbl)s" \
                              % vars(P)
                P.cinp_ano2 = "G  L  G  L"
            else: P.cinp_ano = P.cinp_ano2 = ""
            P.merge_out = "TRUE"

        elif self.mode == "SHARP":
            self.name_ext = ".mtz"
            P.mode_out = "CCP4"

            if P.friedel_out != "FALSE":
                print "Error, no Anomalous information found from XDS."
                sys.exit()
            P.merge_out = "TRUE"

        P.ident = ".".join(P.file_name_in.split(".")[:-1])
        #P.file_name_out = P.ident+"_"+P.ha_name+self.name_ext
        P.file_name_out = P.ident + self.name_ext
        if self.mode == "CCP4" or self.mode == "CCP4F" or \
           self.mode == "CRANK" or self.mode == "SHARP" or self.mode == "PHASER" :
            self.last_name = P.file_name_out
            P.file_name_out = "F2MTZ.HKL.TMP"
            P.last_name = self.last_name

        if self.mode != "SOLVE":
            opWriteCl("XDSCONV.INP", xdsconv_script % vars(P) + P.free_out)


    def run(self):
        ls = os.listdir(".")
        if not self.dir_mode[:-1] in ls:
            os.makedirs(self.dir_mode)

        if self.mode not in ("SOLVE",):
            toexec = os.path.join(xupy.XDSHOME,"xdsconv")
            xupy.exec_prog(toexec, stdout="xdsconv.log")
            if self.mode == "CCP4F":
                os.system("mv XDSCONV2.INP XDSCONV.INP")
                os.system("mv f2mtz.inp2 %s" % self.dir_mode)
                toexec = os.path.join(xupy.XDSHOME,"xdsconv")
                xupy.exec_prog(toexec, stdout="xdsconv_dano.log")


    def post_run(self, P):
        if self.mode == "SHELX":
            # run xprepx
            P.cns_sg = cns_sg_lib[int(XC.spgn_in)]
            opWriteCl("%s/xprep.inp" % P.dir_mode, xprep_script % vars(P))
            opWriteCl("%s/shelxc.inp" % P.dir_mode, shelxc_script % vars(P))
            opWriteCl("%s/run_shelx.sh" % P.dir_mode, shelxall_script % vars(P))
            os.system("chmod a+x %s/run_shelx.sh" % P.dir_mode)
            opWriteCl("%s/run_sitcom.sh" % P.dir_mode, sitcom_script % vars(P))
            os.system("chmod a+x %s/run_sitcom.sh" % P.dir_mode)
            #os.system("cd %s;xprep<xprep.inp>xprep.log" % P.dir_mode)
            os.system("cd %s;shelxc XX1 <shelxc.inp>shelxc.log" % P.dir_mode)

        elif self.mode == "CNS":
            #file_name_out = P.dirmode+dirname+".cv"
            P.cns_sg = cns_sg_lib[int(XC.spgn_in)]
            P.cns_cell = 6*"%.2f, " % tuple(P.cell)
            opWriteCl("%s/cns_xtal.par" % P.dir_mode, cns_par % vars(P))

        elif self.mode == "CCP4":
            opWriteCl("%s/f2mtz.inp" % P.dir_mode, f2mtz_script % vars(P))
            opWriteCl("%s/cad.inp" % P.dir_mode, cad_script % vars(P))
            os.system("cd %s;f2mtz hklout TMP.MTZ<f2mtz.inp >f2mtz.log" % P.dir_mode)
            os.system("cd %s;cad hklin1 TMP.MTZ hklout %s <cad.inp >cad.log"
                          % (P.dir_mode, self.last_name))
            os.system("rm -f F2MTZ.INP ccp4/TMP.MTZ ccp4/F2MTZ.HKL.TMP")

        elif self.mode == "CCP4F":
           opWriteCl("%s/f2mtz.inp" % P.dir_mode, f2mtz_script % vars(P))
           opWriteCl("%s/cad.inp" % P.dir_mode, cad2_script % vars(P))
           os.system("cd %s;f2mtz hklout TMP.MTZ<f2mtz.inp >f2mtz.log" % P.dir_mode)
           os.system("cd %s;f2mtz hklout TMP2.MTZ<f2mtz.inp2 >f2mtz2.log" % P.dir_mode)
           os.system("cd %s;cad hklin1 TMP.MTZ hklin2 TMP2.MTZ hklout %s <cad.inp >cad.log"
                         % (P.dir_mode, self.last_name))
           os.system("rm -f F2MTZ.INP ccp4f/TMP*.MTZ ccp4f/F2MTZ.HKL.TMP*")

        elif self.mode == "CRANK":
            #opWriteCl("%s/f2mtz.inp" % P.dir_mode, f2mtz_script % vars(P))
            opWriteCl("cad.inp", cad_crank_script % vars(P))
            os.system("f2mtz HKLOUT temp.mtz< F2MTZ.INP>f2mtz.log")
            os.system("bash cad.inp >cad.log")
            os.system("cd %s; mv ../output_file_name.mtz %s" % \
                                           (P.dir_mode, self.last_name))
            os.system("rm -f F2MTZ.INP temp.mtz crank/F2MTZ.HKL.TMP")

        elif self.mode == "PHASER":
            opWriteCl("%s/f2mtz.inp" % P.dir_mode, f2mtz_phaser_script % vars(P))
            opWriteCl("%s/cad.inp" % P.dir_mode, cad_script % vars(P))
            os.system("cd %s;f2mtz hklout TMP.MTZ<f2mtz.inp >f2mtz.log" % P.dir_mode)
            os.system("cd %s;cad hklin1 TMP.MTZ hklout %s <cad.inp >cad.log"
                          % (P.dir_mode, self.last_name))
            os.system("rm -f F2MTZ.INP phaser/TMP.MTZ phaser/F2MTZ.HKL.TMP")
            opWriteCl("%s/run_phaser.sh" % P.dir_mode, phaser_script % vars(P))
            os.chmod("%s/run_phaser.sh" % P.dir_mode, 0755)

        elif self.mode == "SHARP":
            opWriteCl("%s/f2mtz.inp" % P.dir_mode, f2mtz_sharp_script % vars(P))
            opWriteCl("%s/cad.inp" % P.dir_mode, cad_script)
            os.system("cd %s;f2mtz hklout TMP.MTZ<f2mtz.inp >f2mtz.log" % P.dir_mode)
            os.system("cd %s;cad hklin1 TMP.MTZ hklout %s <cad.inp >cad.log"
                          % (P.dir_mode, self.last_name))
            os.system("rm -f F2MTZ.INP ccp4/TMP.MTZ ccp4/F2MTZ.HKL.TMP")

        #elif self.mode == "CRANK":
            #opWriteCl("%s/mtzutils.inp" % P.dir_mode, mtzutils_script)
            #os.system("sed -e  's/J *Q *J *Q$/K M K M/' F2MTZ.INP> F2MTZS.INP")
            #os.system("f2mtz hklout %s/TMP.MTZ < F2MTZS.INP > %s/f2mtz.log" % \
                                                        #(P.dir_mode, P.dir_mode))
            #os.system("cd %s;mtzutils hklin TMP.MTZ hklout %s <mtzutils.inp >mtzutils.log"
                          #% (P.dir_mode, "crank.mtz"))
            #os.system("rm -f F2MTZ.INP crank/TMP.MTZ ccp4/F2MTZ.HKL.TMP")

        elif self.mode == "SOLVE":
            os.system("grep -v \! %s > solve/%s" % \
                                      (P.file_name_in, P.file_name_out))
            P.spg_name = (xupy.SPGlib[int(P.spgn_in)][0]).lower()
            opWriteCl("%s/run_solve.sh" % P.dir_mode, solve_script % vars(P))
            os.chmod("%s/run_solve.sh" % P.dir_mode, 0755)

        elif self.mode == "EPMR":
            opWriteCl("%s/cell" % P.dir_mode, "%(cell_str)s %(spgn_in)s\n" % vars(P))

        elif self.mode == "REPLACE":
            t = """awk '{gsub(",","");print}' """
            t += """%(dir_mode)s/%(file_name_out)s >> %(dir_mode)s/data.hkl"""
            os.system(t % vars(P))
            opWriteCl("%s/run_glrf_self.sh" % P.dir_mode, replace_script % vars(P))
            os.system("chmod a+x %(dir_mode)s/run_glrf_self.sh" % vars(P))
            os.system("rm -f %(dir_mode)s/%(file_name_out)s" % vars(P))

        elif self.mode == "AMORE":
            os.system('echo "CELL %(cell_str)s" > format.dat'%vars(P))
            os.system('echo "FORMAT (3I6,2E10.3)" >> format.dat')
            os.system("cat format.dat > %(dir_mode)s/hkl.d" % vars(P))
            t = """awk '{gsub(",","");print}' """
            t += """%(dir_mode)s/%(file_name_out)s >> %(dir_mode)s/hkl.d"""
            os.system(t % vars(P))
            os.system("rm -f format.dat %(dir_mode)s/%(file_name_out)s" % vars(P))
            afmt = " * xds *\n%(cell_str)s\n%(symop)s 0\n95. 0.\n15 3.5\n1 1\n"
            P.symop = amore_symops[int(P.spgn_in)][1]
            opWriteCl("%s/data.d" % P.dir_mode, afmt % vars(P))

if __name__=='__main__':

    import xupy

    # Default options
    __format_out = [] #"CCP4"
    __atom_name = ""
    __num_sites = 0
    __xds_input = ""
    __free_refl_input = ""
    __free_refl_type = ""
    __force_anom = False
    __force_norm = False
    __force_merge = False
    __force_unmerge = False
    __force_free = False
    __force_no_free = False
    __label = ""

    argp_fmt = "<<< %-24s %s"
    progname = sys.argv[0].split("/")[-1]
    if (len(sys.argv) == 1) or ("-h" in sys.argv) \
        or  ("--help" in sys.argv) or ("-help" in sys.argv):
        print usage % (progname, "|".join(options))
        sys.exit(1)

    args = sys.argv[1:]
    print "\n<== OPTIONS:"
    for arg in args:
        nonnum = [i for i in arg if i not in "1234567890"]
        if nonnum == []:
            try:
                n = int(arg)
                __num_sites = n
                print argp_fmt %("nSites:", n)
            except:
                pass
        elif arg.count("-l="):
            __label = "_"+arg[3:]
        elif arg.count("-l"):
            __label = "_"+str(args[args.index("-l") + 1])
        elif arg == "-f":
            __force_free = True
        elif arg == "-nf":
            __force_no_free = True
        elif arg == "-a":
            __force_anom = True
        elif arg == "-n":
            __force_norm = True
        elif arg == "-u":
            __force_unmerge = True
        elif arg == "-m":
            __force_merge = True
        # Geting output format
        elif arg.upper() in options:
            __format_out.append(arg.upper())
            print argp_fmt %("Export mode:", arg.upper())
        # Geting Atom type
        elif arg.title() in atom_names:
            print argp_fmt % ("atomType:", arg)
            __atom_name = arg.title()
        # Identifying file type
        elif os.path.isfile(arg):
            try:
                f = open(arg)
                s = f.read(65)
                f.close()
                ss = []
                try:
                    ss = s.split()
                    ss = map(float,ss)
                except:
                    pass
                if "!FORMAT=XDS_ASCII" in s:
                    __xds_input = arg
                elif (("NREFlection=" in s) and ("ANOMalous=" in s)):
                    __free_refl_input = arg
                    __free_refl_type = "CNS"
                elif (len(ss)>= 12) and (abs(ss[5]) == 1):
                    __free_refl_input = arg
                    __free_refl_type = "SHELX"
                elif (s[:3] == "MTZ") and (arg[-4:].lower() == ".mtz"):
                    __free_refl_input = arg
                    __free_refl_type = "CCP4"
                else:
                    print "\nWarning: Can't define the file type",
                    print "for argument '%s'\n" % arg
            except:
                pass

    # If input file for free reflection set is CCP4 it needs to be converted
    # to shelx format for input in xdsconv.
    if __free_refl_input:
        print argp_fmt % ("free_hkl_to_inherit:", __free_refl_input),
        print "[%s format]." % ( __free_refl_type)
        if __free_refl_type == "CCP4":
            print "   --> Converting CCP4 free reflections to SHELX format."
            script = open("mtz2shelx_free.sh", "w")
            script.write(mtz2various_script)
            script.close()
            os.system("sh mtz2shelx_free.sh %s" % __free_refl_input)
            __free_refl_input = "free_refl_shelx_F3.hkl"
            __free_refl_type = "SHELX"
    #
    if __force_anom and __force_norm:
        print "Warning: Umbiguous options specification (-a and -n), keeping -a."
        __force_norm == False

    # test input file type of inherited reflections
    if __free_refl_input:
        #__force_free = False
        xdsconv_script += "INHERIT_TEST_REFLECTIONS_FROM_FILE=%s %s\n" % \
                          (__free_refl_input, __free_refl_type)

    #### Default values
    XC = Dumy()
    XC.file_type = "XDS_ASCII"
    XC.friedel_out = ""
    XC.free_out = "GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.05\n"
    XC.free_lbl = "FreeR_flag"
    XC.free_code = "X"
    XC.merge_out = ""
    XC.cell = ""
    XC.lbl = ""
    XC.ha_name = __atom_name
    XC.num_sites = __num_sites
    if not __format_out: __format_out = ["CCP4"]
    XC.modes = __format_out
    #XC.wavelength = 1.5418111

    if __force_anom: XC.friedel_out = "FALSE"
    if __force_norm: XC.friedel_out = "TRUE"
    if __force_merge: XC.merge_out = "TRUE"
    if __force_unmerge: XC.merge_out = "FALSE"
    if __force_free:
        XC.free_out = "GENERATE_FRACTION_OF_TEST_REFLECTIONS=0.05\n"
        XC.free_lbl = "FreeR_flag"
        XC.free_code = "X"
    if __force_no_free:
        XC.free_out = ""
        XC.free_lbl = ""
        XC.free_code = ""
    if __label: XC.lbl = __label

    XC.dirname = os.path.split(os.getcwd())[1]
    XC.file_name_in = sys.argv[1]
    sys.argv.remove(XC.file_name_in)

    H = xupy.read_xdsascii_head(XC.file_name_in)
    XC.spgn_in = H["sym"]
    XC.cell = map(float,H["cell"].split())
    XC.cell_str = 6*"%.2f  " % tuple(XC.cell)
    XC.friedel_in = H["friedels_law"]
    XC.merge_in = H["merge"]
    if H["template_name"]:
        XC.ID = H["template_name"]
    else:
        XC.ID = "XDSa"

    if not H["wavelength"]:
    # Try to catch the wavelength from 
        if H["inputfile_name"] \
           and os.path.exists(H["inputfile_name"]) \
           and H["inputfile_type"] == "XDS_ASCII":
               _h = xupy.read_xdsascii_head(H["inputfile_name"])
               if _h["wavelength"] : XC.wavelength = _h["wavelength"]
    else: XC.wavelength = H["wavelength"]
    if H["friedels_law"] == "TRUE" and XC.friedel_out == "FALSE":
        print "\n>>> WARNING!  The input file does not contain Friedel's mate."
    if H["merge"] == "TRUE" and XC.merge_out == "FALSE":
        print "\n>>> WARNING!  The input file does not unmerged reflections."

    if not XC.friedel_out: XC.friedel_out = H["friedels_law"]
    if not XC.merge_out: XC.merge_out = H["merge"]

    # Depending on the chosen mode, predefine suffix_name, merge_out, friedel_out...
    print fmt_inp % vars(XC)

    #-----------------------------
    modes = XC.modes[:]
    for mode in modes:
        XC.mode = mode
        E = DoMode(XC)
        print fmt_outp % vars(XC)
        if XC.friedel_in == "FALSE":
            print fmt_outp_ha % vars(XC)
        E.run()
        E.post_run(XC)

