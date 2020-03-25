#!/bin/bash 

if [ $# -eq 0 ] ; then
  HKLIN=XDS_ASCII.HKL
else
  HKLIN=$1  
fi


function get_ascii_type() {
  head -5 $1 | grep Generated | awk '{print $3}'
}

function get_prefix {
   b=$(head  $1 | grep NAME_TEMPLATE_OF_DATA_FRAMES  | awk '{print $1}')
   test -z $b && b="$1"
   c=$(basename $b)
   echo ${c%_*.*}
}

function run_staraniso() {
PREFIX=$(get_prefix $1)
ANOMALOUS=OFF
prefix=$(basename $1 .HKL)
pointless -copy XDSIN $1 \
          HKLOUT ${PREFIX}_pointless.mtz > ${PREFIX}_pointless.log
aimless hklin ${PREFIX}_pointless.mtz hklout ${PREFIX}_aimless.mtz \
      scales     AIMLESS.scales \
      rogueplot  AIMLESS.rogueplot \
      ROGUES AIMLESS.rogue \
      correlplot AIMLESS.correlplot \
      normplot   AIMLESS.norm \
      anomplot   AIMLESS.anom \
      xmlout   ${PREFIX}_aimless.xml \
      > ${PREFIX}_aimless.log << eof
cycles 0
scales constant    # batch scaling is generally poorer than smoothed
bins 20
ANOM $ANOMALOUS
end
eof

awk '/SUMMARY_BEGIN/{while(getline && !/SUMMARY_END/){print}}' ${PREFIX}_aimless.log > ${PREFIX}_aimless.SUMM.log

echo | staraniso HKLIN  ${PREFIX}_aimless.mtz HKLOUT  ${PREFIX}_staraniso.mtz >  ${PREFIX}_staraniso.log 
grep -A5 "Diffraction limits & principal axes of ellipsoid" ${PREFIX}_staraniso.log
}

#get_ascii_type $1

run_staraniso $HKLIN

