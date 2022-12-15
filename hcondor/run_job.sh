#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt/setup.sh

if [[ $4 == SV ]]
then
  python /nfs/cms/vazqueze/analisisWW/channelSV/selection_both_SV.py --process="$1" --year="$2" --type="$3" --ssos
else
  python /nfs/cms/vazqueze/analisisWW/selection_both.py --process="$1" --year="$2" --type="$3" --ssos 
fi
