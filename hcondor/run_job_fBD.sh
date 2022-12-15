#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt/setup.sh

if [[ $2 == SV ]]
then
  python /nfs/cms/vazqueze/analisisWW/channelSV/qcd_estimation/f_BD_SV.py --year="$1" --ssos 
else
  python /nfs/cms/vazqueze/analisisWW/qcd_estimation/f_BD.py  --year="$1" --ssos
fi
