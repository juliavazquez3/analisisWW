import ROOT, os, sys
from ROOT import *

import numpy as np
import awkward as ak
import uproot
import pandas as pd
import pyarrow as pa
import urllib.request


#Events_set=events.arrays(brlist)

#print(Events_set)

#events_pd=ak.to_pandas(Events_set)

#print(events_pd.head())

#events_pd.to_csv('myWW.csv',header=True)

#print(events_pd[events_pd["typeWW"]==2].groupby(level=[0,1])["GenPart_pdgId"])

#file = uproot.open("/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/260000/062085CD-8DF4-1D40-8042-63998B6A3A95.root")
file = uproot.open("/nfs/cms/vazqueze/analisisWW/old_files/WW.root")
print(file.keys())
print(file.classnames())

events = file['Events']

pdgID = events["GenPart_pdgId"].array()
motherID = events["GenPart_genPartIdxMother"].array()
isHardP = events["ishard"].array()
fromHardP = events["fromhard"].array()
firstC = events["firstC"].array()
isWplusc = events["isWplusc"].array()

ak.to_pandas(pdgID)
ak.to_pandas(motherID)
ak.to_pandas(isHardP)
ak.to_pandas(fromHardP)
ak.to_pandas(firstC)
ak.to_pandas(isWplusc)


print('--------------------------------')
print(pdgID[1271][0:22])
print(motherID[1271][0:22])
print(isHardP[1271][0:22])
print(fromHardP[1271][0:22])
print(firstC[1271][0:22])
print(isWplusc[1271])

print('--------------------------------')
print(pdgID[1581][0:22])
print(motherID[1581][0:22])
print(isHardP[1581][0:22])
print(fromHardP[1581][0:22])
print(firstC[1581][0:22])
print(isWplusc[1581])

print('--------------------------------')
print(pdgID[112701][0:22])
print(motherID[112701][0:22])
print(isHardP[112701][0:22])
print(fromHardP[112701][0:22])
print(firstC[112701][0:22])
print(isWplusc[112701])

print('--------------------------------')
print(pdgID[21803][0:22])
print(motherID[21803][0:22])
print(isHardP[21803][0:22])
print(fromHardP[21803][0:22])
print(firstC[21803][0:22])
print(isWplusc[21803])

print('--------------------------------')
print(pdgID[210704][0:22])
print(motherID[210704][0:22])
print(isHardP[210704][0:22])
print(fromHardP[210704][0:22])
print(firstC[21803][0:22])
print(isWplusc[210704])

print('--------------------------------')
print(pdgID[108011][0:18])
print(motherID[108011][0:18])
print(isHardP[108011][0:18])
print(fromHardP[108011][0:18])
print(firstC[108011][0:18])
print(isWplusc[108011])

print('--------------------------------')
print(pdgID[15811][0:22])
print(motherID[15811][0:22])
print(isHardP[15811][0:22])
print(fromHardP[15811][0:22])
print(firstC[15811][0:22])
print(isWplusc[1581])

print('--------------------------------')
print(pdgID[112702][0:22])
print(motherID[112702][0:22])
print(isHardP[112702][0:22])
print(fromHardP[112702][0:22])
print(firstC[112702][0:22])
print(isWplusc[112702])

print('--------------------------------')
print(pdgID[21804][0:22])
print(motherID[21804][0:22])
print(isHardP[21804][0:22])
print(fromHardP[21804][0:22])
print(firstC[21804][0:22])
print(isWplusc[21804])

print('--------------------------------')
print(pdgID[21071][0:22])
print(motherID[21071][0:22])
print(isHardP[21071][0:22])
print(fromHardP[21071][0:22])
print(firstC[21071][0:22])
print(isWplusc[21071])

print('--------------------------------')
print(pdgID[108011][0:18])
print(motherID[108011][0:18])
print(isHardP[108011][0:18])
print(fromHardP[108011][0:18])
print(firstC[108011][0:18])
print(isWplusc[108011])

