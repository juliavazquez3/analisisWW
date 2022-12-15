import ROOT, os, sys
from ROOT import *
import json
import argparse
import numpy as np

#if not sys.flags.interactive: ROOT.EnableImplicitMT()

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--jets", action="store_true", default=False,
                    help="Wjets and DY HT samples or binned in jets")
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Are these ssos plots?")
parser.add_argument("--sv", action="store_true", default=False,
                    help="SV channel")
parser.add_argument("--qcd", action="store_true", default=False,
                    help="include qcd samples")

# Use like:
# python arg.py --data="No"
# python hist_plot.py --data="No" --stack --ratio

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

## Open hists files

if args.sv:
  if args.ssos: filePath = "/nfs/cms/vazqueze/analisisWW/hists/SV/ssos/"
  else: filePath = "/nfs/cms/vazqueze/analisisWW/hists/SV/"
else:
  if args.ssos: filePath = "/nfs/cms/vazqueze/analisisWW/hists/ssos/"
  else: filePath = "/nfs/cms/vazqueze/analisisWW/hists/"

if args.ssos: ssos_add = "SSOS"
else: ssos_add = "" 

if args.sv: term = "sv_v1v3v5"
else: term = "hists_v1v2v3v4"

datayears = ["2016","2016B","2017","2018"]

samplesHT = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","Wjets1_light","Wjets1_bottom","Wjets1_charm","Wjets1_doublecharm","Wjets2_light","Wjets2_bottom","Wjets2_charm","Wjets2_doublecharm",
	"Wjets3_light","Wjets3_bottom","Wjets3_charm","Wjets3_doublecharm","Wjets4_light","Wjets4_bottom","Wjets4_charm","Wjets4_doublecharm",
	"Wjets5_light","Wjets5_bottom","Wjets5_charm","Wjets5_doublecharm","Wjets6_light","Wjets6_bottom","Wjets6_charm","Wjets6_doublecharm",
	"Wjets7_light","Wjets7_bottom","Wjets7_charm","Wjets7_doublecharm","Wjets8_light","Wjets8_bottom","Wjets8_charm","Wjets8_doublecharm",
	"ttbar_charm","ttbar_nocharm","ttbarlep_charm","ttbarlep_nocharm","ttbarhad_charm","ttbarhad_nocharm","DY1",
        "DY2","DY3","DY4","DY5","DY6","DY7","DY8","WZ","ZZ","ST1","ST2","ST3","ST4","Wjets0J_charm","Wjets0J_doublecharm","Wjets0J_light","Wjets0J_bottom",
	"Wjets1J_charm","Wjets1J_doublecharm","Wjets1J_light","Wjets1J_bottom","Wjets2J_charm","Wjets2J_doublecharm","Wjets2J_light","Wjets2J_bottom",
	"DY0J","DY1J","DY2J"]

## Adding QCD

histFile = {}

for s in samplesHT:
	## mc files
	histFile[s] = {}
	for data_op in datayears:
		if s[0:2] == "WW":
			histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:2]+data_op+s[2:]+".root","READ")
		elif (s[0:5] == "Wjets" and s[6]=="_"):
			histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:6]+data_op+s[6:]+".root","READ")
		elif (s[0:5] == "Wjets" and s[6]=="J"):
			histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:7]+data_op+s[7:]+".root","READ")
		elif s[0:6] == "ttbar_":
			histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:5]+data_op+s[5:]+".root","READ")
		elif s[0:8] == "ttbarlep":
			histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root","READ")
		elif s[0:8] == "ttbarhad":
			histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root","READ")
		else:
			histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s+data_op+".root","READ")

if args.qcd:
  histFile["QCD"] = {}
  for data_op in datayears:
    histFile["QCD"][data_op] = TFile.Open(filePath + term+ssos_add+"_QCD"+data_op+".root","READ")

histFileD = {}

for data_op in datayears:
	histFileD[data_op] = {}
	# data files
	histFileD[data_op]["M"] = TFile.Open(filePath + term+ssos_add+"_"+data_op+"M.root","READ")
	histFileD[data_op]["E"] = TFile.Open(filePath + term+ssos_add+"_"+data_op+"E.root","READ")

samples = ["WW","Wjets1","Wjets2","Wjets3","Wjets4","Wjets5","Wjets6","Wjets7","Wjets8","ttbar","ttbarlep","ttbarhad","DY1",
	"DY2","DY3","DY4","DY5","DY6","DY7","DY8","WZ","ZZ","ST1","ST2","ST3","ST4","QCD","Wjets0J","Wjets1J","Wjets2J","DY0J","DY1J","DY2J"]
samples_d = ["2016","2016B","2017","2018"]
#samples = ["WW","Wjets","ttbar"]

samples_mc = []

for s in samples:
  for y in samples_d:
    samples_mc.append(s+y)

samples_year = []

for y in samples_d:
  samples_year.append(y+"M")
  samples_year.append(y+"E")

## lumi info

lumi = {}
xsecs = {}
nevents = {}

for data_op in datayears:
	files = json.load(open("/nfs/cms/vazqueze/analisisWW/mc_info_v9"+data_op+".json"))
	processes = files.keys()
	lumi[data_op] = {}
	for p in processes:
		num_events = files[p]["events"] # Number of events
		num_files = files[p]["files"] # Number of files
		luminosity = files[p]["lumi"] # Luminosity
		#print(files[p]["type"])
		for s in samples:
			if files[p]["type"]==s+data_op:
				#print(s)
				lumi[data_op][s] = luminosity
#print(lumi)

for data_op in datayears:
	lumi[data_op]["WW_semi_charm"] = lumi[data_op]["WW"]
	lumi[data_op]["WW_semi_nocharm"] = lumi[data_op]["WW"]
	lumi[data_op]["WW_hadronic"] = lumi[data_op]["WW"]
	lumi[data_op]["WW_leptonic"] = lumi[data_op]["WW"]
	lumi[data_op]["ttbar_charm"] = lumi[data_op]["ttbar"]
	lumi[data_op]["ttbar_nocharm"] = lumi[data_op]["ttbar"]
	lumi[data_op]["ttbarlep_charm"] = lumi[data_op]["ttbarlep"]
	lumi[data_op]["ttbarlep_nocharm"] = lumi[data_op]["ttbarlep"]
	lumi[data_op]["ttbarhad_charm"] = lumi[data_op]["ttbarhad"]
	lumi[data_op]["ttbarhad_nocharm"] = lumi[data_op]["ttbarhad"]

	for i in np.arange(8)+1:
		lumi[data_op]["Wjets"+str(i)+"_charm"] = lumi[data_op]["Wjets"+str(i)]
		lumi[data_op]["Wjets"+str(i)+"_doublecharm"] = lumi[data_op]["Wjets"+str(i)]
		lumi[data_op]["Wjets"+str(i)+"_bottom"] = lumi[data_op]["Wjets"+str(i)]
		lumi[data_op]["Wjets"+str(i)+"_light"] = lumi[data_op]["Wjets"+str(i)]

## Correcting lumis for Wjets and DY
lumi["2018"]["DY0J"] = 13.55
lumi["2018"]["DY1J"] = 45.56
lumi["2018"]["DY2J"] = 37.95
lumi["2018"]["Wjets0J"] = 2.53
lumi["2018"]["Wjets1J"] = 9.78
lumi["2018"]["Wjets2J"] = 8.78

lumi["2017"]["DY0J"] = 12.27
lumi["2017"]["DY1J"] = 41.64
lumi["2017"]["DY2J"] = 39.62
lumi["2017"]["Wjets0J"] = 2.48
lumi["2017"]["Wjets1J"] = 9.75
lumi["2017"]["Wjets2J"] = 9.06

lumi["2016"]["DY0J"] = 12.5
lumi["2016"]["DY1J"] = 43.57
lumi["2016"]["DY2J"] = 36.65
lumi["2016"]["Wjets0J"] = 2.46
lumi["2016"]["Wjets1J"] = 10.17
lumi["2016"]["Wjets2J"] = 8.81

lumi["2016B"]["DY0J"] = 11.58
lumi["2016B"]["DY1J"] = 40.61
lumi["2016B"]["DY2J"] = 35.70
lumi["2016B"]["Wjets0J"] = 2.60
lumi["2016B"]["Wjets1J"] = 9.81
lumi["2016B"]["Wjets2J"] = 8.60

for data_op in datayears:
	for i in np.arange(3):
        	lumi[data_op]["Wjets"+str(i)+"J_charm"] = lumi[data_op]["Wjets"+str(i)+"J"]
        	lumi[data_op]["Wjets"+str(i)+"J_doublecharm"] = lumi[data_op]["Wjets"+str(i)+"J"]
        	lumi[data_op]["Wjets"+str(i)+"J_bottom"] = lumi[data_op]["Wjets"+str(i)+"J"]
        	lumi[data_op]["Wjets"+str(i)+"J_light"] = lumi[data_op]["Wjets"+str(i)+"J"]

#print(lumi)

files_d = json.load(open("/nfs/cms/vazqueze/analisisWW/data_info_v9.json"))
processes_d = files_d.keys()

lumi_d = {}
xsecs_d = {}
nevents_d = {}

for p in processes_d:
    num_events = files_d[p]["events"] # Number of events
    num_files = files_d[p]["files"] # Number of files
    luminosity = files_d[p]["lumi"] # Luminosity
    #print(len(list_files))
    for s in samples_d:
      if files_d[p]["type"]==s+"M":
        lumi_d[s] = luminosity
        nevents_d[s] = num_events

print(lumi_d)
print(nevents_d)

#histNames = [h for h in histNames if (h[-1]=="M")]

### Coefficients for fixing W+jets

CoefWJ_M = {}
CoefWJ_E = {}

if args.sv:
  CoefWJ_M["2016"] = 0.8107
  CoefWJ_M["2016B"] = 0.8284
  CoefWJ_M["2017"] = 0.7895
  CoefWJ_M["2018"] = 0.7632

  CoefWJ_E["2016"] = 0.7340
  CoefWJ_E["2016B"] = 0.7408
  CoefWJ_E["2017"] = 0.6895
  CoefWJ_E["2018"] = 0.6824
else:
  CoefWJ_M["2016"] = 0.8159
  CoefWJ_M["2016B"] = 0.8718
  CoefWJ_M["2017"] = 0.8776
  CoefWJ_M["2018"] = 0.8585

  CoefWJ_E["2016"] = 0.7649
  CoefWJ_E["2016B"] = 0.7873
  CoefWJ_E["2017"] = 0.7535
  CoefWJ_E["2018"] = 0.7491

#######################################################
########### Start of plot creation ####################
#######################################################

for name in ["nJetGood_E","nJetGood_M","InvM_2jets_M","InvM_2jets_E"]:

  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_charm","ttbarlep_charm","ttbarhad_charm","ttbar_nocharm","ttbarlep_nocharm","ttbarhad_nocharm","WZ","ZZ","ST1","ST2","ST3","ST4",
	"Wjets1_charm","Wjets1_doublecharm","Wjets1_bottom","Wjets1_light","Wjets2_charm","Wjets2_doublecharm","Wjets2_bottom","Wjets2_light","Wjets3_charm","Wjets3_doublecharm","Wjets3_bottom","Wjets3_light",
	"Wjets4_charm","Wjets4_doublecharm","Wjets4_bottom","Wjets4_light","Wjets5_charm","Wjets5_doublecharm","Wjets5_bottom","Wjets5_light","Wjets6_charm","Wjets6_doublecharm","Wjets6_bottom","Wjets6_light",
	"Wjets7_charm","Wjets7_doublecharm","Wjets7_bottom","Wjets7_light","Wjets8_charm","Wjets8_doublecharm","Wjets8_bottom","Wjets8_light","DY1","DY2","DY3","DY4","DY5","DY6","DY7","DY8",
	"Wjets0J_light","Wjets0J_bottom","Wjets0J_charm","Wjets0J_doublecharm","Wjets1J_light","Wjets1J_bottom","Wjets1J_charm","Wjets1J_doublecharm",
	"Wjets2J_light","Wjets2J_bottom","Wjets2J_charm","Wjets2J_doublecharm","DY0J","DY1J","DY2J"]

  ## HISTS
  histss = {}
  hdataM = {}
  hdataE = {}
  for data_op in datayears:
    histss[data_op] = {}
    for s in samples:
      histss[data_op][s] = histFile[s][data_op].Get(name)
    hdataM[data_op] = histFileD[data_op]["M"].Get(name)
    hdataE[data_op] = histFileD[data_op]["E"].Get(name)
    if (args.qcd and name[-1] == "M") : histss[data_op]["QCD"] = histFile["QCD"][data_op].Get(name)
  
  ## Scaling to lumi
  #print(name) 
  for data_op in datayears:
    lumi_data = lumi_d[data_op]
    for s in samples:
      #print(s)
      if s != "QCD": histss[data_op][s].Scale(lumi_data/lumi[data_op][s])
      ## Fixing single top
    histss[data_op]["ST"] = histss[data_op]["ST1"]
    histss[data_op]["ST"].Add(histss[data_op]["ST2"])
    histss[data_op]["ST"].Add(histss[data_op]["ST3"])
    histss[data_op]["ST"].Add(histss[data_op]["ST4"])
    histss[data_op]["ttbarT_charm"] = histss[data_op]["ttbar_charm"]
    histss[data_op]["ttbarT_charm"].Add(histss[data_op]["ttbarlep_charm"])
    histss[data_op]["ttbarT_charm"].Add(histss[data_op]["ttbarhad_charm"])
    histss[data_op]["ttbarT_nocharm"] = histss[data_op]["ttbar_nocharm"]
    histss[data_op]["ttbarT_nocharm"].Add(histss[data_op]["ttbarlep_nocharm"])
    histss[data_op]["ttbarT_nocharm"].Add(histss[data_op]["ttbarhad_nocharm"])
    histss[data_op]["Wjets_charm"] = histss[data_op]["Wjets1_charm"]
    histss[data_op]["Wjets_charm"].Add(histss[data_op]["Wjets2_charm"])
    histss[data_op]["Wjets_charm"].Add(histss[data_op]["Wjets3_charm"])
    histss[data_op]["Wjets_charm"].Add(histss[data_op]["Wjets4_charm"])
    histss[data_op]["Wjets_charm"].Add(histss[data_op]["Wjets5_charm"])
    histss[data_op]["Wjets_charm"].Add(histss[data_op]["Wjets6_charm"])
    histss[data_op]["Wjets_charm"].Add(histss[data_op]["Wjets7_charm"])
    histss[data_op]["Wjets_charm"].Add(histss[data_op]["Wjets8_charm"])
    histss[data_op]["Wjets_doublecharm"] = histss[data_op]["Wjets1_doublecharm"]
    histss[data_op]["Wjets_doublecharm"].Add(histss[data_op]["Wjets2_doublecharm"])
    histss[data_op]["Wjets_doublecharm"].Add(histss[data_op]["Wjets3_doublecharm"])
    histss[data_op]["Wjets_doublecharm"].Add(histss[data_op]["Wjets4_doublecharm"])
    histss[data_op]["Wjets_doublecharm"].Add(histss[data_op]["Wjets5_doublecharm"])
    histss[data_op]["Wjets_doublecharm"].Add(histss[data_op]["Wjets6_doublecharm"])
    histss[data_op]["Wjets_doublecharm"].Add(histss[data_op]["Wjets7_doublecharm"])
    histss[data_op]["Wjets_doublecharm"].Add(histss[data_op]["Wjets8_doublecharm"])
    histss[data_op]["Wjets_light"] = histss[data_op]["Wjets1_light"]
    histss[data_op]["Wjets_light"].Add(histss[data_op]["Wjets2_light"])
    histss[data_op]["Wjets_light"].Add(histss[data_op]["Wjets3_light"])
    histss[data_op]["Wjets_light"].Add(histss[data_op]["Wjets4_light"])
    histss[data_op]["Wjets_light"].Add(histss[data_op]["Wjets5_light"])
    histss[data_op]["Wjets_light"].Add(histss[data_op]["Wjets6_light"])
    histss[data_op]["Wjets_light"].Add(histss[data_op]["Wjets7_light"])
    histss[data_op]["Wjets_light"].Add(histss[data_op]["Wjets8_light"])
    histss[data_op]["Wjets_bottom"] = histss[data_op]["Wjets1_bottom"]
    histss[data_op]["Wjets_bottom"].Add(histss[data_op]["Wjets2_bottom"])
    histss[data_op]["Wjets_bottom"].Add(histss[data_op]["Wjets3_bottom"])
    histss[data_op]["Wjets_bottom"].Add(histss[data_op]["Wjets4_bottom"])
    histss[data_op]["Wjets_bottom"].Add(histss[data_op]["Wjets5_bottom"])
    histss[data_op]["Wjets_bottom"].Add(histss[data_op]["Wjets6_bottom"])
    histss[data_op]["Wjets_bottom"].Add(histss[data_op]["Wjets7_bottom"])
    histss[data_op]["Wjets_bottom"].Add(histss[data_op]["Wjets8_bottom"])
    histss[data_op]["DY"] = histss[data_op]["DY1"]
    histss[data_op]["DY"].Add(histss[data_op]["DY2"])
    histss[data_op]["DY"].Add(histss[data_op]["DY3"])
    histss[data_op]["DY"].Add(histss[data_op]["DY4"])
    histss[data_op]["DY"].Add(histss[data_op]["DY5"])
    histss[data_op]["DY"].Add(histss[data_op]["DY6"])
    histss[data_op]["DY"].Add(histss[data_op]["DY7"])
    histss[data_op]["DY"].Add(histss[data_op]["DY8"])
    histss[data_op]["WjetsJ_charm"] = histss[data_op]["Wjets0J_charm"]
    histss[data_op]["WjetsJ_charm"].Add(histss[data_op]["Wjets1J_charm"])
    histss[data_op]["WjetsJ_charm"].Add(histss[data_op]["Wjets2J_charm"])
    histss[data_op]["WjetsJ_doublecharm"] = histss[data_op]["Wjets0J_doublecharm"]
    histss[data_op]["WjetsJ_doublecharm"].Add(histss[data_op]["Wjets1J_doublecharm"])
    histss[data_op]["WjetsJ_doublecharm"].Add(histss[data_op]["Wjets2J_doublecharm"])
    histss[data_op]["WjetsJ_light"] = histss[data_op]["Wjets0J_light"]
    histss[data_op]["WjetsJ_light"].Add(histss[data_op]["Wjets1J_light"])
    histss[data_op]["WjetsJ_light"].Add(histss[data_op]["Wjets2J_light"])
    histss[data_op]["WjetsJ_bottom"] = histss[data_op]["Wjets0J_bottom"]
    histss[data_op]["WjetsJ_bottom"].Add(histss[data_op]["Wjets1J_bottom"])
    histss[data_op]["WjetsJ_bottom"].Add(histss[data_op]["Wjets2J_bottom"])
    histss[data_op]["DYJ"] = histss[data_op]["DY0J"]
    histss[data_op]["DYJ"].Add(histss[data_op]["DY1J"])
    histss[data_op]["DYJ"].Add(histss[data_op]["DY2J"])

  #samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_nocharm","ttbar_charm","DY","WZ","ZZ","ST","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light","QCD"]
  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbarT_nocharm","ttbarT_charm","DY","WZ","ZZ","ST","Wjets_charm","Wjets_doublecharm",
    "Wjets_bottom","Wjets_light","WjetsJ_charm","WjetsJ_doublecharm","WjetsJ_light","WjetsJ_bottom","DYJ"]

  ## Correcting W0J
  #for data_op in datayears:
  #  for s in ["WjetsJ_charm","WjetsJ_doublecharm","WjetsJ_light","WjetsJ_bottom"]:
  #    if name[-1] == "M": histss[data_op][s].Scale(CoefWJ_M[data_op])
  #    if name[-1] == "E": histss[data_op][s].Scale(CoefWJ_E[data_op])

  histT = {}
  for s in samples:
    histT[s] = histss["2016"][s]
    histT[s].Add(histss["2016B"][s])
    histT[s].Add(histss["2017"][s])
    histT[s].Add(histss["2018"][s])

  if (args.qcd and name[-1] == "M"):
    histT["QCD"] = histss["2016"]["QCD"]
    histT["QCD"].Add(histss["2016B"]["QCD"])
    histT["QCD"].Add(histss["2017"]["QCD"])
    histT["QCD"].Add(histss["2018"]["QCD"])

  histD = {}
  histD["M"] = hdataM["2016"]
  histD["M"].Add(hdataM["2016B"])
  histD["M"].Add(hdataM["2017"])
  histD["M"].Add(hdataM["2018"])
  histD["E"] = hdataE["2016"]
  histD["E"].Add(hdataE["2016B"])
  histD["E"].Add(hdataE["2017"])
  histD["E"].Add(hdataE["2018"])

  if (name=="nJetGood_E" or name=="nJetGood_M"):
    for s in samples:
      print("Number of events for "+name[-1]+" channel in the sample "+s+" is "+str(histT[s].Integral()))
    if (args.qcd and name[-1] == "M"): print("Number of events for "+name[-1]+" channel in the sample QCD is "+str(histT["QCD"].Integral()))

  if (name=="InvM_2jets_M" or name=="InvM_2jets_E"):
    for s in samples:
      axis = histT[s].GetXaxis()
      bmin = axis.FindBin(60.)
      bmax = axis.FindBin(100.)
      print("Number of events for "+name[-1]+" channel in the sample "+s+" in the (60,100) invariant mass window is "+str(histT[s].Integral(bmin,bmax)))

  if args.jets:
    samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbarT_nocharm","ttbarT_charm","WZ","ZZ","ST","DY","WjetsJ_charm","WjetsJ_doublecharm","WjetsJ_light","WjetsJ_bottom"]
  else:
    samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbarT_nocharm","ttbarT_charm","DY","WZ","ZZ","ST","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light"]

  if (name=="InvM_2jets_M" or name=="InvM_2jets_E"):
    hist_ttbarNC = histT["ttbarT_nocharm"]
    hist_ttbarC = histT["ttbarT_charm"]
    hist_ST = histT["ST"]
    hist_signal = histT[samples[0]]
    hist_backg = histT[samples[1]]
    for s in samples[2:]:
      hist_backg.Add(histT[s])
    axis = hist_backg.GetXaxis()
    bmin = axis.FindBin(60.)
    bmax = axis.FindBin(100.)
    Qsig = hist_signal.Integral(bmin,bmax)
    Qbak = hist_backg.Integral(bmin,bmax)
    QttC = hist_ttbarC.Integral(bmin,bmax)
    QttNC = hist_ttbarNC.Integral(bmin,bmax)
    QST = hist_ST.Integral(bmin,bmax)
    if name[-1]=="E": Qdat = histD["E"].Integral(bmin,bmax)
    if name[-1]=="M": Qdat = histD["M"].Integral(bmin,bmax)
    print("Number of ST mc events for "+name[-1]+" channel in the (60,100) invariant mass window is "+str(QST))
    print("Number of ttbar charm mc events for "+name[-1]+" channel in the (60,100) invariant mass window is "+str(QttC))
    print("Number of ttbar no charm mc events for "+name[-1]+" channel in the (60,100) invariant mass window is "+str(QttNC))
    print("Number of MC signal events for "+name[-1]+" channel in the (60,100) invariant mass window is "+str(Qsig))
    print("Number of MC background events for "+name[-1]+" channel in the (60,100) invariant mass window is "+str(Qbak))
    print("Number of data events for "+name[-1]+" channel in the (60,100) invariant mass window is "+str(Qdat))
    print("Ratio signal over background is "+str(Qsig/Qbak))
    print("Significance is "+str(Qsig/sqrt(Qbak)))

for s in samplesHT:
        for data_op in datayears:
                        histFile[s][data_op].Close()

for data_op in datayears:
	histFileD[data_op]["M"].Close()
	histFileD[data_op]["E"].Close()

