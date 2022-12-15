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
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Are these ssos plots?")
parser.add_argument("--sv", action="store_true", default=False,
                    help="SV channel")

# Use like:
# python arg.py --data="No"
# python hist_plot.py --data="No" --stack --ratio

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

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

datayears = ["2016","2016B","2017","2018"]

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

if args.sv: 
  print('SV channel')
else:
  print('SL channel')

####################################################################################################################################################################

### Obtention of factors of conversion for W+jets

print('W+jets coefficients')

if args.sv:
  if args.ssos:
    filePath = "/nfs/cms/vazqueze/analisisWW/hists/controlR/wjets/SV/ssos/"
  else:
    filePath = "/nfs/cms/vazqueze/analisisWW/hists/controlR/wjets/SV/"
else:
  if args.ssos:
    filePath = "/nfs/cms/vazqueze/analisisWW/hists/controlR/wjets/ssos/"
  else:
    filePath = "/nfs/cms/vazqueze/analisisWW/hists/controlR/wjets/"

if args.ssos: ssos_add = "SSOS"
else: ssos_add = ""

if args.sv: term = "wjetscr_sv_v1v3v5"
else: term = "wjetscr_v1v2v3v4"

datayears = ["2016","2016B","2017","2018"]

samplesHT = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic",
        "ttbar_charm","ttbar_nocharm","ttbarlep_charm","ttbarlep_nocharm","ttbarhad_charm","ttbarhad_nocharm",
        "WZ","ZZ","ST1","ST2","ST3","ST4","Wjets0J_charm","Wjets0J_doublecharm","Wjets0J_light","Wjets0J_bottom",
        "Wjets1J_charm","Wjets1J_doublecharm","Wjets1J_light","Wjets1J_bottom","Wjets2J_charm","Wjets2J_doublecharm","Wjets2J_light","Wjets2J_bottom",
        "DY0J","DY1J","DY2J"]

## Adding QCD

histFileWjets = {}

for s in samplesHT:
        ## mc files
        histFileWjets[s] = {}
        for data_op in datayears:
                if s[0:2] == "WW":
                        histFileWjets[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:2]+data_op+s[2:]+".root","READ")
                elif (s[0:5] == "Wjets" and s[6]=="J"):
                        histFileWjets[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:7]+data_op+s[7:]+".root","READ")
                elif s[0:6] == "ttbar_":
                        histFileWjets[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:5]+data_op+s[5:]+".root","READ")
                elif s[0:8] == "ttbarlep":
                        histFileWjets[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root","READ")
                elif s[0:8] == "ttbarhad":
                        histFileWjets[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root","READ")
                else:
                     	histFileWjets[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s+data_op+".root","READ")

histFileDWjets = {}

for data_op in datayears:
        histFileDWjets[data_op] = {}
        # data files
        histFileDWjets[data_op]["M"] = TFile.Open(filePath + term+ssos_add+"_"+data_op+"M.root","READ")
        histFileDWjets[data_op]["E"] = TFile.Open(filePath + term+ssos_add+"_"+data_op+"E.root","READ")

## Coef dict 

CoefM_WJ = {}
CoefE_WJ = {}

for name in ["transverse_massM","nJetGood_E","nJetGood_M"]:

  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_charm","ttbarlep_charm","ttbarhad_charm","ttbar_nocharm","ttbarlep_nocharm","ttbarhad_nocharm","WZ","ZZ","ST1","ST2","ST3","ST4",
        "Wjets0J_light","Wjets0J_bottom","Wjets0J_charm","Wjets0J_doublecharm","Wjets1J_light","Wjets1J_bottom","Wjets1J_charm","Wjets1J_doublecharm",
        "Wjets2J_light","Wjets2J_bottom","Wjets2J_charm","Wjets2J_doublecharm","DY0J","DY1J","DY2J"]

  ## HISTS
  histssWJ = {}
  hdataMWJ = {}
  hdataEWJ = {}
  for data_op in datayears:
    histssWJ[data_op] = {}
    for s in samples:
      histssWJ[data_op][s] = histFileWjets[s][data_op].Get(name)
    hdataMWJ[data_op] = histFileDWjets[data_op]["M"].Get(name)
    hdataEWJ[data_op] = histFileDWjets[data_op]["E"].Get(name)

  ## Scaling to lumi
  #print(name)
  for data_op in datayears:
    lumi_data = lumi_d[data_op]
    for s in samples:
      #print(s)
      if s != "QCD": histssWJ[data_op][s].Scale(lumi_data/lumi[data_op][s])
      ## Fixing single top
    histssWJ[data_op]["ST"] = histssWJ[data_op]["ST1"]
    histssWJ[data_op]["ST"].Add(histssWJ[data_op]["ST2"])
    histssWJ[data_op]["ST"].Add(histssWJ[data_op]["ST3"])
    histssWJ[data_op]["ST"].Add(histssWJ[data_op]["ST4"])
    histssWJ[data_op]["ttbarT_charm"] = histssWJ[data_op]["ttbar_charm"]
    histssWJ[data_op]["ttbarT_charm"].Add(histssWJ[data_op]["ttbarlep_charm"])
    histssWJ[data_op]["ttbarT_charm"].Add(histssWJ[data_op]["ttbarhad_charm"])
    histssWJ[data_op]["ttbarT_nocharm"] = histssWJ[data_op]["ttbar_nocharm"]
    histssWJ[data_op]["ttbarT_nocharm"].Add(histssWJ[data_op]["ttbarlep_nocharm"])
    histssWJ[data_op]["ttbarT_nocharm"].Add(histssWJ[data_op]["ttbarhad_nocharm"])
    histssWJ[data_op]["WjetsJ_charm"] = histssWJ[data_op]["Wjets0J_charm"]
    histssWJ[data_op]["WjetsJ_charm"].Add(histssWJ[data_op]["Wjets1J_charm"])
    histssWJ[data_op]["WjetsJ_charm"].Add(histssWJ[data_op]["Wjets2J_charm"])
    histssWJ[data_op]["WjetsJ_doublecharm"] = histssWJ[data_op]["Wjets0J_doublecharm"]
    histssWJ[data_op]["WjetsJ_doublecharm"].Add(histssWJ[data_op]["Wjets1J_doublecharm"])
    histssWJ[data_op]["WjetsJ_doublecharm"].Add(histssWJ[data_op]["Wjets2J_doublecharm"])
    histssWJ[data_op]["WjetsJ_light"] = histssWJ[data_op]["Wjets0J_light"]
    histssWJ[data_op]["WjetsJ_light"].Add(histssWJ[data_op]["Wjets1J_light"])
    histssWJ[data_op]["WjetsJ_light"].Add(histssWJ[data_op]["Wjets2J_light"])
    histssWJ[data_op]["WjetsJ_bottom"] = histssWJ[data_op]["Wjets0J_bottom"]
    histssWJ[data_op]["WjetsJ_bottom"].Add(histssWJ[data_op]["Wjets1J_bottom"])
    histssWJ[data_op]["WjetsJ_bottom"].Add(histssWJ[data_op]["Wjets2J_bottom"])
    histssWJ[data_op]["DYJ"] = histssWJ[data_op]["DY0J"]
    histssWJ[data_op]["DYJ"].Add(histssWJ[data_op]["DY1J"])
    histssWJ[data_op]["DYJ"].Add(histssWJ[data_op]["DY2J"])

  #samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_nocharm","ttbar_charm","DY","WZ","ZZ","ST","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light","QCD"]
  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbarT_nocharm","ttbarT_charm","WZ","ZZ","ST","WjetsJ_charm","WjetsJ_doublecharm","WjetsJ_light","WjetsJ_bottom","DYJ"]

  if name=="nJetGood_E":
    for data_op in datayears:
      hTotal = histssWJ[data_op]["WW_semi_charm"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histssWJ[data_op][s])
      hD = hdataEWJ[data_op]
      CoefE_WJ[data_op] = hD.Integral()/hTotal.Integral()
      print("Coefficient E for year"+str(data_op)+" is "+str(CoefE_WJ[data_op]))

  #if name=="transverse_massM":
  #  for data_op in datayears:
  #    hTotal = histssWJ[data_op]["WW_semi_charm"].Clone('hTotal')
  #    for s in samples[1:]:
  #      hTotal.Add(histssWJ[data_op][s])
  #    hD = hdataMWJ[data_op]
  #    print("Events for data are "+str(hD.Integral(30,50))+" and for MC "+str(hTotal.Integral(30,50)))
  #    CoefM_WJ[data_op]	= hD.Integral(30,50)/hTotal.Integral(30,50)
  #    print("Coefficient M for year"+str(data_op)+" is "+str(CoefM_WJ[data_op]))
  #    print("If we used as in E it would be"+str(hD.Integral()/hTotal.Integral()))

  if name=="nJetGood_M":
    for data_op in datayears:
      hTotal = histssWJ[data_op]["WW_semi_charm"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histssWJ[data_op][s])
      hD = hdataMWJ[data_op]
      CoefM_WJ[data_op] = hD.Integral()/hTotal.Integral()
      print("Coefficient M for year"+str(data_op)+" is "+str(CoefM_WJ[data_op]))

print('Coef ended')

#############################################################################################################################################################################         

######## Coefs for top antitop

print('Top antitop coefficients')

if args.sv:
  if args.ssos:
    filePath = "/nfs/cms/vazqueze/analisisWW/hists/controlR/ttbar/SV/ssos/"
  else:
    filePath = "/nfs/cms/vazqueze/analisisWW/hists/controlR/ttbar/SV/"
else:
  if args.ssos:
    filePath = "/nfs/cms/vazqueze/analisisWW/hists/controlR/ttbar/ssos/"
  else:
    filePath = "/nfs/cms/vazqueze/analisisWW/hists/controlR/ttbar/"

if args.ssos: ssos_add = "SSOS"
else: ssos_add = ""

if args.sv: term = "ttbarcr_sv_v1v3v5"
else: term = "ttbarcr_v1v2v3v4"

datayears = ["2016","2016B","2017","2018"]

samplesHT = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic",
        "ttbar_charm","ttbar_nocharm","ttbarlep_charm","ttbarlep_nocharm","ttbarhad_charm","ttbarhad_nocharm",
        "WZ","ZZ","ST1","ST2","ST3","ST4","Wjets0J_charm","Wjets0J_doublecharm","Wjets0J_light","Wjets0J_bottom",
        "Wjets1J_charm","Wjets1J_doublecharm","Wjets1J_light","Wjets1J_bottom","Wjets2J_charm","Wjets2J_doublecharm","Wjets2J_light","Wjets2J_bottom",
        "DY0J","DY1J","DY2J"]

## Adding QCD

histFileTT = {}

for s in samplesHT:
        ## mc files
        histFileTT[s] = {}
        for data_op in datayears:
                if s[0:2] == "WW":
                        histFileTT[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:2]+data_op+s[2:]+".root","READ")
                elif (s[0:5] == "Wjets" and s[6]=="J"):
                        histFileTT[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:7]+data_op+s[7:]+".root","READ")
                elif s[0:6] == "ttbar_":
                        histFileTT[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:5]+data_op+s[5:]+".root","READ")
                elif s[0:8] == "ttbarlep":
                        histFileTT[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root","READ")
                elif s[0:8] == "ttbarhad":
                        histFileTT[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:8]+data_op+s[8:]+".root","READ")
                else:
                     	histFileTT[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s+data_op+".root","READ")

histFileDTT = {}

for data_op in datayears:
        histFileDTT[data_op] = {}
        # data files
        histFileDTT[data_op]["M"] = TFile.Open(filePath + term+ssos_add+"_"+data_op+"M.root","READ")
        histFileDTT[data_op]["E"] = TFile.Open(filePath + term+ssos_add+"_"+data_op+"E.root","READ")

## Coef dict

CoefM_TT = {}
CoefE_TT = {}

for name in ["nJetGood_E","nJetGood_M"]:

  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_charm","ttbarlep_charm","ttbarhad_charm","ttbar_nocharm",
        "ttbarlep_nocharm","ttbarhad_nocharm","WZ","ZZ","ST1","ST2","ST3","ST4",
        "Wjets0J_light","Wjets0J_bottom","Wjets0J_charm","Wjets0J_doublecharm","Wjets1J_light","Wjets1J_bottom","Wjets1J_charm","Wjets1J_doublecharm",
        "Wjets2J_light","Wjets2J_bottom","Wjets2J_charm","Wjets2J_doublecharm","DY0J","DY1J","DY2J"]

  ## HISTS
  histssTT = {}
  hdataMTT = {}
  hdataETT = {}
  for data_op in datayears:
    histssTT[data_op] = {}
    for s in samples:
      histssTT[data_op][s] = histFileTT[s][data_op].Get(name)
    hdataMTT[data_op] = histFileDTT[data_op]["M"].Get(name)
    hdataETT[data_op] = histFileDTT[data_op]["E"].Get(name)

  ## Scaling to lumi
  #print(name)
  for data_op in datayears:
    lumi_data = lumi_d[data_op]
    for s in samples:
      #print(s)
      if s != "QCD": histssTT[data_op][s].Scale(lumi_data/lumi[data_op][s])
      ## Fixing single top
    histssTT[data_op]["ST"] = histssTT[data_op]["ST1"]
    histssTT[data_op]["ST"].Add(histssTT[data_op]["ST2"])
    histssTT[data_op]["ST"].Add(histssTT[data_op]["ST3"])
    histssTT[data_op]["ST"].Add(histssTT[data_op]["ST4"])
    histssTT[data_op]["ttbarT_charm"] = histssTT[data_op]["ttbar_charm"]
    histssTT[data_op]["ttbarT_charm"].Add(histssTT[data_op]["ttbarlep_charm"])
    histssTT[data_op]["ttbarT_charm"].Add(histssTT[data_op]["ttbarhad_charm"])
    histssTT[data_op]["ttbarT_nocharm"] = histssTT[data_op]["ttbar_nocharm"]
    histssTT[data_op]["ttbarT_nocharm"].Add(histssTT[data_op]["ttbarlep_nocharm"])
    histssTT[data_op]["ttbarT_nocharm"].Add(histssTT[data_op]["ttbarhad_nocharm"])
    histssTT[data_op]["WjetsJ_charm"] = histssTT[data_op]["Wjets0J_charm"]
    histssTT[data_op]["WjetsJ_charm"].Add(histssTT[data_op]["Wjets1J_charm"])
    histssTT[data_op]["WjetsJ_charm"].Add(histssTT[data_op]["Wjets2J_charm"])
    histssTT[data_op]["WjetsJ_doublecharm"] = histssTT[data_op]["Wjets0J_doublecharm"]
    histssTT[data_op]["WjetsJ_doublecharm"].Add(histssTT[data_op]["Wjets1J_doublecharm"])
    histssTT[data_op]["WjetsJ_doublecharm"].Add(histssTT[data_op]["Wjets2J_doublecharm"])
    histssTT[data_op]["WjetsJ_light"] = histssTT[data_op]["Wjets0J_light"]
    histssTT[data_op]["WjetsJ_light"].Add(histssTT[data_op]["Wjets1J_light"])
    histssTT[data_op]["WjetsJ_light"].Add(histssTT[data_op]["Wjets2J_light"])
    histssTT[data_op]["WjetsJ_bottom"] = histssTT[data_op]["Wjets0J_bottom"]
    histssTT[data_op]["WjetsJ_bottom"].Add(histssTT[data_op]["Wjets1J_bottom"])
    histssTT[data_op]["WjetsJ_bottom"].Add(histssTT[data_op]["Wjets2J_bottom"])
    histssTT[data_op]["DYJ"] = histssTT[data_op]["DY0J"]
    histssTT[data_op]["DYJ"].Add(histssTT[data_op]["DY1J"])
    histssTT[data_op]["DYJ"].Add(histssTT[data_op]["DY2J"])

  #samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_nocharm","ttbar_charm","DY","WZ","ZZ","ST","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light","QCD"]
  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbarT_nocharm","ttbarT_charm","WZ","ZZ","ST","WjetsJ_charm","WjetsJ_doublecharm","WjetsJ_light","WjetsJ_bottom","DYJ"]

  ## Correcting W0J
  for data_op in datayears:
    for s in ["WjetsJ_charm","WjetsJ_doublecharm","WjetsJ_light","WjetsJ_bottom"]:
      if name[-1] == "M": histssTT[data_op][s].Scale(CoefM_WJ[data_op])
      if name[-1] == "E": histssTT[data_op][s].Scale(CoefE_WJ[data_op])

  if name=="nJetGood_E":
    for data_op in datayears:
      hTotal = histssTT[data_op]["WW_semi_charm"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histssTT[data_op][s])
      hD = hdataETT[data_op]
      CoefE_TT[data_op] = hD.Integral()/hTotal.Integral()
      print("Coefficient E for year"+str(data_op)+" is "+str(CoefE_TT[data_op]))

  if name=="nJetGood_M":
    for data_op in datayears:
      hTotal = histssTT[data_op]["WW_semi_charm"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histssTT[data_op][s])
      hD = hdataMTT[data_op]
      CoefM_TT[data_op] = hD.Integral()/hTotal.Integral()
      print("Coefficient M for year"+str(data_op)+" is "+str(CoefM_TT[data_op]))

print('Coef ended')

################################################################################################################################################################################

for s in samplesHT:
        for data_op in datayears:
                        histFileWjets[s][data_op].Close()

for data_op in datayears:
	histFileDWjets[data_op]["M"].Close()
	histFileDWjets[data_op]["E"].Close()

