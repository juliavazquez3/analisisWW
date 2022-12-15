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
parser.add_argument("--data", type=string, default="No",
                    help="Select type of data and if used or not")
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Are these ssos plots?")
parser.add_argument("--png", action="store_true", default=False,
                    help="png format")

# Use like:
# python arg.py --data="No"
# python hist_plot.py --data="No" --stack --ratio

args = parser.parse_args()

if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
else: raise NameError('Incorrect data option')

if args.ssos: plotdir = '/nfs/cms/vazqueze/analisisWW/plots/SV/ssos/' # this is where we'll save your plots
else: plotdir = '/nfs/cms/vazqueze/analisisWW/plots/SV/'
if args.png: plotdir = '/nfs/cms/vazqueze/analisisWW/plots/png_f/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

## Open hists files

if args.ssos: filePath = "/nfs/cms/vazqueze/analisisWW/hists/SV/ssos/"
else: filePath = "/nfs/cms/vazqueze/analisisWW/hists/SV/"

if args.ssos: ssos_add = "SSOS"
else: ssos_add = ""

## mc files
histFile_signal = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_WW"+data_op+"_semi_charm.root","READ")
histFile_WWnocharm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_WW"+data_op+"_semi_nocharm.root","READ")
histFile_WWhadronic = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_WW"+data_op+"_hadronic.root","READ")
histFile_WWleptonic = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_WW"+data_op+"_leptonic.root","READ")
histFile_Wjet1_charm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets1"+data_op+"_charm.root","READ")
histFile_Wjet1_doublecharm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets1"+data_op+"_doublecharm.root","READ")
histFile_Wjet1_bottom = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets1"+data_op+"_bottom.root","READ")
histFile_Wjet1_light = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets1"+data_op+"_light.root","READ")
histFile_Wjet2_charm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets2"+data_op+"_charm.root","READ")
histFile_Wjet2_doublecharm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets2"+data_op+"_doublecharm.root","READ")
histFile_Wjet2_bottom = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets2"+data_op+"_bottom.root","READ")
histFile_Wjet2_light = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets2"+data_op+"_light.root","READ")
histFile_Wjet3_charm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets3"+data_op+"_charm.root","READ")
histFile_Wjet3_doublecharm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets3"+data_op+"_doublecharm.root","READ")
histFile_Wjet3_bottom = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets3"+data_op+"_bottom.root","READ")
histFile_Wjet3_light = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets3"+data_op+"_light.root","READ")
histFile_Wjet4_charm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets4"+data_op+"_charm.root","READ")
histFile_Wjet4_doublecharm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets4"+data_op+"_doublecharm.root","READ")
histFile_Wjet4_bottom = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets4"+data_op+"_bottom.root","READ")
histFile_Wjet4_light = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets4"+data_op+"_light.root","READ")
histFile_Wjet5_charm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets5"+data_op+"_charm.root","READ")
histFile_Wjet5_doublecharm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets5"+data_op+"_doublecharm.root","READ")
histFile_Wjet5_bottom = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets5"+data_op+"_bottom.root","READ")
histFile_Wjet5_light = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets5"+data_op+"_light.root","READ")
histFile_Wjet6_charm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets6"+data_op+"_charm.root","READ")
histFile_Wjet6_doublecharm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets6"+data_op+"_doublecharm.root","READ")
histFile_Wjet6_bottom = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets6"+data_op+"_bottom.root","READ")
histFile_Wjet6_light = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets6"+data_op+"_light.root","READ")
histFile_Wjet7_charm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets7"+data_op+"_charm.root","READ")
histFile_Wjet7_doublecharm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets7"+data_op+"_doublecharm.root","READ")
histFile_Wjet7_bottom = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets7"+data_op+"_bottom.root","READ")
histFile_Wjet7_light = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets7"+data_op+"_light.root","READ")
histFile_Wjet8_charm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets8"+data_op+"_charm.root","READ")
histFile_Wjet8_doublecharm = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets8"+data_op+"_doublecharm.root","READ")
histFile_Wjet8_bottom = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets8"+data_op+"_bottom.root","READ")
histFile_Wjet8_light = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_Wjets8"+data_op+"_light.root","READ")
histFile_ttbarC = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_ttbar"+data_op+"_charm.root","READ")
histFile_ttbarlepC = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_ttbarlep"+data_op+"_charm.root","READ")
histFile_ttbarhadC = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_ttbarhad"+data_op+"_charm.root","READ")
histFile_ttbarNC = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_ttbar"+data_op+"_nocharm.root","READ")
histFile_ttbarlepNC = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_ttbarlep"+data_op+"_nocharm.root","READ")
histFile_ttbarhadNC = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_ttbarhad"+data_op+"_nocharm.root","READ")
histFile_DY1 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_DY1"+data_op+".root","READ")
histFile_DY2 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_DY2"+data_op+".root","READ")
histFile_DY3 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_DY3"+data_op+".root","READ")
histFile_DY4 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_DY4"+data_op+".root","READ")
histFile_DY5 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_DY5"+data_op+".root","READ")
histFile_DY6 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_DY6"+data_op+".root","READ")
histFile_DY7 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_DY7"+data_op+".root","READ")
histFile_DY8 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_DY8"+data_op+".root","READ")
histFile_ZZ = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_ZZ"+data_op+".root","READ")
histFile_WZ = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_WZ"+data_op+".root","READ")
histFile_ST1 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_ST1"+data_op+".root","READ")
histFile_ST2 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_ST2"+data_op+".root","READ")
histFile_ST3 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_ST3"+data_op+".root","READ")
histFile_ST4 = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_ST4"+data_op+".root","READ")
#histFile_QCD = TFile.Open(filePath + "sv_v1v3v4v5v94"+ssos_add+"_QCD"+data_op+".root","READ")

histNames = []

histNames.append("InvM_2jets_M")
histNames.append("InvM_2jets_E")

not_rebin = ["nJetGood_M","nJetGood_E","SSOS_M","SSOS_E","nSV_E","nSV_M","sv_jet_ntracks_E","sv_jet_ntracks_M","sv_jet_pangle_E","sv_jet_pangle_M","sv_jet_ndof_E","sv_jet_ndof_M",
	"sv_jet_charge_E","sv_jet_charge_M","sv_max_ntracks_E","sv_max_ntracks_M","sv_max_pangle_E","sv_max_pangle_M","sv_max_ndof_E","sv_max_ndof_M","sv_max_charge_E","sv_max_charge_M"]

samples = ["WW","Wjets1","Wjets2","Wjets3","Wjets4","Wjets5","Wjets6","Wjets7","Wjets8","ttbar","ttbarlep","ttbarhad","DY1","DY2","DY3","DY4","DY5","DY6","DY7","DY8","WZ","ZZ","ST1","ST2","ST3","ST4","QCD"]
samples_d = ["2016","2017","2018"]
#samples = ["WW","Wjets","ttbar"]

## lumi info

files = json.load(open("/nfs/cms/vazqueze/analisisWW/mc_info_v92018.json"))
processes = files.keys()

lumi = {}
xsecs = {}
nevents = {}

for p in processes:
    num_events = files[p]["events"] # Number of events
    num_files = files[p]["files"] # Number of files
    luminosity = files[p]["lumi"] # Luminosity
    #print(files[p]["type"])
    for s in samples:
      if files[p]["type"]==s+data_op:
        #print(s)
        lumi[s] = luminosity
        xsecs[s] = files[p]["xsec"]
        nevents[s] = num_events

lumi["WW_semi_charm"] = lumi["WW"]
lumi["WW_semi_nocharm"] = lumi["WW"]
lumi["WW_hadronic"] = lumi["WW"]
lumi["WW_leptonic"] = lumi["WW"]
lumi["ttbar_charm"] = lumi["ttbar"]
lumi["ttbar_nocharm"] = lumi["ttbar"]
lumi["ttbarlep_charm"] = lumi["ttbarlep"]
lumi["ttbarlep_nocharm"] = lumi["ttbarlep"]
lumi["ttbarhad_charm"] = lumi["ttbarhad"]
lumi["ttbarhad_nocharm"] = lumi["ttbarhad"]

for i in np.arange(8)+1:
	lumi["Wjets"+str(i)+"_charm"] = lumi["Wjets"+str(i)]
	lumi["Wjets"+str(i)+"_doublecharm"] = lumi["Wjets"+str(i)]
	lumi["Wjets"+str(i)+"_bottom"] = lumi["Wjets"+str(i)]
	lumi["Wjets"+str(i)+"_light"] = lumi["Wjets"+str(i)]

print(lumi)
print(nevents)

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

## Get histograms from files and draw
for name in histNames:

  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_charm","ttbarlep_charm","ttbarhad_charm","ttbar_nocharm","ttbarlep_nocharm","ttbarhad_nocharm","WZ","ZZ","ST1","ST2","ST3","ST4",
	"Wjets1_charm","Wjets1_doublecharm","Wjets1_bottom","Wjets1_light","Wjets2_charm","Wjets2_doublecharm","Wjets2_bottom","Wjets2_light","Wjets3_charm","Wjets3_doublecharm","Wjets3_bottom","Wjets3_light",
	"Wjets4_charm","Wjets4_doublecharm","Wjets4_bottom","Wjets4_light","Wjets5_charm","Wjets5_doublecharm","Wjets5_bottom","Wjets5_light","Wjets6_charm","Wjets6_doublecharm","Wjets6_bottom","Wjets6_light",
	"Wjets7_charm","Wjets7_doublecharm","Wjets7_bottom","Wjets7_light","Wjets8_charm","Wjets8_doublecharm","Wjets8_bottom","Wjets8_light","DY1","DY2","DY3","DY4","DY5","DY6","DY7","DY8"]

  ## HISTS
  hS = histFile_signal.Get(name)
  hWWnC = histFile_WWnocharm.Get(name)
  hWWH = histFile_WWhadronic.Get(name)
  hWWL = histFile_WWleptonic.Get(name)
  hWJC1 = histFile_Wjet1_charm.Get(name)
  hWJDC1 = histFile_Wjet1_doublecharm.Get(name)
  hWJB1 = histFile_Wjet1_bottom.Get(name)
  hWJL1 = histFile_Wjet1_light.Get(name)
  hWJC2 = histFile_Wjet2_charm.Get(name)
  hWJDC2 = histFile_Wjet2_doublecharm.Get(name)
  hWJB2 = histFile_Wjet2_bottom.Get(name)
  hWJL2 = histFile_Wjet2_light.Get(name)
  hWJC3 = histFile_Wjet3_charm.Get(name)
  hWJDC3 = histFile_Wjet3_doublecharm.Get(name)
  hWJB3 = histFile_Wjet3_bottom.Get(name)
  hWJL3 = histFile_Wjet3_light.Get(name)
  hWJC4 = histFile_Wjet4_charm.Get(name)
  hWJDC4 = histFile_Wjet4_doublecharm.Get(name)
  hWJB4 = histFile_Wjet4_bottom.Get(name)
  hWJL4 = histFile_Wjet4_light.Get(name)
  hWJC5 = histFile_Wjet5_charm.Get(name)
  hWJDC5 = histFile_Wjet5_doublecharm.Get(name)
  hWJB5 = histFile_Wjet5_bottom.Get(name)
  hWJL5 = histFile_Wjet5_light.Get(name)
  hWJC6 = histFile_Wjet6_charm.Get(name)
  hWJDC6 = histFile_Wjet6_doublecharm.Get(name)
  hWJB6 = histFile_Wjet6_bottom.Get(name)
  hWJL6 = histFile_Wjet6_light.Get(name)
  hWJC7 = histFile_Wjet7_charm.Get(name)
  hWJDC7 = histFile_Wjet7_doublecharm.Get(name)
  hWJB7 = histFile_Wjet7_bottom.Get(name)
  hWJL7 = histFile_Wjet7_light.Get(name)
  hWJC8 = histFile_Wjet8_charm.Get(name)
  hWJDC8 = histFile_Wjet8_doublecharm.Get(name)
  hWJB8 = histFile_Wjet8_bottom.Get(name)
  hWJL8 = histFile_Wjet8_light.Get(name)
  hTTC = histFile_ttbarC.Get(name)
  hTTlepC = histFile_ttbarlepC.Get(name)
  hTThadC = histFile_ttbarhadC.Get(name)
  hTTNC = histFile_ttbarNC.Get(name)
  hTTlepNC = histFile_ttbarlepNC.Get(name)
  hTThadNC = histFile_ttbarhadNC.Get(name)
  hDY1 = histFile_DY1.Get(name)
  hDY2 = histFile_DY2.Get(name)
  hDY3 = histFile_DY3.Get(name)
  hDY4 = histFile_DY4.Get(name)
  hDY5 = histFile_DY5.Get(name)
  hDY6 = histFile_DY6.Get(name)
  hDY7 = histFile_DY7.Get(name)
  hDY8 = histFile_DY8.Get(name)
  hZZ = histFile_ZZ.Get(name)
  hWZ = histFile_WZ.Get(name)
  hST1 = histFile_ST1.Get(name)
  hST2 = histFile_ST2.Get(name)
  hST3 = histFile_ST3.Get(name)
  hST4 = histFile_ST4.Get(name)
  #hQCD = histFile_QCD.Get(name)

  histss = {}
  histss["WW_semi_charm"] = hS
  histss["WW_semi_nocharm"] = hWWnC
  histss["WW_hadronic"] = hWWH
  histss["WW_leptonic"] = hWWL
  histss["Wjets1_charm"] = hWJC1
  histss["Wjets1_doublecharm"] = hWJDC1
  histss["Wjets1_bottom"] = hWJB1
  histss["Wjets1_light"] = hWJL1
  histss["Wjets2_charm"] = hWJC2
  histss["Wjets2_doublecharm"] = hWJDC2
  histss["Wjets2_bottom"] = hWJB2
  histss["Wjets2_light"] = hWJL2
  histss["Wjets3_charm"] = hWJC3
  histss["Wjets3_doublecharm"] = hWJDC3
  histss["Wjets3_bottom"] = hWJB3
  histss["Wjets3_light"] = hWJL3
  histss["Wjets4_charm"] = hWJC4
  histss["Wjets4_doublecharm"] = hWJDC4
  histss["Wjets4_bottom"] = hWJB4
  histss["Wjets4_light"] = hWJL4
  histss["Wjets5_charm"] = hWJC5
  histss["Wjets5_doublecharm"] = hWJDC5
  histss["Wjets5_bottom"] = hWJB5
  histss["Wjets5_light"] = hWJL5
  histss["Wjets6_charm"] = hWJC6
  histss["Wjets6_doublecharm"] = hWJDC6
  histss["Wjets6_bottom"] = hWJB6
  histss["Wjets6_light"] = hWJL6
  histss["Wjets7_charm"] = hWJC7
  histss["Wjets7_doublecharm"] = hWJDC7
  histss["Wjets7_bottom"] = hWJB7
  histss["Wjets7_light"] = hWJL7
  histss["Wjets8_charm"] = hWJC8
  histss["Wjets8_doublecharm"] = hWJDC8
  histss["Wjets8_bottom"] = hWJB8
  histss["Wjets8_light"] = hWJL8
  histss["ttbar_charm"] = hTTC
  histss["ttbarlep_charm"] = hTTlepC
  histss["ttbarhad_charm"] = hTThadC
  histss["ttbar_nocharm"] = hTTNC
  histss["ttbarlep_nocharm"] = hTTlepNC
  histss["ttbarhad_nocharm"] = hTThadNC
  histss["DY1"] = hDY1
  histss["DY2"] = hDY2
  histss["DY3"] = hDY3
  histss["DY4"] = hDY4
  histss["DY5"] = hDY5
  histss["DY6"] = hDY6
  histss["DY7"] = hDY7
  histss["DY8"] = hDY8
  histss["ZZ"] = hZZ
  histss["WZ"] = hWZ
  histss["ST1"] = hST1
  histss["ST2"] = hST2
  histss["ST3"] = hST3
  histss["ST4"] = hST4
  #histss["QCD"] = hQCD

  if (name=="nJetGood E" or name=="nJetGood M"):
    for s in samples:
      print("Number of events before scaling for "+name[-1]+" channel in the sample "+s+" is "+str(histss[s].Integral()))

  ## Scaling to lumi
  #print(name) 
  lumi_data = lumi_d[data_op]
  for s in samples:
    #print(s)
    if s != "QCD": histss[s].Scale(lumi_data/lumi[s])
    ## Fixing single top
  hST = hST1
  hST.Add(hST2)
  hST.Add(hST3)
  hST.Add(hST4)
  hTTCtot = hTTC
  hTTCtot.Add(hTTlepC)
  hTTCtot.Add(hTThadC)
  hTTNCtot = hTTNC
  hTTNCtot.Add(hTTlepNC)
  hTTNCtot.Add(hTThadNC)
  hWJC = hWJC1
  hWJC.Add(hWJC2)
  hWJC.Add(hWJC3)
  hWJC.Add(hWJC4)
  hWJC.Add(hWJC5)
  hWJC.Add(hWJC6)
  hWJC.Add(hWJC7)
  hWJC.Add(hWJC8)
  hWJDC = hWJDC1
  hWJDC.Add(hWJDC2)
  hWJDC.Add(hWJDC3)
  hWJDC.Add(hWJDC4)
  hWJDC.Add(hWJDC5)
  hWJDC.Add(hWJDC6)
  hWJDC.Add(hWJDC7)
  hWJDC.Add(hWJDC8)
  hWJB = hWJB1
  hWJB.Add(hWJB2)
  hWJB.Add(hWJB3)
  hWJB.Add(hWJB4)
  hWJB.Add(hWJB5)
  hWJB.Add(hWJB6)
  hWJB.Add(hWJB7)
  hWJB.Add(hWJB8)
  hWJL = hWJL1
  hWJL.Add(hWJL2)
  hWJL.Add(hWJL3)
  hWJL.Add(hWJL4)
  hWJL.Add(hWJL5)
  hWJL.Add(hWJL6)
  hWJL.Add(hWJL7)
  hWJL.Add(hWJL8)
  hDY = hDY1
  hDY.Add(hDY2)
  hDY.Add(hDY3)
  hDY.Add(hDY4)
  hDY.Add(hDY5)
  hDY.Add(hDY6)
  hDY.Add(hDY7)
  hDY.Add(hDY8)
  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_nocharm","ttbar_charm","DY","WZ","ZZ","ST","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light"]
  histss["ST"] = hST
  histss["ttbar_charm"] = hTTCtot
  histss["ttbar_nocharm"] = hTTNCtot
  histss["DY"] = hDY
  histss["Wjets_charm"] = hWJC
  histss["Wjets_doublecharm"] = hWJDC
  histss["Wjets_bottom"] = hWJB
  histss["Wjets_light"] = hWJL

  if (name=="nJetGood_E" or name=="nJetGood_M"):
    for s in samples:
      print("Number of events for "+name[-1]+" channel in the sample "+s+" is "+str(histss[s].Integral()))

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  sign = hS.Clone("sign")
  backT = hWWnC.Clone('backT')
  backT.Add(hWWH)
  backT.Add(hWWL)
  backT.Add(hWJC)
  backT.Add(hWJDC)
  backT.Add(hWJL)
  backT.Add(hWJB)
  backT.Add(hTTCtot)
  backT.Add(hTTNCtot)
  backT.Add(hDY)
  backT.Add(hWZ)
  backT.Add(hZZ)
  backT.Add(hST)
  #backT.Add(hQCD)

  print("Number of signal events in Ivariant mass between 60 and 99 in the channel "+name[-1]+" is "+str(sign.Integral(21,33)))
  print("Number of backgorund events in Ivariant mass between 60 and 99 in the channel "+name[-1]+" is "+str(backT.Integral(21,33)))
  print("The ratio of events is "+str(sign.Integral(21,33)/backT.Integral(21,33))+" (signal/backgorund)")

histFile_signal.Close()
histFile_WWnocharm.Close()
histFile_WWhadronic.Close()
histFile_WWleptonic.Close()
histFile_Wjet1_charm.Close()
histFile_Wjet1_doublecharm.Close()
histFile_Wjet1_bottom.Close()
histFile_Wjet1_light.Close()
histFile_Wjet2_charm.Close()
histFile_Wjet2_doublecharm.Close()
histFile_Wjet2_bottom.Close()
histFile_Wjet2_light.Close()
histFile_Wjet3_charm.Close()
histFile_Wjet3_doublecharm.Close()
histFile_Wjet3_bottom.Close()
histFile_Wjet3_light.Close()
histFile_Wjet4_charm.Close()
histFile_Wjet4_doublecharm.Close()
histFile_Wjet4_bottom.Close()
histFile_Wjet4_light.Close()
histFile_Wjet5_charm.Close()
histFile_Wjet5_doublecharm.Close()
histFile_Wjet5_bottom.Close()
histFile_Wjet5_light.Close()
histFile_Wjet6_charm.Close()
histFile_Wjet6_doublecharm.Close()
histFile_Wjet6_bottom.Close()
histFile_Wjet6_light.Close()
histFile_Wjet7_charm.Close()
histFile_Wjet7_doublecharm.Close()
histFile_Wjet7_bottom.Close()
histFile_Wjet7_light.Close()
histFile_Wjet8_charm.Close()
histFile_Wjet8_doublecharm.Close()
histFile_Wjet8_bottom.Close()
histFile_Wjet8_light.Close()
histFile_ttbarC.Close()
histFile_ttbarlepC.Close()
histFile_ttbarhadC.Close()
histFile_ttbarNC.Close()
histFile_ttbarlepNC.Close()
histFile_ttbarhadNC.Close()
histFile_DY1.Close()
histFile_DY3.Close()
histFile_DY3.Close()
histFile_DY4.Close()
histFile_DY5.Close()
histFile_DY6.Close()
histFile_DY7.Close()
histFile_DY8.Close()
histFile_WZ.Close()
histFile_ZZ.Close()
histFile_ST1.Close()
histFile_ST2.Close()
histFile_ST3.Close()
histFile_ST4.Close()
#histFile_QCD.Close()

