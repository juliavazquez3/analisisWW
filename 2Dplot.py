import ROOT, os, sys
from ROOT import *
import json
import argparse

#if not sys.flags.interactive: ROOT.EnableImplicitMT()

plotdir = '/nfs/cms/vazqueze/analisisWW/plots/2Dplots/' # this is where we'll save your plots
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

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
#parser.add_argument("--stack", action="store_true", default=False,
#                    help="Stack simulation or not")
#parser.add_argument("--ratio", action="store_true", default=False,
#                    help="Plot ratio or not")
#parser.add_argument("--linear", action="store_true", default=False,
#                    help="Plot linearly")

# Use like:
# python arg.py --data="No"
# python hist_plot.py --data="No" --stack --ratio

args = parser.parse_args()

if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
else: raise NameError('Incorrect data option')

## Open hists files

filePath = "/nfs/cms/vazqueze/analisisWW/hists/"

## mc files
histFile_signal = TFile.Open(filePath + "hists_v1v2v3_WW2016_semi_charm.root","READ")
histFile_WWnocharm = TFile.Open(filePath + "hists_v1v2v3_WW2016_semi_nocharm.root","READ")
histFile_WWhadronic = TFile.Open(filePath + "hists_v1v2v3_WW2016_hadronic.root","READ")
histFile_WWleptonic = TFile.Open(filePath + "hists_v1v2v3_WW2016_leptonic.root","READ")
histFile_Wjet_charm = TFile.Open(filePath + "hists_v1v2v3_Wjets2016_charm.root","READ")
histFile_Wjet_bottom = TFile.Open(filePath + "hists_v1v2v3_Wjets2016_bottom.root","READ")
histFile_Wjet_light = TFile.Open(filePath + "hists_v1v2v3_Wjets2016_light.root","READ")
histFile_Wjet_doublecharm = TFile.Open(filePath + "hists_v1v2v3_Wjets2016_doublecharm.root","READ")
histFile_ttbar = TFile.Open(filePath + "hists_v1v2v3_ttbar2016.root","READ")
histFile_ttbarlep = TFile.Open(filePath + "hists_v1v2v3_ttbarlep2016.root","READ")
histFile_ttbarhad = TFile.Open(filePath + "hists_v1v2v3_ttbarhad2016.root","READ")
histFile_DY = TFile.Open(filePath + "hists_v1v2v3_DY2016.root","READ")
histFile_DYjets = TFile.Open(filePath + "hists_v1v2v3_DYjets2016.root","READ")
histFile_ZZ = TFile.Open(filePath + "hists_v1v2v3_ZZ2016.root","READ")
histFile_WZ = TFile.Open(filePath + "hists_v1v2v3_WZ2016.root","READ")
histFile_ST1 = TFile.Open(filePath + "hists_v1v2v3_ST12016.root","READ")
histFile_ST2 = TFile.Open(filePath + "hists_v1v2v3_ST22016.root","READ")
histFile_ST3 = TFile.Open(filePath + "hists_v1v2v3_ST32016.root","READ")
histFile_ST4 = TFile.Open(filePath + "hists_v1v2v3_ST42016.root","READ")
histFile_QCD = TFile.Open(filePath + "hists_v1v2v3_QCD2016.root","READ")

## data files

if data_op == "2016":
  histFile_2016_1 = TFile.Open(filePath + "hists_v1v2v3_2016M.root","READ")
  histFile_2016_2 = TFile.Open(filePath + "hists_v1v2v3_2016E.root","READ")

if data_op == "2017":
  histFile_2017_1 = TFile.Open(filePath + "hists_2017M_range_0_200.root","READ")
  histFile_2017_2 = TFile.Open(filePath + "hists_2017M_range_200_400.root","READ")
  histFile_2017_3 = TFile.Open(filePath + "hists_2017M_range_400_600.root","READ")
  histFile_2017_4 = TFile.Open(filePath + "hists_2017M_range_600_790.root","READ")

if data_op == "2018":
  histFile_2018_1 = TFile.Open(filePath + "hists_2018M_range_0_200.root","READ")
  histFile_2018_2 = TFile.Open(filePath + "hists_2018M_range_200_456.root","READ")

histNames = []

#histNames.append("pT 2jets M")
#histNames.append("pT 2jets E")
histNames.append("deltaeta InvM M")
histNames.append("deltaeta InvM E")
histNames.append("deltaeta qgl M")
histNames.append("deltaeta qgl E")
histNames.append("deltaR InvM M")
histNames.append("deltaR InvM E")
histNames.append("qgl InvM M")
histNames.append("qgl InvM E")

samples = ["WW","Wjets","ttbar","ttbarlep","ttbarhad","DY","DYjets","WZ","ZZ","ST1","ST2","ST3","ST4","QCD"]
samples_d = ["2016","2017","2018"]
#samples = ["WW","Wjets","ttbar"]

## lumi info

files = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mc_info.json")))
processes = files.keys()

lumi = {}
xsecs = {}
nevents = {}

for p in processes:
    num_events = files[p]["events"] # Number of events
    num_files = files[p]["files"] # Number of files
    luminosity = files[p]["lumi"] # Luminosity
    #print(len(list_files))
    for s in samples:
      if files[p]["type"]==s+"2016":
        lumi[s] = luminosity
        xsecs[s] = files[p]["xsec"]
        nevents[s] = num_events

lumi["WW_semi_charm"] = lumi["WW"]
lumi["WW_semi_nocharm"] = lumi["WW"]
lumi["WW_hadronic"] = lumi["WW"]
lumi["WW_leptonic"] = lumi["WW"]
lumi["Wjets_charm"] = lumi["Wjets"]
lumi["Wjets_doublecharm"] = lumi["Wjets"]
lumi["Wjets_bottom"] = lumi["Wjets"]
lumi["Wjets_light"] = lumi["Wjets"]

print(lumi)
print(nevents)

files_d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_info.json")))
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
  #print(name)
  samples = ["WW_semi_charm","Wjets_charm","Wjets_light"]
  ## HISTS
  hS = histFile_signal.Get(name)
  hWWnC = histFile_WWnocharm.Get(name)
  hWWH = histFile_WWhadronic.Get(name)
  hWWL = histFile_WWleptonic.Get(name)
  hWJC = histFile_Wjet_charm.Get(name)
  hWJDC = histFile_Wjet_doublecharm.Get(name)
  hWJB = histFile_Wjet_bottom.Get(name)
  hWJL = histFile_Wjet_light.Get(name)
  hTT = histFile_ttbar.Get(name)
  hTTlep = histFile_ttbarlep.Get(name)
  hTThad = histFile_ttbarhad.Get(name)
  hDY = histFile_DY.Get(name)
  hDYjets = histFile_DYjets.Get(name)
  hZZ = histFile_ZZ.Get(name)
  hWZ = histFile_WZ.Get(name)
  hST1 = histFile_ST1.Get(name)
  hST2 = histFile_ST2.Get(name)
  hST3 = histFile_ST3.Get(name)
  hST4 = histFile_ST4.Get(name)
  hQCD = histFile_QCD.Get(name)

  ## data
  if data_op == "2016":
    hdataM = histFile_2016_1.Get(name)
    hdataE = histFile_2016_2.Get(name)

  if data_op == "2017":
    hdataM = histFile_2017_1.Get(name)
    hdataM.Add(histFile_2017_2.Get(name))
    hdataM.Add(histFile_2017_3.Get(name))
    hdataM.Add(histFile_2017_4.Get(name))

  if data_op == "2018":
    hdataM = histFile_2018_1.Get(name)
    hdataM.Add(histFile_2018_2.Get(name))

  histss = {}
  histss["WW_semi_charm"] = hS
  histss["WW_semi_nocharm"] = hWWnC
  histss["WW_hadronic"] = hWWH
  histss["WW_leptonic"] = hWWL
  histss["Wjets_charm"] = hWJC
  histss["Wjets_doublecharm"] = hWJDC
  histss["Wjets_bottom"] = hWJB
  histss["Wjets_light"] = hWJL
  histss["ttbar"] = hTT
  histss["ttbarlep"] = hTTlep
  histss["ttbarhad"] = hTThad
  histss["DY"] = hDY
  histss["DYjets"] = hDYjets
  histss["ZZ"] = hZZ
  histss["WZ"] = hWZ
  histss["ST1"] = hST1
  histss["ST2"] = hST2
  histss["ST3"] = hST3
  histss["ST4"] = hST4
  histss["QCD"] = hQCD

  hST = hST1
  hST.Add(hST2)
  hST.Add(hST3)
  hST.Add(hST4)
  hTTtot = hTT
  hTTtot.Add(hTTlep)
  hTTtot.Add(hTThad)
  samples = ["WW_semi_charm","Wjets_charm","Wjets_light"]
  histss["ST"] = hST
  histss["ttbar"] = hTTtot

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  for s in samples:
    c1 = TCanvas("c1","",600,400)

    histss[s].Draw("COLZ")
    print("The correlation factor for the "+name+" hist of the sample "+s+" is "+str(histss[s].GetCorrelationFactor()))

    c1.Print(plotdir+'2D_'+s+'_'+ name + ".pdf")

  c1 = TCanvas("c1","",600,400)

  #if name[-1]=='E':
    #hdataE.Draw("COLZ")
    #c1.Print(plotdir+'2D_'+'dataE_'+ name + ".pdf")

  #if name[-1]=='M':
    #hdataM.Draw("COLZ")
    #c1.Print(plotdir+'2D_'+'dataM_'+ name + ".pdf")


histFile_signal.Close()
histFile_WWnocharm.Close()
histFile_WWhadronic.Close()
histFile_WWleptonic.Close()
histFile_Wjet_charm.Close()
histFile_Wjet_doublecharm.Close()
histFile_Wjet_bottom.Close()
histFile_Wjet_light.Close()
histFile_ttbar.Close()
histFile_ttbarlep.Close()
histFile_ttbarhad.Close()
histFile_DY.Close()
histFile_DYjets.Close()
histFile_WZ.Close()
histFile_ZZ.Close()
histFile_ST1.Close()
histFile_ST2.Close()
histFile_ST3.Close()
histFile_ST4.Close()
histFile_QCD.Close()

if data_op == "2016":
  histFile_2016_1.Close()
  histFile_2016_2.Close()

if data_op == "2017":
  histFile_2017_1.Close()
  histFile_2017_2.Close()
  histFile_2017_3.Close()
  histFile_2017_4.Close()

if data_op == "2018":
  histFile_2018_1.Close()
  histFile_2018_2.Close()

