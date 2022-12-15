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
parser.add_argument("--stack", action="store_true", default=False,
                    help="Stack simulation or not")
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Are these ssos plots?")
parser.add_argument("--ratio", action="store_true", default=False,
                    help="Plot ratio or not")
parser.add_argument("--linear", action="store_true", default=False,
                    help="Plot linearly")
parser.add_argument("--png", action="store_true", default=False,
                    help="png format")

# Use like:
# python arg.py --data="No"
# python hist_plot.py --data="No" --stack --ratio

args = parser.parse_args()

if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
else: raise NameError('Incorrect data option')

if args.ssos: plotdir = '/nfs/cms/vazqueze/analisisWW/plots/controlR/SV/ssos/' # this is where we'll save your plots
else: plotdir = '/nfs/cms/vazqueze/analisisWW/plots/controlR/SV/'
if args.png:plotdir = '/nfs/cms/vazqueze/analisisWW/plots/png_f/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

## Open hists files

if args.ssos: filePath = "/nfs/cms/vazqueze/analisisWW/hists/controlR/SV/ssos/"
else: filePath = "/nfs/cms/vazqueze/analisisWW/hists/controlR/SV/"

if args.ssos: ssos_add = "SSOS"
else: ssos_add = "" 

## mc files
histFile_signal = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_WW"+data_op+"_semi_charm.root","READ")
histFile_WWnocharm = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_WW"+data_op+"_semi_nocharm.root","READ")
histFile_WWhadronic = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_WW"+data_op+"_hadronic.root","READ")
histFile_WWleptonic = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_WW"+data_op+"_leptonic.root","READ")
histFile_Wjet0J_charm = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets0J"+data_op+"_charm.root","READ")
histFile_Wjet0J_doublecharm = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets0J"+data_op+"_doublecharm.root","READ")
histFile_Wjet0J_bottom = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets0J"+data_op+"_bottom.root","READ")
histFile_Wjet0J_light = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets0J"+data_op+"_light.root","READ")
histFile_Wjet1J_charm = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets1J"+data_op+"_charm.root","READ")
histFile_Wjet1J_doublecharm = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets1J"+data_op+"_doublecharm.root","READ")
histFile_Wjet1J_bottom = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets1J"+data_op+"_bottom.root","READ")
histFile_Wjet1J_light = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets1J"+data_op+"_light.root","READ")
histFile_Wjet2J_charm = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets2J"+data_op+"_charm.root","READ")
histFile_Wjet2J_doublecharm = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets2J"+data_op+"_doublecharm.root","READ")
histFile_Wjet2J_bottom = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets2J"+data_op+"_bottom.root","READ")
histFile_Wjet2J_light = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_Wjets2J"+data_op+"_light.root","READ")
histFile_ttbarC = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_ttbar"+data_op+"_charm.root","READ")
histFile_ttbarlepC = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_ttbarlep"+data_op+"_charm.root","READ")
histFile_ttbarhadC = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_ttbarhad"+data_op+"_charm.root","READ")
histFile_ttbarNC = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_ttbar"+data_op+"_nocharm.root","READ")
histFile_ttbarlepNC = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_ttbarlep"+data_op+"_nocharm.root","READ")
histFile_ttbarhadNC = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_ttbarhad"+data_op+"_nocharm.root","READ")
histFile_DY0J = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_DY0J"+data_op+".root","READ")
histFile_DY1J = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_DY1J"+data_op+".root","READ")
histFile_DY2J = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_DY2J"+data_op+".root","READ")
histFile_ZZ = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_ZZ"+data_op+".root","READ")
histFile_WZ = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_WZ"+data_op+".root","READ")
histFile_ST1 = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_ST1"+data_op+".root","READ")
histFile_ST2 = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_ST2"+data_op+".root","READ")
histFile_ST3 = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_ST3"+data_op+".root","READ")
histFile_ST4 = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_ST4"+data_op+".root","READ")
#histFile_QCD = TFile.Open(filePath + "ver9_v1v3v4v5v6v635pTjets"+ssos_add+"_QCD"+data_op+".root","READ")

## data files

if data_op == "2016":
  histFile_2016_1 = TFile.Open(filePath + "sv_njet_v1v3btagMedTight"+ssos_add+"_2016M.root","READ")
  histFile_2016_2 = TFile.Open(filePath + "sv_njet_v1v3btagMedTight"+ssos_add+"_2016E.root","READ")

if data_op == "2017":
  histFile_2017_1 = TFile.Open(filePath + "sv_njet_2017M_range_0_200.root","READ")
  histFile_2017_2 = TFile.Open(filePath + "sv_njet_2017M_range_200_400.root","READ")
  histFile_2017_3 = TFile.Open(filePath + "sv_njet_2017M_range_400_600.root","READ")
  histFile_2017_4 = TFile.Open(filePath + "sv_njet_2017M_range_600_790.root","READ")

if data_op == "2018":
  histFile_2018_1 = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_2018M.root","READ")
  histFile_2018_2 = TFile.Open(filePath + "sv_njet_v1v3v4v5"+ssos_add+"_2018E.root","READ")

histNames = []

histNames.append("genW1_M")
histNames.append("genW1_E")
histNames.append("genW2_M")
histNames.append("genW2_E")

histNames = ["genW1_M","genW1_E","genW2_M","genW2_E"]

not_rebin = ["nJetGood_M","nJetGood_E","nMuoninJet_M","nMuoninJet_E","jet_muon_nmu_M","jet_muon_nmu_E","SSOS_M","SSOS_E", "jet_sv_nmu_M","jet_sv_nmu_E","sv_jet_chi_M","sv_jet_chi_E"
	"sv_max_chi_M","sv_max_chi_E","sv_jet_ndof_M","sv_jet_ndof_E","sv_max_ndof_M","sv_max_ndof_E","sv_jet_pangle_M","sv_jet_pangle_E","sv_max_pangle_M","sv_max_pangle_E",
	"sv_jet_ntracks_M","sv_jet_ntracks_E","sv_max_ntracks_M","sv_max_ntracks_E","genW1_M","genW1_E","genW2_M","genW2_E"]

not_data = ["genW1_M","genW1_E","genW2_M","genW2_E"]

samples = ["WW","Wjets1J","Wjets2J","Wjets0J","ttbar","ttbarlep","ttbarhad","DY0J","DY1J","DY2J","WZ","ZZ","ST1","ST2","ST3","ST4"]
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

## Correcting lumis for Wjets

lumi["Wjets0J"] = 2.59
lumi["Wjets1J"] = 10.03
lumi["Wjets2J"] = 8.88

for i in np.arange(3):
	lumi["Wjets"+str(i)+"J_charm"] = lumi["Wjets"+str(i)+"J"]
	lumi["Wjets"+str(i)+"J_doublecharm"] = lumi["Wjets"+str(i)+"J"]
	lumi["Wjets"+str(i)+"J_bottom"] = lumi["Wjets"+str(i)+"J"]
	lumi["Wjets"+str(i)+"J_light"] = lumi["Wjets"+str(i)+"J"]

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

  if args.png: c1 = TCanvas("c1","",1200,800)
  else: c1 = TCanvas("c1","",600,400)

  if not args.linear: c1.SetLogy()

  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_charm","ttbarlep_charm","ttbarhad_charm","ttbar_nocharm","ttbarlep_nocharm","ttbarhad_nocharm","WZ","ZZ","ST1","ST2","ST3","ST4",
	"Wjets1J_charm","Wjets1J_doublecharm","Wjets1J_bottom","Wjets1J_light","Wjets2J_charm","Wjets2J_doublecharm","Wjets2J_bottom","Wjets2J_light","Wjets0J_charm","Wjets0J_doublecharm","Wjets0J_bottom","Wjets0J_light","DY0J","DY1J","DY2J"]

  ## HISTS
  hS = histFile_signal.Get(name)
  hWWnC = histFile_WWnocharm.Get(name)
  hWWH = histFile_WWhadronic.Get(name)
  hWWL = histFile_WWleptonic.Get(name)
  hWJC1 = histFile_Wjet1J_charm.Get(name)
  hWJDC1 = histFile_Wjet1J_doublecharm.Get(name)
  hWJB1 = histFile_Wjet1J_bottom.Get(name)
  hWJL1 = histFile_Wjet1J_light.Get(name)
  hWJC2 = histFile_Wjet2J_charm.Get(name)
  hWJDC2 = histFile_Wjet2J_doublecharm.Get(name)
  hWJB2 = histFile_Wjet2J_bottom.Get(name)
  hWJL2 = histFile_Wjet2J_light.Get(name)
  hWJC3 = histFile_Wjet0J_charm.Get(name)
  hWJDC3 = histFile_Wjet0J_doublecharm.Get(name)
  hWJB3 = histFile_Wjet0J_bottom.Get(name)
  hWJL3 = histFile_Wjet0J_light.Get(name)
  hTTC = histFile_ttbarC.Get(name)
  hTTlepC = histFile_ttbarlepC.Get(name)
  hTThadC = histFile_ttbarhadC.Get(name)
  hTTNC = histFile_ttbarNC.Get(name)
  hTTlepNC = histFile_ttbarlepNC.Get(name)
  hTThadNC = histFile_ttbarhadNC.Get(name)
  hDY0J = histFile_DY0J.Get(name)
  hDY1J = histFile_DY1J.Get(name)
  hDY2J = histFile_DY2J.Get(name)
  hZZ = histFile_ZZ.Get(name)
  hWZ = histFile_WZ.Get(name)
  hST1 = histFile_ST1.Get(name)
  hST2 = histFile_ST2.Get(name)
  hST3 = histFile_ST3.Get(name)
  hST4 = histFile_ST4.Get(name)
  #hQCD = histFile_QCD.Get(name)

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
    hdataE = histFile_2018_2.Get(name)

  histss = {}
  histss["WW_semi_charm"] = hS
  histss["WW_semi_nocharm"] = hWWnC
  histss["WW_hadronic"] = hWWH
  histss["WW_leptonic"] = hWWL
  histss["Wjets1J_charm"] = hWJC1
  histss["Wjets1J_doublecharm"] = hWJDC1
  histss["Wjets1J_bottom"] = hWJB1
  histss["Wjets1J_light"] = hWJL1
  histss["Wjets2J_charm"] = hWJC2
  histss["Wjets2J_doublecharm"] = hWJDC2
  histss["Wjets2J_bottom"] = hWJB2
  histss["Wjets2J_light"] = hWJL2
  histss["Wjets0J_charm"] = hWJC3
  histss["Wjets0J_doublecharm"] = hWJDC3
  histss["Wjets0J_bottom"] = hWJB3
  histss["Wjets0J_light"] = hWJL3
  histss["ttbar_charm"] = hTTC
  histss["ttbarlep_charm"] = hTTlepC
  histss["ttbarhad_charm"] = hTThadC
  histss["ttbar_nocharm"] = hTTNC
  histss["ttbarlep_nocharm"] = hTTlepNC
  histss["ttbarhad_nocharm"] = hTThadNC
  histss["DY0J"] = hDY0J
  histss["DY1J"] = hDY1J
  histss["DY2J"] = hDY2J
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
  if not args.stack: 
    lumi_data = 1
    ## Fixing single top and top antitop
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
    hWJDC = hWJDC1                         
    hWJDC.Add(hWJDC2)
    hWJDC.Add(hWJDC3)
    hWJB = hWJB1                         
    hWJB.Add(hWJB2)
    hWJB.Add(hWJB3)
    hWJL = hWJL1                         
    hWJL.Add(hWJL2)
    hWJL.Add(hWJL3)
    hDY = hDY0J
    hDY.Add(hDY1J)
    hDY.Add(hDY2J)
    samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_charm","ttbar_nocharm","DY","WZ","ZZ","ST","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light"]
    histss["ST"] = hST
    histss["ttbar_charm"] = hTTCtot
    histss["ttbar_nocharm"] = hTTNCtot
    histss["DY"] = hDY
    histss["Wjets_charm"] = hWJC
    histss["Wjets_doublecharm"] = hWJDC
    histss["Wjets_bottom"] = hWJB
    histss["Wjets_light"] = hWJL
    for s in samples:
      if histss[s].Integral()!= 0: histss[s].Scale(lumi_data/histss[s].Integral())
  else:
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
    hWJDC = hWJDC1
    hWJDC.Add(hWJDC2)
    hWJDC.Add(hWJDC3)
    hWJB = hWJB1
    hWJB.Add(hWJB2)
    hWJB.Add(hWJB3)
    hWJL = hWJL1
    hWJL.Add(hWJL2)
    hWJL.Add(hWJL3)
    hDY = hDY0J
    hDY.Add(hDY1J)
    hDY.Add(hDY2J)
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

  colors = {}
  colors["WW_semi_charm"] = (222,90,106)
  colors["WW_semi_nocharm"] = (246,165,42)
  colors["WW_hadronic"] = (183,2,2)
  colors["WW_leptonic"] = (153,76,0)
  colors["Wjets_doublecharm"] = (51,51,255)
  colors["Wjets_charm"] = (155,152,204)
  colors["Wjets_bottom"] = (255,0,127)
  colors["Wjets_light"] = (255,153,204)
  colors["ttbar_charm"] = (204,255,153)
  colors["ttbar_nocharm"] = (120,154,86)
  colors["DY"] = (153,255,255)
  colors["DYjets"] = (0,153,153)
  colors["ZZ"] = (255,255,102)
  colors["WZ"] = (153,153,0)
  colors["ST"] = (153,51,255)
  colors["QCD"] = (0,153,76)

  if args.stack:
    ymax = 0
    for s in samples:
      histss[s].SetLineWidth(1)
      histss[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      ## Axis
      histss[s].GetYaxis().SetTitle("Number of events")
      histss[s].GetXaxis().SetTitle(name)
      if args.ssos: 
         if name not in not_rebin: histss[s].Rebin(5)
      histss[s].GetXaxis().SetRange(histss[s].GetXaxis().GetFirst(),histss[s].GetXaxis().GetLast()+1)

      y = histss[s].GetMaximum()
      if y>ymax: ymax=y

      hS.SetMinimum(1.)
      hS.SetMaximum(5*ymax)
      if args.linear: hS.SetMaximum(1.1*ymax)

    ## Stack creation
    stack = ROOT.THStack()
    for s in samples:
      stack.Add(histss[s])
    stack.SetMinimum(1.)
    stack.SetMaximum(2*ymax)
    stack.Draw("HIST")

  else:
    for s in samples:
      histss[s].SetLineColor(ROOT.TColor.GetColor(*colors[s]))
      histss[s].SetLineWidth(1)
      ## Axis
      histss[s].GetYaxis().SetTitle("Number of events")
      histss[s].GetXaxis().SetTitle(name)
      if name not in not_rebin: histss[s].Rebin(5)
      hS.SetMaximum(0.7)
    
    ## Drawing not stacked hists
    hS.Draw("HIST C")
    #hWWnC.Draw("HIST SAME C")
    #hWWH.Draw("HIST SAME C")
    #hWWL.Draw("HIST SAME C")
    #hTTCtot.Draw("HIST SAME C")
    #hTTNCtot.Draw("HIST SAME C")
    #hDY.Draw("HIST SAME C")
    #hZZ.Draw("HIST SAME C")
    #hWZ.Draw("HIST SAME C")
    #hST.Draw("HIST SAME C")
    hWJC.Draw("HIST SAME C")
    #hWJDC.Draw("HIST SAME C")
    #hWJB.Draw("HIST SAME C")
    #hWJL.Draw("HIST SAME C")
    #hQCD.Draw("HIST SAME C")

  ## Legends
  leg = TLegend(0.77,0.7,0.89,0.89)
  leg.SetBorderSize(0)
  leg.AddEntry(hS,"Signal WW semi charm","f")
  leg.AddEntry(hWWnC,"WW no charm","f")
  leg.AddEntry(hWWH,"WW hadronic","f")
  leg.AddEntry(hWWL,"WW leptonic","f")
  leg.AddEntry(hWJDC,"W plus double charm","f")
  leg.AddEntry(hWJC,"W plus c","f")
  leg.AddEntry(hWJB,"W plus bottom","f")
  leg.AddEntry(hWJL,"W plus light","f")
  leg.AddEntry(hTTCtot,"top antitop charm","f")
  leg.AddEntry(hTTNCtot,"top antitop no charm","f")
  leg.AddEntry(hDY,"Drell Yan","f")
  leg.AddEntry(hWZ,"WZ","f")
  leg.AddEntry(hZZ,"ZZ","f")
  leg.AddEntry(hST,"single top","f")
  #leg.AddEntry(hQCD,"QCD","f")
  leg.Draw()

  ## Plot settings

  ## print image to a gif file
  if args.ratio: 
    notation = "CR_sv_v1v2v3v4v5_ratio_"
    if args.linear:
      notation = "CR_sv_v1v2v3v4v5_linratio_"
  else: 
    notation = "CR_sv_v1v2v3v4v5_normed_"

  if args.ssos: ssos_add="ssos_"
  else: ssos_add=""

  if args.png: c1.Print(plotdir+ssos_add+notation+data_op + name + ".png")
  else: c1.Print(plotdir+ssos_add+notation+data_op + name + ".pdf")

histFile_signal.Close()
histFile_WWnocharm.Close()
histFile_WWhadronic.Close()
histFile_WWleptonic.Close()
histFile_Wjet0J_charm.Close()
histFile_Wjet0J_doublecharm.Close()
histFile_Wjet0J_bottom.Close()
histFile_Wjet0J_light.Close()
histFile_Wjet1J_charm.Close()
histFile_Wjet1J_doublecharm.Close()
histFile_Wjet1J_bottom.Close()
histFile_Wjet1J_light.Close()
histFile_Wjet2J_charm.Close()
histFile_Wjet2J_doublecharm.Close()
histFile_Wjet2J_bottom.Close()
histFile_Wjet2J_light.Close()
histFile_ttbarC.Close()
histFile_ttbarlepC.Close()
histFile_ttbarhadC.Close()
histFile_ttbarNC.Close()
histFile_ttbarlepNC.Close()
histFile_ttbarhadNC.Close()
histFile_DY0J.Close()
histFile_DY1J.Close()
histFile_DY2J.Close()
histFile_WZ.Close()
histFile_ZZ.Close()
histFile_ST1.Close()
histFile_ST2.Close()
histFile_ST3.Close()
histFile_ST4.Close()
#histFile_QCD.Close()

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

