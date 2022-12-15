import ROOT, os, sys
from ROOT import *
import json
import argparse

#if not sys.flags.interactive: ROOT.EnableImplicitMT()

plotdir = '/nfs/cms/vazqueze/analisisWW/plots/' # this is where we'll save your plots
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
parser.add_argument("--stack", action="store_true", default=False,
                    help="Stack simulation or not")
parser.add_argument("--linear", action="store_true", default=False,
                    help="Plot linearly")

# Use like:
# python arg.py --data="No"
# python hist_plot.py --data="No" --stack --ratio

args = parser.parse_args()

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

histNames = []

histNames.append("nJetGood M")
histNames.append("nJetGood E")
#histNames.append("nMuoninJet M")
#histNames.append("nMuoninJet E")
#histNames.append("jet muon pt M")
#histNames.append("jet muon nmu M")
#histNames.append("jet not muon pt M")
#histNames.append("jet muon eta M")
#histNames.append("jet not muon eta M")
histNames.append("jet notmuon qgl M")
#histNames.append("jet notmuon nmu M")
#histNames.append("jet notmuon mass M")
#histNames.append("jet muon pt E")
#histNames.append("jet muon nmu E")
#histNames.append("jet not muon pt E")
#histNames.append("jet muon eta E")
#histNames.append("jet not muon eta E")
histNames.append("jet notmuon qgl E")
#histNames.append("jet notmuon nmu E")
#histNames.append("jet notmuon mass E")
#histNames.append("lepton pt M")
#histNames.append("lepton eta M")
#histNames.append("lepton pt E")
#histNames.append("lepton eta E")
histNames.append("muon jet pt M")
histNames.append("muon jet eta M")
histNames.append("muon jet pt E")
histNames.append("muon jet eta E")
#histNames.append("InvM 2jets M")
#histNames.append("InvM 2jets E")
#histNames.append("InvM jetM lepM")
#histNames.append("InvM jetM lepE")
histNames.append("deltaR jetM lepM")
histNames.append("deltaR jetM lepE")
histNames.append("deltaR jetM jetNM M")
histNames.append("deltaR jetM jetNM E")
histNames.append("deltapt jetM jetNM M")
histNames.append("deltapt jetM jetNM E")
histNames.append("deltaeta jetM jetNM M")
histNames.append("deltaeta jetM jetNM E")
histNames.append("deltaphi jetM jetNM M")
histNames.append("deltaphi jetM jetNM E")
histNames.append("muon jet relpt M")
histNames.append("muon jet relpt E")
histNames.append("muon jet sigxy M")
histNames.append("muon jet sigxy E")
histNames.append("muon jet sigz M")
histNames.append("muon jet sigz E")
histNames.append("muon jet sigr M")
histNames.append("muon jet sigr E")
histNames.append("pT sum M")
histNames.append("pT sum E")
histNames.append("pT product M")
histNames.append("pT product E")
histNames.append("pT sum 2J M")
histNames.append("pT sum 2J E")
histNames.append("pT proy M")
histNames.append("pT proy E")
histNames.append("deltaR lep 2jets M")
histNames.append("deltaR lep 2jets E")
histNames.append("eta 2jets M")
histNames.append("eta 2jets E")
histNames.append("pt 2jets M")
histNames.append("pt 2jets E")
histNames.append("deltaphi MET 2jets M")
histNames.append("deltaphi MET 2jets E")
histNames.append("deltaphi lephad M")
histNames.append("deltaphi lephad E")
histNames.append("deltaR lephad M")
histNames.append("deltaR lephad E")
histNames.append("deltaphi lep 2jets M")
histNames.append("deltaphi lep 2jets E")
histNames.append("deltaeta lephad M")
histNames.append("deltaeta lephad E")
histNames.append("deltaeta lep 2jets M")
histNames.append("deltaeta lep 2jets E")
histNames.append("deltapt lephad M")
histNames.append("deltapt lephad E")
histNames.append("deltapt lep 2jets M")
histNames.append("deltapt lep 2jets E")
#histNames.append("jet muon btag M")
#histNames.append("jet muon btag E")
#histNames.append("jet notmuon btag M")
#histNames.append("jet notmuon btag E")

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


## Get histograms from files and draw
for name in histNames:

  c1 = TCanvas("c1","",600,400)
  if args.stack and (not args.linear): c1.SetLogy()

  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light","ttbar_semi","ttbar_lep","ttbar_had"]
  ## HISTS
  hS = histFile_signal.Get(name)
  hWWnC = histFile_WWnocharm.Get(name)
  hWWH = histFile_WWhadronic.Get(name)
  hWWL = histFile_WWleptonic.Get(name)
  hWJC = histFile_Wjet_charm.Get(name)
  hWJDC = histFile_Wjet_doublecharm.Get(name)
  hWJB = histFile_Wjet_bottom.Get(name)
  hWJL = histFile_Wjet_light.Get(name)
  hTTS = histFile_ttbar.Get(name)
  hTTL = histFile_ttbarlep.Get(name)
  hTTH = histFile_ttbarhad.Get(name)

  histss = {}
  histss["WW_semi_charm"] = hS
  histss["WW_semi_nocharm"] = hWWnC
  histss["WW_hadronic"] = hWWH
  histss["WW_leptonic"] = hWWL
  histss["Wjets_charm"] = hWJC
  histss["Wjets_doublecharm"] = hWJDC
  histss["Wjets_bottom"] = hWJB
  histss["Wjets_light"] = hWJL
  histss["ttbar_semi"] = hTTS
  histss["ttbar_lep"] = hTTL
  histss["ttbar_had"] = hTTH

  ## Scaling to lumi
  if args.stack:
    #print(name)
    lumi_data = 100
    for s in samples:
      #print(s)
      if s != "QCD": histss[s].Scale(lumi_data/lumi[s])
    hTT = hTTS
    hTT.Add(hTTL)
    hTT.Add(hTTH)
    samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light","ttbar"]
    histss["ttbar"]=hTT
  else: 
    lumi_data = 1
    hTT = hTTS
    hTT.Add(hTTL)
    hTT.Add(hTTH)
    samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light","ttbar"]
    histss["ttbar"]=hTT
    for s in samples:
      if histss[s].Integral()!= 0: histss[s].Scale(lumi_data/histss[s].Integral())

  if (name=="nJetGood E" or name=="nJetGood M"):
    for s in samples:
      print("Number of events for "+name[-1]+" channel in the sample "+s+" is "+str(histss[s].Integral()))

  ##get maximum to plot (to avoid cutting histograms when superimposing)
  #print(name)
  ymax = 0
  for s in samples:
   #print(s)
   y = histss[s].GetMaximum()
   if y>ymax:
     ymax=y

  if args.stack:
    hS.SetMinimum(1.)
    hS.SetMaximum(5*ymax)
    if args.linear: hS.SetMaximum(1.1*ymax)
  else:
    hS.SetMaximum(1.5*ymax)

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
  colors["ttbar"] = (204,255,153)

  if args.stack:
    for s in samples:
      histss[s].SetLineWidth(1)
      histss[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      ## Axis
      histss[s].GetYaxis().SetTitle("Number of events")
      histss[s].GetXaxis().SetTitle(name)

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
      histss[s].SetLineWidth(2)
      ## Axis
      histss[s].GetYaxis().SetTitle("Number of events")
      histss[s].GetXaxis().SetTitle(name)
      histss[s].Rebin(5)
    
    ## Drawing not stacked hists
    hS.Draw("HIST C")
    #hWWnC.Draw("HIST SAME C")
    #hWWH.Draw("HIST SAME C")
    #hWWL.Draw("HIST SAME C")
    hWJC.Draw("HIST SAME C")
    #hWJDC.Draw("HIST SAME C")
    #hWJB.Draw("HIST SAME C")
    #hWJL.Draw("HIST SAME C")
    hTT.Draw("HIST SAME C")

  ## Legends
  leg = TLegend(0.77,0.7,0.89,0.89)
  leg.SetBorderSize(0)
  leg.AddEntry(hS,"Signal WW semi charm","f")
  leg.AddEntry(hWWnC,"WW no charm","f")
  leg.AddEntry(hWWH,"WW hadronic","f")
  leg.AddEntry(hWWL,"WW leptonic","f")
  leg.AddEntry(hWJDC,"W plus double c","f")
  leg.AddEntry(hWJC,"W plus c","f")
  leg.AddEntry(hWJB,"W plus bottom","f")
  leg.AddEntry(hWJL,"W plus light","f")
  leg.AddEntry(hTT,"top antitop","f")
  leg.Draw()

  ## Plot settings

  ## print image to a gif file

  if args.linear:
      notation = "wplusc_v1v2v3_linratio_"
  elif args.stack:
      notation = "wplusc_v1v2v3_ratio_"
  else: 
      notation = "wplusc_v1v2v3_normed_"
  c1.Print(plotdir+notation+ name + ".pdf")


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

