import ROOT, os, sys
from ROOT import *
import json
import numpy as np
#import awkward as ak
import uproot


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

### signal

filePath = "/nfs/cms/vazqueze/analisisWW/hists/"

histFile_signal = TFile.Open(filePath + "hists_WW.root","READ")
histFile_Wjet = TFile.Open(filePath + "hists_Wjets.root","READ")
histFile_ttbar1 = TFile.Open(filePath + "hists_ttbar_range_0_100.root","READ")
histFile_ttbar2 = TFile.Open(filePath + "hists_ttbar_range_101_265.root","READ")

histNames = []

histNames.append("nJetGood")
histNames.append("jet muon pt")
histNames.append("jet not muon pt")
histNames.append("jet muon eta")
histNames.append("jet not muon eta")
histNames.append("leptonM pt")
histNames.append("leptonM eta")
histNames.append("leptonE pt")
histNames.append("leptonE eta")
histNames.append("muon jet pt")
histNames.append("muon jet eta")
histNames.append("InvM 2jets")
histNames.append("InvM jetM lepM")
histNames.append("InvM jetM lepE")
histNames.append("MET pt")
histNames.append("transverse massM")
histNames.append("tranverse massE")
histNames.append("tracks jetM")
histNames.append("tracks jetNM")

samples = ["WW","Wjets","ttbar"]

def drawHistos(hist_signal,hist_Wjet, hist_ttbar):

  ##get maximum to plot (to avoid cutting histograms when superimposing)
  ymax = 0
  y = hist_signal.GetMaximum()
  if y>ymax:
     ymax=y
  y = hist_Wjet.GetMaximum()
  if y>ymax: 
     ymax=y
  y = hist_ttbar.GetMaximum()
  if y>ymax: 
     ymax=y
  y = hist_ttbar.GetMaximum()
  if y>ymax:
     ymax=y

  hist_signal.SetMaximum(1.1*ymax);
  hist_signal.SetMinimum(0.0001);

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos
  hist_signal.SetLineColor(kRed)
  hist_signal.SetLineWidth(1)
  hist_signal.Draw()
  hist_Wjet.SetLineColor(kBlue)
  hist_Wjet.SetLineWidth(1)
  hist_Wjet.Draw("SAME")
  hist_ttbar.SetLineColor(kGreen)
  hist_ttbar.SetLineWidth(1)
  hist_ttbar.Draw("SAME")


## lumi info

files = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_info_xrootd.json")))
processes = files.keys()

lumi = {}
xsecs = {}
nevents = {}

for p in processes:
    # Construct the dataframes
    num_events = files[p]["events"] # Number of events
    cross_sec = files[p]["events"] # Cross section
    luminosity = files[p]["events"] # Luminosity
    if files[p]["type"]=="WW":
        lumi["WW"] = luminosity
        xsecs["WW"] = cross_sec
        nevents["WW"] = num_events
    if files[p]["type"]=="Wjets":
        lumi["Wjets"] = luminosity
        xsecs["Wjets"] = cross_sec
        nevents["Wjets"] = num_events
    if files[p]["type"]=="ttbar":
        lumi["ttbar"] = luminosity
        xsecs["ttbar"] = cross_sec
        nevents["ttbar"] = num_events


## Get histograms from files and draw
for name in histNames:

  c1 = TCanvas("c1","",600,400)
  c1.SetLogy()

  hS = histFile_signal.Get(name)
  hWJ = histFile_Wjet.Get(name)
  hTT = histFile_ttbar1.Get(name)

  ## Adding in case hists are in separate files
  hTT.Add(histFile_ttbar2.Get(name))

  histss = {}
  histss["WW"] = hS
  histss["Wjets"] = hWJ
  histss["ttbar"] = hTT

  ## Correction of luminosity and scaling
  nevents_cor = {}
  lumi_cor = {}
  scales = {}
  for s in samples:
    nevents_cor[s] = histss[s].GetEntries()
    lumi_cor[s] = nevents_cor[s]/(xsecs[s]*1000) 
    scales[s] = 100/lumi_cor[s]

  ## Scaling to lumi 
  hS.Scale(scales["WW"])
  hWJ.Scale(scales["Wjets"])
  hTT.Scale(scales["ttbar"])

  ## Axis
  hS.GetYaxis().SetTitle("Number of events")
  hS.GetXaxis().SetTitle(name)
  hWJ.GetYaxis().SetTitle("Number of events")
  hWJ.GetXaxis().SetTitle(name)
  hTT.GetYaxis().SetTitle("Number of events")
  hTT.GetXaxis().SetTitle(name)

  drawHistos(hS, hWJ, hTT)

  ## Legends
  leg = TLegend(0.77,0.7,0.89,0.89)
  leg.SetBorderSize(0)
  leg.AddEntry(hS,"Signal WW","l")
  leg.AddEntry(hWJ,"W plus jets","l")
  leg.AddEntry(hTT,"top antitop","l")
  leg.Draw()

  ## Plot settings

  ## print image to a gif file
  c1.Print(plotdir+"test_together_" + name + ".pdf")


histFile_signal.Close()
histFile_Wjet.Close()
histFile_ttbar1.Close()
histFile_ttbar2.Close()
