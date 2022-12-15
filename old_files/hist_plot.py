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
parser.add_argument("--data", type=string, default="No",
                    help="Select type of data and if used or not")
parser.add_argument("--stack", action="store_true", default=False,
                    help="Stack simulation or not")
parser.add_argument("--ratio", action="store_true", default=False,
                    help="Plot ratio or not")

# Use like:
# python arg.py --data="No"

args = parser.parse_args()

if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
else: raise NameError('Incorrect data option')

## Open hists files

filePath = "/nfs/cms/vazqueze/analisisWW/hists/"

histFile_signal = TFile.Open(filePath + "hists_WW_semi_charm.root","READ")
histFile_WWnocharm = TFile.Open(filePath + "hists_WW_semi_nocharm.root","READ")
histFile_WWhadronic = TFile.Open(filePath + "hists_WW_hadronic.root","READ")
histFile_WWleptonic = TFile.Open(filePath + "hists_WW_leptonic.root","READ")
histFile_Wjet = TFile.Open(filePath + "hists_Wjets.root","READ")
histFile_ttbar1 = TFile.Open(filePath + "hists_ttbar_range_0_100.root","READ")
histFile_ttbar2 = TFile.Open(filePath + "hists_ttbar_range_101_265.root","READ")
histFile_2016_1 = TFile.Open(filePath + "hists_2016_range_0_200.root","READ")
histFile_2016_2 = TFile.Open(filePath + "hists_2016_range_201_330.root","READ")
histFile_2017_1 = TFile.Open(filePath + "hists_2017_range_0_200.root","READ")
histFile_2017_2 = TFile.Open(filePath + "hists_2017_range_201_400.root","READ")
histFile_2017_3 = TFile.Open(filePath + "hists_2017_range_401_600.root","READ")
histFile_2017_4 = TFile.Open(filePath + "hists_2017_range_601_789.root","READ")
histFile_2018_1 = TFile.Open(filePath + "hists_2018_range_0_200.root","READ")
histFile_2018_2 = TFile.Open(filePath + "hists_2018_range_201_455.root","READ")

histNames = []

histNames.append("nJetGood M")
histNames.append("nJetGood E")
histNames.append("jet muon pt M")
histNames.append("jet not muon pt M")
histNames.append("jet muon eta M")
histNames.append("jet not muon eta M")
histNames.append("jet muon pt E")
histNames.append("jet not muon pt E")
histNames.append("jet muon eta E")
histNames.append("jet not muon eta E")
histNames.append("leptonM pt")
histNames.append("leptonM eta")
histNames.append("leptonE pt")
histNames.append("leptonE eta")
histNames.append("muon jet pt M")
histNames.append("muon jet eta M")
histNames.append("muon jet pt E")
histNames.append("muon jet eta E")
histNames.append("InvM 2jets M")
histNames.append("InvM 2jets E")
histNames.append("InvM jetM lepM")
histNames.append("InvM jetM lepE")
histNames.append("deltaR jetM lepM")
histNames.append("deltaR jetM lepE")
histNames.append("deltaR jetM jetNM M")
histNames.append("deltaR jetM jetNM E")
histNames.append("MET pt M")
histNames.append("MET pt E")
histNames.append("transverse massM")
histNames.append("tranverse massE")
histNames.append("tracks jetM M")
histNames.append("tracks jetNM M")
histNames.append("tracks jetM E")
histNames.append("tracks jetNM E")

samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","Wjets","ttbar","2016","2017","2018"]
#samples = ["WW","Wjets","ttbar"]

## lumi info

files = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_info.json")))
processes = files.keys()

lumi = {}
xsecs = {}
nevents = {}

for p in processes:
    # Construct the dataframes
    num_events = files[p]["events"] # Number of events
    luminosity = files[p]["lumi"] # Luminosity
    if files[p]["type"]=="WW":
        lumi["WW"] = luminosity
        xsecs["WW"] = files[p]["xsec"]
        nevents["WW"] = num_events
    if files[p]["type"]=="Wjets":
        lumi["Wjets"] = luminosity
        xsecs["Wjets"] = files[p]["xsec"]
        nevents["Wjets"] = num_events
    if files[p]["type"]=="ttbar":
        lumi["ttbar"] = luminosity
        xsecs["ttbar"] = files[p]["xsec"]
        nevents["ttbar"] = num_events
    if files[p]["type"]=="2016":
        lumi["2016"] = luminosity
        nevents["2016"] = files[p]["events_total"]
    if files[p]["type"]=="2017":
        lumi["2017"] = luminosity
        nevents["2017"] = files[p]["events_total"]
    if files[p]["type"]=="2018":
        lumi["2018"] = luminosity
        nevents["2018"] = files[p]["events_total"]

print(lumi)
print(nevents)

## Get histograms from files and draw
for name in histNames:

  c1 = TCanvas("c1","",600,400)

  if args.ratio:
    ## In case of ratio plot
    upper_pad = ROOT.TPad("upper_pad", "", 0, 0.35, 1, 1)
    lower_pad = ROOT.TPad("lower_pad", "", 0, 0, 1, 0.35)
    for p in [upper_pad, lower_pad]:
        p.SetLeftMargin(0.14)
        p.SetRightMargin(0.05)
        p.SetTickx(False)
        p.SetTicky(False)
    upper_pad.SetBottomMargin(0)
    lower_pad.SetTopMargin(0)
    lower_pad.SetBottomMargin(0.3)
 
    upper_pad.SetLogy()
    upper_pad.Draw()
    lower_pad.Draw()

  else:
    c1.SetLogy()

  ## HISTS
  hS = histFile_signal.Get(name)
  hWWnC = histFile_WWnocharm.Get(name)
  hWWH = histFile_WWhadronic.Get(name)
  hWWL = histFile_WWleptonic.Get(name)
  hWJ = histFile_Wjet.Get(name)
  hTT = histFile_ttbar1.Get(name)
  h2016 = histFile_2016_1.Get(name)
  h2017 = histFile_2017_1.Get(name)
  h2018 = histFile_2018_1.Get(name)
  
  ## Adding in case hists are in separate files
  hTT.Add(histFile_ttbar2.Get(name))
  h2016.Add(histFile_2016_2.Get(name))
  h2017.Add(histFile_2017_2.Get(name))
  h2017.Add(histFile_2017_3.Get(name))
  h2017.Add(histFile_2017_4.Get(name))
  h2018.Add(histFile_2018_2.Get(name))

  histss = {}
  histss["WW_semi_charm"] = hS
  histss["WW_semi_nocharm"] = hWWnC
  histss["WW_hadronic"] = hWWH
  histss["WW_leptonic"] = hWWL
  histss["Wjets"] = hWJ
  histss["ttbar"] = hTT
  histss["2016"] = h2016
  histss["2017"] = h2017
  histss["2018"] = h2018
  
  ## Scaling to lumi
  if data_op == "No": 
    lumi_data = 1
    hS.Scale(lumi_data/hS.Integral())
    hWWnC.Scale(lumi_data/hWWnC.Integral())
    hWWH.Scale(lumi_data/hWWH.Integral())
    hWWL.Scale(lumi_data/hWWL.Integral())
    hWJ.Scale(lumi_data/hWJ.Integral())
    hTT.Scale(lumi_data/hTT.Integral())
  else: 
    lumi_data = lumi[data_op]
    hS.Scale(lumi_data/lumi["WW"])
    hWWnC.Scale(lumi_data/lumi["WW"])
    hWWH.Scale(lumi_data/lumi["WW"])
    hWWL.Scale(lumi_data/lumi["WW"])
    hWJ.Scale(lumi_data/lumi["Wjets"])
    hTT.Scale(lumi_data/lumi["ttbar"])

  ##get maximum to plot (to avoid cutting histograms when superimposing)
  ymax = 0
  y = hS.GetMaximum()
  if y>ymax:
     ymax=y
  y = hWWnC.GetMaximum()
  if y>ymax:
     ymax=y
  y = hWWH.GetMaximum()
  if y>ymax:
     ymax=y
  y = hWWL.GetMaximum()
  if y>ymax:
     ymax=y
  y = hWJ.GetMaximum()
  if y>ymax:
     ymax=y
  y = hTT.GetMaximum()
  if y>ymax:
     ymax=y

  hS.SetMinimum(1.)
  hS.SetMaximum(2*ymax)

  gStyle.SetOptStat(kFALSE);  ## remove statistics box in histos

  if args.stack:
    hS.SetLineColor(1)
    hS.SetLineWidth(1)
    hS.SetFillColor(ROOT.TColor.GetColor(*(222,90,106)))
    hWWnC.SetLineColor(1)
    hWWnC.SetLineWidth(1)
    hWWnC.SetFillColor(ROOT.TColor.GetColor(*(246,165,42)))
    hWWH.SetLineColor(1)
    hWWH.SetLineWidth(1)
    hWWH.SetFillColor(ROOT.TColor.GetColor(*(246,239,42)))
    hWWL.SetLineColor(1)
    hWWL.SetLineWidth(1)
    hWWL.SetFillColor(ROOT.TColor.GetColor(*(153,76,0)))
    hWJ.SetLineColor(1)
    hWJ.SetLineWidth(1)
    hWJ.SetFillColor(ROOT.TColor.GetColor(*(155,152,204)))
    hTT.SetLineColor(1)
    hTT.SetLineWidth(1)
    hTT.SetFillColor(ROOT.TColor.GetColor(*(208,240,193)))

    ## Axis
    hS.GetYaxis().SetTitle("Number of events")
    hS.GetXaxis().SetTitle(name)
    hWWnC.GetYaxis().SetTitle("Number of events")
    hWWnC.GetXaxis().SetTitle(name)
    hWWH.GetYaxis().SetTitle("Number of events")
    hWWH.GetXaxis().SetTitle(name)
    hWWL.GetYaxis().SetTitle("Number of events")
    hWWL.GetXaxis().SetTitle(name)
    hWJ.GetYaxis().SetTitle("Number of events")
    hWJ.GetXaxis().SetTitle(name)
    hTT.GetYaxis().SetTitle("Number of events")
    hTT.GetXaxis().SetTitle(name)

    ## Stack creation
    if args.ratio: upper_pad.cd()
    stack = ROOT.THStack()
    stack.Add(hS)
    stack.Add(hWWnC)
    stack.Add(hWWH)
    stack.Add(hWWL)
    stack.Add(hWJ)
    stack.Add(hTT)
    stack.SetMinimum(1.)
    stack.SetMaximum(2*ymax)
    stack.Draw("HIST")

  else:
    hS.SetLineColor(2)
    hS.SetLineWidth(1)
    hWWnC.SetLineColor(6)
    hWWnC.SetLineWidth(1)
    hWWH.SetLineColor(807)
    hWWH.SetLineWidth(1)
    hWWL.SetLineColor(794)
    hWWL.SetLineWidth(1)
    hWJ.SetLineColor(4)
    hWJ.SetLineWidth(1)
    hTT.SetLineColor(3)
    hTT.SetLineWidth(1)

    ## Axis
    hS.GetYaxis().SetTitle("Number of events")
    hS.GetXaxis().SetTitle(name)
    hWWnC.GetYaxis().SetTitle("Number of events")
    hWWnC.GetXaxis().SetTitle(name)
    hWWH.GetYaxis().SetTitle("Number of events")
    hWWH.GetXaxis().SetTitle(name)
    hWWL.GetYaxis().SetTitle("Number of events")
    hWWL.GetXaxis().SetTitle(name)
    hWJ.GetYaxis().SetTitle("Number of events")
    hWJ.GetXaxis().SetTitle(name)
    hTT.GetYaxis().SetTitle("Number of events")
    hTT.GetXaxis().SetTitle(name)
    
    ## Drawing not stacked hists
    if args.ratio: upper_pad.cd()
    hS.Draw("HIST C")
    hWWnC.Draw("HIST SAME C")
    #hWWH.Draw("HIST SAME C")
    hWWL.Draw("HIST SAME C")
    hWJ.Draw("HIST SAME C")
    hTT.Draw("HIST SAME C")

  if data_op != "No":
    # Draw data
    data = histss[data_op]
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.3)
    data.SetLineWidth(1)
    data.SetLineColor(ROOT.kBlack)
    data.Draw("E SAME")

    if args.ratio:
      lower_pad.cd()
      ratio = data.Clone("ratio")
      ratio.SetLineColor(kBlack)
      ratio.SetMarkerStyle(21)
      ratio.SetTitle("")
      ratio.SetMinimum(0)
      ratio.SetMaximum(2)
      ratio.GetYaxis().SetTitle("Data/MC")
      ratio.GetXaxis().SetTitle(name)
      ratio.GetXaxis().SetLabelSize(0.08)
      ratio.GetXaxis().SetTitleSize(0.12)
      ratio.GetXaxis().SetTitleOffset(1.0)
      ratio.GetYaxis().SetLabelSize(0.05)
      ratio.GetYaxis().SetTitleSize(0.09)
      ratio.GetYaxis().CenterTitle()
      ratio.GetYaxis().SetTitleOffset(0.5)
      # Set up plot for markers and errors
      ratio.Sumw2()
      ratio.SetStats(0)
      hTotal = hS.Clone('hTotal')
      hTotal.Add(hWWnC)
      hTotal.Add(hWWH)
      hTotal.Add(hWWL)
      hTotal.Add(hWJ)
      hTotal.Add(hTT)
      ratio.Divide(hTotal)
      ratio.Draw("ep")

  ## Legends
  if args.ratio: upper_pad.cd()
  leg = TLegend(0.77,0.7,0.89,0.89)
  leg.SetBorderSize(0)
  leg.AddEntry(hS,"Signal WW semi charm","f")
  leg.AddEntry(hWWnC,"WW no charm","f")
  leg.AddEntry(hWWH,"WW hadronic","f")
  leg.AddEntry(hWWL,"WW leptonic","f")
  leg.AddEntry(hWJ,"W plus jets","f")
  leg.AddEntry(hTT,"top antitop","f")
  if data_op != "No": leg.AddEntry(data, "Data "+data_op  ,"lep")
  leg.Draw()

  ## Plot settings

  ## print image to a gif file
  if args.ratio: notation = "version1_ratio_"
  else: notation = "version1_"
  if data_op == "No": data_name = ""
  else: data_name = data_op
  c1.Print(plotdir+notation+data_name + name + ".pdf")
  #c1.Print(plotdir + 'test2_together_'  + name + ".pdf")

histFile_signal.Close()
histFile_WWnocharm.Close()
histFile_WWhadronic.Close()
histFile_WWleptonic.Close()
histFile_Wjet.Close()
histFile_ttbar1.Close()
histFile_ttbar2.Close()
histFile_2016_1.Close()
histFile_2016_2.Close()
histFile_2017_1.Close()
histFile_2017_2.Close()
histFile_2017_3.Close()
histFile_2017_4.Close()
histFile_2018_1.Close()
histFile_2018_2.Close()

