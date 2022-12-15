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

if args.ssos: plotdir = '/nfs/cms/vazqueze/analisisWW/plots/SV/ssos/' # this is where we'll save your plots
else: plotdir = '/nfs/cms/vazqueze/analisisWW/plots/SV/'
if args.png: plotdir = '/nfs/cms/vazqueze/analisisWW/plots/png_f/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

## Open hists files

if args.ssos: filePath = "/nfs/cms/vazqueze/analisisWW/hists/SV/ssos/cuts/"
else: filePath = "/nfs/cms/vazqueze/analisisWW/hists/SV/cuts/"

if args.ssos: ssos_add = "SSOS"
else: ssos_add = ""

term = "sv_v1v3v5v6"

## mc files
histFile_signal = TFile.Open(filePath + term+ssos_add+"_WW"+data_op+"_semi_charm.root","READ")
histFile_WWnocharm = TFile.Open(filePath + term+ssos_add+"_WW"+data_op+"_semi_nocharm.root","READ")
histFile_WWhadronic = TFile.Open(filePath + term+ssos_add+"_WW"+data_op+"_hadronic.root","READ")
histFile_WWleptonic = TFile.Open(filePath + term+ssos_add+"_WW"+data_op+"_leptonic.root","READ")
histFile_Wjet1_charm = TFile.Open(filePath + term+ssos_add+"_Wjets1"+data_op+"_charm.root","READ")
histFile_Wjet1_doublecharm = TFile.Open(filePath + term+ssos_add+"_Wjets1"+data_op+"_doublecharm.root","READ")
histFile_Wjet1_bottom = TFile.Open(filePath + term+ssos_add+"_Wjets1"+data_op+"_bottom.root","READ")
histFile_Wjet1_light = TFile.Open(filePath + term+ssos_add+"_Wjets1"+data_op+"_light.root","READ")
histFile_Wjet2_charm = TFile.Open(filePath + term+ssos_add+"_Wjets2"+data_op+"_charm.root","READ")
histFile_Wjet2_doublecharm = TFile.Open(filePath + term+ssos_add+"_Wjets2"+data_op+"_doublecharm.root","READ")
histFile_Wjet2_bottom = TFile.Open(filePath + term+ssos_add+"_Wjets2"+data_op+"_bottom.root","READ")
histFile_Wjet2_light = TFile.Open(filePath + term+ssos_add+"_Wjets2"+data_op+"_light.root","READ")
histFile_Wjet3_charm = TFile.Open(filePath + term+ssos_add+"_Wjets3"+data_op+"_charm.root","READ")
histFile_Wjet3_doublecharm = TFile.Open(filePath + term+ssos_add+"_Wjets3"+data_op+"_doublecharm.root","READ")
histFile_Wjet3_bottom = TFile.Open(filePath + term+ssos_add+"_Wjets3"+data_op+"_bottom.root","READ")
histFile_Wjet3_light = TFile.Open(filePath + term+ssos_add+"_Wjets3"+data_op+"_light.root","READ")
histFile_Wjet4_charm = TFile.Open(filePath + term+ssos_add+"_Wjets4"+data_op+"_charm.root","READ")
histFile_Wjet4_doublecharm = TFile.Open(filePath + term+ssos_add+"_Wjets4"+data_op+"_doublecharm.root","READ")
histFile_Wjet4_bottom = TFile.Open(filePath + term+ssos_add+"_Wjets4"+data_op+"_bottom.root","READ")
histFile_Wjet4_light = TFile.Open(filePath + term+ssos_add+"_Wjets4"+data_op+"_light.root","READ")
histFile_Wjet5_charm = TFile.Open(filePath + term+ssos_add+"_Wjets5"+data_op+"_charm.root","READ")
histFile_Wjet5_doublecharm = TFile.Open(filePath + term+ssos_add+"_Wjets5"+data_op+"_doublecharm.root","READ")
histFile_Wjet5_bottom = TFile.Open(filePath + term+ssos_add+"_Wjets5"+data_op+"_bottom.root","READ")
histFile_Wjet5_light = TFile.Open(filePath + term+ssos_add+"_Wjets5"+data_op+"_light.root","READ")
histFile_Wjet6_charm = TFile.Open(filePath + term+ssos_add+"_Wjets6"+data_op+"_charm.root","READ")
histFile_Wjet6_doublecharm = TFile.Open(filePath + term+ssos_add+"_Wjets6"+data_op+"_doublecharm.root","READ")
histFile_Wjet6_bottom = TFile.Open(filePath + term+ssos_add+"_Wjets6"+data_op+"_bottom.root","READ")
histFile_Wjet6_light = TFile.Open(filePath + term+ssos_add+"_Wjets6"+data_op+"_light.root","READ")
histFile_Wjet7_charm = TFile.Open(filePath + term+ssos_add+"_Wjets7"+data_op+"_charm.root","READ")
histFile_Wjet7_doublecharm = TFile.Open(filePath + term+ssos_add+"_Wjets7"+data_op+"_doublecharm.root","READ")
histFile_Wjet7_bottom = TFile.Open(filePath + term+ssos_add+"_Wjets7"+data_op+"_bottom.root","READ")
histFile_Wjet7_light = TFile.Open(filePath + term+ssos_add+"_Wjets7"+data_op+"_light.root","READ")
histFile_Wjet8_charm = TFile.Open(filePath + term+ssos_add+"_Wjets8"+data_op+"_charm.root","READ")
histFile_Wjet8_doublecharm = TFile.Open(filePath + term+ssos_add+"_Wjets8"+data_op+"_doublecharm.root","READ")
histFile_Wjet8_bottom = TFile.Open(filePath + term+ssos_add+"_Wjets8"+data_op+"_bottom.root","READ")
histFile_Wjet8_light = TFile.Open(filePath + term+ssos_add+"_Wjets8"+data_op+"_light.root","READ")
histFile_ttbarC = TFile.Open(filePath + term+ssos_add+"_ttbar"+data_op+"_charm.root","READ")
histFile_ttbarlepC = TFile.Open(filePath + term+ssos_add+"_ttbarlep"+data_op+"_charm.root","READ")
histFile_ttbarhadC = TFile.Open(filePath + term+ssos_add+"_ttbarhad"+data_op+"_charm.root","READ")
histFile_ttbarNC = TFile.Open(filePath + term+ssos_add+"_ttbar"+data_op+"_nocharm.root","READ")
histFile_ttbarlepNC = TFile.Open(filePath + term+ssos_add+"_ttbarlep"+data_op+"_nocharm.root","READ")
histFile_ttbarhadNC = TFile.Open(filePath + term+ssos_add+"_ttbarhad"+data_op+"_nocharm.root","READ")
histFile_DY1 = TFile.Open(filePath + term+ssos_add+"_DY1"+data_op+".root","READ")
histFile_DY2 = TFile.Open(filePath + term+ssos_add+"_DY2"+data_op+".root","READ")
histFile_DY3 = TFile.Open(filePath + term+ssos_add+"_DY3"+data_op+".root","READ")
histFile_DY4 = TFile.Open(filePath + term+ssos_add+"_DY4"+data_op+".root","READ")
histFile_DY5 = TFile.Open(filePath + term+ssos_add+"_DY5"+data_op+".root","READ")
histFile_DY6 = TFile.Open(filePath + term+ssos_add+"_DY6"+data_op+".root","READ")
histFile_DY7 = TFile.Open(filePath + term+ssos_add+"_DY7"+data_op+".root","READ")
histFile_DY8 = TFile.Open(filePath + term+ssos_add+"_DY8"+data_op+".root","READ")
histFile_ZZ = TFile.Open(filePath + term+ssos_add+"_ZZ"+data_op+".root","READ")
histFile_WZ = TFile.Open(filePath + term+ssos_add+"_WZ"+data_op+".root","READ")
histFile_ST1 = TFile.Open(filePath + term+ssos_add+"_ST1"+data_op+".root","READ")
histFile_ST2 = TFile.Open(filePath + term+ssos_add+"_ST2"+data_op+".root","READ")
histFile_ST3 = TFile.Open(filePath + term+ssos_add+"_ST3"+data_op+".root","READ")
histFile_ST4 = TFile.Open(filePath + term+ssos_add+"_ST4"+data_op+".root","READ")
histFile_QCD = TFile.Open(filePath + term+ssos_add+"_QCD"+data_op+".root","READ")

## data files

if data_op == "2016":
  histFile_2016_1 = TFile.Open(filePath + "hists_v1v2v3btagMedTight"+ssos_add+"_2016M.root","READ")
  histFile_2016_2 = TFile.Open(filePath + "hists_v1v2v3btagMedTight"+ssos_add+"_2016E.root","READ")

if data_op == "2017":
  histFile_2017_1 = TFile.Open(filePath + "hists_2017M_range_0_200.root","READ")
  histFile_2017_2 = TFile.Open(filePath + "hists_2017M_range_200_400.root","READ")
  histFile_2017_3 = TFile.Open(filePath + "hists_2017M_range_400_600.root","READ")
  histFile_2017_4 = TFile.Open(filePath + "hists_2017M_range_600_790.root","READ")

if data_op == "2018":
  histFile_2018_1 = TFile.Open(filePath + term+ssos_add+"_2018M.root","READ")
  histFile_2018_2 = TFile.Open(filePath + term+ssos_add+"_2018E.root","READ")

histNames = []

histNames.append("nJetGood_M")
histNames.append("nJetGood_E")
histNames.append("nSV_M")
histNames.append("nSV_E")
#histNames.append("nLooseLepton_M")
#histNames.append("nLooseLepton_E")
histNames.append("jet_sv_pt_M")
histNames.append("jet_not_sv_pt_M")
histNames.append("jet_sv_eta_M")
histNames.append("jet_not_sv_eta_M")
histNames.append("jet_sv_pt_E")
histNames.append("jet_not_sv_pt_E")
histNames.append("jet_sv_eta_E")
histNames.append("jet_not_sv_eta_E")
histNames.append("lepton_pt_M")
histNames.append("lepton_eta_M")
histNames.append("lepton_pt_E")
histNames.append("lepton_eta_E")
histNames.append("sv_jet_pt_M")
histNames.append("sv_jet_eta_M")
histNames.append("sv_jet_pt_E")
histNames.append("sv_jet_eta_E")
histNames.append("sv_jet_ntracks_M")
histNames.append("sv_jet_ntracks_E")
histNames.append("sv_jet_mass_M")
histNames.append("sv_jet_pangle_M")
histNames.append("sv_jet_mass_E")
histNames.append("sv_jet_pangle_E")
histNames.append("sv_jet_chi_M")
histNames.append("sv_jet_chi_E")
histNames.append("sv_jet_ndof_M")
histNames.append("sv_jet_ndof_E")
histNames.append("sv_jet_charge_M")
histNames.append("sv_jet_charge_E")
histNames.append("sv_jet_xy_M")
histNames.append("sv_jet_xy_E")
histNames.append("sv_jet_r_M")
histNames.append("sv_jet_r_E")
histNames.append("sv_jet_sigxy_M")
histNames.append("sv_jet_sigxy_E")
histNames.append("sv_jet_sigr_M")
histNames.append("sv_jet_sigr_E")
histNames.append("sv_jet_relpt_M")
histNames.append("sv_jet_relpt_E")
histNames.append("sv_max_pt_M")
histNames.append("sv_max_eta_M")
histNames.append("sv_max_pt_E")
histNames.append("sv_max_eta_E")
histNames.append("sv_max_ntracks_M")
histNames.append("sv_max_ntracks_E")
histNames.append("sv_max_mass_M")
histNames.append("sv_max_pangle_M")
histNames.append("sv_max_mass_E")
histNames.append("sv_max_pangle_E")
histNames.append("sv_max_chi_M")
histNames.append("sv_max_chi_E")
histNames.append("sv_max_ndof_M")
histNames.append("sv_max_ndof_E")
histNames.append("sv_max_charge_M")
histNames.append("sv_max_charge_E")
#histNames.append("sv_max_xy_M")
#histNames.append("sv_max_xy_E")
#histNames.append("sv_max_r_M")
#histNames.append("sv_max_r_E")
#histNames.append("sv_max_sigxy_M")
#histNames.append("sv_max_sigxy_E")
#histNames.append("sv_max_sigr_M")
#histNames.append("sv_max_sigr_E")
histNames.append("InvM_2jets_M")
histNames.append("InvM_2jets_E")
histNames.append("InvM_jetM_lepM")
histNames.append("InvM_jetM_lepE")
histNames.append("deltaR_jetM_lepM")
histNames.append("deltaR_jetM_lepE")
histNames.append("deltaR_jetM_jetNM_M")
histNames.append("deltaR_jetM_jetNM_E")
histNames.append("MET_pt_M")
histNames.append("MET_pt_E")
histNames.append("transverse_massM")
histNames.append("transverse_massE")
histNames.append("tracks_jetM_M")
histNames.append("tracks_jetNM_M")
histNames.append("tracks_jetM_E")
histNames.append("tracks_jetNM_E")
histNames.append("EMN_jetM_M")
histNames.append("EMC_jetM_M")
histNames.append("EMN_jetM_E")
histNames.append("EMC_jetM_E")
histNames.append("EMtotal_jetM_M")
histNames.append("EMtotal_jetM_E")
histNames.append("SSOS_M")
histNames.append("SSOS_E")
histNames.append("jet_sv_btag_M")
histNames.append("jet_sv_btag_E")
histNames.append("jet_notsv_btag_M")
histNames.append("jet_notsv_btag_E")

not_rebin = ["nJetGood_M","nJetGood_E","SSOS_M","SSOS_E","nSV_E","nSV_M","sv_jet_ntracks_E","sv_jet_ntracks_M","sv_jet_pangle_E","sv_jet_pangle_M","sv_jet_ndof_E","sv_jet_ndof_M",
	"sv_jet_charge_E","sv_jet_charge_M","sv_max_ntracks_E","sv_max_ntracks_M","sv_max_pangle_E","sv_max_pangle_M","sv_max_ndof_E","sv_max_ndof_M","sv_max_charge_E",
	"sv_max_charge_M","nLooseLepton_M","nLooseLepton_E"]

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

  if args.png: c1 = TCanvas("c1","",1200,800)
  else: c1 = TCanvas("c1","",600,400)

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
 
    if not args.linear: upper_pad.SetLogy()
    upper_pad.Draw()
    lower_pad.Draw()

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
    hdataE = histFile_2018_2.Get(name)

  #print(hdataM.Integral())

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
  histss["QCD"] = hQCD

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
    #samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_charm","ttbar_nocharm","DY","WZ","ZZ","ST","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light"]
    samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_charm","ttbar_nocharm","DY","WZ","ZZ","ST","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light","QCD"]
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
    #samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_nocharm","ttbar_charm","DY","WZ","ZZ","ST","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light"]
    samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_nocharm","ttbar_charm","DY","WZ","ZZ","ST","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light","QCD"]
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
    for s in samples:
      histss[s].SetLineWidth(1)
      histss[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      ## Axis
      histss[s].GetYaxis().SetTitle("Number of events")
      histss[s].GetXaxis().SetTitle(name)
      if args.ssos:
         if name not in not_rebin: histss[s].Rebin(5)
      histss[s].GetXaxis().SetRange(histss[s].GetXaxis().GetFirst(),histss[s].GetXaxis().GetLast()+1)

      ymax = 0
      y = histss[s].GetMaximum()
      if y>ymax: ymax=y
      y = hdataM.GetMaximum()
      if y>ymax: ymax=y
      y = hdataE.GetMaximum()
      if y>ymax: ymax=y

      hS.SetMinimum(1.)
      hS.SetMaximum(5*ymax)
      if args.linear: hS.SetMaximum(1.1*ymax)

    ## Stack creation
    if args.ratio: upper_pad.cd()
    stack = ROOT.THStack()
    for s in samples:
      stack.Add(histss[s])
    y = stack.GetMaximum()
    if y>ymax: ymax=y
    stack.SetMinimum(1.)
    stack.SetMaximum(5*ymax)
    if args.linear: stack.SetMaximum(1.2*ymax)
    stack.Draw("HIST")
   
  else:
    for s in samples:
      histss[s].SetLineColor(ROOT.TColor.GetColor(*colors[s]))
      histss[s].SetLineWidth(1)
      ## Axis
      histss[s].GetYaxis().SetTitle("Number of events")
      histss[s].GetXaxis().SetTitle(name)
      if (name[0]!='n'and name != 'sv jet ntracks M' and name != 'sv jet ntracks E' ): histss[s].Rebin(5)
      hS.SetMaximum(0.7)
    
    ## Drawing not stacked hists
    if args.ratio: upper_pad.cd()
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

  if args.stack:
    # Draw data
    if name[-1] == "M": data = hdataM
    if name[-1] == "E": data = hdataE
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.3)
    data.SetLineWidth(1)
    data.SetLineColor(ROOT.kBlack)
    if args.ssos:
         if name not in not_rebin: data.Rebin(5)
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
      hTotal.Add(hWJC)
      hTotal.Add(hWJDC)
      hTotal.Add(hWJL)
      hTotal.Add(hWJB)
      hTotal.Add(hTTCtot)
      hTotal.Add(hTTNCtot)
      hTotal.Add(hDY)
      hTotal.Add(hWZ)
      hTotal.Add(hZZ)
      hTotal.Add(hST)
      hTotal.Add(hQCD)
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
  leg.AddEntry(hQCD,"QCD","f")
  if args.stack: leg.AddEntry(data, "Data "+data_op  ,"lep")
  leg.Draw()

  ## Plot settings

  ## print image to a gif file
  if args.ratio: 
    notation = "sv_v1v3v5v6_ratio_"
    if args.linear:
      notation = "sv_v1v3v5v6_linratio_"
  else: 
    notation = "sv_v1v3v5v6_normed_"

  if args.ssos: ssos_add="ssos_"
  else: ssos_add=""

  if args.png: c1.Print(plotdir+ssos_add+notation+data_op + name + ".png")
  else: c1.Print(plotdir+ssos_add+notation+data_op + name + ".pdf")

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

