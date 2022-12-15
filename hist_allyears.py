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
parser.add_argument("--stack", action="store_true", default=False,
                    help="Stack simulation or not")
parser.add_argument("--jets", action="store_true", default=False,
                    help="Wjets and DY HT samples or binned in jets")
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Are these ssos plots?")
parser.add_argument("--ratio", action="store_true", default=False,
                    help="Plot ratio or not")
parser.add_argument("--linear", action="store_true", default=False,
                    help="Plot linearly")
parser.add_argument("--png", action="store_true", default=False,
                    help="png format")
parser.add_argument("--qcd", action="store_true", default=False,
                    help="include qcd samples")

# Use like:
# python arg.py --data="No"
# python hist_plot.py --data="No" --stack --ratio

args = parser.parse_args()

#if (args.data == "No" or args.data == "2016" or args.data == "2017" or args.data == "2018"): data_op = str(args.data)
#else: raise NameError('Incorrect data option')

if args.ssos: plotdir = '/nfs/cms/vazqueze/analisisWW/plots/ssos/' # this is where we'll save your plots
else: plotdir = '/nfs/cms/vazqueze/analisisWW/plots/'
if args.png:plotdir = '/nfs/cms/vazqueze/analisisWW/plots/png_f/'

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

## Open hists files

if args.ssos: filePath = "/nfs/cms/vazqueze/analisisWW/hists/ssos/"
else: filePath = "/nfs/cms/vazqueze/analisisWW/hists/"

if args.ssos: ssos_add = "SSOS"
else: ssos_add = "" 

term = "hists_v1v2v3v4"

datayears = ["2016","2016B","2017","2018"]

samplesHT = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","Wjets1_light","Wjets1_bottom","Wjets1_charm","Wjets1_doublecharm","Wjets2_light","Wjets2_bottom","Wjets2_charm","Wjets2_doublecharm",
	"Wjets3_light","Wjets3_bottom","Wjets3_charm","Wjets3_doublecharm","Wjets4_light","Wjets4_bottom","Wjets4_charm","Wjets4_doublecharm",
	"Wjets5_light","Wjets5_bottom","Wjets5_charm","Wjets5_doublecharm","Wjets6_light","Wjets6_bottom","Wjets6_charm","Wjets6_doublecharm",
	"Wjets7_light","Wjets7_bottom","Wjets7_charm","Wjets7_doublecharm","Wjets8_light","Wjets8_bottom","Wjets8_charm","Wjets8_doublecharm",
	"ttbar_charm","ttbar_nocharm","ttbarlep_charm","ttbarlep_nocharm","ttbarhad_charm","ttbarhad_nocharm","DY1",
        "DY2","DY3","DY4","DY5","DY6","DY7","DY8","WZ","ZZ","ST1_charm","ST2_charm","ST3_charm","ST4_charm","ST1_nocharm","ST2_nocharm","ST3_nocharm","ST4_nocharm",
	"Wjets0J_charm","Wjets0J_doublecharm","Wjets0J_light","Wjets0J_bottom",
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
		elif s[0:2] == "ST":
			histFile[s][data_op] = TFile.Open(filePath + term+ssos_add+"_"+s[0:3]+data_op+s[3:]+".root","READ")
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

histNames = []

histNames.append("nJetGood_M")
histNames.append("nJetGood_E")
histNames.append("nMuoninJet_M")
histNames.append("nMuoninJet_E")
#histNames.append("nLooseLepton_M")
#histNames.append("nLooseLepton_E")
histNames.append("jet_muon_pt_M")
histNames.append("jet_muon_nmu_M")
histNames.append("jet_not_muon_pt_M")
histNames.append("jet_muon_eta_M")
histNames.append("jet_not_muon_eta_M")
histNames.append("jet_notmuon_qgl_M")
histNames.append("jet_notmuon_nmu_M")
histNames.append("jet_notmuon_mass_M")
histNames.append("jet_muon_pt_E")
histNames.append("jet_muon_nmu_E")
histNames.append("jet_not_muon_pt_E")
histNames.append("jet_muon_eta_E")
histNames.append("jet_not_muon_eta_E")
histNames.append("jet_notmuon_qgl_E")
histNames.append("jet_notmuon_nmu_E")
histNames.append("jet_notmuon_mass_E")
histNames.append("lepton_pt_M")
histNames.append("lepton_eta_M")
histNames.append("lepton_pt_E")
histNames.append("lepton_eta_E")
histNames.append("muon_jet_pt_M")
histNames.append("muon_jet_eta_M")
histNames.append("muon_jet_pt_E")
histNames.append("muon_jet_eta_E")
histNames.append("InvM_2jets_M")
histNames.append("InvM_2jets_E")
histNames.append("InvM_jetM_lepM")
histNames.append("InvM_jetM_lepE")
histNames.append("deltaR_jetM_lepM")
histNames.append("deltaR_jetM_lepE")
histNames.append("deltaR_jetM_jetNM_M")
histNames.append("deltaR_jetM_jetNM_E")
histNames.append("deltapt_jetM_jetNM_M")
histNames.append("deltapt_jetM_jetNM_E")
histNames.append("deltaeta_jetM_jetNM_M")
histNames.append("deltaeta_jetM_jetNM_E")
histNames.append("deltaphi_jetM_jetNM_M")
histNames.append("deltaphi_jetM_jetNM_E")
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
histNames.append("InvM_muon_jet_M")
histNames.append("muon_jet_mva_M")
histNames.append("muon_jet_mva_E")
histNames.append("muon_jet_tight_M")
histNames.append("muon_jet_tight_E")
histNames.append("muon_jet_relpt_M")
histNames.append("muon_jet_relpt_E")
histNames.append("muon_jet_sigxy_M")
histNames.append("muon_jet_sigxy_E")
histNames.append("muon_jet_sigz_M")
histNames.append("muon_jet_sigz_E")
histNames.append("muon_jet_sigr_M")
histNames.append("muon_jet_sigr_E")
histNames.append("SSOS_M")
histNames.append("SSOS_E")
histNames.append("lepton_iso_E")
histNames.append("lepton_mva_E")
histNames.append("pT_sum_M")
histNames.append("pT_sum_E")
histNames.append("pT_product_M")
histNames.append("pT_product_E")
histNames.append("deltaR_lep_2jets_M")
histNames.append("deltaR_lep_2jets_E")
histNames.append("eta_2jets_M")
histNames.append("eta_2jets_E")
histNames.append("pt_2jets_M")
histNames.append("pt_2jets_E")
histNames.append("pt_Wlep_M")
histNames.append("pt_Wlep_E")
histNames.append("pT_proy_M")
histNames.append("pT_proy_E")
histNames.append("deltaphi_MET_2jets_M")
histNames.append("deltaphi_MET_2jets_E")
histNames.append("deltaphi_lephad_M")
histNames.append("deltaphi_lephad_E")
histNames.append("deltaR_lep_jet1_M")
histNames.append("deltaR_lep_jet2_M")
histNames.append("deltaPhi_lep_jet1_M")
histNames.append("deltaPhi_lep_jet2_M")
histNames.append("deltaEta_lep_jet1_M")
histNames.append("deltaEta_lep_jet2_M")
histNames.append("deltaR_lep_jet1_E")
histNames.append("deltaR_lep_jet2_E")
histNames.append("deltaPhi_lep_jet1_E")
histNames.append("deltaPhi_lep_jet2_E")
histNames.append("deltaEta_lep_jet1_E")
histNames.append("deltaEta_lep_jet2_E")
histNames.append("jet_muon_btag_M")
histNames.append("jet_muon_btag_E")
histNames.append("jet_notmuon_btag_M")
histNames.append("jet_notmuon_btag_E")
histNames.append("jet_notmuon_deeptagG_M")
histNames.append("jet_notmuon_deeptagG_E")
histNames.append("jet_notmuon_deeptagC_M")
histNames.append("jet_notmuon_deeptagC_E")

#histNames = ["jet_notmuon_deeptagG_M","jet_notmuon_deeptagG_E","jet_notmuon_deeptagC_M","jet_notmuon_deeptagC_E"]

not_rebin = ["nJetGood_M","nJetGood_E","nMuoninJet_M","nMuoninJet_E","jet_muon_nmu_M","jet_muon_nmu_E","SSOS_M","SSOS_E","nLooseLepton_M","nLooseLepton_E", "muon_jet_tight_M","muon_jet_tight_E"]

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

	for i in np.arange(4)+1:
		lumi[data_op]["ST"+str(i)+"_charm"] = lumi[data_op]["ST"+str(i)]
		lumi[data_op]["ST"+str(i)+"_nocharm"] = lumi[data_op]["ST"+str(i)]

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

### Coefficients for fixing W+jets and ttbar

CoefWJ_M = {}
CoefWJ_E = {}

CoefWJ_M["2016"] = 0.8524
CoefWJ_M["2016B"] = 0.8678
CoefWJ_M["2017"] = 0.9159
CoefWJ_M["2018"] = 0.8935

CoefWJ_E["2016"] = 0.7649
CoefWJ_E["2016B"] = 0.7873
CoefWJ_E["2017"] = 0.7535
CoefWJ_E["2018"] = 0.7479

CoefTT_M = {}
CoefTT_E = {}

CoefTT_M["2016"] = 0.9096
CoefTT_M["2016B"] = 0.9527
CoefTT_M["2017"] = 1.1457
CoefTT_M["2018"] = 1.0136

CoefTT_E["2016"] = 0.9653
CoefTT_E["2016B"] = 0.8520
CoefTT_E["2017"] = 0.8075
CoefTT_E["2018"] = 0.9789


#######################################################
########### Start of plot creation ####################
#######################################################

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

  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_charm","ttbarlep_charm","ttbarhad_charm","ttbar_nocharm","ttbarlep_nocharm","ttbarhad_nocharm",
	"WZ","ZZ","ST1_charm","ST2_charm","ST3_charm","ST4_charm","ST1_nocharm","ST2_nocharm","ST3_nocharm","ST4_nocharm",
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
    histss[data_op]["ST_charm"] = histss[data_op]["ST1_charm"]
    histss[data_op]["ST_charm"].Add(histss[data_op]["ST2_charm"])
    histss[data_op]["ST_charm"].Add(histss[data_op]["ST3_charm"])
    histss[data_op]["ST_charm"].Add(histss[data_op]["ST4_charm"])
    histss[data_op]["ST_nocharm"] = histss[data_op]["ST1_nocharm"]
    histss[data_op]["ST_nocharm"].Add(histss[data_op]["ST2_nocharm"])
    histss[data_op]["ST_nocharm"].Add(histss[data_op]["ST3_nocharm"])
    histss[data_op]["ST_nocharm"].Add(histss[data_op]["ST4_nocharm"])
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

  #samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbar_nocharm","ttbar_charm","DY","WZ","ZZ","ST_charm","ST_nocharm","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light","QCD"]
  samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbarT_nocharm","ttbarT_charm","DY","WZ","ZZ","ST_charm","ST_nocharm","Wjets_charm","Wjets_doublecharm",
    "Wjets_bottom","Wjets_light","WjetsJ_charm","WjetsJ_doublecharm","WjetsJ_light","WjetsJ_bottom","DYJ"]

  # Correcting W0J
  #for data_op in datayears:
  #  for s in ["WjetsJ_charm","WjetsJ_doublecharm","WjetsJ_light","WjetsJ_bottom"]:
  #    if name[-1] == "M": histss[data_op][s].Scale(CoefWJ_M[data_op])
  #    if name[-1] == "E": histss[data_op][s].Scale(CoefWJ_E[data_op])

  # Correcting ttbar
  #for data_op in datayears:
  #  for s in ["ttbarT_charm","ttbarT_nocharm"]:
  #    if name[-1] == "M": histss[data_op][s].Scale(CoefTT_M[data_op])
  #    if name[-1] == "E": histss[data_op][s].Scale(CoefTT_E[data_op])

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
  colors["WjetsJ_doublecharm"] = (51,51,255)
  colors["WjetsJ_charm"] = (155,152,204)
  colors["WjetsJ_bottom"] = (255,0,127)
  colors["WjetsJ_light"] = (255,153,204)
  colors["ttbarT_charm"] = (204,255,153)
  colors["ttbarT_nocharm"] = (120,154,86)
  colors["DY"] = (153,255,255)
  colors["DYJ"] = (153,255,255)
  colors["DYjets"] = (0,153,153)
  colors["ZZ"] = (255,255,102)
  colors["WZ"] = (153,153,0)
  colors["ST_charm"] = (153,51,255)
  colors["ST_nocharm"] = (190,153,228)
  colors["QCD"] = (0,153,76)

  if args.jets:
    samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbarT_nocharm","ttbarT_charm","WZ","ZZ","ST_charm","ST_nocharm","DYJ","WjetsJ_charm","WjetsJ_doublecharm","WjetsJ_light","WjetsJ_bottom"]
    samplesB = ["WW_hadronic","WW_leptonic","ttbarT_nocharm","WZ","ZZ","ST_nocharm","DYJ","WjetsJ_charm","WjetsJ_doublecharm","WjetsJ_light","WjetsJ_bottom"]
  else:
    samples = ["WW_semi_charm","WW_semi_nocharm","WW_hadronic","WW_leptonic","ttbarT_nocharm","ttbarT_charm","DY","WZ","ZZ","ST_charm","ST_nocharm","Wjets_charm","Wjets_doublecharm","Wjets_bottom","Wjets_light"]
    samplesB = ["WW_hadronic","WW_leptonic","ttbarT_nocharm","WZ","ZZ","ST_nocharm","DY","Wjets_charm","Wjets_doublecharm","Wjets_light","Wjets_bottom"]

  if args.stack:
    for s in samples:
      histT[s].SetLineWidth(1)
      histT[s].SetFillColor(ROOT.TColor.GetColor(*colors[s]))
      ## Axis
      histT[s].GetYaxis().SetTitle("Number of events")
      histT[s].GetXaxis().SetTitle(name)
      if args.ssos:
         if name not in not_rebin: histT[s].Rebin(2)
      histT[s].GetXaxis().SetRange(histT[s].GetXaxis().GetFirst(),histT[s].GetXaxis().GetLast()+1)

      ymax = 0
      y = histT[s].GetMaximum()
      if y>ymax: ymax=y
      if name[-1] == "M": y = histD["M"].GetMaximum()
      if name[-1] == "E": y = histD["E"].GetMaximum()
      if y>ymax: ymax=y

      histT["WW_semi_charm"].SetMinimum(1.)
      histT["WW_semi_charm"].SetMaximum(5*ymax)
      if args.linear: histT["WW_semi_charm"].SetMaximum(1.3*ymax)

    if (args.qcd and name[-1] == "M"): 
      histT["QCD"].SetLineWidth(1)
      histT["QCD"].SetFillColor(ROOT.TColor.GetColor(*colors["QCD"]))
      ## Axis
      histT["QCD"].GetYaxis().SetTitle("Number of events")
      histT["QCD"].GetXaxis().SetTitle(name)
      if args.ssos:
         if name not in not_rebin: histT["QCD"].Rebin(2)
      histT["QCD"].GetXaxis().SetRange(histT["QCD"].GetXaxis().GetFirst(),histT["QCD"].GetXaxis().GetLast()+1)

    ## Stack creation
    if args.ratio: upper_pad.cd()
    stack = ROOT.THStack()
    for s in samples:
      stack.Add(histT[s])
    if (args.qcd and name[-1] == "M"): stack.Add(histT["QCD"])
    y = stack.GetMaximum()
    if y>ymax: ymax=y
    stack.SetMinimum(1.)
    stack.SetMaximum(5*ymax)
    if args.linear: stack.SetMaximum(1.3*ymax)
    stack.Draw("HIST")

  if args.stack:
    # Draw data
    if name[-1] == "M": data = histD["M"]
    if name[-1] == "E": data = histD["E"]
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.3)
    data.SetLineWidth(1)
    data.SetLineColor(ROOT.kBlack)
    if args.ssos:
         if name not in not_rebin: data.Rebin(2)
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
      hTotal = histT["WW_semi_charm"].Clone('hTotal')
      for s in samples[1:]:
        hTotal.Add(histT[s])
      if (args.qcd and name[-1] == "M"): hTotal.Add(histT["QCD"])
      ratio.Divide(hTotal)
      ratio.Draw("ep")

      ratioB = data.Clone("ratio")
      ratioB.SetLineColor(kRed)
      ratioB.SetMarkerStyle(21)
      ratioB.SetMinimum(0)
      ratioB.SetMaximum(2)
      ratioB.GetYaxis().SetTitle("Data/MC backg only")
      ratioB.GetYaxis().SetLabelSize(0.05)
      ratioB.GetYaxis().SetTitleSize(0.09)
      ratioB.GetYaxis().CenterTitle()
      ratioB.GetYaxis().SetTitleOffset(0.5)
      # Set up plot for markers and errors
      ratioB.Sumw2()
      ratioB.SetStats(0)
      hTotalB = histT["WW_semi_nocharm"].Clone('hTotalB')
      for s in samplesB:
        hTotalB.Add(histT[s])
      if (args.qcd and name[-1] == "M"): hTotalB.Add(histT["QCD"])
      ratioB.Divide(hTotalB)
      ratioB.Draw("ep SAME")

  ## Legends
  if args.ratio: upper_pad.cd()
  leg = TLegend(0.77,0.7,0.89,0.89)
  leg.SetBorderSize(0)
  for s in samples:
    leg.AddEntry(histT[s],s,"f")
  if (args.qcd and name[-1] == "M"): leg.AddEntry(histT["QCD"],"QCD","f")
  if args.stack: leg.AddEntry(data, "Data" ,"lep")
  leg.Draw()

  ## Plot settings

  ## print image to a gif file
  if args.jets: 
    term= "total0Jnocor_v1v2v3v4"
  else:
    term= "totalHT_v1v2v3v4"

  if args.ratio: 
    notation = "_ratio_"
    if args.linear:
      notation = "_linratio_"
  else: 
    notation = "_normed_"

  if args.ssos: ssos_add="ssos_"
  else: ssos_add=""

  if args.png: c1.Print(plotdir+ssos_add+term+notation+ name + ".png")
  else: c1.Print(plotdir+ssos_add+term+notation + name + ".pdf")

for s in samplesHT:
        for data_op in datayears:
                        histFile[s][data_op].Close()

for data_op in datayears:
	histFileD[data_op]["M"].Close()
	histFileD[data_op]["E"].Close()

