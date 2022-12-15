###################################                                        
###################################                                        
#######     QCD f_BD est    #######
###################################                                        
###################################                                        

## Meant to be run locally, not sent to condor 

## Script to calculate the constant f_BD necessary to compute a QCD background estimation

import ROOT, os, sys
from ROOT import *
from os import listdir
from os.path import isfile, join

import json
import argparse

# Some defaults
gROOT.SetStyle("Plain")
gStyle.SetOptStat(1111111)
gStyle.SetPadGridX(True)
gStyle.SetPadGridY(True)
gStyle.SetGridStyle(3)
gStyle.SetCanvasDefW(1600)
gStyle.SetCanvasDefH(800)

ROOT.EnableImplicitMT()

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--year", type=string, default="2018",
                    help="Select data year to run")
parser.add_argument("--notfull", action="store_true", default=False,
                    help="Use just a range of the sample")
parser.add_argument('-l','--list', nargs='+', help='range of sample to use')
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Perform ssos substraction")
# Use like:
# python arg.py -l 1234 2345 3456 4567

args = parser.parse_args()

if (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): years = str(args.year)
else: raise NameError('Incorrect year')

samples_data = []

proc = ["M","E"]

for p in proc:
    samples_data.append(years+p)

print("the samples treated are",samples_data)

#proc_mc = ["WW","Wjets1","DY1","Wjets2","DY2","Wjets3","DY3","Wjets4","DY4","Wjets5","DY5","Wjets6","DY6","Wjets7","DY7","Wjets8","DY8",
#        "ttbar","ZZ","ttbarlep","WZ","ttbarhad","ST1","ST2","ST3","ST4"]
proc_mc = ["WW","Wjets0J","Wjets1J","Wjets2J","DY0J","DY1J","DY2J",
        "ttbar","ZZ","ttbarlep","WZ","ttbarhad","ST1","ST2","ST3","ST4"]
samples_mc = []

for p in proc_mc:
    samples_mc.append(p+years)

print("the mc samples treated are",samples_mc)

if args.notfull:
	if len(args.list) != 2: raise NameError('List has to have 2 elements')
	print("the range of files is",args.list)
 


# Create a ROOT dataframe for each dataset
# Note that we load the filenames from the external json file placed in the same folder than this script.
# Example "python analisisWW/selection_v2.py --process="WW" --notfull -l 0 50"

files_data = json.load(open("/nfs/cms/vazqueze/analisisWW/data_info_v9.json"))
processes = files_data.keys()

df_data = {}
archives_data = {}
lumi_data = {}

for s in samples_data:
  archives_data[s]=[]

for p in processes:
    # Construct the dataframes
    folder = files_data[p]["folder_name"] # Folder name
    filename = files_data[p]["sample"] # Sample name
    num_events = files_data[p]["events"] # Number of events
    num_files = files_data[p]["files"] # Number of files
    luminosity = files_data[p]["lumi"] # Luminosity
    list_files = [f for f in listdir(folder) if isfile(join(folder, f))] # Lista de archivos
    #print(len(list_files))
    for s in samples_data:
      if files_data[p]["type"]==s:
          lumi_data[s] = luminosity
          if (num_files == len(list_files)):
              for f in list_files:
                  archives_data[s].append(join(folder,f))

for s in samples_data:
  if args.notfull: archives_data[s]=archives_data[s][int(args.list[0]):int(args.list[1])]
  df_data[s] = ROOT.RDataFrame("Events",set(archives_data[s]))
  print("Number of files for",s,len(archives_data[s]))

files_mc = json.load(open("/nfs/cms/vazqueze/analisisWW/mc_info_v9"+years+".json"))
processes_mc = files_mc.keys()

special_weights = ["Wjets0J2016","Wjets1J2016","Wjets2J2016","DY0J2016","DY1J2016","DY2J2016",
        "Wjets0J2016B","Wjets1J2016B","Wjets2J2016B","DY0J2016B","DY1J2016B","DY2J2016B",
        "Wjets0J2017","Wjets1J2017","Wjets2J2017","DY0J2017","DY1J2017","DY2J2017",
        "Wjets0J2018","Wjets1J2018","Wjets2J2018","DY0J2018","DY1J2018","DY2J2018"]

df_mc = {}
archives_mc = {}
lumi_mc = {}

for s in samples_mc:
  archives_mc[s]=[]

for p in processes_mc:
    # Construct the dataframes
    folder = files_mc[p]["folder_name"] # Folder name
    filename = files_mc[p]["sample"] # Sample name
    num_events = files_mc[p]["events"] # Number of events
    num_files = files_mc[p]["files"] # Number of files
    luminosity = files_mc[p]["lumi"] # Luminosity
    list_files = [f for f in listdir(folder) if isfile(join(folder, f))] # Lista de archivos
    #print(len(list_files))
    for s in samples_mc:
      if files_mc[p]["type"]==s:
          lumi_mc[s] = luminosity
          if (num_files == len(list_files)):
              for f in list_files:
                  archives_mc[s].append(join(folder,f))

for s in samples_mc:
  if args.notfull: archives_mc[s]=archives_mc[s][int(args.list[0]):int(args.list[1])]
  df_mc[s] = ROOT.RDataFrame("Events",set(archives_mc[s]))
  print("Number of files for",s,len(archives_mc[s]))

## Correcting lumis for Wjets and DY per year

if years == "2016":
  lumi_mc["DY0J"] = 12.5
  lumi_mc["DY1J"] = 43.57
  lumi_mc["DY2J"] = 36.65
  lumi_mc["Wjets0J"] = 2.46
  lumi_mc["Wjets1J"] = 10.17
  lumi_mc["Wjets2J"] = 8.81

if years == "2016B":
  lumi_mc["DY0J"] = 11.58
  lumi_mc["DY1J"] = 40.61
  lumi_mc["DY2J"] = 35.70
  lumi_mc["Wjets0J"] = 2.60
  lumi_mc["Wjets1J"] = 9.81
  lumi_mc["Wjets2J"] = 8.60

if years == "2017":
  lumi_mc["DY0J"] = 12.27
  lumi_mc["DY1J"] = 41.64
  lumi_mc["DY2J"] = 39.62
  lumi_mc["Wjets0J"] = 2.48
  lumi_mc["Wjets1J"] = 9.75
  lumi_mc["Wjets2J"] = 9.06

if years == "2018":
  lumi_mc["DY0J"] = 13.55
  lumi_mc["DY1J"] = 45.56
  lumi_mc["DY2J"] = 37.95
  lumi_mc["Wjets0J"] = 2.53
  lumi_mc["Wjets1J"] = 9.78
  lumi_mc["Wjets2J"] = 8.78

## Cuts per year

cuts_btag = {}

cuts_btag["2016"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2016B"]=[0.0480, 0.2489, 0.6377]
cuts_btag["2017"]=[0.0532, 0.3040, 0.7476]
cuts_btag["2018"]=[0.0490,0.2783,0.7100]

## Muon and Electron pT cuts per year

muon_pt = {}

muon_pt["2016"]=26
muon_pt["2016B"]=26
muon_pt["2017"]=29
muon_pt["2018"]=26

el_pt = {}

el_pt["2016"]=29
el_pt["2016B"]=29
el_pt["2017"]=34
el_pt["2018"]=34

## Triggers per year

muon_trig = {}

muon_trig["2016"]="HLT_IsoMu24 || HLT_IsoTkMu24"
muon_trig["2016B"]="HLT_IsoMu24 || HLT_IsoTkMu24"
muon_trig["2017"]="HLT_IsoMu27"
muon_trig["2018"]="HLT_IsoMu24"

el_trig = {}

el_trig["2016"]="HLT_Ele27_WPTight_Gsf"
el_trig["2016B"]="HLT_Ele27_WPTight_Gsf"
el_trig["2017"]="HLT_Ele32_WPTight_Gsf"
el_trig["2018"]="HLT_Ele32_WPTight_Gsf"

### Coefficients for fixing W+jets

CoefWJ_M = {}
CoefWJ_E = {}

CoefWJ_M["2016"] = 0.8159 
CoefWJ_M["2016B"] = 0.8718
CoefWJ_M["2017"] = 0.8776
CoefWJ_M["2018"] = 0.8585

CoefWJ_E["2016"] = 0.7649
CoefWJ_E["2016B"] = 0.7873
CoefWJ_E["2017"] = 0.7535
CoefWJ_E["2018"] = 0.7491

#######################################
########      GEN_PART      ###########
#######################################

# Attempt to classify with GenPart
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;
 
// Function to print the
// index of an element
// 1: Full hadronic 2: One W hadronic, the other leptonic (SEMI) 3: Full leptonic
      auto getIndex(UInt_t nGenPart,Vint v, int K){
            int ind = 2;
            while (v[ind]!=K) {
                ++ind;
            }
            return ind;
      }
      auto typeWW(UInt_t nGenPart, Vint partId, Vint motherId) {
	    int type = -1;
	    auto c = (partId == 24)||(partId == -24);
	    std::vector<int> v1(size(partId));
	    std::iota(std::begin(v1),std::end(v1),0);
	    ROOT::VecOps::RVec<int> myRVec(v1.data(), v1.size());
	    int v2 = -1;
	    auto Windexes = ROOT::VecOps::Where(c,myRVec,v2);
	    int idWlast1 = ROOT::VecOps::Max(Windexes);
	    int idWlast2 = idWlast1-1;
	    int idpart1 = getIndex(nGenPart,motherId, idWlast1);
            int idpart2 = getIndex(nGenPart,motherId, idWlast2);
	    
	    if (motherId[idpart1]==motherId[idpart1+1] && motherId[idpart2]==motherId[idpart2+1]){
		if (fabs(partId[idpart1])< 9 && fabs(partId[idpart2])< 9) {
  			type = 1;
		} else if (fabs(partId[idpart1])< 9 || fabs(partId[idpart2])< 9) {
  			type = 2;
		} else {
  			type = 3;
		}
	    }
    	    return type;
      };
""")


# Attempt to classify with GenPart
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;

// Function to print the
// index of an element
// 1: Charm -1: No Charm
      auto typeC(UInt_t nGenPart,Vint partId, Vint motherId){
	    int type = -1;
	    auto c = (partId == 24)||(partId == -24);
            std::vector<int> v1(size(partId));
            std::iota(std::begin(v1),std::end(v1),0);
            ROOT::VecOps::RVec<int> myRVec(v1.data(), v1.size());
            int v2 = -1;
            auto Windexes = ROOT::VecOps::Where(c,myRVec,v2);
            int idWlast1 = ROOT::VecOps::Max(Windexes);
            int idWlast2 = idWlast1-1;
            int idpart1 = getIndex(nGenPart,motherId, idWlast1);
            int idpart2 = getIndex(nGenPart,motherId, idWlast2);

            if (motherId[idpart1]==motherId[idpart1+1] && motherId[idpart2]==motherId[idpart2+1]){
                if ((fabs(partId[idpart1])==4 || fabs(partId[idpart2])==4) || (fabs(partId[idpart1+1])==4 || fabs(partId[idpart2+1])==4)) {
                        type = 1;
                }
            }
            return type;
      };
""")

#######################    main selection    #########################

###################################################
################   DEFINITIONS   ##################
###################################################

### LEPTONS

## Funciones para seleccionar muones y electrones
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vbool tID, Float_t cutpt) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
		if (pt[i]>cutpt && fabs(eta[i])<2.4 && tID[i]){
                	vb.push_back(i);
		}
            }
            return vb;
      };
      auto elInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vint cutB, Vbool mva80 , Vbool mvaL, Float_t cutpt) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
                if (pt[i]>cutpt && fabs(eta[i])<2.5 && mva80[i]){
                	vb.push_back(i);
                }
            }
            return vb;
      };
""")

#######   JETS   #######

## Funciones para seleccionar JETS
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetInd(UInt_t njet, Vfloat pt, Vfloat eta, Vfloat phi, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vint el_good, Vfloat el_eta, Vfloat el_phi, Vint puID, Vint jetID) {
            vector<int> vb;
	    bool cond = false;
	    bool cond1 = false;
            for (unsigned int i=0; i<njet; ++i) {
		if(mu_good.size()>0){
			pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
			cond = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],eta[i],mu_phi[mu_good[0]],phi[i]) > 0.4;
                        if (pt[i]>30. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                                vb.push_back(i);
                        }
		}
                if(el_good.size()>0){
			pt[i]<50. ? cond1 = puID[i]>=4 : cond1 = true;
                        cond = ROOT::VecOps::DeltaR(el_eta[el_good[0]],eta[i],el_phi[el_good[0]],phi[i]) > 0.4;
                        if (pt[i]>30. && fabs(eta[i])<2.4 && cond && cond1 && jetID[i]>1){
                                vb.push_back(i);
                        }
                }
            }
            return vb;
      };

      auto InvariantM(const float pt, const float eta, const float phi, const float mass, const float pt1, const float eta1, const float phi1, const float mass1) {
            auto x = pt*std::cos(phi);
	    auto x1 = pt1*std::cos(phi1);
            auto y = pt*std::sin(phi);
            auto y1 = pt1*std::sin(phi1);
            auto z = pt*std::sinh(eta);
	    auto z1 = pt1*std::sinh(eta1);
            auto e = std::sqrt(x*x+y*y+z*z+mass*mass);
	    auto e1 = std::sqrt(x1*x1+y1*y1+z1*z1+mass1*mass1);

            auto mJet = std::sqrt((e+e1)*(e+e1)-(x+x1)*(x+x1)-(y+y1)*(y+y1)-(z+z1)*(z+z1));
            return mJet;
      };
""")

## Muon dentro del jet

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetMuonIndJet(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_pt, Vfloat mu_iso, Vint el_good, Vint mu_jetid) {
            vector<int> vb;
            bool cond = false;
            bool cond1 = true;
            int ind = -1;
            float ptM{-10.};
            float ptJ{-10.};
            if (good.size() == 1){
                for (unsigned int i=0; i<nmu; ++i){
                        cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[0]],mu_phi[i],phi[good[0]]) < 0.4;
                        if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                        if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0] && mu_iso[i] > 0.2){
                                ind = good[0];
                                ptM = mu_pt[i];
                        }
                }
            }
            if (good.size() > 1){
                for (unsigned int j=0; j<good.size(); ++j){
                        for (unsigned int i=0; i<nmu; ++i){
                                cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[j]],mu_phi[i],phi[good[j]]) < 0.4;
                                if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                                if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && mu_iso[i] > 0.2){
                                        ind = good[j];
                                        ptM = mu_pt[i];
                                        ptJ = pt[good[j]];
                                }
                        }
                }
            }
            if (ind>-1) {
                vb.push_back(ind);
            }
            return vb;
      };
      auto JetMuonInd(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_pt, Vfloat mu_iso, Vint el_good, Vint mu_jetid) {
            vector<int> vb;
            bool cond = false;
            bool cond1 = true;
            int ind = -1;
            float ptM{-10.};
            float ptJ{-10.};
            if (good.size() == 1){
                for (unsigned int i=0; i<nmu; ++i){
                        cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[0]],mu_phi[i],phi[good[0]]) < 0.4;
                        if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                        if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[0] && mu_iso[i] > 0.2){
                                ind = i;
                                ptM = mu_pt[i];
                        }
                }
            }
            if (good.size() > 1){
                for (unsigned int j=0; j<good.size(); ++j){
                        for (unsigned int i=0; i<nmu; ++i){
                                cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[j]],mu_phi[i],phi[good[j]]) < 0.4;
                                if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                                if(cond && mu_pt[i] > ptM && cond1 && mu_jetid[i]==good[j] && mu_iso[i] > 0.2){
                                        ind = i;
                                        ptM = mu_pt[i];
                                        ptJ = pt[good[j]];
                                }
                        }
                }
            }
            if (ind>-1) {
                vb.push_back(ind);
            }
            return vb;
      };
""")

######### pedimos algunas condiciones al muon seleccionado

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonCond( Vint mu_jet,Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_iso, Vbool mu_softid) {
            bool cond = false;
            if (mu_jet.size()>0){
                if (mu_pt[mu_jet[0]]<25. && fabs(mu_eta[mu_jet[0]])<2.4 && mu_softid[mu_jet[0]]) {
                        cond = true;
                }
            }
            return cond;
      };
""")

####   SSOS   ####

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto SSOS(Vint mu_good, Vint el_good, Vint mu_jet,Vint el_charge, Vint mu_charge) {
	    int ind = 0;
	    if(mu_jet.size()>0){
            	if(mu_good.size()>0){
			ind = mu_charge[mu_good[0]]*mu_charge[mu_jet[0]];
            	}
            	if(el_good.size()>0){
                	ind = el_charge[el_good[0]]*mu_charge[mu_jet[0]];
            	}
	    }
            return ind;
      };
""")

## Trigger function for 2017 electron data

gInterpreter.Declare("""
      #include <iomanip>
      #include <math.h>
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto calDeltaR(float eta1, float phi1, float eta2, float phi2) {
          float dphi = phi1-phi2;
          float deta = eta1-eta2;
          if (dphi <= - M_PI) {
             dphi = dphi + M_PI;
          }
          if (dphi >= M_PI) {
             dphi = dphi - M_PI;
          }
          float prod = deta*deta + dphi*dphi;
          return prod;
      };
      auto triggeremulator(UInt_t nel, UInt_t ntrig, Vfloat el_eta, Vfloat el_deltaeta, Vfloat el_phi, Vfloat trig_eta, Vfloat trig_phi, Vint trig_bits, Vint trig_id) {
          bool cond = false;
          bool cond1 = false;
          bool cond2 = false;
          for (unsigned int j=0; j<nel; ++j) {
            for (unsigned int i=0; i<ntrig; ++i) {
                cond1 = calDeltaR(el_eta[j]+el_deltaeta[j],el_phi[j], trig_eta[i], trig_phi[i]) < 0.01;
                cond2 = trig_bits[i] & (0x1 << 10);
                if( trig_id[i]==11 && cond1 && cond2) {
                    cond = true;
                }
            }
          }
          return cond;
      };
""")

#############################################################
######################## Definitions ########################
#############################################################

for s in samples_mc:
      if s in special_weights:
             df_mc[s] = df_mc[s].Define('weight_aux','fabs(genWeight) > 0 ? genWeight/fabs(genWeight) : 0')
      else:
             df_mc[s] = df_mc[s].Define('weight_aux','1')


df_data_muon = {}
df_data_electron = {}

for s in samples_data:
	df_data[s] = df_data[s].Define('MuonGoodInd','muonInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all, Muon_tightId,'+str(muon_pt[years])+')')
	df_data[s] = df_data[s].Define('ElectronGoodInd','elInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_cutBased, Electron_mvaFall17V2noIso_WP80, Electron_mvaFall17V2noIso_WP90,'+str(el_pt[years])+')')
	df_data[s] = df_data[s].Define('nMuonGood','MuonGoodInd.size()')
	df_data[s] = df_data[s].Define('nElectronGood','ElectronGoodInd.size()')
	df_data[s] = df_data[s].Define('JetGoodInd','JetInd(nJet, Jet_pt, Jet_eta, Jet_phi, MuonGoodInd, Muon_eta, Muon_phi, ElectronGoodInd, Electron_eta, Electron_phi, Jet_puId, Jet_jetId)')
	df_data[s] = df_data[s].Define('nJetGood','JetGoodInd.size()')
	df_data[s] = df_data[s].Define('MuonJetInd','JetMuonInd(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx)')
	df_data[s] = df_data[s].Define('JetMuonInd','JetMuonIndJet(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx)')
	df_data[s] = df_data[s].Define('MuonJetGood','muonCond( MuonJetInd, Muon_pt, Muon_eta, Muon_pfRelIso04_all, Muon_softId)')
	df_data[s] = df_data[s].Define('MuonLepSign','SSOS(MuonGoodInd, ElectronGoodInd, MuonJetInd, Electron_charge, Muon_charge)')
	if args.ssos: df_data[s] = df_data[s].Define('weightSSOS','-1*MuonLepSign')
	else: df_data[s] = df_data[s].Define('weightSSOS','1')

df_mc_muon = {}
df_mc_electron = {}

for s in samples_mc:
        df_mc[s] = df_mc[s].Define('MuonGoodInd','muonInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all, Muon_tightId,'+str(muon_pt[years])+')')
        df_mc[s] = df_mc[s].Define('ElectronGoodInd','elInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_cutBased, Electron_mvaFall17V2noIso_WP80, Electron_mvaFall17V2noIso_WP90,'+str(el_pt[years])+')')
        df_mc[s] = df_mc[s].Define('nMuonGood','MuonGoodInd.size()')
        df_mc[s] = df_mc[s].Define('nElectronGood','ElectronGoodInd.size()')
        df_mc[s] = df_mc[s].Define('JetGoodInd','JetInd(nJet, Jet_pt, Jet_eta, Jet_phi, MuonGoodInd, Muon_eta, Muon_phi, ElectronGoodInd, Electron_eta, Electron_phi, Jet_puId, Jet_jetId)')
        df_mc[s] = df_mc[s].Define('nJetGood','JetGoodInd.size()')
        df_mc[s] = df_mc[s].Define('MuonJetInd','JetMuonInd(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx)')
        df_mc[s] = df_mc[s].Define('JetMuonInd','JetMuonIndJet(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx)')
        df_mc[s] = df_mc[s].Define('MuonJetGood','muonCond( MuonJetInd, Muon_pt, Muon_eta, Muon_pfRelIso04_all, Muon_softId)')
        df_mc[s] = df_mc[s].Define('MuonLepSign','SSOS(MuonGoodInd, ElectronGoodInd, MuonJetInd, Electron_charge, Muon_charge)')
        if args.ssos: df_mc[s] = df_mc[s].Define('weightSSOS','-1*MuonLepSign*weight_aux')
        else: df_mc[s] = df_mc[s].Define('weightSSOS','weight_aux')

#################     FILTERS     #######################

##### New cuts, compared to verison 0 and 1
##### Exactly 2 jets, not more or less
##### ETA restrictions: 2.4 for jets and muons and 2.5 for electrons

###########    Extra definitions

for s in samples_data:
	df_data[s] = df_data[s].Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood==2').Filter('MuonJetInd.size() >= 1').Filter('MuonJetGood')
	df_data[s] = df_data[s].Define('JetnotMuonInd','nJetGood>1 ? (JetMuonInd[0] == JetGoodInd[0] ? JetGoodInd[1] : JetGoodInd[0]) : -1')
	### hists definitions
	df_data[s] = df_data[s].Define('jet_muon_pt','Jet_pt[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('jet_notmuon_pt','nJetGood>1 ? Jet_pt[JetnotMuonInd] : 0')
	df_data[s] = df_data[s].Define('jet_muon_eta','Jet_eta[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('jet_notmuon_eta','nJetGood>1 ? Jet_eta[JetnotMuonInd] : 0')
	df_data[s] = df_data[s].Define('muon_jet_pt','Muon_pt[MuonJetInd[0]]')
	df_data[s] = df_data[s].Define('muon_jet_eta','Muon_eta[MuonJetInd[0]]')
	df_data[s] = df_data[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Jet_pt[JetnotMuonInd],Jet_eta[JetnotMuonInd],Jet_phi[JetnotMuonInd],Jet_mass[JetnotMuonInd]) : 0')
	df_data[s] = df_data[s].Define('deltaR_jetM_jetNM','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotMuonInd], Jet_eta[JetMuonInd[0]] , Jet_phi[JetnotMuonInd], Jet_phi[JetMuonInd[0]])  : 10')
	df_data[s] = df_data[s].Define('tracks_jetM','Jet_nConstituents[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('tracks_jetNM','nJetGood>1 ? Jet_nConstituents[JetnotMuonInd] : 0')
	df_data[s] = df_data[s].Define('EMN_jetM','Jet_neEmEF[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('EMC_jetM','Jet_chEmEF[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('EMtotal_jetM','Jet_chEmEF[JetMuonInd[0]]+Jet_neEmEF[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('muon_jet_mva','Muon_softMva[MuonJetInd[0]]')
	df_data[s] = df_data[s].Define('muon_jet_relpt','Muon_pt[MuonJetInd[0]]/Jet_pt[JetMuonInd[0]]')                                                                                        
	df_data[s] = df_data[s].Define('jet_muon_btag','Jet_btagDeepFlavB[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('jet_notmuon_btag','Jet_btagDeepFlavB[JetnotMuonInd]')

for s in samples_mc:
        df_mc[s] = df_mc[s].Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood==2').Filter('MuonJetInd.size() >= 1').Filter('MuonJetGood')
        df_mc[s] = df_mc[s].Define('JetnotMuonInd','nJetGood>1 ? (JetMuonInd[0] == JetGoodInd[0] ? JetGoodInd[1] : JetGoodInd[0]) : -1')
        ### hists definitions
        df_mc[s] = df_mc[s].Define('jet_muon_pt','Jet_pt[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('jet_notmuon_pt','nJetGood>1 ? Jet_pt[JetnotMuonInd] : 0')
        df_mc[s] = df_mc[s].Define('jet_muon_eta','Jet_eta[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('jet_notmuon_eta','nJetGood>1 ? Jet_eta[JetnotMuonInd] : 0')
        df_mc[s] = df_mc[s].Define('muon_jet_pt','Muon_pt[MuonJetInd[0]]')
        df_mc[s] = df_mc[s].Define('muon_jet_eta','Muon_eta[MuonJetInd[0]]')
        df_mc[s] = df_mc[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Jet_pt[JetnotMuonInd],Jet_eta[JetnotMuonInd],Jet_phi[JetnotMuonInd],Jet_mass[JetnotMuonInd]) : 0')
        df_mc[s] = df_mc[s].Define('deltaR_jetM_jetNM','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotMuonInd], Jet_eta[JetMuonInd[0]] , Jet_phi[JetnotMuonInd], Jet_phi[JetMuonInd[0]])  : 10')
        df_mc[s] = df_mc[s].Define('tracks_jetM','Jet_nConstituents[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('tracks_jetNM','nJetGood>1 ? Jet_nConstituents[JetnotMuonInd] : 0')
        df_mc[s] = df_mc[s].Define('EMN_jetM','Jet_neEmEF[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('EMC_jetM','Jet_chEmEF[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('EMtotal_jetM','Jet_chEmEF[JetMuonInd[0]]+Jet_neEmEF[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('muon_jet_mva','Muon_softMva[MuonJetInd[0]]')
        df_mc[s] = df_mc[s].Define('muon_jet_relpt','Muon_pt[MuonJetInd[0]]/Jet_pt[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('jet_muon_btag','Jet_btagDeepFlavB[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('jet_notmuon_btag','Jet_btagDeepFlavB[JetnotMuonInd]')

## Differentiated definitions
for s in samples_data:
	df_data_muon[s] = df_data[s].Filter('nMuonGood>0')
	df_data_electron[s] = df_data[s].Filter('nElectronGood >0')
	if s[-1] == "M":
		df_data_muon[s] = df_data_muon[s].Filter(muon_trig[years])
		df_data_electron[s] = df_data_electron[s].Filter(muon_trig[years])
	if s[-1] == "E":
		if years == "2017":
			df_data_muon[s] = df_data_muon[s].Filter('triggeremulator(nElectron, nTrigObj, Electron_eta, Electron_deltaEtaSC, Electron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)')
			df_data_electron[s] = df_data_electron[s].Filter('triggeremulator(nElectron, nTrigObj, Electron_eta, Electron_deltaEtaSC, Electron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)')
		else:
			df_data_muon[s] = df_data_muon[s].Filter(el_trig[years])
			df_data_electron[s] = df_data_electron[s].Filter(el_trig[years])
	df_data_muon[s] = df_data_muon[s].Define('lepton_pt','Muon_pt[MuonGoodInd[0]]')
	df_data_muon[s] = df_data_muon[s].Define('lepton_eta','Muon_eta[MuonGoodInd[0]]')
	df_data_muon[s] = df_data_muon[s].Define('lepton_iso','Muon_pfRelIso04_all[MuonGoodInd[0]]')
	df_data_electron[s] = df_data_electron[s].Define('lepton_pt','Electron_pt[ElectronGoodInd[0]]')
	df_data_electron[s] = df_data_electron[s].Define('lepton_eta','Electron_eta[ElectronGoodInd[0]]')
	df_data_electron[s] = df_data_electron[s].Define('lepton_iso','Electron_pfRelIso03_all[ElectronGoodInd[0]]')
	df_data_muon[s] = df_data_muon[s].Define('deltaR_jetM_lepM','ROOT::VecOps::DeltaR(Muon_eta[MuonGoodInd[0]],Jet_eta[JetMuonInd[0]] , Muon_phi[MuonGoodInd[0]], Jet_phi[JetMuonInd[0]])')
	df_data_electron[s] = df_data_electron[s].Define('deltaR_jetM_lepE','ROOT::VecOps::DeltaR(Electron_eta[ElectronGoodInd[0]],Jet_eta[JetMuonInd[0]] , Electron_phi[ElectronGoodInd[0]], Jet_phi[JetMuonInd[0]])')
	df_data_muon[s] = df_data_muon[s].Define('InvM_jetM_lep', 'InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]])')
	df_data_electron[s] = df_data_electron[s].Define('InvM_jetM_lep', 'InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Electron_pt[ElectronGoodInd[0]],Electron_eta[ElectronGoodInd[0]],Electron_phi[ElectronGoodInd[0]], Electron_mass[ElectronGoodInd[0]])')
	df_data_muon[s] = df_data_muon[s].Define('transverse_mass','std::sqrt(2*Muon_pt[MuonGoodInd[0]]*MET_pt*(1-std::cos(Muon_phi[MuonGoodInd[0]]-MET_phi)))')
	df_data_electron[s] = df_data_electron[s].Define('transverse_mass','std::sqrt(2*Electron_pt[ElectronGoodInd[0]]*MET_pt*(1-std::cos(Electron_phi[ElectronGoodInd[0]]-MET_phi)))')
	df_data_muon[s] = df_data_muon[s].Define('InvM_muon_jet','InvariantM(Muon_pt[MuonJetInd[0]],Muon_eta[MuonJetInd[0]],Muon_phi[MuonJetInd[0]],Muon_mass[MuonJetInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]])')
	## New cuts
	df_data_muon[s] = df_data_muon[s].Filter('InvM_muon_jet >12').Filter('InvM_muon_jet > 110 || InvM_muon_jet < 70')
	#df_data_electron[s] = df_data_electron[s].Filter('transverse_mass > 50')
	#df_data_muon[s] = df_data_muon[s].Filter('transverse_mass > 50')
	df_data_electron[s] = df_data_electron[s].Filter('muon_jet_relpt<0.5')
	df_data_muon[s] = df_data_muon[s].Filter('muon_jet_relpt<0.5')
	df_data_electron[s] = df_data_electron[s].Filter('EMtotal_jetM<0.4')
	df_data_muon[s] = df_data_muon[s].Filter('EMtotal_jetM<0.4')
	#df_data_electron[s] = df_data_electron[s].Filter('InvM_2jets<100').Filter('InvM_2jets>60')
	#df_data_muon[s] = df_data_muon[s].Filter('InvM_2jets<100').Filter('InvM_2jets>60')
	df_data_electron[s] = df_data_electron[s].Filter('jet_notmuon_btag<'+str(cuts_btag[years][1]))
	df_data_muon[s] = df_data_muon[s].Filter('jet_notmuon_btag<'+str(cuts_btag[years][1]))
	df_data_electron[s] = df_data_electron[s].Filter('jet_muon_btag<'+str(cuts_btag[years][2]))
	df_data_muon[s] = df_data_muon[s].Filter('jet_muon_btag<'+str(cuts_btag[years][2]))

## Differentiated definitions
for s in samples_mc:
        df_mc_muon[s] = df_mc[s].Filter('nMuonGood>0')
        df_mc_electron[s] = df_mc[s].Filter('nElectronGood >0')
        df_mc_muon[s] = df_mc_muon[s].Filter(muon_trig[years])
        df_mc_electron[s] = df_mc_electron[s].Filter(el_trig[years])
        df_mc_muon[s] = df_mc_muon[s].Define('lepton_pt','Muon_pt[MuonGoodInd[0]]')
        df_mc_muon[s] = df_mc_muon[s].Define('lepton_eta','Muon_eta[MuonGoodInd[0]]')
        df_mc_muon[s] = df_mc_muon[s].Define('lepton_iso','Muon_pfRelIso04_all[MuonGoodInd[0]]')
        df_mc_electron[s] = df_mc_electron[s].Define('lepton_pt','Electron_pt[ElectronGoodInd[0]]')
        df_mc_electron[s] = df_mc_electron[s].Define('lepton_eta','Electron_eta[ElectronGoodInd[0]]')
        df_mc_electron[s] = df_mc_electron[s].Define('lepton_iso','Electron_pfRelIso03_all[ElectronGoodInd[0]]')
        df_mc_muon[s] = df_mc_muon[s].Define('deltaR_jetM_lepM','ROOT::VecOps::DeltaR(Muon_eta[MuonGoodInd[0]],Jet_eta[JetMuonInd[0]] , Muon_phi[MuonGoodInd[0]], Jet_phi[JetMuonInd[0]])')
        df_mc_electron[s] = df_mc_electron[s].Define('deltaR_jetM_lepE','ROOT::VecOps::DeltaR(Electron_eta[ElectronGoodInd[0]],Jet_eta[JetMuonInd[0]] , Electron_phi[ElectronGoodInd[0]], Jet_phi[JetMuonInd[0]])')
        df_mc_muon[s] = df_mc_muon[s].Define('InvM_jetM_lep', 'InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]])')
        df_mc_electron[s] = df_mc_electron[s].Define('InvM_jetM_lep', 'InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Electron_pt[ElectronGoodInd[0]],Electron_eta[ElectronGoodInd[0]],Electron_phi[ElectronGoodInd[0]], Electron_mass[ElectronGoodInd[0]])')
        df_mc_muon[s] = df_mc_muon[s].Define('transverse_mass','std::sqrt(2*Muon_pt[MuonGoodInd[0]]*MET_pt*(1-std::cos(Muon_phi[MuonGoodInd[0]]-MET_phi)))')
        df_mc_electron[s] = df_mc_electron[s].Define('transverse_mass','std::sqrt(2*Electron_pt[ElectronGoodInd[0]]*MET_pt*(1-std::cos(Electron_phi[ElectronGoodInd[0]]-MET_phi)))')
        df_mc_muon[s] = df_mc_muon[s].Define('InvM_muon_jet','InvariantM(Muon_pt[MuonJetInd[0]],Muon_eta[MuonJetInd[0]],Muon_phi[MuonJetInd[0]],Muon_mass[MuonJetInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]])')
        ## New cuts
        df_mc_muon[s] = df_mc_muon[s].Filter('InvM_muon_jet >12').Filter('InvM_muon_jet > 110 || InvM_muon_jet < 70')
        #df_mc_electron[s] = df_mc_electron[s].Filter('transverse_mass > 50')
        #df_mc_muon[s] = df_mc_muon[s].Filter('transverse_mass > 50')
        df_mc_electron[s] = df_mc_electron[s].Filter('muon_jet_relpt<0.5')
        df_mc_muon[s] = df_mc_muon[s].Filter('muon_jet_relpt<0.5')
        df_mc_electron[s] = df_mc_electron[s].Filter('EMtotal_jetM<0.4')
        df_mc_muon[s] = df_mc_muon[s].Filter('EMtotal_jetM<0.4')
        #df_mc_electron[s] = df_mc_electron[s].Filter('InvM_2jets<100').Filter('InvM_2jets>60')
        #df_mc_muon[s] = df_mc_muon[s].Filter('InvM_2jets<100').Filter('InvM_2jets>60')
        df_mc_electron[s] = df_mc_electron[s].Filter('jet_notmuon_btag<'+str(cuts_btag[years][1]))
        df_mc_muon[s] = df_mc_muon[s].Filter('jet_notmuon_btag<'+str(cuts_btag[years][1]))
        df_mc_electron[s] = df_mc_electron[s].Filter('jet_muon_btag<'+str(cuts_btag[years][2]))
        df_mc_muon[s] = df_mc_muon[s].Filter('jet_muon_btag<'+str(cuts_btag[years][2]))

######################################
######### B and D definition #########
######################################

df_dataB_muon = {}
df_dataB_electron = {}
df_dataD_muon = {}
df_dataD_electron = {}

for s in samples_data:
	df_dataB_muon[s] = df_data_muon[s].Filter('lepton_iso<0.15').Filter('transverse_mass<50')
	df_dataD_muon[s] = df_data_muon[s].Filter('lepton_iso>0.2').Filter('transverse_mass<50')
	df_dataB_electron[s] = df_data_electron[s].Filter('lepton_iso<0.15').Filter('transverse_mass<50')
	df_dataD_electron[s] = df_data_electron[s].Filter('lepton_iso>0.2').Filter('transverse_mass<50')

df_mcB_muon = {}
df_mcB_electron = {}
df_mcD_muon = {}
df_mcD_electron = {}

for s in samples_mc:
        df_mcB_muon[s] = df_mc_muon[s].Filter('lepton_iso<0.15').Filter('transverse_mass<50')
        df_mcD_muon[s] = df_mc_muon[s].Filter('lepton_iso>0.2').Filter('transverse_mass<50')
        df_mcB_electron[s] = df_mc_electron[s].Filter('lepton_iso<0.15').Filter('transverse_mass<50')
        df_mcD_electron[s] = df_mc_electron[s].Filter('lepton_iso>0.2').Filter('transverse_mass<50')

############################################################
####################     HISTS    ##########################
############################################################

hist_dataB_nJetGood_M = {}
hist_dataD_nJetGood_M = {}
hist_dataB_nJetGood_E = {}
hist_dataD_nJetGood_E = {}

                
for s in samples_data:
	hist_dataB_nJetGood_M[s] = df_dataB_muon[s].Histo1D(("nJetGood M","",10,0,10),"nJetGood","weightSSOS")
	hist_dataD_nJetGood_M[s] = df_dataD_muon[s].Histo1D(("nJetGood M","",10,0,10),"nJetGood","weightSSOS")
	hist_dataB_nJetGood_E[s] = df_dataB_electron[s].Histo1D(("nJetGood E","",10,0,10),"nJetGood","weightSSOS")
	hist_dataD_nJetGood_E[s] = df_dataD_electron[s].Histo1D(("nJetGood E","",10,0,10),"nJetGood","weightSSOS")
	# For testing
	if s[-1] == "E":
		hist_dataB_lepton_iso_E = df_dataB_electron[s].Histo1D(("lepton iso E","",50,0,2),"lepton_iso","weightSSOS")
		hist_dataD_lepton_iso_E = df_dataD_electron[s].Histo1D(("lepton iso E","",50,0,2),"lepton_iso","weightSSOS")
		hist_dataB_transverse_mass_E = df_dataB_electron[s].Histo1D(("transverse mass E","",50,0,200),"transverse_mass","weightSSOS")
		hist_dataD_transverse_mass_E = df_dataD_electron[s].Histo1D(("transverse mass E","",50,0,200),"transverse_mass","weightSSOS")

hist_mcB_nJetGood_M = {}
hist_mcD_nJetGood_M = {}
hist_mcB_nJetGood_E = {}
hist_mcD_nJetGood_E = {}

hist_mcB_lepton_iso_E = {}
hist_mcD_lepton_iso_E = {}
hist_mcB_transverse_mass_E = {}
hist_mcD_transverse_mass_E = {}


for s in samples_mc:
        hist_mcB_nJetGood_M[s] = df_mcB_muon[s].Histo1D(("nJetGood M","",10,0,10),"nJetGood","weightSSOS")
        hist_mcD_nJetGood_M[s] = df_mcD_muon[s].Histo1D(("nJetGood M","",10,0,10),"nJetGood","weightSSOS")
        hist_mcB_nJetGood_E[s] = df_mcB_electron[s].Histo1D(("nJetGood E","",10,0,10),"nJetGood","weightSSOS")
        hist_mcD_nJetGood_E[s] = df_mcD_electron[s].Histo1D(("nJetGood E","",10,0,10),"nJetGood","weightSSOS")
        # For testing
        hist_mcB_lepton_iso_E[s] = df_mcB_electron[s].Histo1D(("lepton iso M","",50,0,2),"lepton_iso","weightSSOS")
        hist_mcD_lepton_iso_E[s] = df_mcD_electron[s].Histo1D(("lepton iso M","",50,0,2),"lepton_iso","weightSSOS")
        hist_mcB_transverse_mass_E[s] = df_mcB_electron[s].Histo1D(("transverse mass M","",50,0,200),"transverse_mass","weightSSOS")
        hist_mcD_transverse_mass_E[s] = df_mcD_electron[s].Histo1D(("transverse mass M","",50,0,200),"transverse_mass","weightSSOS")


####### Normalization

for s in samples_data:
  lumiD = lumi_data[s]

for s in samples_mc:
  hist_mcB_nJetGood_M[s].Scale(lumiD/lumi_mc[s])
  hist_mcD_nJetGood_M[s].Scale(lumiD/lumi_mc[s])
  hist_mcB_nJetGood_E[s].Scale(lumiD/lumi_mc[s])
  hist_mcD_nJetGood_E[s].Scale(lumiD/lumi_mc[s])
  hist_mcB_lepton_iso_E[s].Scale(lumiD/lumi_mc[s])
  hist_mcD_lepton_iso_E[s].Scale(lumiD/lumi_mc[s])
  hist_mcB_transverse_mass_E[s].Scale(lumiD/lumi_mc[s])
  hist_mcD_transverse_mass_E[s].Scale(lumiD/lumi_mc[s])

## Special W+jet normalization (control region)

for s in samples_mc:
  if s[0]+s[1]+s[2]+s[3]+s[4] == "Wjets":
    hist_mcB_nJetGood_M[s].Scale(CoefWJ_M[years])
    hist_mcD_nJetGood_M[s].Scale(CoefWJ_M[years])
    hist_mcB_nJetGood_E[s].Scale(CoefWJ_E[years])
    hist_mcD_nJetGood_E[s].Scale(CoefWJ_E[years])
    hist_mcB_lepton_iso_E[s].Scale(CoefWJ_E[years])
    hist_mcD_lepton_iso_E[s].Scale(CoefWJ_E[years])
    hist_mcB_transverse_mass_E[s].Scale(CoefWJ_E[years])
    hist_mcD_transverse_mass_E[s].Scale(CoefWJ_E[years])

######## Joining mc 

hMCB_M = hist_mcB_nJetGood_M[samples_mc[0]]
hMCD_M = hist_mcD_nJetGood_M[samples_mc[0]]
hMCB_E = hist_mcB_nJetGood_E[samples_mc[0]]
hMCD_E = hist_mcD_nJetGood_E[samples_mc[0]]
hMCB_lepiso_E = hist_mcB_lepton_iso_E[samples_mc[0]]
hMCD_lepiso_E = hist_mcD_lepton_iso_E[samples_mc[0]]
hMCB_transM_E = hist_mcB_transverse_mass_E[samples_mc[0]]
hMCD_transM_E = hist_mcD_transverse_mass_E[samples_mc[0]]

print(type(hMCB_M))
print(type(hist_mcB_nJetGood_M[samples_mc[1]]))

for s in samples_mc:
  if s != samples_mc[0]:
    hMCB_M.Add(hist_mcB_nJetGood_M[s].GetPtr())
    hMCD_M.Add(hist_mcD_nJetGood_M[s].GetPtr())
    hMCB_E.Add(hist_mcB_nJetGood_E[s].GetPtr())
    hMCD_E.Add(hist_mcD_nJetGood_E[s].GetPtr())
    hMCB_lepiso_E.Add(hist_mcB_lepton_iso_E[s].GetPtr())
    hMCD_lepiso_E.Add(hist_mcD_lepton_iso_E[s].GetPtr())
    hMCB_transM_E.Add(hist_mcB_transverse_mass_E[s].GetPtr())
    hMCD_transM_E.Add(hist_mcD_transverse_mass_E[s].GetPtr())


Nb_mc_M = hMCB_M.Integral()
Nd_mc_M = hMCD_M.Integral()
Nb_mc_E = hMCB_E.Integral()
Nd_mc_E = hMCD_E.Integral()

for s in samples_data:
  if s[-1] == "M": 
    Nb_data_M = hist_dataB_nJetGood_M[s].Integral()
    Nd_data_M = hist_dataD_nJetGood_M[s].Integral()
  if s[-1] == "E":
    Nb_data_E = hist_dataB_nJetGood_E[s].Integral()
    Nd_data_E = hist_dataD_nJetGood_E[s].Integral()

print('SL up to v4')

print('Number of data events for the muon channel in the B set is',str(Nb_data_M))
print('Number of data events for the muon channel in the D set is',str(Nd_data_M))
print('Number of mc events for the muon channel in the B set is',str(Nb_mc_M))
print('Number of mc events for the muon channel in the D set is',str(Nd_mc_M))

print('Number of data events for the electron channel in the B set is',str(Nb_data_E))
print('Number of data events for the electron channel in the D set is',str(Nd_data_E))
print('Number of mc events for the electron channel in the B set is',str(Nb_mc_E))
print('Number of mc events for the electron channel in the D set is',str(Nd_mc_E))

print('SL channel')

print('f_BD factor for the muon channel is', str((Nb_data_M - Nb_mc_M)/(Nd_data_M - Nd_mc_M)))
print('f_BD factor for the electron channel is', str((Nb_data_E - Nb_mc_E)/(Nd_data_M - Nd_mc_E)))


print('Ended successfully')
