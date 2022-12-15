
###################################                                        
###################################                                        
#######    QCD estimation   #######
###################################                                        
###################################                                        

## Meant to be run locally, not sent to condor 

## QCD background shape estimation, obtention of a dataset simulating QCD background
## Only for muon channel, we checked electron channel doesn't need a QCD sample

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
parser.add_argument("--year", type=string, default="2016",
                    help="Select data year to run")
parser.add_argument("--notfull", action="store_true", default=False,
                    help="Use just a range of the sample")
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Perform ssos substraction")
parser.add_argument('-l','--list', nargs='+', help='range of sample to use')
# Use like:
# python arg.py -l 1234 2345 3456 4567

args = parser.parse_args()

if (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): years = str(args.year)
else: raise NameError('Incorrect year')

samples_data = []

special_weights = ["Wjets0J2016","Wjets1J2016","Wjets2J2016","DY0J2016","DY1J2016","DY2J2016",
        "Wjets0J2016B","Wjets1J2016B","Wjets2J2016B","DY0J2016B","DY1J2016B","DY2J2016B",
        "Wjets0J2017","Wjets1J2017","Wjets2J2017","DY0J2017","DY1J2017","DY2J2017",
        "Wjets0J2018","Wjets1J2018","Wjets2J2018","DY0J2018","DY1J2018","DY2J2018"]

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

## f_BD factors

f_BD = {}

#f_BD["2016"] = 1.1873
#f_BD["2016B"] = 1.8003
#f_BD["2017"] = 1.2303
#f_BD["2018"] = 1.1872

### temporary

f_BD["2016"] = 2.1451
f_BD["2016B"] = 3.2609
f_BD["2017"] = 2.3450
f_BD["2018"] = 2.1792

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

print('QCD muon channel samples for SL')

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

## Número de muones dentro de un jet
gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muoninjet(UInt_t nmu, Vint mu_jetid, Vint mu_good) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
                if (mu_good.size()>0){
                        if (i!=mu_good[0] && mu_jetid[i]!=-1){
                                vb.push_back(i);
                        }
                } else {
                        if (mu_jetid[i]!=-1){
                                vb.push_back(i);
                        }
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

## pT component calculations

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto pTsum(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            auto ptot = plep1+plep2+phad1+phad2;
            auto plepX = plep.Px(); auto plepY = plep.Py(); auto phadX = phad.Px(); auto phadY = phad.Py(); auto ptotX = ptot.Px(); auto ptotY = ptot.Py();
            float plepmod; float phadmod; float ptotmod;
            plepmod = std::sqrt(plepX*plepX + plepY*plepY);
            phadmod = std::sqrt(phadX*phadX + phadY*phadY);
            ptotmod = std::sqrt(ptotX*ptotX + ptotY*ptotY);
            return ptotmod/(plepmod+phadmod);
      };
      auto pTprod(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            auto plepX = plep.Px(); auto plepY = plep.Py(); auto phadX = phad.Px(); auto phadY = phad.Py();
            float plepmod; float phadmod; float ptotmod;
            plepmod = std::sqrt(plepX*plepX + plepY*plepY);
            phadmod = std::sqrt(phadX*phadX + phadY*phadY);
            ptotmod = plepX*phadX + plepY*phadY;
            return ptotmod/(plepmod*phadmod);
      };
      auto variousSUM(Vint mu_good, Vint el_good, Vint jet_muon, UInt_t jet_notmuon, Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_mass,Vfloat el_pt, Vfloat el_eta, Vfloat el_phi, Vfloat el_mass, Float_t met_pt, Float_t met_phi, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_mass) {
            vector<float> vb;
            ROOT::Math::PtEtaPhiMVector plep1;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
            }
            ROOT::Math::PtEtaPhiMVector plep2(met_pt, 0., met_phi, 0.);
            ROOT::Math::PtEtaPhiMVector phad1(jet_pt[jet_muon[0]], jet_eta[jet_muon[0]], jet_phi[jet_muon[0]], jet_mass[jet_muon[0]]);
            ROOT::Math::PtEtaPhiMVector phad2(jet_pt[jet_notmuon], jet_eta[jet_notmuon], jet_phi[jet_notmuon], jet_mass[jet_notmuon]);
            auto plep = plep1+plep2;
            auto phad = phad1+phad2;
            vb.push_back(ROOT::VecOps::DeltaR(plep1.Eta(),phad.Eta(),plep1.Phi(),phad.Phi()));
            vb.push_back(fabs(plep2.Phi()-phad.Phi()));
            vb.push_back(phad.Eta());
            vb.push_back(phad.Pt());
            vb.push_back(fabs(plep.Phi()-phad.Phi()));
            vb.push_back(ROOT::VecOps::DeltaR(plep.Eta(),phad.Eta(),plep.Phi(),phad.Phi()));
            vb.push_back(fabs(plep1.Phi()-phad.Phi()));
            vb.push_back(fabs(plep.Eta()-phad.Eta()));
            vb.push_back(fabs(plep1.Eta()-phad.Eta()));
            vb.push_back(fabs(plep.Pt()-phad.Pt()));
            vb.push_back(fabs(plep1.Pt()-phad.Pt()));
            vb.push_back(phad2.Pt()*std::sin(phad1.Phi()-phad2.Phi()));
            vb.push_back(phad.Pt()/(phad1.Pt()+phad2.Pt()));
            vb.push_back(plep.Pt());
            return vb;
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
      if s in special_weights:
             df_mc[s] = df_mc[s].Define('weight_aux','fabs(genWeight) > 0 ? genWeight/fabs(genWeight) : 0')
      else:
             df_mc[s] = df_mc[s].Define('weight_aux','1')

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
	df_data[s] = df_data[s].Define('jet_muon_nmu','Jet_nMuons[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('jet_notmuon_pt','nJetGood>1 ? Jet_pt[JetnotMuonInd] : 0')
	df_data[s] = df_data[s].Define('jet_muon_eta','Jet_eta[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('jet_notmuon_eta','nJetGood>1 ? Jet_eta[JetnotMuonInd] : 0')
	df_data[s] = df_data[s].Define('jet_notmuon_mass','nJetGood>1 ? Jet_mass[JetnotMuonInd] : 0')
	df_data[s] = df_data[s].Define('jet_notmuon_qgl','nJetGood>1 ? Jet_qgl[JetnotMuonInd] : 0')
	df_data[s] = df_data[s].Define('jet_notmuon_nmu','nJetGood>1 ? Jet_nMuons[JetnotMuonInd] : 0')
	df_data[s] = df_data[s].Define('muon_jet_pt','Muon_pt[MuonJetInd[0]]')
	df_data[s] = df_data[s].Define('muon_jet_eta','Muon_eta[MuonJetInd[0]]')
	df_data[s] = df_data[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Jet_pt[JetnotMuonInd],Jet_eta[JetnotMuonInd],Jet_phi[JetnotMuonInd],Jet_mass[JetnotMuonInd]) : 0')
	df_data[s] = df_data[s].Define('deltaR_jetM_jetNM','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotMuonInd], Jet_eta[JetMuonInd[0]] , Jet_phi[JetnotMuonInd], Jet_phi[JetMuonInd[0]])  : 10')
	df_data[s] = df_data[s].Define('deltaphi_jetM_jetNM','fabs(Jet_phi[JetMuonInd[0]]-Jet_phi[JetnotMuonInd])')
	df_data[s] = df_data[s].Define('deltaeta_jetM_jetNM','fabs(Jet_eta[JetMuonInd[0]]-Jet_eta[JetnotMuonInd])')
	df_data[s] = df_data[s].Define('deltapt_jetM_jetNM','fabs(Jet_pt[JetMuonInd[0]]-Jet_pt[JetnotMuonInd])')
	df_data[s] = df_data[s].Define('tracks_jetM','Jet_nConstituents[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('tracks_jetNM','nJetGood>1 ? Jet_nConstituents[JetnotMuonInd] : 0')
	df_data[s] = df_data[s].Define('EMN_jetM','Jet_neEmEF[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('EMC_jetM','Jet_chEmEF[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('EMtotal_jetM','Jet_chEmEF[JetMuonInd[0]]+Jet_neEmEF[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('muon_jet_mva','Muon_softMva[MuonJetInd[0]]')
	df_data[s] = df_data[s].Define('muon_jet_relpt','Muon_pt[MuonJetInd[0]]/Jet_pt[JetMuonInd[0]]')                                                                                        
	df_data[s] = df_data[s].Define('MuoninJetAux','muoninjet(nMuon, Muon_jetIdx, MuonGoodInd)')
	df_data[s] = df_data[s].Define('nMuoninJet','MuoninJetAux.size()')
	df_data[s] = df_data[s].Define('muon_jet_sigxy','Muon_dxyErr[MuonJetInd[0]]>0 ? Muon_dxy[MuonJetInd[0]]/Muon_dxyErr[MuonJetInd[0]] : -10')
	df_data[s] = df_data[s].Define('muon_jet_sigz','Muon_dzErr[MuonJetInd[0]]>0 ? Muon_dz[MuonJetInd[0]]/Muon_dzErr[MuonJetInd[0]] : -10')
	df_data[s] = df_data[s].Define('muon_jet_r','pow(pow(Muon_dz[MuonJetInd[0]],2) + pow(Muon_dxy[MuonJetInd[0]],2),0.5)')
	df_data[s] = df_data[s].Define('muon_jet_Err','muon_jet_r > 0 ? pow(pow(Muon_dz[MuonJetInd[0]]*Muon_dzErr[MuonJetInd[0]],2)+pow(Muon_dxy[MuonJetInd[0]]*Muon_dxyErr[MuonJetInd[0]],2),0.5)/muon_jet_r : -10')
	df_data[s] = df_data[s].Define('muon_jet_sigr','(muon_jet_Err>0) ? muon_jet_r/muon_jet_Err : -10')
	df_data[s] = df_data[s].Define('pT_sum','pTsum(MuonGoodInd, ElectronGoodInd, JetMuonInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
	df_data[s] = df_data[s].Define('pT_product','pTprod(MuonGoodInd, ElectronGoodInd, JetMuonInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
	df_data[s] = df_data[s].Define('aux_various','variousSUM(MuonGoodInd, ElectronGoodInd, JetMuonInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
	df_data[s] = df_data[s].Define('deltaR_lep_2jets','aux_various[0]')
	df_data[s] = df_data[s].Define('deltaphi_MET_2jets','aux_various[1]')
	df_data[s] = df_data[s].Define('eta_2jets','aux_various[2]')
	df_data[s] = df_data[s].Define('pt_2jets','aux_various[3]')
	df_data[s] = df_data[s].Define('deltaphi_lephad','aux_various[4]')
	df_data[s] = df_data[s].Define('deltaR_lephad','aux_various[5]')
	df_data[s] = df_data[s].Define('deltaphi_lep_2jets','aux_various[6]')
	df_data[s] = df_data[s].Define('deltaeta_lephad','aux_various[7]')
	df_data[s] = df_data[s].Define('deltaeta_lep_2jets','aux_various[8]')
	df_data[s] = df_data[s].Define('deltapt_lephad','aux_various[9]')
	df_data[s] = df_data[s].Define('deltapt_lep_2jets','aux_various[10]')
	df_data[s] = df_data[s].Define('pT_proy','aux_various[11]')
	df_data[s] = df_data[s].Define('pT_sum_2J','aux_various[12]')
	df_data[s] = df_data[s].Define('pT_Wlep','aux_various[13]')
	df_data[s] = df_data[s].Define('jet_muon_btag','Jet_btagDeepFlavB[JetMuonInd[0]]')
	df_data[s] = df_data[s].Define('jet_notmuon_btag','Jet_btagDeepFlavB[JetnotMuonInd]')
	df_data[s] = df_data[s].Define('jet_notmuon_deeptagG','Jet_btagDeepFlavQG[JetnotMuonInd]')
	df_data[s] = df_data[s].Define('jet_notmuon_deeptagC','Jet_btagDeepFlavCvL[JetnotMuonInd]')

for s in samples_mc:
        df_mc[s] = df_mc[s].Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood==2').Filter('MuonJetInd.size() >= 1').Filter('MuonJetGood')
        df_mc[s] = df_mc[s].Define('JetnotMuonInd','nJetGood>1 ? (JetMuonInd[0] == JetGoodInd[0] ? JetGoodInd[1] : JetGoodInd[0]) : -1')
        ### hists definitions
        df_mc[s] = df_mc[s].Define('jet_muon_pt','Jet_pt[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('jet_muon_nmu','Jet_nMuons[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('jet_notmuon_pt','nJetGood>1 ? Jet_pt[JetnotMuonInd] : 0')
        df_mc[s] = df_mc[s].Define('jet_muon_eta','Jet_eta[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('jet_notmuon_eta','nJetGood>1 ? Jet_eta[JetnotMuonInd] : 0')
        df_mc[s] = df_mc[s].Define('jet_notmuon_mass','nJetGood>1 ? Jet_mass[JetnotMuonInd] : 0')
        df_mc[s] = df_mc[s].Define('jet_notmuon_qgl','nJetGood>1 ? Jet_qgl[JetnotMuonInd] : 0')
        df_mc[s] = df_mc[s].Define('jet_notmuon_nmu','nJetGood>1 ? Jet_nMuons[JetnotMuonInd] : 0')
        df_mc[s] = df_mc[s].Define('muon_jet_pt','Muon_pt[MuonJetInd[0]]')
        df_mc[s] = df_mc[s].Define('muon_jet_eta','Muon_eta[MuonJetInd[0]]')
        df_mc[s] = df_mc[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Jet_pt[JetnotMuonInd],Jet_eta[JetnotMuonInd],Jet_phi[JetnotMuonInd],Jet_mass[JetnotMuonInd]) : 0')
        df_mc[s] = df_mc[s].Define('deltaR_jetM_jetNM','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotMuonInd], Jet_eta[JetMuonInd[0]] , Jet_phi[JetnotMuonInd], Jet_phi[JetMuonInd[0]])  : 10')
        df_mc[s] = df_mc[s].Define('deltaphi_jetM_jetNM','fabs(Jet_phi[JetMuonInd[0]]-Jet_phi[JetnotMuonInd])')
        df_mc[s] = df_mc[s].Define('deltaeta_jetM_jetNM','fabs(Jet_eta[JetMuonInd[0]]-Jet_eta[JetnotMuonInd])')
        df_mc[s] = df_mc[s].Define('deltapt_jetM_jetNM','fabs(Jet_pt[JetMuonInd[0]]-Jet_pt[JetnotMuonInd])')
        df_mc[s] = df_mc[s].Define('tracks_jetM','Jet_nConstituents[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('tracks_jetNM','nJetGood>1 ? Jet_nConstituents[JetnotMuonInd] : 0')
        df_mc[s] = df_mc[s].Define('EMN_jetM','Jet_neEmEF[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('EMC_jetM','Jet_chEmEF[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('EMtotal_jetM','Jet_chEmEF[JetMuonInd[0]]+Jet_neEmEF[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('muon_jet_mva','Muon_softMva[MuonJetInd[0]]')
        df_mc[s] = df_mc[s].Define('muon_jet_relpt','Muon_pt[MuonJetInd[0]]/Jet_pt[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('MuoninJetAux','muoninjet(nMuon, Muon_jetIdx, MuonGoodInd)')
        df_mc[s] = df_mc[s].Define('nMuoninJet','MuoninJetAux.size()')
        df_mc[s] = df_mc[s].Define('muon_jet_sigxy','Muon_dxyErr[MuonJetInd[0]]>0 ? Muon_dxy[MuonJetInd[0]]/Muon_dxyErr[MuonJetInd[0]] : -10')
        df_mc[s] = df_mc[s].Define('muon_jet_sigz','Muon_dzErr[MuonJetInd[0]]>0 ? Muon_dz[MuonJetInd[0]]/Muon_dzErr[MuonJetInd[0]] : -10')
        df_mc[s] = df_mc[s].Define('muon_jet_r','pow(pow(Muon_dz[MuonJetInd[0]],2) + pow(Muon_dxy[MuonJetInd[0]],2),0.5)')
        df_mc[s] = df_mc[s].Define('muon_jet_Err','muon_jet_r > 0 ? pow(pow(Muon_dz[MuonJetInd[0]]*Muon_dzErr[MuonJetInd[0]],2)+pow(Muon_dxy[MuonJetInd[0]]*Muon_dxyErr[MuonJetInd[0]],2),0.5)/muon_jet_r : -10')
        df_mc[s] = df_mc[s].Define('muon_jet_sigr','(muon_jet_Err>0) ? muon_jet_r/muon_jet_Err : -10')
        df_mc[s] = df_mc[s].Define('pT_sum','pTsum(MuonGoodInd, ElectronGoodInd, JetMuonInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
        df_mc[s] = df_mc[s].Define('pT_product','pTprod(MuonGoodInd, ElectronGoodInd, JetMuonInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
        df_mc[s] = df_mc[s].Define('aux_various','variousSUM(MuonGoodInd, ElectronGoodInd, JetMuonInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
        df_mc[s] = df_mc[s].Define('deltaR_lep_2jets','aux_various[0]')
        df_mc[s] = df_mc[s].Define('deltaphi_MET_2jets','aux_various[1]')
        df_mc[s] = df_mc[s].Define('eta_2jets','aux_various[2]')
        df_mc[s] = df_mc[s].Define('pt_2jets','aux_various[3]')
        df_mc[s] = df_mc[s].Define('deltaphi_lephad','aux_various[4]')
        df_mc[s] = df_mc[s].Define('deltaR_lephad','aux_various[5]')
        df_mc[s] = df_mc[s].Define('deltaphi_lep_2jets','aux_various[6]')
        df_mc[s] = df_mc[s].Define('deltaeta_lephad','aux_various[7]')
        df_mc[s] = df_mc[s].Define('deltaeta_lep_2jets','aux_various[8]')
        df_mc[s] = df_mc[s].Define('deltapt_lephad','aux_various[9]')
        df_mc[s] = df_mc[s].Define('deltapt_lep_2jets','aux_various[10]')
        df_mc[s] = df_mc[s].Define('pT_proy','aux_various[11]')
        df_mc[s] = df_mc[s].Define('pT_sum_2J','aux_various[12]')
        df_mc[s] = df_mc[s].Define('pT_Wlep','aux_various[13]')
        df_mc[s] = df_mc[s].Define('jet_muon_btag','Jet_btagDeepFlavB[JetMuonInd[0]]')
        df_mc[s] = df_mc[s].Define('jet_notmuon_btag','Jet_btagDeepFlavB[JetnotMuonInd]')
        df_mc[s] = df_mc[s].Define('jet_notmuon_deeptagG','Jet_btagDeepFlavQG[JetnotMuonInd]')
        df_mc[s] = df_mc[s].Define('jet_notmuon_deeptagC','Jet_btagDeepFlavCvL[JetnotMuonInd]')

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
	df_data_electron[s] = df_data_electron[s].Define('lepton_mva', 'Electron_mvaFall17V2noIso[ElectronGoodInd[0]]')
	## New cuts
	df_data_muon[s] = df_data_muon[s].Filter('InvM_muon_jet >12').Filter('InvM_muon_jet > 110 || InvM_muon_jet < 70')
	#df_data_electron[s] = df_data_electron[s].Filter('transverse_mass > 50')
	#df_data_muon[s] = df_data_muon[s].Filter('transverse_mass > 50')
	df_data_electron[s] = df_data_electron[s].Filter('muon_jet_relpt<0.5')
	df_data_muon[s] = df_data_muon[s].Filter('muon_jet_relpt<0.5')
	df_data_electron[s] = df_data_electron[s].Filter('EMtotal_jetM<0.4')
	df_data_muon[s] = df_data_muon[s].Filter('EMtotal_jetM<0.4')
	df_data_electron[s] = df_data_electron[s].Filter('jet_notmuon_btag<'+str(cuts_btag[years][1]))
	df_data_muon[s] = df_data_muon[s].Filter('jet_notmuon_btag<'+str(cuts_btag[years][1]))
	df_data_electron[s] = df_data_electron[s].Filter('jet_muon_btag<'+str(cuts_btag[years][2]))
	df_data_muon[s] = df_data_muon[s].Filter('jet_muon_btag<'+str(cuts_btag[years][2]))
	#df_data_electron[s] = df_data_electron[s].Filter('InvM_2jets<100').Filter('InvM_2jets>60')
	#df_data_muon[s] = df_data_muon[s].Filter('InvM_2jets<100').Filter('InvM_2jets>60')

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
        df_mc_electron[s] = df_mc_electron[s].Define('lepton_mva', 'Electron_mvaFall17V2noIso[ElectronGoodInd[0]]')
        ## New cuts
        df_mc_muon[s] = df_mc_muon[s].Filter('InvM_muon_jet >12').Filter('InvM_muon_jet > 110 || InvM_muon_jet < 70')
        #df_mc_electron[s] = df_mc_electron[s].Filter('transverse_mass > 50')
        #df_mc_muon[s] = df_mc_muon[s].Filter('transverse_mass > 50')
        df_mc_electron[s] = df_mc_electron[s].Filter('muon_jet_relpt<0.5')
        df_mc_muon[s] = df_mc_muon[s].Filter('muon_jet_relpt<0.5')
        df_mc_electron[s] = df_mc_electron[s].Filter('EMtotal_jetM<0.4')
        df_mc_muon[s] = df_mc_muon[s].Filter('EMtotal_jetM<0.4')
        df_mc_electron[s] = df_mc_electron[s].Filter('jet_notmuon_btag<'+str(cuts_btag[years][1]))
        df_mc_muon[s] = df_mc_muon[s].Filter('jet_notmuon_btag<'+str(cuts_btag[years][1]))
        df_mc_electron[s] = df_mc_electron[s].Filter('jet_muon_btag<'+str(cuts_btag[years][2]))
        df_mc_muon[s] = df_mc_muon[s].Filter('jet_muon_btag<'+str(cuts_btag[years][2]))
        #df_mc_electron[s] = df_mc_electron[s].Filter('InvM_2jets<100').Filter('InvM_2jets>60')
        #df_mc_muon[s] = df_mc_muon[s].Filter('InvM_2jets<100').Filter('InvM_2jets>60')

######################################
######### B and D definition #########
######################################

df_dataC_muon = {}
df_dataC_electron = {}

for s in samples_data:
	df_dataC_muon[s] = df_data_muon[s].Filter('lepton_iso>0.2').Filter('transverse_mass>50')
	df_dataC_electron[s] = df_data_electron[s].Filter('lepton_iso>0.2').Filter('transverse_mass>50')

df_mcC_muon = {}
df_mcC_electron = {}

for s in samples_mc:
        df_mcC_muon[s] = df_mc_muon[s].Filter('lepton_iso>0.2').Filter('transverse_mass>50')
        df_mcC_electron[s] = df_mc_electron[s].Filter('lepton_iso>0.2').Filter('transverse_mass>50')

############################################################
####################     HISTS    ##########################
############################################################

for s in samples_data:
   if s[-1] == "M":
        hist_dataC_nJetGood_M = df_dataC_muon[s].Histo1D(("nJetGood_M","",10,0,10),"nJetGood","weightSSOS")
        hist_dataC_nMuoninJet_M = df_dataC_muon[s].Histo1D(("nMuoninJet_M","",10,0,10),"nMuoninJet","weightSSOS")
        hist_dataC_jet_muon_pt_M = df_dataC_muon[s].Histo1D(("jet_muon_pt_M","",50,20,120),"jet_muon_pt","weightSSOS")
        hist_dataC_jet_muon_nmu_M = df_dataC_muon[s].Histo1D(("jet_muon_nmu_M","",10,0,10),"jet_muon_nmu","weightSSOS")
        hist_dataC_jet_notmuon_pt_M = df_dataC_muon[s].Histo1D(("jet_not_muon_pt_M","",50,20,120),"jet_notmuon_pt","weightSSOS")
        hist_dataC_jet_muon_eta_M = df_dataC_muon[s].Histo1D(("jet_muon_eta_M","",80,-4,4),"jet_muon_eta","weightSSOS")
        hist_dataC_jet_notmuon_eta_M = df_dataC_muon[s].Histo1D(("jet_not_muon_eta_M","",80,-4,4),"jet_notmuon_eta","weightSSOS")
        hist_dataC_jet_notmuon_nmu_M = df_dataC_muon[s].Histo1D(("jet_notmuon_nmu_M","",10,0,10),"jet_notmuon_nmu","weightSSOS")
        hist_dataC_jet_notmuon_mass_M = df_dataC_muon[s].Histo1D(("jet_notmuon_mass_M","",40,0,40),"jet_notmuon_mass","weightSSOS")
        hist_dataC_jet_notmuon_qgl_M = df_dataC_muon[s].Histo1D(("jet_notmuon_qgl_M","",100,0,1),"jet_notmuon_qgl","weightSSOS")
        hist_dataC_lepton_pt_M = df_dataC_muon[s].Histo1D(("lepton_pt_M","",50,20,120),"lepton_pt","weightSSOS")
        hist_dataC_lepton_eta_M = df_dataC_muon[s].Histo1D(("lepton_eta_M","",80,-4,4),"lepton_eta","weightSSOS")
        hist_dataC_muon_jet_pt_M = df_dataC_muon[s].Histo1D(("muon_jet_pt_M","",50,0,25),"muon_jet_pt","weightSSOS")
        hist_dataC_muon_jet_eta_M = df_dataC_muon[s].Histo1D(("muon_jet_eta_M","",80,-4,4),"muon_jet_eta","weightSSOS")
        hist_dataC_InvM_2jets_M = df_dataC_muon[s].Histo1D(("InvM_2jets_M","",100,0,300),"InvM_2jets","weightSSOS")
        hist_dataC_InvM_jetM_lepM = df_dataC_muon[s].Histo1D(("InvM_jetM_lepM","",100,0,300),"InvM_jetM_lep","weightSSOS")
        hist_dataC_MET_M = df_dataC_muon[s].Histo1D(("MET_pt_M","",100,0,150),"MET_pt","weightSSOS")
        hist_dataC_deltaR_jetM_lepM = df_dataC_muon[s].Histo1D(("deltaR_jetM_lepM","",100,0,5),"deltaR_jetM_lepM","weightSSOS")
        hist_dataC_deltaR_jetM_jetNM_M = df_dataC_muon[s].Histo1D(("deltaR_jetM_jetNM_M","",100,0,5),"deltaR_jetM_jetNM","weightSSOS")
        hist_dataC_deltaphi_jetM_jetNM_M = df_dataC_muon[s].Histo1D(("deltaphi_jetM_jetNM_M","",100,0,5),"deltaphi_jetM_jetNM","weightSSOS")
        hist_dataC_deltaeta_jetM_jetNM_M = df_dataC_muon[s].Histo1D(("deltaeta_jetM_jetNM_M","",100,0,5),"deltaeta_jetM_jetNM","weightSSOS")
        hist_dataC_deltapt_jetM_jetNM_M = df_dataC_muon[s].Histo1D(("deltapt_jetM_jetNM_M","",100,0,50),"deltapt_jetM_jetNM","weightSSOS")
        hist_dataC_tranverse_massM = df_dataC_muon[s].Histo1D(("transverse_massM","",50,0,150),"transverse_mass","weightSSOS")
        hist_dataC_tracks_jetM_M = df_dataC_muon[s].Histo1D(("tracks_jetM_M","",60,0,60),"tracks_jetM","weightSSOS")
        hist_dataC_tracks_jetNM_M = df_dataC_muon[s].Histo1D(("tracks_jetNM_M","",60,0,60),"tracks_jetNM","weightSSOS")
        hist_dataC_EMN_jetM_M = df_dataC_muon[s].Histo1D(("EMN_jetM_M","",60,0,1),"EMN_jetM","weightSSOS")
        hist_dataC_EMC_jetM_M = df_dataC_muon[s].Histo1D(("EMC_jetM_M","",60,0,1),"EMC_jetM","weightSSOS")
        hist_dataC_EMtotal_jetM_M = df_dataC_muon[s].Histo1D(("EMtotal_jetM_M","",60,0,1),"EMtotal_jetM","weightSSOS")
        hist_dataC_InvM_muon_jet_M = df_dataC_muon[s].Histo1D(("InvM_muon_jet_M","",50,0,200),"InvM_muon_jet","weightSSOS")
        hist_dataC_muon_jet_mva_M = df_dataC_muon[s].Histo1D(("muon_jet_mva_M","",50,0,1),"muon_jet_mva","weightSSOS")
        hist_dataC_muon_jet_relpt_M = df_dataC_muon[s].Histo1D(("muon_jet_relpt_M","",50,0,1),"muon_jet_relpt","weightSSOS")
        hist_dataC_muon_jet_sigxy_M = df_dataC_muon[s].Histo1D(("muon_jet_sigxy_M","",100,-4,4),"muon_jet_sigxy","weightSSOS")
        hist_dataC_muon_jet_sigz_M = df_dataC_muon[s].Histo1D(("muon_jet_sigz_M","",100,-4,4),"muon_jet_sigz","weightSSOS")
        hist_dataC_muon_jet_sigr_M = df_dataC_muon[s].Histo1D(("muon_jet_sigr_M","",100,0,6),"muon_jet_sigr","weightSSOS")
        hist_dataC_SSOS_M = df_dataC_muon[s].Histo1D(("SSOS_M","",4,-2,2),"MuonLepSign","weightSSOS")
        hist_dataC_pT_sum_M = df_dataC_muon[s].Histo1D(("pT_sum_M","",100,0,1),"pT_sum","weightSSOS")
        hist_dataC_pT_product_M = df_dataC_muon[s].Histo1D(("pT_product_M","",100,-1,1),"pT_product","weightSSOS")
        hist_dataC_deltaR_lep_2jets_M = df_dataC_muon[s].Histo1D(("deltaR_lep_2jets_M","",100,0,5),"deltaR_lep_2jets","weightSSOS")
        hist_dataC_deltaphi_MET_2jets_M = df_dataC_muon[s].Histo1D(("deltaphi_MET_2jets_M","",100,0,5),"deltaphi_MET_2jets","weightSSOS")
        hist_dataC_deltaphi_lephad_M = df_dataC_muon[s].Histo1D(("deltaphi_lephad_M","",100,0,5),"deltaphi_lephad","weightSSOS")
        hist_dataC_eta_2jets_M = df_dataC_muon[s].Histo1D(("eta_2jets_M","",50,-5,5),"eta_2jets","weightSSOS")
        hist_dataC_pt_2jets_M = df_dataC_muon[s].Histo1D(("pt_2jets_M","",70,0,200),"pt_2jets","weightSSOS")
        hist_dataC_pt_Wlep_M = df_dataC_muon[s].Histo1D(("pt_Wlep_M","",100,0,200),"pT_Wlep","weightSSOS")
        hist_dataC_deltaR_lephad_M = df_dataC_muon[s].Histo1D(("deltaR_lephad_M","",100,0,5),"deltaR_lephad","weightSSOS")
        hist_dataC_deltaphi_lep_2jets_M = df_dataC_muon[s].Histo1D(("deltaphi_lep_2jets_M","",100,0,5),"deltaphi_lep_2jets","weightSSOS")
        hist_dataC_deltaeta_lephad_M = df_dataC_muon[s].Histo1D(("deltaeta_lephad_M","",100,0,5),"deltaeta_lephad","weightSSOS")
        hist_dataC_deltaeta_lep_2jets_M = df_dataC_muon[s].Histo1D(("deltaeta_lep_2jets_M","",100,0,5),"deltaeta_lep_2jets","weightSSOS")
        hist_dataC_deltapt_lephad_M = df_dataC_muon[s].Histo1D(("deltapt_lephad_M","",100,0,5),"deltapt_lephad","weightSSOS")
        hist_dataC_deltapt_lep_2jets_M = df_dataC_muon[s].Histo1D(("deltapt_lep_2jets_M","",100,0,5),"deltapt_lep_2jets","weightSSOS")
        hist_dataC_jet_muon_btag_M = df_dataC_muon[s].Histo1D(("jet_muon_btag_M","",50,0,1),"jet_muon_btag","weightSSOS")
        hist_dataC_jet_notmuon_btag_M = df_dataC_muon[s].Histo1D(("jet_notmuon_btag_M","",50,0,1),"jet_notmuon_btag","weightSSOS")
        hist_dataC_jet_notmuon_deeptagG_M = df_dataC_muon[s].Histo1D(("jet_notmuon_deeptagG_M","",50,0,1),"jet_notmuon_deeptagG","weightSSOS")
        hist_dataC_jet_notmuon_deeptagC_M = df_dataC_muon[s].Histo1D(("jet_notmuon_deeptagC_M","",50,0,1),"jet_notmuon_deeptagC","weightSSOS")
        hist_dataC_pT_proy_M = df_dataC_muon[s].Histo1D(("pT_proy_M","",100,-100,100),"pT_proy","weightSSOS")
        hist_dataC_pT_sum_2J_M = df_dataC_muon[s].Histo1D(("pT_sum_2J_M","",100,0,1),"pT_sum_2J","weightSSOS")

hist_mcC_nJetGood_M = {}
hist_mcC_nMuoninJet_M = {}
hist_mcC_jet_muon_pt_M = {}
hist_mcC_jet_muon_nmu_M = {}
hist_mcC_jet_notmuon_pt_M = {}
hist_mcC_jet_muon_eta_M = {}
hist_mcC_jet_notmuon_eta_M = {}
hist_mcC_jet_notmuon_nmu_M = {}
hist_mcC_jet_notmuon_mass_M = {}
hist_mcC_jet_notmuon_qgl_M = {}
hist_mcC_lepton_pt_M = {}
hist_mcC_lepton_eta_M = {}
hist_mcC_muon_jet_pt_M = {}
hist_mcC_muon_jet_eta_M = {}
hist_mcC_InvM_2jets_M = {}
hist_mcC_InvM_jetM_lepM = {}
hist_mcC_deltaR_jetM_lepM = {}
hist_mcC_deltaR_jetM_jetNM_M = {}
hist_mcC_deltaphi_jetM_jetNM_M = {}
hist_mcC_deltaeta_jetM_jetNM_M = {}
hist_mcC_deltapt_jetM_jetNM_M = {}
hist_mcC_MET_M = {}
hist_mcC_tranverse_massM = {}
hist_mcC_tracks_jetM_M = {}
hist_mcC_tracks_jetNM_M = {}
hist_mcC_EMN_jetM_M = {}
hist_mcC_EMC_jetM_M = {}
hist_mcC_EMtotal_jetM_M = {}
hist_mcC_InvM_muon_jet_M = {}
hist_mcC_muon_jet_mva_M = {}
hist_mcC_muon_jet_relpt_M = {}
hist_mcC_SSOS_M = {}
hist_mcC_muon_jet_sigxy_M = {}
hist_mcC_muon_jet_sigz_M = {}
hist_mcC_muon_jet_sigr_M = {}
hist_mcC_pT_sum_M = {}
hist_mcC_pT_product_M = {}
hist_mcC_deltaR_lep_2jets_M = {}
hist_mcC_deltaphi_MET_2jets_M = {}
hist_mcC_deltaphi_lephad_M = {}
hist_mcC_eta_2jets_M = {}
hist_mcC_pt_2jets_M = {}
hist_mcC_pt_Wlep_M = {}
hist_mcC_jet_muon_btag_M = {}
hist_mcC_jet_notmuon_btag_M = {}
hist_mcC_jet_notmuon_deeptagG_M = {}
hist_mcC_jet_notmuon_deeptagC_M = {}
hist_mcC_deltaR_lephad_M = {}
hist_mcC_deltaphi_lep_2jets_M = {}
hist_mcC_deltaeta_lephad_M = {}
hist_mcC_deltaeta_lep_2jets_M = {}
hist_mcC_deltapt_lephad_M = {}
hist_mcC_deltapt_lep_2jets_M = {}
hist_mcC_pT_proy_M = {}
hist_mcC_pT_sum_2J_M = {}

for s in samples_mc:
        hist_mcC_nJetGood_M[s] = df_mcC_muon[s].Histo1D(("nJetGood_M","",10,0,10),"nJetGood","weightSSOS")

        hist_mcC_nMuoninJet_M[s] = df_mcC_muon[s].Histo1D(("nMuoninJet_M","",10,0,10),"nMuoninJet","weightSSOS")

        hist_mcC_jet_muon_pt_M[s] = df_mcC_muon[s].Histo1D(("jet_muon_pt_M","",50,20,120),"jet_muon_pt","weightSSOS")
        hist_mcC_jet_muon_nmu_M[s] = df_mcC_muon[s].Histo1D(("jet_muon_nmu_M","",10,0,10),"jet_muon_nmu","weightSSOS")
        hist_mcC_jet_notmuon_pt_M[s] = df_mcC_muon[s].Histo1D(("jet_not_muon_pt_M","",50,20,120),"jet_notmuon_pt","weightSSOS")
        hist_mcC_jet_muon_eta_M[s] = df_mcC_muon[s].Histo1D(("jet_muon_eta_M","",80,-4,4),"jet_muon_eta","weightSSOS")
        hist_mcC_jet_notmuon_eta_M[s] = df_mcC_muon[s].Histo1D(("jet_not_muon_eta_M","",80,-4,4),"jet_notmuon_eta","weightSSOS")
        hist_mcC_jet_notmuon_nmu_M[s] = df_mcC_muon[s].Histo1D(("jet_notmuon_nmu_M","",10,0,10),"jet_notmuon_nmu","weightSSOS")
        hist_mcC_jet_notmuon_mass_M[s] = df_mcC_muon[s].Histo1D(("jet_notmuon_mass_M","",40,0,40),"jet_notmuon_mass","weightSSOS")
        hist_mcC_jet_notmuon_qgl_M[s] = df_mcC_muon[s].Histo1D(("jet_notmuon_qgl_M","",100,0,1),"jet_notmuon_qgl","weightSSOS")

        hist_mcC_lepton_pt_M[s] = df_mcC_muon[s].Histo1D(("lepton_pt_M","",50,20,120),"lepton_pt","weightSSOS")
        hist_mcC_lepton_eta_M[s] = df_mcC_muon[s].Histo1D(("lepton_eta_M","",80,-4,4),"lepton_eta","weightSSOS")

        hist_mcC_muon_jet_pt_M[s] = df_mcC_muon[s].Histo1D(("muon_jet_pt_M","",50,0,25),"muon_jet_pt","weightSSOS")
        hist_mcC_muon_jet_eta_M[s] = df_mcC_muon[s].Histo1D(("muon_jet_eta_M","",80,-4,4),"muon_jet_eta","weightSSOS")

        hist_mcC_InvM_2jets_M[s] = df_mcC_muon[s].Histo1D(("InvM_2jets_M","",100,0,300),"InvM_2jets","weightSSOS")
        hist_mcC_InvM_jetM_lepM[s] = df_mcC_muon[s].Histo1D(("InvM_jetM_lepM","",100,0,300),"InvM_jetM_lep","weightSSOS")

        hist_mcC_MET_M[s] = df_mcC_muon[s].Histo1D(("MET_pt_M","",100,0,150),"MET_pt","weightSSOS")

        hist_mcC_deltaR_jetM_lepM[s] = df_mcC_muon[s].Histo1D(("deltaR_jetM_lepM","",100,0,5),"deltaR_jetM_lepM","weightSSOS")
        hist_mcC_deltaR_jetM_jetNM_M[s] = df_mcC_muon[s].Histo1D(("deltaR_jetM_jetNM_M","",100,0,5),"deltaR_jetM_jetNM","weightSSOS")

        hist_mcC_deltaphi_jetM_jetNM_M[s] = df_mcC_muon[s].Histo1D(("deltaphi_jetM_jetNM_M","",100,0,5),"deltaphi_jetM_jetNM","weightSSOS")
        hist_mcC_deltaeta_jetM_jetNM_M[s] = df_mcC_muon[s].Histo1D(("deltaeta_jetM_jetNM_M","",100,0,5),"deltaeta_jetM_jetNM","weightSSOS")
        hist_mcC_deltapt_jetM_jetNM_M[s] = df_mcC_muon[s].Histo1D(("deltapt_jetM_jetNM_M","",100,0,50),"deltapt_jetM_jetNM","weightSSOS")

        hist_mcC_tranverse_massM[s] = df_mcC_muon[s].Histo1D(("transverse_massM","",50,0,150),"transverse_mass","weightSSOS")

        hist_mcC_tracks_jetM_M[s] = df_mcC_muon[s].Histo1D(("tracks_jetM_M","",60,0,60),"tracks_jetM","weightSSOS")
        hist_mcC_tracks_jetNM_M[s] = df_mcC_muon[s].Histo1D(("tracks_jetNM_M","",60,0,60),"tracks_jetNM","weightSSOS")

        hist_mcC_EMN_jetM_M[s] = df_mcC_muon[s].Histo1D(("EMN_jetM_M","",60,0,1),"EMN_jetM","weightSSOS")
        hist_mcC_EMC_jetM_M[s] = df_mcC_muon[s].Histo1D(("EMC_jetM_M","",60,0,1),"EMC_jetM","weightSSOS")

        hist_mcC_EMtotal_jetM_M[s] = df_mcC_muon[s].Histo1D(("EMtotal_jetM_M","",60,0,1),"EMtotal_jetM","weightSSOS")

        hist_mcC_InvM_muon_jet_M[s] = df_mcC_muon[s].Histo1D(("InvM_muon_jet_M","",50,0,200),"InvM_muon_jet","weightSSOS")

        hist_mcC_muon_jet_mva_M[s] = df_mcC_muon[s].Histo1D(("muon_jet_mva_M","",50,0,1),"muon_jet_mva","weightSSOS")
        hist_mcC_muon_jet_relpt_M[s] = df_mcC_muon[s].Histo1D(("muon_jet_relpt_M","",50,0,1),"muon_jet_relpt","weightSSOS")
        hist_mcC_muon_jet_sigxy_M[s] = df_mcC_muon[s].Histo1D(("muon_jet_sigxy_M","",100,-4,4),"muon_jet_sigxy","weightSSOS")
        hist_mcC_muon_jet_sigz_M[s] = df_mcC_muon[s].Histo1D(("muon_jet_sigz_M","",100,-4,4),"muon_jet_sigz","weightSSOS")
        hist_mcC_muon_jet_sigr_M[s] = df_mcC_muon[s].Histo1D(("muon_jet_sigr_M","",100,0,6),"muon_jet_sigr","weightSSOS")

        hist_mcC_SSOS_M[s] = df_mcC_muon[s].Histo1D(("SSOS_M","",4,-2,2),"MuonLepSign","weightSSOS")
        hist_mcC_pT_sum_M[s] = df_mcC_muon[s].Histo1D(("pT_sum_M","",100,0,1),"pT_sum","weightSSOS")
        hist_mcC_pT_product_M[s] = df_mcC_muon[s].Histo1D(("pT_product_M","",100,-1,1),"pT_product","weightSSOS")

        hist_mcC_deltaR_lep_2jets_M[s] = df_mcC_muon[s].Histo1D(("deltaR_lep_2jets_M","",100,0,5),"deltaR_lep_2jets","weightSSOS")
        hist_mcC_deltaphi_MET_2jets_M[s] = df_mcC_muon[s].Histo1D(("deltaphi_MET_2jets_M","",100,0,5),"deltaphi_MET_2jets","weightSSOS")
        hist_mcC_deltaphi_lephad_M[s] = df_mcC_muon[s].Histo1D(("deltaphi_lephad_M","",100,0,5),"deltaphi_lephad","weightSSOS")

        hist_mcC_eta_2jets_M[s] = df_mcC_muon[s].Histo1D(("eta_2jets_M","",50,-5,5),"eta_2jets","weightSSOS")
        hist_mcC_pt_2jets_M[s] = df_mcC_muon[s].Histo1D(("pt_2jets_M","",70,0,200),"pt_2jets","weightSSOS")
        hist_mcC_pt_Wlep_M[s] = df_mcC_muon[s].Histo1D(("pt_Wlep_M","",100,0,200),"pT_Wlep","weightSSOS")

        hist_mcC_deltaR_lephad_M[s] = df_mcC_muon[s].Histo1D(("deltaR_lephad_M","",100,0,5),"deltaR_lephad","weightSSOS")
        hist_mcC_deltaphi_lep_2jets_M[s] = df_mcC_muon[s].Histo1D(("deltaphi_lep_2jets_M","",100,0,5),"deltaphi_lep_2jets","weightSSOS")
        hist_mcC_deltaeta_lephad_M[s] = df_mcC_muon[s].Histo1D(("deltaeta_lephad_M","",100,0,5),"deltaeta_lephad","weightSSOS")
        hist_mcC_deltaeta_lep_2jets_M[s] = df_mcC_muon[s].Histo1D(("deltaeta_lep_2jets_M","",100,0,5),"deltaeta_lep_2jets","weightSSOS")
        hist_mcC_deltapt_lephad_M[s] = df_mcC_muon[s].Histo1D(("deltapt_lephad_M","",100,0,5),"deltapt_lephad","weightSSOS")
        hist_mcC_deltapt_lep_2jets_M[s] = df_mcC_muon[s].Histo1D(("deltapt_lep_2jets_M","",100,0,5),"deltapt_lep_2jets","weightSSOS")

        hist_mcC_jet_muon_btag_M[s] = df_mcC_muon[s].Histo1D(("jet_muon_btag_M","",50,0,1),"jet_muon_btag","weightSSOS")
        hist_mcC_jet_notmuon_btag_M[s] = df_mcC_muon[s].Histo1D(("jet_notmuon_btag_M","",50,0,1),"jet_notmuon_btag","weightSSOS")

        hist_mcC_jet_notmuon_deeptagG_M[s] = df_mcC_muon[s].Histo1D(("jet_notmuon_deeptagG_M","",50,0,1),"jet_notmuon_deeptagG","weightSSOS")
        hist_mcC_jet_notmuon_deeptagC_M[s] = df_mcC_muon[s].Histo1D(("jet_notmuon_deeptagC_M","",50,0,1),"jet_notmuon_deeptagG","weightSSOS")

        hist_mcC_pT_proy_M[s] = df_mcC_muon[s].Histo1D(("pT_proy_M","",100,-100,100),"pT_proy","weightSSOS")

        hist_mcC_pT_sum_2J_M[s] = df_mcC_muon[s].Histo1D(("pT_sum_2J_M","",100,0,1),"pT_sum_2J","weightSSOS")

####### Normalization

for s in samples_data:
  lumiD = lumi_data[s]

for s in samples_mc:
        hist_mcC_nJetGood_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_nMuoninJet_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_muon_pt_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_muon_nmu_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_notmuon_pt_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_muon_eta_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_notmuon_eta_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_notmuon_nmu_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_notmuon_qgl_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_notmuon_mass_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_lepton_pt_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_lepton_eta_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_muon_jet_pt_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_muon_jet_eta_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_InvM_2jets_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_InvM_jetM_lepM[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_MET_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltaR_jetM_lepM[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltaR_jetM_jetNM_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltaeta_jetM_jetNM_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltapt_jetM_jetNM_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltaphi_jetM_jetNM_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_tranverse_massM[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_tracks_jetM_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_tracks_jetNM_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_EMN_jetM_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_EMC_jetM_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_EMtotal_jetM_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_InvM_muon_jet_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_muon_jet_mva_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_muon_jet_relpt_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_muon_jet_sigxy_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_muon_jet_sigz_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_muon_jet_sigr_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_SSOS_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_pT_sum_M[s].Scale(lumiD/lumi_mc[s])
       	hist_mcC_pT_product_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltaR_lep_2jets_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltaphi_MET_2jets_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltaphi_lephad_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_eta_2jets_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_pt_2jets_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_pt_Wlep_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltaR_lephad_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltaphi_lep_2jets_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltaeta_lephad_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltaeta_lep_2jets_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltapt_lephad_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_deltapt_lep_2jets_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_muon_btag_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_notmuon_btag_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_notmuon_deeptagG_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_jet_notmuon_deeptagC_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_pT_proy_M[s].Scale(lumiD/lumi_mc[s])
        hist_mcC_pT_sum_2J_M[s].Scale(lumiD/lumi_mc[s])

######## Special scaling for W+jets samples

for s in samples_mc:
   if s[0]+s[1]+s[2]+s[3]+s[4] == "Wjets":
        hist_mcC_nJetGood_M[s].Scale(CoefWJ_M[years])
        hist_mcC_nMuoninJet_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_muon_pt_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_muon_nmu_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_notmuon_pt_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_muon_eta_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_notmuon_eta_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_notmuon_nmu_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_notmuon_qgl_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_notmuon_mass_M[s].Scale(CoefWJ_M[years])
        hist_mcC_lepton_pt_M[s].Scale(CoefWJ_M[years])
        hist_mcC_lepton_eta_M[s].Scale(CoefWJ_M[years])
        hist_mcC_muon_jet_pt_M[s].Scale(CoefWJ_M[years])
        hist_mcC_muon_jet_eta_M[s].Scale(CoefWJ_M[years])
        hist_mcC_InvM_2jets_M[s].Scale(CoefWJ_M[years])
        hist_mcC_InvM_jetM_lepM[s].Scale(CoefWJ_M[years])
        hist_mcC_MET_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltaR_jetM_lepM[s].Scale(CoefWJ_M[years])
        hist_mcC_deltaR_jetM_jetNM_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltaeta_jetM_jetNM_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltapt_jetM_jetNM_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltaphi_jetM_jetNM_M[s].Scale(CoefWJ_M[years])
        hist_mcC_tranverse_massM[s].Scale(CoefWJ_M[years])
        hist_mcC_tracks_jetM_M[s].Scale(CoefWJ_M[years])
        hist_mcC_tracks_jetNM_M[s].Scale(CoefWJ_M[years])
        hist_mcC_EMN_jetM_M[s].Scale(CoefWJ_M[years])
        hist_mcC_EMC_jetM_M[s].Scale(CoefWJ_M[years])
        hist_mcC_EMtotal_jetM_M[s].Scale(CoefWJ_M[years])
        hist_mcC_InvM_muon_jet_M[s].Scale(CoefWJ_M[years])
        hist_mcC_muon_jet_mva_M[s].Scale(CoefWJ_M[years])
        hist_mcC_muon_jet_relpt_M[s].Scale(CoefWJ_M[years])
        hist_mcC_muon_jet_sigxy_M[s].Scale(CoefWJ_M[years])
        hist_mcC_muon_jet_sigz_M[s].Scale(CoefWJ_M[years])
        hist_mcC_muon_jet_sigr_M[s].Scale(CoefWJ_M[years])
        hist_mcC_SSOS_M[s].Scale(CoefWJ_M[years])
        hist_mcC_pT_sum_M[s].Scale(CoefWJ_M[years])
        hist_mcC_pT_product_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltaR_lep_2jets_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltaphi_MET_2jets_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltaphi_lephad_M[s].Scale(CoefWJ_M[years])
        hist_mcC_eta_2jets_M[s].Scale(CoefWJ_M[years])
        hist_mcC_pt_2jets_M[s].Scale(CoefWJ_M[years])
        hist_mcC_pt_Wlep_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltaR_lephad_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltaphi_lep_2jets_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltaeta_lephad_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltaeta_lep_2jets_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltapt_lephad_M[s].Scale(CoefWJ_M[years])
        hist_mcC_deltapt_lep_2jets_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_muon_btag_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_notmuon_btag_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_notmuon_deeptagG_M[s].Scale(CoefWJ_M[years])
        hist_mcC_jet_notmuon_deeptagC_M[s].Scale(CoefWJ_M[years])
        hist_mcC_pT_proy_M[s].Scale(CoefWJ_M[years])
        hist_mcC_pT_sum_2J_M[s].Scale(CoefWJ_M[years])

######## Joining mc 

h_nJetGood_MCC_M = hist_mcC_nJetGood_M[samples_mc[0]]
h_nMuoninJet_MCC_M = hist_mcC_nMuoninJet_M[samples_mc[0]]
h_jet_muon_pt_MCC_M = hist_mcC_jet_muon_pt_M[samples_mc[0]]
h_jet_muon_nmu_MCC_M = hist_mcC_jet_muon_nmu_M[samples_mc[0]]
h_jet_notmuon_pt_MCC_M = hist_mcC_jet_notmuon_pt_M[samples_mc[0]]
h_jet_muon_eta_MCC_M = hist_mcC_jet_muon_eta_M[samples_mc[0]]
h_jet_notmuon_eta_MCC_M = hist_mcC_jet_notmuon_eta_M[samples_mc[0]]
h_jet_notmuon_nmu_MCC_M = hist_mcC_jet_notmuon_nmu_M[samples_mc[0]]
h_jet_notmuon_mass_MCC_M = hist_mcC_jet_notmuon_mass_M[samples_mc[0]]
h_jet_notmuon_qgl_MCC_M = hist_mcC_jet_notmuon_qgl_M[samples_mc[0]]
h_lepton_pt_MCC_M = hist_mcC_lepton_pt_M[samples_mc[0]]
h_lepton_eta_MCC_M = hist_mcC_lepton_eta_M[samples_mc[0]]
h_muon_jet_pt_MCC_M = hist_mcC_muon_jet_pt_M[samples_mc[0]]
h_muon_jet_eta_MCC_M = hist_mcC_muon_jet_eta_M[samples_mc[0]]
h_InvM_2jets_MCC_M = hist_mcC_InvM_2jets_M[samples_mc[0]]
h_InvM_jetM_lepMCC_M = hist_mcC_InvM_jetM_lepM[samples_mc[0]]
h_MET_MCC_M = hist_mcC_MET_M[samples_mc[0]]
h_deltaR_jetM_lepMCC_M = hist_mcC_deltaR_jetM_lepM[samples_mc[0]]
h_deltaR_jetM_jetNM_MCC_M = hist_mcC_deltaR_jetM_jetNM_M[samples_mc[0]]
h_deltaeta_jetM_jetNM_MCC_M = hist_mcC_deltaeta_jetM_jetNM_M[samples_mc[0]]
h_deltapt_jetM_jetNM_MCC_M = hist_mcC_deltapt_jetM_jetNM_M[samples_mc[0]]
h_deltaphi_jetM_jetNM_MCC_M = hist_mcC_deltaphi_jetM_jetNM_M[samples_mc[0]]
h_tranverse_massMCC_M = hist_mcC_tranverse_massM[samples_mc[0]]
h_tracks_jetM_MCC_M = hist_mcC_tracks_jetM_M[samples_mc[0]]
h_tracks_jetNM_MCC_M = hist_mcC_tracks_jetNM_M[samples_mc[0]]
h_EMN_jetM_MCC_M = hist_mcC_EMN_jetM_M[samples_mc[0]]
h_EMC_jetM_MCC_M = hist_mcC_EMC_jetM_M[samples_mc[0]]
h_EMtotal_jetM_MCC_M = hist_mcC_EMtotal_jetM_M[samples_mc[0]]
h_InvM_muon_jet_MCC_M = hist_mcC_InvM_muon_jet_M[samples_mc[0]]
h_muon_jet_mva_MCC_M = hist_mcC_muon_jet_mva_M[samples_mc[0]]
h_muon_jet_relpt_MCC_M = hist_mcC_muon_jet_relpt_M[samples_mc[0]]
h_muon_jet_sigxy_MCC_M = hist_mcC_muon_jet_sigxy_M[samples_mc[0]]
h_muon_jet_sigz_MCC_M = hist_mcC_muon_jet_sigz_M[samples_mc[0]]
h_muon_jet_sigr_MCC_M = hist_mcC_muon_jet_sigr_M[samples_mc[0]]
h_SSOS_MCC_M = hist_mcC_SSOS_M[samples_mc[0]]
h_pT_sum_MCC_M = hist_mcC_pT_sum_M[samples_mc[0]]
h_pT_product_MCC_M = hist_mcC_pT_product_M[samples_mc[0]]
h_deltaR_lep_2jets_MCC_M = hist_mcC_deltaR_lep_2jets_M[samples_mc[0]]
h_deltaphi_MET_2jets_MCC_M = hist_mcC_deltaphi_MET_2jets_M[samples_mc[0]]
h_deltaphi_lephad_MCC_M = hist_mcC_deltaphi_lephad_M[samples_mc[0]]
h_eta_2jets_MCC_M = hist_mcC_eta_2jets_M[samples_mc[0]]
h_pt_2jets_MCC_M = hist_mcC_pt_2jets_M[samples_mc[0]]
h_pt_Wlep_MCC_M = hist_mcC_pt_Wlep_M[samples_mc[0]]
h_deltaR_lephad_MCC_M = hist_mcC_deltaR_lephad_M[samples_mc[0]]
h_deltaphi_lep_2jets_MCC_M = hist_mcC_deltaphi_lep_2jets_M[samples_mc[0]]
h_deltaeta_lephad_MCC_M = hist_mcC_deltaeta_lephad_M[samples_mc[0]]
h_deltaeta_lep_2jets_MCC_M = hist_mcC_deltaeta_lep_2jets_M[samples_mc[0]]
h_deltapt_lephad_MCC_M = hist_mcC_deltapt_lephad_M[samples_mc[0]]
h_deltapt_lep_2jets_MCC_M = hist_mcC_deltapt_lep_2jets_M[samples_mc[0]]
h_jet_muon_btag_MCC_M = hist_mcC_jet_muon_btag_M[samples_mc[0]]
h_jet_notmuon_btag_MCC_M = hist_mcC_jet_notmuon_btag_M[samples_mc[0]]
h_jet_notmuon_deeptagG_MCC_M = hist_mcC_jet_notmuon_deeptagG_M[samples_mc[0]]
h_jet_notmuon_deeptagC_MCC_M = hist_mcC_jet_notmuon_deeptagC_M[samples_mc[0]]
h_pT_proy_MCC_M = hist_mcC_pT_proy_M[samples_mc[0]]
h_pT_sum_2J_MCC_M = hist_mcC_pT_sum_2J_M[samples_mc[0]]

for s in samples_mc:
  if s != samples_mc[0]:
       h_nJetGood_MCC_M.Add(hist_mcC_nJetGood_M[s].GetPtr())
       h_nMuoninJet_MCC_M.Add(hist_mcC_nMuoninJet_M[s].GetPtr())
       h_jet_muon_pt_MCC_M.Add(hist_mcC_jet_muon_pt_M[s].GetPtr())
       h_jet_muon_nmu_MCC_M.Add(hist_mcC_jet_muon_nmu_M[s].GetPtr())
       h_jet_notmuon_pt_MCC_M.Add(hist_mcC_jet_notmuon_pt_M[s].GetPtr())
       h_jet_muon_eta_MCC_M.Add(hist_mcC_jet_muon_eta_M[s].GetPtr())
       h_jet_notmuon_eta_MCC_M.Add(hist_mcC_jet_notmuon_eta_M[s].GetPtr())
       h_jet_notmuon_nmu_MCC_M.Add(hist_mcC_jet_notmuon_nmu_M[s].GetPtr())
       h_jet_notmuon_qgl_MCC_M.Add(hist_mcC_jet_notmuon_qgl_M[s].GetPtr())
       h_jet_notmuon_mass_MCC_M.Add(hist_mcC_jet_notmuon_mass_M[s].GetPtr())
       h_lepton_pt_MCC_M.Add(hist_mcC_lepton_pt_M[s].GetPtr())
       h_lepton_eta_MCC_M.Add(hist_mcC_lepton_eta_M[s].GetPtr())
       h_muon_jet_pt_MCC_M.Add(hist_mcC_muon_jet_pt_M[s].GetPtr())
       h_muon_jet_eta_MCC_M.Add(hist_mcC_muon_jet_eta_M[s].GetPtr())
       h_InvM_2jets_MCC_M.Add(hist_mcC_InvM_2jets_M[s].GetPtr())
       h_InvM_jetM_lepMCC_M.Add(hist_mcC_InvM_jetM_lepM[s].GetPtr())
       h_MET_MCC_M.Add(hist_mcC_MET_M[s].GetPtr())
       h_deltaR_jetM_lepMCC_M.Add(hist_mcC_deltaR_jetM_lepM[s].GetPtr())
       h_deltaR_jetM_jetNM_MCC_M.Add(hist_mcC_deltaR_jetM_jetNM_M[s].GetPtr())
       h_deltaphi_jetM_jetNM_MCC_M.Add(hist_mcC_deltaphi_jetM_jetNM_M[s].GetPtr())
       h_deltapt_jetM_jetNM_MCC_M.Add(hist_mcC_deltapt_jetM_jetNM_M[s].GetPtr())
       h_deltaeta_jetM_jetNM_MCC_M.Add(hist_mcC_deltaeta_jetM_jetNM_M[s].GetPtr())
       h_tranverse_massMCC_M.Add(hist_mcC_tranverse_massM[s].GetPtr())
       h_tracks_jetM_MCC_M.Add(hist_mcC_tracks_jetM_M[s].GetPtr())
       h_tracks_jetNM_MCC_M.Add(hist_mcC_tracks_jetNM_M[s].GetPtr())
       h_EMN_jetM_MCC_M.Add(hist_mcC_EMN_jetM_M[s].GetPtr())
       h_EMC_jetM_MCC_M.Add(hist_mcC_EMC_jetM_M[s].GetPtr())
       h_EMtotal_jetM_MCC_M.Add(hist_mcC_EMtotal_jetM_M[s].GetPtr())
       h_InvM_muon_jet_MCC_M.Add(hist_mcC_InvM_muon_jet_M[s].GetPtr())
       h_muon_jet_mva_MCC_M.Add(hist_mcC_muon_jet_mva_M[s].GetPtr())
       h_muon_jet_relpt_MCC_M.Add(hist_mcC_muon_jet_relpt_M[s].GetPtr())
       h_muon_jet_sigxy_MCC_M.Add(hist_mcC_muon_jet_sigxy_M[s].GetPtr())
       h_muon_jet_sigz_MCC_M.Add(hist_mcC_muon_jet_sigz_M[s].GetPtr())
       h_muon_jet_sigr_MCC_M.Add(hist_mcC_muon_jet_sigr_M[s].GetPtr())
       h_SSOS_MCC_M.Add(hist_mcC_SSOS_M[s].GetPtr())
       h_pT_sum_MCC_M.Add(hist_mcC_pT_sum_M[s].GetPtr())
       h_pT_product_MCC_M.Add(hist_mcC_pT_product_M[s].GetPtr())
       h_deltaR_lep_2jets_MCC_M.Add(hist_mcC_deltaR_lep_2jets_M[s].GetPtr())
       h_deltaphi_MET_2jets_MCC_M.Add(hist_mcC_deltaphi_MET_2jets_M[s].GetPtr())
       h_deltaphi_lephad_MCC_M.Add(hist_mcC_deltaphi_lephad_M[s].GetPtr())
       h_eta_2jets_MCC_M.Add(hist_mcC_eta_2jets_M[s].GetPtr())
       h_pt_2jets_MCC_M.Add(hist_mcC_pt_2jets_M[s].GetPtr())
       h_pt_Wlep_MCC_M.Add(hist_mcC_pt_Wlep_M[s].GetPtr())
       h_jet_muon_btag_MCC_M.Add(hist_mcC_jet_muon_btag_M[s].GetPtr())
       h_jet_notmuon_btag_MCC_M.Add(hist_mcC_jet_notmuon_btag_M[s].GetPtr())
       h_jet_notmuon_deeptagG_MCC_M.Add(hist_mcC_jet_notmuon_deeptagG_M[s].GetPtr())
       h_jet_notmuon_deeptagC_MCC_M.Add(hist_mcC_jet_notmuon_deeptagC_M[s].GetPtr())
       h_deltaR_lephad_MCC_M.Add(hist_mcC_deltaR_lephad_M[s].GetPtr())
       h_deltaphi_lep_2jets_MCC_M.Add(hist_mcC_deltaphi_lep_2jets_M[s].GetPtr())
       h_deltaeta_lephad_MCC_M.Add(hist_mcC_deltaeta_lephad_M[s].GetPtr())
       h_deltaeta_lep_2jets_MCC_M.Add(hist_mcC_deltaeta_lep_2jets_M[s].GetPtr())
       h_deltapt_lephad_MCC_M.Add(hist_mcC_deltapt_lephad_M[s].GetPtr())
       h_deltapt_lep_2jets_MCC_M.Add(hist_mcC_deltapt_lep_2jets_M[s].GetPtr())
       h_pT_proy_MCC_M.Add(hist_mcC_pT_proy_M[s].GetPtr())
       h_pT_sum_2J_MCC_M.Add(hist_mcC_pT_sum_2J_M[s].GetPtr())

hist_dataC_nJetGood_M.Add(h_nJetGood_MCC_M.GetPtr(),-1)
hist_dataC_nMuoninJet_M.Add(h_nMuoninJet_MCC_M.GetPtr(),-1)
hist_dataC_jet_muon_pt_M.Add(h_jet_muon_pt_MCC_M.GetPtr(),-1)
hist_dataC_jet_muon_nmu_M.Add(h_jet_muon_nmu_MCC_M.GetPtr(),-1)
hist_dataC_jet_notmuon_pt_M.Add(h_jet_notmuon_pt_MCC_M.GetPtr(),-1)
hist_dataC_jet_muon_eta_M.Add(h_jet_muon_eta_MCC_M.GetPtr(),-1)
hist_dataC_jet_notmuon_eta_M.Add(h_jet_notmuon_eta_MCC_M.GetPtr(),-1)
hist_dataC_jet_notmuon_nmu_M.Add(h_jet_notmuon_nmu_MCC_M.GetPtr(),-1)
hist_dataC_jet_notmuon_qgl_M.Add(h_jet_notmuon_qgl_MCC_M.GetPtr(),-1)
hist_dataC_jet_notmuon_mass_M.Add(h_jet_notmuon_mass_MCC_M.GetPtr(),-1)
hist_dataC_lepton_pt_M.Add(h_lepton_pt_MCC_M.GetPtr(),-1)
hist_dataC_lepton_eta_M.Add(h_lepton_eta_MCC_M.GetPtr(),-1)
hist_dataC_muon_jet_pt_M.Add(h_muon_jet_pt_MCC_M.GetPtr(),-1)
hist_dataC_muon_jet_eta_M.Add(h_muon_jet_eta_MCC_M.GetPtr(),-1)
hist_dataC_InvM_2jets_M.Add(h_InvM_2jets_MCC_M.GetPtr(),-1)
hist_dataC_InvM_jetM_lepM.Add(h_InvM_jetM_lepMCC_M.GetPtr(),-1)
hist_dataC_InvM_muon_jet_M.Add(h_InvM_muon_jet_MCC_M.GetPtr(),-1)
hist_dataC_MET_M.Add(h_MET_MCC_M.GetPtr(),-1)
hist_dataC_deltaR_jetM_lepM.Add(h_deltaR_jetM_lepMCC_M.GetPtr(),-1)
hist_dataC_deltaR_jetM_jetNM_M.Add(h_deltaR_jetM_jetNM_MCC_M.GetPtr(),-1)
hist_dataC_deltapt_jetM_jetNM_M.Add(h_deltapt_jetM_jetNM_MCC_M.GetPtr(),-1)
hist_dataC_deltaeta_jetM_jetNM_M.Add(h_deltaeta_jetM_jetNM_MCC_M.GetPtr(),-1)
hist_dataC_deltaphi_jetM_jetNM_M.Add(h_deltaphi_jetM_jetNM_MCC_M.GetPtr(),-1)
hist_dataC_tranverse_massM.Add(h_tranverse_massMCC_M.GetPtr(),-1)
hist_dataC_tracks_jetM_M.Add(h_tracks_jetM_MCC_M.GetPtr(),-1)
hist_dataC_tracks_jetNM_M.Add(h_tracks_jetNM_MCC_M.GetPtr(),-1)
hist_dataC_EMN_jetM_M.Add(h_EMN_jetM_MCC_M.GetPtr(),-1)
hist_dataC_EMC_jetM_M.Add(h_EMC_jetM_MCC_M.GetPtr(),-1)
hist_dataC_EMtotal_jetM_M.Add(h_EMtotal_jetM_MCC_M.GetPtr(),-1)
hist_dataC_muon_jet_mva_M.Add(h_muon_jet_mva_MCC_M.GetPtr(),-1)
hist_dataC_muon_jet_relpt_M.Add(h_muon_jet_relpt_MCC_M.GetPtr(),-1)
hist_dataC_muon_jet_sigxy_M.Add(h_muon_jet_sigxy_MCC_M.GetPtr(),-1)
hist_dataC_muon_jet_sigz_M.Add(h_muon_jet_sigz_MCC_M.GetPtr(),-1)
hist_dataC_muon_jet_sigr_M.Add(h_muon_jet_sigr_MCC_M.GetPtr(),-1)
hist_dataC_SSOS_M.Add(h_SSOS_MCC_M.GetPtr(),-1)
hist_dataC_pT_sum_M.Add(h_pT_sum_MCC_M.GetPtr(),-1)
hist_dataC_pT_product_M.Add(h_pT_product_MCC_M.GetPtr(),-1)
hist_dataC_deltaR_lep_2jets_M.Add(h_deltaR_lep_2jets_MCC_M.GetPtr(),-1)
hist_dataC_deltaphi_MET_2jets_M.Add(h_deltaphi_MET_2jets_MCC_M.GetPtr(),-1)
hist_dataC_deltaphi_lephad_M.Add(h_deltaphi_lephad_MCC_M.GetPtr(),-1)
hist_dataC_eta_2jets_M.Add(h_eta_2jets_MCC_M.GetPtr(),-1)
hist_dataC_pt_2jets_M.Add(h_pt_2jets_MCC_M.GetPtr(),-1)
hist_dataC_pt_Wlep_M.Add(h_pt_Wlep_MCC_M.GetPtr(),-1)
hist_dataC_jet_muon_btag_M.Add(h_jet_muon_btag_MCC_M.GetPtr(),-1)
hist_dataC_jet_notmuon_btag_M.Add(h_jet_notmuon_btag_MCC_M.GetPtr(),-1)
hist_dataC_jet_notmuon_deeptagG_M.Add(h_jet_notmuon_deeptagG_MCC_M.GetPtr(),-1)
hist_dataC_jet_notmuon_deeptagC_M.Add(h_jet_notmuon_deeptagC_MCC_M.GetPtr(),-1)
hist_dataC_deltaR_lephad_M.Add(h_deltaR_lephad_MCC_M.GetPtr(),-1)
hist_dataC_deltaphi_lep_2jets_M.Add(h_deltaphi_lep_2jets_MCC_M.GetPtr(),-1)
hist_dataC_deltaeta_lephad_M.Add(h_deltaeta_lephad_MCC_M.GetPtr(),-1)
hist_dataC_deltaeta_lep_2jets_M.Add(h_deltaeta_lep_2jets_MCC_M.GetPtr(),-1)
hist_dataC_deltapt_lephad_M.Add(h_deltapt_lephad_MCC_M.GetPtr(),-1)
hist_dataC_deltapt_lep_2jets_M.Add(h_deltapt_lep_2jets_MCC_M.GetPtr(),-1)
hist_dataC_pT_proy_M.Add(h_pT_proy_MCC_M.GetPtr(),-1)
hist_dataC_pT_sum_2J_M.Add(h_pT_sum_2J_MCC_M.GetPtr(),-1)

hist_dataC_nJetGood_M.Scale(f_BD[years])
hist_dataC_nMuoninJet_M.Scale(f_BD[years])
hist_dataC_jet_muon_pt_M.Scale(f_BD[years])
hist_dataC_jet_muon_nmu_M.Scale(f_BD[years])
hist_dataC_jet_notmuon_pt_M.Scale(f_BD[years])
hist_dataC_jet_muon_eta_M.Scale(f_BD[years])
hist_dataC_jet_notmuon_eta_M.Scale(f_BD[years])
hist_dataC_jet_notmuon_mass_M.Scale(f_BD[years])
hist_dataC_jet_notmuon_qgl_M.Scale(f_BD[years])
hist_dataC_jet_notmuon_nmu_M.Scale(f_BD[years])
hist_dataC_lepton_pt_M.Scale(f_BD[years])
hist_dataC_lepton_eta_M.Scale(f_BD[years])
hist_dataC_muon_jet_pt_M.Scale(f_BD[years])
hist_dataC_muon_jet_eta_M.Scale(f_BD[years])
hist_dataC_InvM_2jets_M.Scale(f_BD[years])
hist_dataC_InvM_jetM_lepM.Scale(f_BD[years])
hist_dataC_InvM_muon_jet_M.Scale(f_BD[years])
hist_dataC_MET_M.Scale(f_BD[years])
hist_dataC_deltaR_jetM_lepM.Scale(f_BD[years])
hist_dataC_deltaR_jetM_jetNM_M.Scale(f_BD[years])
hist_dataC_deltaeta_jetM_jetNM_M.Scale(f_BD[years])
hist_dataC_deltapt_jetM_jetNM_M.Scale(f_BD[years])
hist_dataC_deltaphi_jetM_jetNM_M.Scale(f_BD[years])
hist_dataC_tranverse_massM.Scale(f_BD[years])
hist_dataC_tracks_jetM_M.Scale(f_BD[years])
hist_dataC_tracks_jetNM_M.Scale(f_BD[years])
hist_dataC_EMN_jetM_M.Scale(f_BD[years])
hist_dataC_EMC_jetM_M.Scale(f_BD[years])
hist_dataC_EMtotal_jetM_M.Scale(f_BD[years])
hist_dataC_muon_jet_mva_M.Scale(f_BD[years])
hist_dataC_muon_jet_relpt_M.Scale(f_BD[years])
hist_dataC_muon_jet_sigxy_M.Scale(f_BD[years])
hist_dataC_muon_jet_sigz_M.Scale(f_BD[years])
hist_dataC_muon_jet_sigr_M.Scale(f_BD[years])
hist_dataC_SSOS_M.Scale(f_BD[years])
hist_dataC_pT_sum_M.Scale(f_BD[years])
hist_dataC_pT_product_M.Scale(f_BD[years])
hist_dataC_deltaR_lep_2jets_M.Scale(f_BD[years])
hist_dataC_deltaphi_MET_2jets_M.Scale(f_BD[years])
hist_dataC_deltaphi_lephad_M.Scale(f_BD[years])
hist_dataC_eta_2jets_M.Scale(f_BD[years])
hist_dataC_pt_2jets_M.Scale(f_BD[years])
hist_dataC_pt_Wlep_M.Scale(f_BD[years])
hist_dataC_jet_muon_btag_M.Scale(f_BD[years])
hist_dataC_jet_notmuon_btag_M.Scale(f_BD[years])
hist_dataC_jet_notmuon_deeptagG_M.Scale(f_BD[years])
hist_dataC_jet_notmuon_deeptagC_M.Scale(f_BD[years])
hist_dataC_deltaR_lephad_M.Scale(f_BD[years])
hist_dataC_deltaphi_lep_2jets_M.Scale(f_BD[years])
hist_dataC_deltaeta_lephad_M.Scale(f_BD[years])
hist_dataC_deltaeta_lep_2jets_M.Scale(f_BD[years])
hist_dataC_deltapt_lephad_M.Scale(f_BD[years])
hist_dataC_deltapt_lep_2jets_M.Scale(f_BD[years])
hist_dataC_pT_proy_M.Scale(f_BD[years])
hist_dataC_pT_sum_2J_M.Scale(f_BD[years])

if args.ssos: path_hist = '/nfs/cms/vazqueze/analisisWW/hists/ssos/hists_v1v2v3v4SSOS_QCD'+years+'.root'
else: path_hist = '/nfs/cms/vazqueze/analisisWW/hists/hists_v1v2v3v4_QCD'+years+'.root'
myfile = TFile( path_hist, 'RECREATE' )

hist_dataC_nJetGood_M.Write()
hist_dataC_nMuoninJet_M.Write()
hist_dataC_jet_muon_pt_M.Write()
hist_dataC_jet_muon_nmu_M.Write()
hist_dataC_jet_notmuon_pt_M.Write()
hist_dataC_jet_muon_eta_M.Write()
hist_dataC_jet_notmuon_eta_M.Write()
hist_dataC_jet_notmuon_qgl_M.Write()
hist_dataC_jet_notmuon_nmu_M.Write()
hist_dataC_jet_notmuon_mass_M.Write()
hist_dataC_lepton_pt_M.Write()
hist_dataC_lepton_eta_M.Write()
hist_dataC_muon_jet_pt_M.Write()
hist_dataC_muon_jet_eta_M.Write()
hist_dataC_InvM_2jets_M.Write()
hist_dataC_InvM_jetM_lepM.Write()
hist_dataC_InvM_muon_jet_M.Write()
hist_dataC_MET_M.Write()
hist_dataC_deltaR_jetM_lepM.Write()
hist_dataC_deltaR_jetM_jetNM_M.Write()
hist_dataC_deltapt_jetM_jetNM_M.Write()
hist_dataC_deltaeta_jetM_jetNM_M.Write()
hist_dataC_deltaphi_jetM_jetNM_M.Write()
hist_dataC_tranverse_massM.Write()
hist_dataC_tracks_jetM_M.Write()
hist_dataC_tracks_jetNM_M.Write()
hist_dataC_EMN_jetM_M.Write()
hist_dataC_EMC_jetM_M.Write()
hist_dataC_EMtotal_jetM_M.Write()
hist_dataC_muon_jet_mva_M.Write()
hist_dataC_muon_jet_relpt_M.Write()
hist_dataC_muon_jet_sigxy_M.Write()
hist_dataC_muon_jet_sigz_M.Write()
hist_dataC_muon_jet_sigr_M.Write()
hist_dataC_SSOS_M.Write()
hist_dataC_pT_sum_M.Write()
hist_dataC_pT_product_M.Write()
hist_dataC_deltaR_lep_2jets_M.Write()                          
hist_dataC_deltaphi_MET_2jets_M.Write()                            
hist_dataC_deltaphi_lephad_M.Write()                         
hist_dataC_eta_2jets_M.Write()                   
hist_dataC_pt_2jets_M.Write()
hist_dataC_pt_Wlep_M.Write()                 
hist_dataC_jet_muon_btag_M.Write()                       
hist_dataC_jet_notmuon_btag_M.Write()
hist_dataC_jet_notmuon_deeptagG_M.Write()
hist_dataC_jet_notmuon_deeptagC_M.Write()
hist_dataC_deltaR_lephad_M.Write()
hist_dataC_deltaphi_lep_2jets_M.Write()
hist_dataC_deltaeta_lephad_M.Write()
hist_dataC_deltaeta_lep_2jets_M.Write()
hist_dataC_deltapt_lephad_M.Write()
hist_dataC_deltapt_lep_2jets_M.Write()
hist_dataC_pT_proy_M.Write()
hist_dataC_pT_sum_2J_M.Write()

myfile.Close()

if args.ssos: print('SSOS version')

print('Ended successfully')
