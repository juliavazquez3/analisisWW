######################################                                        
######################################                                        
#######      GOOD  VERSION     #######
######################################                                     
######################################                                     

print('SL CHANNEL')

## Modified version, meant to be run locally, not sent to condor 

## Selection for noth MC and data samples, created to distinguish between years
## Different histogram files are produced for each situation 
## It includes an option for SSOS substraction

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
parser.add_argument("--process", type=string, default="WW",
                    help="Select type of process to run")
parser.add_argument("--year", type=string, default="2016",
                    help="Select year of process to run")
parser.add_argument("--type", type=string, default="data",
                    help="Selec type of data to run")
parser.add_argument("--notfull", action="store_true", default=False,
                    help="Use just a range of the sample")
parser.add_argument("--ssos", action="store_true", default=False,
                    help="Perform ssos substraction")
parser.add_argument('-l','--list', nargs='+', help='range of sample to use')
# Use like:
# python arg.py -l 1234 2345 3456 4567

args = parser.parse_args()

if args.process == "allMC": proc = ["WW","Wjets1","Wjets2","Wjets3","Wjets4","Wjets5","Wjets6","Wjets7","Wjets8","DY1","DY2","DY3","DY4","DY5","DY6","DY7","DY8",
        "ttbar","ttbarlep","ttbarhad","ZZ","WZ","ST1","ST2","ST3","ST4","Wjets0J","Wjets1J","Wjets2J","DY0J","DY1J","DY2J"]
elif (args.process == "WW" or args.process == "Wjets" or args.process == "ttbar" or args.process == "DY" or args.process == "WZ"or args.process == "ZZ" or args.process == "ST1"
        or args.process == "ST2" or args.process == "ST3" or args.process == "ST4"  or args.process == "DYjets"  or args.process == "ttbarlep"
        or args.process == "ttbarhad" or args.process == "M" or args.process == "E"): proc = [str(args.process)]
else: raise NameError('Incorrect process name')

if args.year == "all": years = ["2016","2016B","2017","2018"]
elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): years = [str(args.year)]
else: raise NameError('Incorrect year')

if args.type == "data": mode = "data"
elif args.type == "mc": mode = "mc"
else: raise NameError('Incorrect type')

#proc = ["Wjets0J","Wjets1J","Wjets2J"]

special_weights = ["Wjets0J2016","Wjets1J2016","Wjets2J2016","DY0J2016","DY1J2016","DY2J2016",
	"Wjets0J2016B","Wjets1J2016B","Wjets2J2016B","DY0J2016B","DY1J2016B","DY2J2016B",
	"Wjets0J2017","Wjets1J2017","Wjets2J2017","DY0J2017","DY1J2017","DY2J2017",
	"Wjets0J2018","Wjets1J2018","Wjets2J2018","DY0J2018","DY1J2018","DY2J2018"]

samples = []

for p in proc:
  for y in years:
    if mode == "mc": samples.append(p+y)
    else: samples.append(y+p)

samples_test = samples[:]

print("the samples treated are",samples)

if args.notfull:
	if len(args.list) != 2: raise NameError('List has to have 2 elements')
	print("the range of files is",args.list)
 

# Create a ROOT dataframe for each dataset
# Note that we load the filenames from the external json file placed in the same folder than this script.
# Example "python analisisWW/selection_v2.py --process="WW" --notfull -l 0 50"

if mode == "mc": files = json.load(open("/nfs/cms/vazqueze/analisisWW/mc_info_v9"+years[0]+".json"))
else: files = json.load(open("/nfs/cms/vazqueze/analisisWW/data_info_v9.json"))

processes = files.keys()

df = {}
xsecs = {}
sumws = {}
archives = {}
event_test = {}

for s in samples:
  archives[s]=[]

for p in processes:
    # Construct the dataframes
    folder = files[p]["folder_name"] # Folder name
    filename = files[p]["sample"] # Sample name
    if mode == "mc": num_events = files[p]["events"] # Number of events
    else: num_events = files[p]["events_total"] # Number of events
    num_files = files[p]["files"] # Number of files
    list_files = [f for f in listdir(folder) if isfile(join(folder, f))] # Lista de archivos
    #print(len(list_files))
    for s in samples:
      if files[p]["type"]==s:
          event_test[s] = num_events
          if (num_files == len(list_files)):
              for f in list_files:
                  archives[s].append(join(folder,f))

for s in samples:
  if args.notfull: archives[s]=archives[s][int(args.list[0]):int(args.list[1])]
  df[s] = ROOT.RDataFrame("Events",set(archives[s]))
  print("Number of files for",s,len(archives[s]))

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

######## Gen identification for W plus jets

gInterpreter.Declare("""
      #include <bitset>
      #include <string>
      #include <iostream>
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using Vbool = const ROOT::RVec<bool>&;
      struct ASCII
      {
                std::string toBinary(int n)
                {
                        std::string r;
                        while(n!=0) {r=(n%2==0 ?"0":"1")+r; n/=2;}
                        return r;
                }

      };
      auto testSF(UInt_t part_ind, UInt_t status, UInt_t n) {
            char ind;
            bool hardP = false;
            std::bitset<15> b(status);
            auto statusflags_string = b.to_string();
            for(unsigned int j=statusflags_string.length()-1; j>=statusflags_string.length()-(n+1); j--)
            {
                   //std::cout << "statusflags bit " << j << " " << statusflags_string[j] <<std::endl;
                   ind = statusflags_string.at(j);
            }
            if(ind=='1') hardP = true;
            return hardP;
      };
      auto vectorHP(UInt_t nPart, Vint status, Vint pdg, UInt_t n) {
            vector<bool> vb;
            for (unsigned int i=0; i<nPart; ++i) {
                vb.push_back(testSF(i,status[i],n));
            }
            return vb;
      };
      auto wpluscbool(UInt_t nPart, Vint status, Vint pdg, Vbool hardP, Vbool firstC) {
            int typeC = 0;
            int indC = 0;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==5 && hardP[i]) {
                          typeC = 3;
                } else {
                     if (fabs(pdg[i])==4 && hardP[i] && firstC[i]) indC++; 
                }
            }
            if (typeC!=3 && (indC % 2 != 0)) typeC = 2;
            if (typeC!=3 && (indC % 2 == 0) && (indC > 0)) typeC = 1; 
            return typeC;
      };
      auto ttbarcharm(UInt_t nPart, Vint status, Vint pdg, Vint mother) {
            bool typeC = false;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==4 && fabs(pdg[mother[i]])==24) {
                          typeC = true;
                }
            }
            return typeC;
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
		if (pt[i]>cutpt && fabs(eta[i])<2.4 && iso[i]<0.15 && tID[i]){
                	vb.push_back(i);
		}
            }
            return vb;
      };
      auto elInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vint cutB, Vbool mva80 , Vbool mva90, Float_t cutpt) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
                if (pt[i]>cutpt && fabs(eta[i])<2.5 && iso[i]<0.15 && mva80[i]){
                	vb.push_back(i);
                }
            }
            return vb;
      };
""")

## Numero de muones dentro de un jet
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
            int ind=-1;
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
            int ind=-1;
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
                //if (mu_pt[mu_jet[0]]<25. && fabs(mu_eta[mu_jet[0]])<2.4) {
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
            float deltaRlepjet1;
       	    float deltaRlepjet2;
       	    float deltaPhilepjet1;
            float deltaPhilepjet2;
            float deltaEtalepjet1;
            float deltaEtalepjet2;
            if(mu_good.size()>0){
                  plep1.SetPt(mu_pt[mu_good[0]]); plep1.SetEta(mu_eta[mu_good[0]]); plep1.SetPhi(mu_phi[mu_good[0]]); plep1.SetM(mu_mass[mu_good[0]]);
                  deltaRlepjet1 = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_muon[0]],mu_phi[mu_good[0]],jet_phi[jet_muon[0]]);
       	       	  deltaRlepjet2	= ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],jet_eta[jet_notmuon],mu_phi[mu_good[0]],jet_phi[jet_notmuon]);
                  deltaPhilepjet1 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_muon[0]]);
       	       	  deltaPhilepjet2 = fabs(mu_phi[mu_good[0]]-jet_phi[jet_notmuon]);
       	       	  deltaEtalepjet1 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(mu_eta[mu_good[0]]-jet_eta[jet_notmuon]);
            } else if(el_good.size()>0){
                  plep1.SetPt(el_pt[el_good[0]]); plep1.SetEta(el_eta[el_good[0]]); plep1.SetPhi(el_phi[el_good[0]]); plep1.SetM(el_mass[el_good[0]]);
       	       	  deltaRlepjet1	= ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_muon[0]],el_phi[el_good[0]],jet_phi[jet_muon[0]]);
                  deltaRlepjet2 = ROOT::VecOps::DeltaR(el_eta[el_good[0]],jet_eta[jet_notmuon],el_phi[el_good[0]],jet_phi[jet_notmuon]);
       	       	  deltaPhilepjet1 = fabs(el_phi[el_good[0]]-jet_phi[jet_muon[0]]);
                  deltaPhilepjet2 = fabs(el_phi[el_good[0]]-jet_phi[jet_notmuon]);
                  deltaEtalepjet1 = fabs(el_eta[el_good[0]]-jet_eta[jet_muon[0]]);
                  deltaEtalepjet2 = fabs(el_eta[el_good[0]]-jet_eta[jet_notmuon]);
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
            vb.push_back(deltaRlepjet1);
            vb.push_back(deltaRlepjet2);
            vb.push_back(deltaPhilepjet1);
            vb.push_back(deltaPhilepjet2);
            vb.push_back(deltaEtalepjet1);
            vb.push_back(deltaEtalepjet2);
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

df_muon = {}
df_electron = {}
df_test = {}

for s in samples:
      if s in special_weights:
             df[s] = df[s].Define('weight_aux','fabs(genWeight) > 0 ? genWeight/fabs(genWeight) : 0')
      else:
             df[s] = df[s].Define('weight_aux','1')

for s in samples:
	df[s] = df[s].Define('MuonGoodInd','muonInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all, Muon_tightId,'+str(muon_pt[years[0]])+')')
	df[s] = df[s].Define('ElectronGoodInd','elInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_cutBased, Electron_mvaFall17V2noIso_WP80, Electron_mvaFall17V2noIso_WP90,'+str(el_pt[years[0]])+')')
	df[s] = df[s].Define('nMuonGood','MuonGoodInd.size()')
	df[s] = df[s].Define('nElectronGood','ElectronGoodInd.size()')
	df[s] = df[s].Define('JetGoodInd','JetInd(nJet, Jet_pt, Jet_eta, Jet_phi, MuonGoodInd, Muon_eta, Muon_phi, ElectronGoodInd, Electron_eta, Electron_phi, Jet_puId, Jet_jetId)')
	df[s] = df[s].Define('nJetGood','JetGoodInd.size()')
	df[s] = df[s].Define('MuonJetInd','JetMuonInd(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx)')
	df[s] = df[s].Define('JetMuonInd','JetMuonIndJet(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx)')
	df[s] = df[s].Define('MuonJetGood','muonCond( MuonJetInd, Muon_pt, Muon_eta, Muon_pfRelIso04_all,Muon_softId)')
	df[s] = df[s].Define('MuonLepSign','SSOS(MuonGoodInd, ElectronGoodInd, MuonJetInd, Electron_charge, Muon_charge)')
	if args.ssos: df[s] = df[s].Define('weightSSOS','-1*MuonLepSign*weight_aux')
	else: df[s] = df[s].Define('weightSSOS','weight_aux')
	df_test[s] = df[s]

#################     FILTERS     #######################

##### New cuts, compared to verison 0 and 1
##### Exactly 2 jets, not more or less
##### ETA restrictions: 2.4 for jets and muons and 2.5 for electrons

###########    Extra definitions

for s in samples:
	df[s] = df[s].Filter('nMuonGood<2 && nElectronGood<2').Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood==2').Filter('MuonJetInd.size() >= 1').Filter('MuonJetGood')
	df[s] = df[s].Define('JetnotMuonInd','nJetGood>1 ? (JetMuonInd[0] == JetGoodInd[0] ? JetGoodInd[1] : JetGoodInd[0]) : -1')
	### hists definitions
	df[s] = df[s].Define('jet_muon_pt','Jet_pt[JetMuonInd[0]]')
	df[s] = df[s].Define('jet_muon_nmu','Jet_nMuons[JetMuonInd[0]]')
	df[s] = df[s].Define('jet_notmuon_pt','nJetGood>1 ? Jet_pt[JetnotMuonInd] : 0')
	df[s] = df[s].Define('jet_muon_eta','Jet_eta[JetMuonInd[0]]')
	df[s] = df[s].Define('jet_notmuon_eta','nJetGood>1 ? Jet_eta[JetnotMuonInd] : 0')
	df[s] = df[s].Define('jet_notmuon_mass','nJetGood>1 ? Jet_mass[JetnotMuonInd] : 0')
	df[s] = df[s].Define('jet_notmuon_qgl','nJetGood>1 ? Jet_qgl[JetnotMuonInd] : 0')
	df[s] = df[s].Define('jet_notmuon_nmu','nJetGood>1 ? Jet_nMuons[JetnotMuonInd] : 0')
	df[s] = df[s].Define('muon_jet_pt','Muon_pt[MuonJetInd[0]]')
	df[s] = df[s].Define('muon_jet_eta','Muon_eta[MuonJetInd[0]]')
	df[s] = df[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Jet_pt[JetnotMuonInd],Jet_eta[JetnotMuonInd],Jet_phi[JetnotMuonInd],Jet_mass[JetnotMuonInd]) : 0')
	df[s] = df[s].Define('deltaR_jetM_jetNM','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotMuonInd], Jet_eta[JetMuonInd[0]] , Jet_phi[JetnotMuonInd], Jet_phi[JetMuonInd[0]])  : 10')
	df[s] = df[s].Define('deltaphi_jetM_jetNM','fabs(Jet_phi[JetMuonInd[0]]-Jet_phi[JetnotMuonInd])')
	df[s] = df[s].Define('deltaeta_jetM_jetNM','fabs(Jet_eta[JetMuonInd[0]]-Jet_eta[JetnotMuonInd])')
	df[s] = df[s].Define('deltapt_jetM_jetNM','fabs(Jet_pt[JetMuonInd[0]]-Jet_pt[JetnotMuonInd])')
	df[s] = df[s].Define('tracks_jetM','Jet_nConstituents[JetMuonInd[0]]')
	df[s] = df[s].Define('tracks_jetNM','nJetGood>1 ? Jet_nConstituents[JetnotMuonInd] : 0')
	df[s] = df[s].Define('EMN_jetM','Jet_neEmEF[JetMuonInd[0]]')
	df[s] = df[s].Define('EMC_jetM','Jet_chEmEF[JetMuonInd[0]]')
	df[s] = df[s].Define('EMtotal_jetM','Jet_chEmEF[JetMuonInd[0]]+Jet_neEmEF[JetMuonInd[0]]')
	df[s] = df[s].Define('muon_jet_mva','Muon_softMva[MuonJetInd[0]]')
	df[s] = df[s].Define('muon_jet_tight','Muon_tightId[MuonJetInd[0]]')
	df[s] = df[s].Define('muon_jet_relpt','Muon_pt[MuonJetInd[0]]/Jet_pt[JetMuonInd[0]]')
	df[s] = df[s].Define('MuoninJetAux','muoninjet(nMuon, Muon_jetIdx, MuonGoodInd)')
	df[s] = df[s].Define('nMuoninJet','MuoninJetAux.size()')
	df[s] = df[s].Define('muon_jet_sigxy','Muon_dxyErr[MuonJetInd[0]]>0 ? Muon_dxy[MuonJetInd[0]]/Muon_dxyErr[MuonJetInd[0]] : -10')
	df[s] = df[s].Define('muon_jet_sigz','Muon_dzErr[MuonJetInd[0]]>0 ? Muon_dz[MuonJetInd[0]]/Muon_dzErr[MuonJetInd[0]] : -10')
	df[s] = df[s].Define('muon_jet_r','pow(pow(Muon_dz[MuonJetInd[0]],2) + pow(Muon_dxy[MuonJetInd[0]],2),0.5)')
	df[s] = df[s].Define('muon_jet_Err','muon_jet_r > 0 ? pow(pow(Muon_dz[MuonJetInd[0]]*Muon_dzErr[MuonJetInd[0]],2)+pow(Muon_dxy[MuonJetInd[0]]*Muon_dxyErr[MuonJetInd[0]],2),0.5)/muon_jet_r : -10')
	df[s] = df[s].Define('muon_jet_sigr','(muon_jet_Err>0) ? muon_jet_r/muon_jet_Err : -10')
	df[s] = df[s].Define('pT_sum','pTsum(MuonGoodInd, ElectronGoodInd, JetMuonInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
	df[s] = df[s].Define('pT_product','pTprod(MuonGoodInd, ElectronGoodInd, JetMuonInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
	df[s] = df[s].Define('aux_various','variousSUM(MuonGoodInd, ElectronGoodInd, JetMuonInd, JetnotMuonInd, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Electron_pt, Electron_eta, Electron_phi, Electron_mass, MET_pt, MET_phi, Jet_pt, Jet_eta, Jet_phi, Jet_mass)')
	df[s] = df[s].Define('deltaR_lep_2jets','aux_various[0]')
	df[s] = df[s].Define('deltaphi_MET_2jets','aux_various[1]')
	df[s] = df[s].Define('eta_2jets','aux_various[2]')
	df[s] = df[s].Define('pt_2jets','aux_various[3]')
	df[s] = df[s].Define('deltaphi_lephad','aux_various[4]')
	df[s] = df[s].Define('deltaR_lephad','aux_various[5]')
	df[s] = df[s].Define('deltaphi_lep_2jets','aux_various[6]')
	df[s] = df[s].Define('deltaeta_lephad','aux_various[7]')
	df[s] = df[s].Define('deltaeta_lep_2jets','aux_various[8]')
	df[s] = df[s].Define('deltapt_lephad','aux_various[9]')
	df[s] = df[s].Define('deltapt_lep_2jets','aux_various[10]')
	df[s] = df[s].Define('pT_proy','aux_various[11]')
	df[s] = df[s].Define('pT_sum_2J','aux_various[12]')
	df[s] = df[s].Define('pT_Wlep','aux_various[13]')
	df[s] = df[s].Define('deltaR_lep_jet1','aux_various[14]')
	df[s] = df[s].Define('deltaR_lep_jet2','aux_various[15]')
	df[s] = df[s].Define('deltaPhi_lep_jet1','aux_various[16]')
	df[s] = df[s].Define('deltaPhi_lep_jet2','aux_various[17]')
	df[s] = df[s].Define('deltaEta_lep_jet1','aux_various[18]')
	df[s] = df[s].Define('deltaEta_lep_jet2','aux_various[19]')
	df[s] = df[s].Define('jet_muon_btag','Jet_btagDeepFlavB[JetMuonInd[0]]')
	df[s] = df[s].Define('jet_notmuon_btag','Jet_btagDeepFlavB[JetnotMuonInd]')
	df[s] = df[s].Define('jet_notmuon_deeptagG','Jet_btagDeepFlavQG[JetnotMuonInd]')
	df[s] = df[s].Define('jet_notmuon_deeptagC','Jet_btagDeepFlavCvL[JetnotMuonInd]')
	df[s] = df[s].Define('nLooseLepton','nMuon+nElectron-1')

#################################
#######    GEN FILTER    ########
#################################

######### WW sectioning in events of interest (semi charm) and not 

if mode == "mc":
	for s in samples:
		if (s[0]+s[1] == "WW" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
			df[s] = df[s].Define('typeWW','typeWW(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother)')
			df[s] = df[s].Define('typeC','typeC(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother)')
			df[s+"_hadronic"] = df[s].Filter('typeWW == 1')
			df[s+"_leptonic"] = df[s].Filter('typeWW == 3')
			df[s+"_semi_charm"] = df[s].Filter('typeWW == 2 && typeC == 1')
			df[s+"_semi_nocharm"] = df[s].Filter('typeWW == 2 && typeC != 1')
			## Samples correction
			samples.append(s+"_hadronic")
			samples.append(s+"_leptonic")
			samples.append(s+"_semi_charm")
			samples.append(s+"_semi_nocharm")

samples = [s for s in samples if not (s[0]+s[1] == "WW" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

########## Wjets sectioning for W plus c discrimination

if mode	== "mc":
	for s in samples:
        	if (s[0]+s[1]+s[2]+s[3]+s[4] == "Wjets" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                	df[s] = df[s].Define('ishard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,7)')
                	df[s] = df[s].Define('first_copy','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,12)')
                	df[s] = df[s].Define('isWplusc','wpluscbool(nGenPart,GenPart_statusFlags,GenPart_pdgId,ishard,first_copy)')
                	df[s+"_charm"] = df[s].Filter('isWplusc == 2')
                	df[s+"_bottom"] = df[s].Filter('isWplusc == 3')
                	df[s+"_doublecharm"] = df[s].Filter('isWplusc == 1')
                	df[s+"_light"] = df[s].Filter('isWplusc == 0')
                	## Samples correction
                	samples.append(s+"_charm")
                	samples.append(s+"_bottom")
                	samples.append(s+"_doublecharm")
                	samples.append(s+"_light")

samples = [s for s in samples if not (s[0]+s[1]+s[2]+s[3]+s[4] == "Wjets" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

########## ttbar sectioning for charm discrimination

if mode	== "mc":
	for s in samples:
        	if (s[0]+s[1]+s[2]+s[3]+s[4] == "ttbar" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                	df[s] = df[s].Define('isttbarC','ttbarcharm(nGenPart,GenPart_statusFlags,GenPart_pdgId,GenPart_genPartIdxMother)')
                	df[s+"_charm"] = df[s].Filter('isttbarC')
                	df[s+"_nocharm"] = df[s].Filter('!isttbarC')
                	## Samples correction
                	samples.append(s+"_charm")
                	samples.append(s+"_nocharm")

samples = [s for s in samples if not (s[0]+s[1]+s[2]+s[3]+s[4] == "ttbar"  and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

########## ST sectioning for charm discrimination

if mode == "mc":
        for s in samples:
                if (s[0]+s[1] == "ST" and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B")):
                        df[s] = df[s].Define('isSTC','ttbarcharm(nGenPart,GenPart_statusFlags,GenPart_pdgId,GenPart_genPartIdxMother)')
                        df[s+"_charm"] = df[s].Filter('isSTC')
                        df[s+"_nocharm"] = df[s].Filter('!isSTC')
                        ## Samples correction
                        samples.append(s+"_charm")
                        samples.append(s+"_nocharm")

samples = [s for s in samples if not (s[0]+s[1] == "ST"  and (s[-1]=="6" or s[-1]=="7" or s[-1]=="8" or s[-1]=="B"))]

print(samples)

## Differentiated definitions
for s in samples:
	df_muon[s] = df[s].Filter('nMuonGood>0')
	df_electron[s] = df[s].Filter('nElectronGood >0')
	if args.type == "mc":
		df_muon[s] = df_muon[s].Filter(muon_trig[years[0]])
		df_electron[s] = df_electron[s].Filter(el_trig[years[0]])
	else:
		if args.process == "M": 
			df_muon[s] = df_muon[s].Filter(muon_trig[years[0]])
			df_electron[s] = df_electron[s].Filter(muon_trig[years[0]])
		if args.process == "E":
			if years[0] == "2017":
				df_muon[s] = df_muon[s].Filter('triggeremulator(nElectron, nTrigObj, Electron_eta, Electron_deltaEtaSC, Electron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)')
				df_electron[s] = df_electron[s].Filter('triggeremulator(nElectron, nTrigObj, Electron_eta, Electron_deltaEtaSC, Electron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)')
			else:
				df_muon[s] = df_muon[s].Filter(el_trig[years[0]])
				df_electron[s] = df_electron[s].Filter(el_trig[years[0]])
	df_muon[s] = df_muon[s].Define('lepton_pt','Muon_pt[MuonGoodInd[0]]')
	df_muon[s] = df_muon[s].Define('lepton_eta','Muon_eta[MuonGoodInd[0]]')
	df_electron[s] = df_electron[s].Define('lepton_pt','Electron_pt[ElectronGoodInd[0]]')
	df_electron[s] = df_electron[s].Define('lepton_eta','Electron_eta[ElectronGoodInd[0]]')
	df_muon[s] = df_muon[s].Define('deltaR_jetM_lepM','ROOT::VecOps::DeltaR(Muon_eta[MuonGoodInd[0]],Jet_eta[JetMuonInd[0]] , Muon_phi[MuonGoodInd[0]], Jet_phi[JetMuonInd[0]])')
	df_electron[s] = df_electron[s].Define('deltaR_jetM_lepE','ROOT::VecOps::DeltaR(Electron_eta[ElectronGoodInd[0]],Jet_eta[JetMuonInd[0]] , Electron_phi[ElectronGoodInd[0]], Jet_phi[JetMuonInd[0]])')
	df_muon[s] = df_muon[s].Define('InvM_jetM_lep', 'InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]])')
	df_electron[s] = df_electron[s].Define('InvM_jetM_lep', 'InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Electron_pt[ElectronGoodInd[0]],Electron_eta[ElectronGoodInd[0]],Electron_phi[ElectronGoodInd[0]], Electron_mass[ElectronGoodInd[0]])')
	df_muon[s] = df_muon[s].Define('transverse_mass','std::sqrt(2*Muon_pt[MuonGoodInd[0]]*MET_pt*(1-std::cos(Muon_phi[MuonGoodInd[0]]-MET_phi)))')
	df_electron[s] = df_electron[s].Define('transverse_mass','std::sqrt(2*Electron_pt[ElectronGoodInd[0]]*MET_pt*(1-std::cos(Electron_phi[ElectronGoodInd[0]]-MET_phi)))')
	df_muon[s] = df_muon[s].Define('InvM_muon_jet','InvariantM(Muon_pt[MuonJetInd[0]],Muon_eta[MuonJetInd[0]],Muon_phi[MuonJetInd[0]],Muon_mass[MuonJetInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]])')
	df_electron[s] = df_electron[s].Define('lepton_iso', 'Electron_pfRelIso03_all[ElectronGoodInd[0]]')
	df_electron[s] = df_electron[s].Define('lepton_mva', 'Electron_mvaFall17V2noIso[ElectronGoodInd[0]]')
	## New cuts
	df_muon[s] = df_muon[s].Filter('InvM_muon_jet >12').Filter('InvM_muon_jet > 110 || InvM_muon_jet < 70')
	df_electron[s] = df_electron[s].Filter('transverse_mass > 50')
	df_muon[s] = df_muon[s].Filter('transverse_mass > 50')
	df_electron[s] = df_electron[s].Filter('muon_jet_relpt<0.5')
	df_muon[s] = df_muon[s].Filter('muon_jet_relpt<0.5')
	df_electron[s] = df_electron[s].Filter('EMtotal_jetM<0.4')
	df_muon[s] = df_muon[s].Filter('EMtotal_jetM<0.4')
	df_electron[s] = df_electron[s].Filter('jet_notmuon_btag<'+str(cuts_btag[years[0]][1]))
	df_muon[s] = df_muon[s].Filter('jet_notmuon_btag<'+str(cuts_btag[years[0]][1]))
	df_electron[s] = df_electron[s].Filter('jet_muon_btag<'+str(cuts_btag[years[0]][2]))
	df_muon[s] = df_muon[s].Filter('jet_muon_btag<'+str(cuts_btag[years[0]][2]))
	#df_electron[s] = df_electron[s].Filter('InvM_2jets<100').Filter('InvM_2jets>60')
	#df_muon[s] = df_muon[s].Filter('InvM_2jets<100').Filter('InvM_2jets>60')
	## jet pt
	#df_electron[s] = df_electron[s].Filter('jet_muon_pt + jet_notmuon_pt > 90.')
	#df_muon[s] = df_muon[s].Filter('jet_muon_pt + jet_notmuon_pt > 90.')
	#df_electron[s] = df_electron[s].Filter('deltaR_jetM_jetNM>1.5')
	#df_muon[s] = df_muon[s].Filter('deltaR_jetM_jetNM>1.5')
	#df_electron[s] = df_electron[s].Filter('muon_jet_sigr>2')
	#df_muon[s] = df_muon[s].Filter('muon_jet_sigr>2')

############################################################
####################     HISTS    ##########################
############################################################

hist_nJetGood_M = {}
hist_nJetGood_E = {}
hist_nMuoninJet_M = {}
hist_nMuoninJet_E = {}
hist_nLooseLepton_M = {}
hist_nLooseLepton_E = {}
hist_jet_muon_pt_M = {}
hist_jet_muon_nmu_M = {}
hist_jet_notmuon_pt_M = {}
hist_jet_muon_eta_M = {}
hist_jet_notmuon_eta_M = {}
hist_jet_notmuon_nmu_M = {}
hist_jet_notmuon_qgl_M = {}
hist_jet_notmuon_mass_M = {}
hist_jet_muon_pt_E = {}
hist_jet_muon_nmu_E = {}
hist_jet_notmuon_pt_E = {}
hist_jet_muon_eta_E = {}
hist_jet_notmuon_eta_E = {}
hist_jet_notmuon_nmu_E = {}
hist_jet_notmuon_qgl_E = {}
hist_jet_notmuon_mass_E = {}
hist_lepton_pt_M = {}
hist_lepton_eta_M = {}
hist_lepton_pt_E = {}
hist_lepton_eta_E = {}
hist_muon_jet_pt_M = {}
hist_muon_jet_eta_M = {}
hist_muon_jet_pt_E = {}
hist_muon_jet_eta_E = {}
hist_InvM_2jets_M = {}
hist_InvM_2jets_E = {}
hist_InvM_jetM_lepM = {}
hist_InvM_jetM_lepE = {}
hist_deltaR_jetM_lepM = {}
hist_deltaR_jetM_lepE =	{}
hist_deltaR_jetM_jetNM_M = {}
hist_deltaR_jetM_jetNM_E = {}
hist_deltaphi_jetM_jetNM_M = {}
hist_deltaphi_jetM_jetNM_E = {}
hist_deltaeta_jetM_jetNM_M = {}
hist_deltaeta_jetM_jetNM_E = {}
hist_deltapt_jetM_jetNM_M = {}
hist_deltapt_jetM_jetNM_E = {}
hist_MET_M = {}
hist_MET_E = {}
hist_tranverse_massM = {}
hist_tranverse_massE = {}
hist_tracks_jetM_M = {}
hist_tracks_jetNM_M = {}
hist_tracks_jetM_E = {}
hist_tracks_jetNM_E = {}
hist_EMN_jetM_M = {}
hist_EMC_jetM_M = {}
hist_EMN_jetM_E = {}
hist_EMC_jetM_E = {}
hist_EMtotal_jetM_M = {}
hist_EMtotal_jetM_E = {}
hist_InvM_muon_jet_M = {}
hist_muon_jet_mva_M = {}
hist_muon_jet_mva_E = {}
hist_muon_jet_tight_M = {}
hist_muon_jet_tight_E = {}
hist_muon_jet_relpt_M = {}
hist_muon_jet_relpt_E = {}
hist_SSOS_M = {}
hist_SSOS_E = {}
hist_lepton_iso_E = {}
hist_lepton_mva_E = {}
hist_muon_jet_sigxy_M = {}
hist_muon_jet_sigxy_E = {}
hist_muon_jet_sigz_M = {}
hist_muon_jet_sigz_E = {}
hist_muon_jet_sigr_M = {}
hist_muon_jet_sigr_E = {}
hist_pT_sum_M = {}
hist_pT_sum_E = {}
hist_pT_product_M = {}
hist_pT_product_E = {}
hist_deltaR_lep_2jets_M = {}
hist_deltaR_lep_2jets_E	= {}
hist_deltaphi_MET_2jets_M = {}
hist_deltaphi_MET_2jets_E = {}
hist_deltaphi_lephad_M = {}
hist_deltaphi_lephad_E = {}
hist_eta_2jets_M = {}
hist_eta_2jets_E = {}
hist_pt_2jets_M = {}
hist_pt_2jets_E = {}
hist_pt_Wlep_M = {}
hist_pt_Wlep_E = {}
hist_deltaR_lephad_M = {}
hist_deltaR_lephad_E = {}
hist_deltaphi_lep_2jets_M = {}
hist_deltaphi_lep_2jets_E = {}
hist_deltaeta_lephad_M = {}
hist_deltaeta_lephad_E = {}
hist_deltaeta_lep_2jets_M = {}
hist_deltaeta_lep_2jets_E = {}
hist_deltapt_lephad_M = {}
hist_deltapt_lephad_E = {}
hist_deltapt_lep_2jets_M = {}
hist_deltapt_lep_2jets_E = {}
hist_jet_muon_btag_M = {}
hist_jet_muon_btag_E = {}
hist_jet_notmuon_btag_M = {}
hist_jet_notmuon_btag_E = {}
hist_jet_notmuon_deeptagG_M = {}
hist_jet_notmuon_deeptagG_E = {}
hist_jet_notmuon_deeptagC_M = {}
hist_jet_notmuon_deeptagC_E = {}
hist_deltaR_lep_jet1_M = {}
hist_deltaR_lep_jet2_M = {}
hist_deltaPhi_lep_jet1_M = {}
hist_deltaPhi_lep_jet2_M = {}
hist_deltaEta_lep_jet1_M = {}
hist_deltaEta_lep_jet2_M = {}
hist_deltaR_lep_jet1_E = {}
hist_deltaR_lep_jet2_E = {}
hist_deltaPhi_lep_jet1_E = {}
hist_deltaPhi_lep_jet2_E = {}
hist_deltaEta_lep_jet1_E = {}
hist_deltaEta_lep_jet2_E = {}
hist_pT_proy_M = {} 
hist_pT_proy_E = {}
hist_pT_sum_2J_M = {}
hist_pT_sum_2J_E = {}
hist_Deta_InvM_M = {}
hist_Deta_InvM_E = {}
hist_Deta_qgl_M = {}
hist_Deta_qgl_E = {}
hist_DR_InvM_M = {}
hist_DR_InvM_E = {}
hist_qgl_InvM_M = {}
hist_qgl_InvM_E = {}

for s in samples:
        hist_nJetGood_M[s] = df_muon[s].Histo1D(("nJetGood_M","",10,0,10),"nJetGood","weightSSOS")
        hist_nJetGood_E[s] = df_electron[s].Histo1D(("nJetGood_E","",10,0,10),"nJetGood","weightSSOS")

        hist_nMuoninJet_M[s] = df_muon[s].Histo1D(("nMuoninJet_M","",10,0,10),"nMuoninJet","weightSSOS")
        hist_nMuoninJet_E[s] = df_electron[s].Histo1D(("nMuoninJet_E","",10,0,10),"nMuoninJet","weightSSOS")

        hist_nLooseLepton_M[s] = df_muon[s].Histo1D(("nLooseLepton_M","",10,0,10),"nLooseLepton","weightSSOS")
        hist_nLooseLepton_E[s] = df_electron[s].Histo1D(("nLooseLepton_E","",10,0,10),"nLooseLepton","weightSSOS")

        hist_jet_muon_pt_M[s] = df_muon[s].Histo1D(("jet_muon_pt_M","",50,20,120),"jet_muon_pt","weightSSOS")
        hist_jet_muon_nmu_M[s] = df_muon[s].Histo1D(("jet_muon_nmu_M","",10,0,10),"jet_muon_nmu","weightSSOS")
        hist_jet_notmuon_pt_M[s] = df_muon[s].Histo1D(("jet_not_muon_pt_M","",50,20,120),"jet_notmuon_pt","weightSSOS")
        hist_jet_muon_eta_M[s] = df_muon[s].Histo1D(("jet_muon_eta_M","",80,-4,4),"jet_muon_eta","weightSSOS")
        hist_jet_notmuon_eta_M[s] = df_muon[s].Histo1D(("jet_not_muon_eta_M","",80,-4,4),"jet_notmuon_eta","weightSSOS")
        hist_jet_notmuon_nmu_M[s] = df_muon[s].Histo1D(("jet_notmuon_nmu_M","",10,0,10),"jet_notmuon_nmu","weightSSOS")
        hist_jet_notmuon_mass_M[s] = df_muon[s].Histo1D(("jet_notmuon_mass_M","",40,0,40),"jet_notmuon_mass","weightSSOS")
        hist_jet_notmuon_qgl_M[s] = df_muon[s].Histo1D(("jet_notmuon_qgl_M","",100,0,1),"jet_notmuon_qgl","weightSSOS")

        hist_jet_muon_pt_E[s] = df_electron[s].Histo1D(("jet_muon_pt_E","",50,20,120),"jet_muon_pt","weightSSOS")
        hist_jet_muon_nmu_E[s] = df_electron[s].Histo1D(("jet_muon_nmu_E","",10,0,10),"jet_muon_nmu","weightSSOS")
        hist_jet_notmuon_pt_E[s] = df_electron[s].Histo1D(("jet_not_muon_pt_E","",50,20,120),"jet_notmuon_pt","weightSSOS")
        hist_jet_muon_eta_E[s] = df_electron[s].Histo1D(("jet_muon_eta_E","",80,-4,4),"jet_muon_eta","weightSSOS")
        hist_jet_notmuon_eta_E[s] = df_electron[s].Histo1D(("jet_not_muon_eta_E","",80,-4,4),"jet_notmuon_eta","weightSSOS")
        hist_jet_notmuon_nmu_E[s] = df_electron[s].Histo1D(("jet_notmuon_nmu_E","",10,0,10),"jet_notmuon_nmu","weightSSOS")
        hist_jet_notmuon_mass_E[s] = df_electron[s].Histo1D(("jet_notmuon_mass_E","",40,0,40),"jet_notmuon_mass","weightSSOS")
        hist_jet_notmuon_qgl_E[s] = df_electron[s].Histo1D(("jet_notmuon_qgl_E","",100,0,1),"jet_notmuon_qgl","weightSSOS")

        hist_lepton_pt_M[s] = df_muon[s].Histo1D(("lepton_pt_M","",50,20,120),"lepton_pt","weightSSOS")
        hist_lepton_eta_M[s] = df_muon[s].Histo1D(("lepton_eta_M","",80,-4,4),"lepton_eta","weightSSOS")

        hist_lepton_pt_E[s] = df_electron[s].Histo1D(("lepton_pt_E","",50,20,120),"lepton_pt","weightSSOS")
        hist_lepton_eta_E[s] = df_electron[s].Histo1D(("lepton_eta_E","",80,-4,4),"lepton_eta","weightSSOS")

        hist_muon_jet_pt_M[s] = df_muon[s].Histo1D(("muon_jet_pt_M","",50,0,25),"muon_jet_pt","weightSSOS")
        hist_muon_jet_eta_M[s] = df_muon[s].Histo1D(("muon_jet_eta_M","",80,-4,4),"muon_jet_eta","weightSSOS")

        hist_muon_jet_pt_E[s] = df_electron[s].Histo1D(("muon_jet_pt_E","",50,0,25),"muon_jet_pt","weightSSOS")
        hist_muon_jet_eta_E[s] = df_electron[s].Histo1D(("muon_jet_eta_E","",80,-4,4),"muon_jet_eta","weightSSOS")

        hist_InvM_2jets_M[s] = df_muon[s].Histo1D(("InvM_2jets_M","",100,0,300),"InvM_2jets","weightSSOS")
        hist_InvM_2jets_E[s] = df_electron[s].Histo1D(("InvM_2jets_E","",100,0,300),"InvM_2jets","weightSSOS")
        hist_InvM_jetM_lepM[s] = df_muon[s].Histo1D(("InvM_jetM_lepM","",100,0,300),"InvM_jetM_lep","weightSSOS")
        hist_InvM_jetM_lepE[s] = df_electron[s].Histo1D(("InvM_jetM_lepE","",100,0,300),"InvM_jetM_lep","weightSSOS")

        hist_MET_M[s] = df_muon[s].Histo1D(("MET_pt_M","",100,0,150),"MET_pt","weightSSOS")
        hist_MET_E[s] = df_electron[s].Histo1D(("MET_pt_E","",100,0,150),"MET_pt","weightSSOS")

        hist_deltaR_jetM_lepM[s] = df_muon[s].Histo1D(("deltaR_jetM_lepM","",100,0,5),"deltaR_jetM_lepM","weightSSOS")
        hist_deltaR_jetM_lepE[s] = df_electron[s].Histo1D(("deltaR_jetM_lepE","",100,0,5),"deltaR_jetM_lepE","weightSSOS")
        hist_deltaR_jetM_jetNM_M[s] = df_muon[s].Histo1D(("deltaR_jetM_jetNM_M","",100,0,5),"deltaR_jetM_jetNM","weightSSOS")
        hist_deltaR_jetM_jetNM_E[s] = df_electron[s].Histo1D(("deltaR_jetM_jetNM_E","",100,0,5),"deltaR_jetM_jetNM","weightSSOS")

        hist_deltaphi_jetM_jetNM_M[s] = df_muon[s].Histo1D(("deltaphi_jetM_jetNM_M","",100,0,5),"deltaphi_jetM_jetNM","weightSSOS")
        hist_deltaphi_jetM_jetNM_E[s] = df_electron[s].Histo1D(("deltaphi_jetM_jetNM_E","",100,0,5),"deltaphi_jetM_jetNM","weightSSOS")

        hist_deltaeta_jetM_jetNM_M[s] = df_muon[s].Histo1D(("deltaeta_jetM_jetNM_M","",100,0,5),"deltaeta_jetM_jetNM","weightSSOS")
        hist_deltaeta_jetM_jetNM_E[s] = df_electron[s].Histo1D(("deltaeta_jetM_jetNM_E","",100,0,5),"deltaeta_jetM_jetNM","weightSSOS")

        hist_deltapt_jetM_jetNM_M[s] = df_muon[s].Histo1D(("deltapt_jetM_jetNM_M","",100,0,50),"deltapt_jetM_jetNM","weightSSOS")
        hist_deltapt_jetM_jetNM_E[s] = df_electron[s].Histo1D(("deltapt_jetM_jetNM_E","",100,0,50),"deltapt_jetM_jetNM","weightSSOS")

        hist_tranverse_massM[s] = df_muon[s].Histo1D(("transverse_massM","",50,0,150),"transverse_mass","weightSSOS")
        hist_tranverse_massE[s] = df_electron[s].Histo1D(("transverse_massE","",50,0,150),"transverse_mass","weightSSOS")

        hist_tracks_jetM_M[s] = df_muon[s].Histo1D(("tracks_jetM_M","",60,0,60),"tracks_jetM","weightSSOS")
        hist_tracks_jetNM_M[s] = df_muon[s].Histo1D(("tracks_jetNM_M","",60,0,60),"tracks_jetNM","weightSSOS")

        hist_tracks_jetM_E[s] = df_electron[s].Histo1D(("tracks_jetM_E","",60,0,60),"tracks_jetM","weightSSOS")
        hist_tracks_jetNM_E[s] = df_electron[s].Histo1D(("tracks_jetNM_E","",60,0,60),"tracks_jetNM","weightSSOS")

        hist_EMN_jetM_M[s] = df_muon[s].Histo1D(("EMN_jetM_M","",60,0,1),"EMN_jetM","weightSSOS")
        hist_EMC_jetM_M[s] = df_muon[s].Histo1D(("EMC_jetM_M","",60,0,1),"EMC_jetM","weightSSOS")
        hist_EMN_jetM_E[s] = df_electron[s].Histo1D(("EMN_jetM_E","",60,0,1),"EMN_jetM","weightSSOS")
        hist_EMC_jetM_E[s] = df_electron[s].Histo1D(("EMC_jetM_E","",60,0,1),"EMC_jetM","weightSSOS")

        hist_EMtotal_jetM_M[s] = df_muon[s].Histo1D(("EMtotal_jetM_M","",60,0,1),"EMtotal_jetM","weightSSOS")
        hist_EMtotal_jetM_E[s] = df_electron[s].Histo1D(("EMtotal_jetM_E","",60,0,1),"EMtotal_jetM","weightSSOS")

        hist_InvM_muon_jet_M[s] = df_muon[s].Histo1D(("InvM_muon_jet_M","",50,0,200),"InvM_muon_jet","weightSSOS")

        hist_muon_jet_mva_M[s] = df_muon[s].Histo1D(("muon_jet_mva_M","",50,0,1),"muon_jet_mva","weightSSOS")
        hist_muon_jet_mva_E[s] = df_electron[s].Histo1D(("muon_jet_mva_E","",50,0,1),"muon_jet_mva","weightSSOS")

        hist_muon_jet_tight_M[s] = df_muon[s].Histo1D(("muon_jet_tight_M","",2,0,2),"muon_jet_tight","weightSSOS")
        hist_muon_jet_tight_E[s] = df_electron[s].Histo1D(("muon_jet_tight_E","",2,0,2),"muon_jet_tight","weightSSOS")

        hist_muon_jet_relpt_M[s] = df_muon[s].Histo1D(("muon_jet_relpt_M","",50,0,1),"muon_jet_relpt","weightSSOS")
        hist_muon_jet_relpt_E[s] = df_electron[s].Histo1D(("muon_jet_relpt_E","",50,0,1),"muon_jet_relpt","weightSSOS")

        hist_muon_jet_sigxy_M[s] = df_muon[s].Histo1D(("muon_jet_sigxy_M","",100,-4,4),"muon_jet_sigxy","weightSSOS")
        hist_muon_jet_sigxy_E[s] = df_electron[s].Histo1D(("muon_jet_sigxy_E","",100,-4,4),"muon_jet_sigxy","weightSSOS")

        hist_muon_jet_sigz_M[s] = df_muon[s].Histo1D(("muon_jet_sigz_M","",100,-4,4),"muon_jet_sigz","weightSSOS")
        hist_muon_jet_sigz_E[s] = df_electron[s].Histo1D(("muon_jet_sigz_E","",100,-4,4),"muon_jet_sigz","weightSSOS")

        hist_muon_jet_sigr_M[s] = df_muon[s].Histo1D(("muon_jet_sigr_M","",100,0,30),"muon_jet_sigr","weightSSOS")
        hist_muon_jet_sigr_E[s] = df_electron[s].Histo1D(("muon_jet_sigr_E","",100,0,30),"muon_jet_sigr","weightSSOS")

        hist_SSOS_M[s] = df_muon[s].Histo1D(("SSOS_M","",4,-2,2),"MuonLepSign","weightSSOS")
        hist_SSOS_E[s] = df_electron[s].Histo1D(("SSOS_E","",4,-2,2),"MuonLepSign","weightSSOS")

        hist_lepton_iso_E[s] = df_electron[s].Histo1D(("lepton_iso_E","",100,0,0.4),"lepton_iso","weightSSOS")
        hist_lepton_mva_E[s] = df_electron[s].Histo1D(("lepton_mva_E","",80,0,1),"lepton_mva","weightSSOS")

        hist_pT_sum_M[s] = df_muon[s].Histo1D(("pT_sum_M","",100,0,1),"pT_sum","weightSSOS")
        hist_pT_sum_E[s] = df_electron[s].Histo1D(("pT_sum_E","",100,0,1),"pT_sum","weightSSOS")

        hist_pT_product_M[s] = df_muon[s].Histo1D(("pT_product_M","",100,-1,1),"pT_product","weightSSOS")
        hist_pT_product_E[s] = df_electron[s].Histo1D(("pT_product_E","",100,-1,1),"pT_product","weightSSOS")

        hist_deltaR_lep_2jets_M[s] = df_muon[s].Histo1D(("deltaR_lep_2jets_M","",100,0,5),"deltaR_lep_2jets","weightSSOS")
        hist_deltaR_lep_2jets_E[s] = df_electron[s].Histo1D(("deltaR_lep_2jets_E","",100,0,5),"deltaR_lep_2jets","weightSSOS")
        hist_deltaphi_MET_2jets_M[s] = df_muon[s].Histo1D(("deltaphi_MET_2jets_M","",100,0,5),"deltaphi_MET_2jets","weightSSOS")
        hist_deltaphi_MET_2jets_E[s] = df_electron[s].Histo1D(("deltaphi_MET_2jets_E","",100,0,5),"deltaphi_MET_2jets","weightSSOS")
        hist_deltaphi_lephad_M[s] = df_muon[s].Histo1D(("deltaphi_lephad_M","",100,0,5),"deltaphi_lephad","weightSSOS")
        hist_deltaphi_lephad_E[s] = df_electron[s].Histo1D(("deltaphi_lephad_E","",100,0,5),"deltaphi_lephad","weightSSOS")

        hist_eta_2jets_M[s] = df_muon[s].Histo1D(("eta_2jets_M","",50,-5,5),"eta_2jets","weightSSOS")
        hist_eta_2jets_E[s] = df_electron[s].Histo1D(("eta_2jets_E","",50,-5,5),"eta_2jets","weightSSOS")
        hist_pt_2jets_M[s] = df_muon[s].Histo1D(("pt_2jets_M","",100,0,200),"pt_2jets","weightSSOS")
        hist_pt_2jets_E[s] = df_electron[s].Histo1D(("pt_2jets_E","",100,0,200),"pt_2jets","weightSSOS")

        hist_pt_Wlep_M[s] = df_muon[s].Histo1D(("pt_Wlep_M","",100,0,200),"pT_Wlep","weightSSOS")
        hist_pt_Wlep_E[s] = df_electron[s].Histo1D(("pt_Wlep_E","",100,0,200),"pT_Wlep","weightSSOS")

        hist_deltaR_lephad_M[s] = df_muon[s].Histo1D(("deltaR_lephad_M","",100,0,5),"deltaR_lephad","weightSSOS")
        hist_deltaR_lephad_E[s] = df_electron[s].Histo1D(("deltaR_lephad_E","",100,0,5),"deltaR_lephad","weightSSOS")
        hist_deltaphi_lep_2jets_M[s] = df_muon[s].Histo1D(("deltaphi_lep_2jets_M","",100,0,5),"deltaphi_lep_2jets","weightSSOS")
        hist_deltaphi_lep_2jets_E[s] = df_electron[s].Histo1D(("deltaphi_lep_2jets_E","",100,0,5),"deltaphi_lep_2jets","weightSSOS")
        hist_deltaeta_lephad_M[s] = df_muon[s].Histo1D(("deltaeta_lephad_M","",100,0,5),"deltaeta_lephad","weightSSOS")
        hist_deltaeta_lephad_E[s] = df_electron[s].Histo1D(("deltaeta_lephad_E","",100,0,5),"deltaeta_lephad","weightSSOS")
        hist_deltaeta_lep_2jets_M[s] = df_muon[s].Histo1D(("deltaeta_lep_2jets_M","",100,0,5),"deltaeta_lep_2jets","weightSSOS")
        hist_deltaeta_lep_2jets_E[s] = df_electron[s].Histo1D(("deltaeta_lep_2jets_E","",100,0,5),"deltaeta_lep_2jets","weightSSOS")
        hist_deltapt_lephad_M[s] = df_muon[s].Histo1D(("deltapt_lephad_M","",100,0,5),"deltapt_lephad","weightSSOS")
        hist_deltapt_lephad_E[s] = df_electron[s].Histo1D(("deltapt_lephad_E","",100,0,5),"deltapt_lephad","weightSSOS")
        hist_deltapt_lep_2jets_M[s] = df_muon[s].Histo1D(("deltapt_lep_2jets_M","",100,0,5),"deltapt_lep_2jets","weightSSOS")
        hist_deltapt_lep_2jets_E[s] = df_electron[s].Histo1D(("deltapt_lep_2jets_E","",100,0,5),"deltapt_lep_2jets","weightSSOS")

        hist_deltaR_lep_jet1_M[s] = df_muon[s].Histo1D(("deltaR_lep_jet1_M","",100,0,5),"deltaR_lep_jet1","weightSSOS")
        hist_deltaR_lep_jet2_M[s] = df_muon[s].Histo1D(("deltaR_lep_jet2_M","",100,0,5),"deltaR_lep_jet2","weightSSOS")
        hist_deltaPhi_lep_jet1_M[s] = df_muon[s].Histo1D(("deltaPhi_lep_jet1_M","",100,0,5),"deltaPhi_lep_jet1","weightSSOS")
        hist_deltaPhi_lep_jet2_M[s] = df_muon[s].Histo1D(("deltaPhi_lep_jet2_M","",100,0,5),"deltaPhi_lep_jet2","weightSSOS")
        hist_deltaEta_lep_jet1_M[s] = df_muon[s].Histo1D(("deltaEta_lep_jet1_M","",100,0,5),"deltaEta_lep_jet1","weightSSOS")
        hist_deltaEta_lep_jet2_M[s] = df_muon[s].Histo1D(("deltaEta_lep_jet2_M","",100,0,5),"deltaEta_lep_jet2","weightSSOS")

        hist_deltaR_lep_jet1_E[s] = df_electron[s].Histo1D(("deltaR_lep_jet1_E","",100,0,5),"deltaR_lep_jet1","weightSSOS")
        hist_deltaR_lep_jet2_E[s] = df_electron[s].Histo1D(("deltaR_lep_jet2_E","",100,0,5),"deltaR_lep_jet2","weightSSOS")
        hist_deltaPhi_lep_jet1_E[s] = df_electron[s].Histo1D(("deltaPhi_lep_jet1_E","",100,0,5),"deltaPhi_lep_jet1","weightSSOS")
        hist_deltaPhi_lep_jet2_E[s] = df_electron[s].Histo1D(("deltaPhi_lep_jet2_E","",100,0,5),"deltaPhi_lep_jet2","weightSSOS")
        hist_deltaEta_lep_jet1_E[s] = df_electron[s].Histo1D(("deltaEta_lep_jet1_E","",100,0,5),"deltaEta_lep_jet1","weightSSOS")
        hist_deltaEta_lep_jet2_E[s] = df_electron[s].Histo1D(("deltaEta_lep_jet2_E","",100,0,5),"deltaEta_lep_jet2","weightSSOS")

        hist_jet_muon_btag_M[s] = df_muon[s].Histo1D(("jet_muon_btag_M","",50,0,1),"jet_muon_btag","weightSSOS")
        hist_jet_muon_btag_E[s] = df_electron[s].Histo1D(("jet_muon_btag_E","",50,0,1),"jet_muon_btag","weightSSOS")
        hist_jet_notmuon_btag_M[s] = df_muon[s].Histo1D(("jet_notmuon_btag_M","",50,0,1),"jet_notmuon_btag","weightSSOS")
        hist_jet_notmuon_btag_E[s] = df_electron[s].Histo1D(("jet_notmuon_btag_E","",50,0,1),"jet_notmuon_btag","weightSSOS")

        hist_jet_notmuon_deeptagG_M[s] = df_muon[s].Histo1D(("jet_notmuon_deeptagG_M","",50,0,1),"jet_notmuon_deeptagG","weightSSOS")
        hist_jet_notmuon_deeptagG_E[s] = df_electron[s].Histo1D(("jet_notmuon_deeptagG_E","",50,0,1),"jet_notmuon_deeptagG","weightSSOS")
        hist_jet_notmuon_deeptagC_M[s] = df_muon[s].Histo1D(("jet_notmuon_deeptagC_M","",50,0,1),"jet_notmuon_deeptagC","weightSSOS")
        hist_jet_notmuon_deeptagC_E[s] = df_electron[s].Histo1D(("jet_notmuon_deeptagC_E","",50,0,1),"jet_notmuon_deeptagC","weightSSOS")

        hist_pT_proy_M[s] = df_muon[s].Histo1D(("pT_proy_M","",100,-100,100),"pT_proy","weightSSOS")
        hist_pT_proy_E[s] = df_electron[s].Histo1D(("pT_proy_E","",100,-100,100),"pT_proy","weightSSOS")

        hist_pT_sum_2J_M[s] = df_muon[s].Histo1D(("pT_sum_2J_M","",100,0,1),"pT_sum_2J","weightSSOS")
        hist_pT_sum_2J_E[s] = df_electron[s].Histo1D(("pT_sum_2J_E","",100,0,1),"pT_sum_2J","weightSSOS")

        ## 2D hists, all for deltas between the 2 jets, the only discriminant variables
        hist_Deta_InvM_M[s] = df_muon[s].Histo2D(("deltaeta_InvM_M","",100,0,5,100,0,300),"deltaeta_jetM_jetNM","InvM_2jets","weightSSOS")
        hist_Deta_InvM_E[s] = df_electron[s].Histo2D(("deltaeta_InvM_E","",100,0,5,100,0,300),"deltaeta_jetM_jetNM","InvM_2jets","weightSSOS")
        hist_Deta_qgl_M[s] = df_muon[s].Histo2D(("deltaeta_qgl_M","",100,0,5,100,0,1),"deltaeta_jetM_jetNM","jet_notmuon_qgl","weightSSOS")
        hist_Deta_qgl_E[s] = df_electron[s].Histo2D(("deltaeta_qgl_E","",100,0,5,100,0,1),"deltaeta_jetM_jetNM","jet_notmuon_qgl","weightSSOS")
        hist_DR_InvM_M[s] = df_muon[s].Histo2D(("deltaR_InvM_M","",100,0,5,100,0,300),"deltaR_jetM_jetNM","InvM_2jets","weightSSOS")
        hist_DR_InvM_E[s] = df_electron[s].Histo2D(("deltaR_InvM_E","",100,0,5,100,0,300),"deltaR_jetM_jetNM","InvM_2jets","weightSSOS")
        hist_qgl_InvM_M[s] = df_muon[s].Histo2D(("qgl_InvM_M","",100,0,1,100,0,300),"jet_notmuon_qgl","InvM_2jets","weightSSOS")
        hist_qgl_InvM_E[s] = df_electron[s].Histo2D(("qgl_InvM_E","",100,0,1,100,0,300),"jet_notmuon_qgl","InvM_2jets","weightSSOS")

#############################
####     DATA SAVING     ####
#############################

if args.notfull:
	for s in samples:
		if args.ssos : path_hist = '/nfs/cms/vazqueze/analisisWW/hists/ssos/hists_v1v2v3v4SSOS_'+s+'_range_'+str(args.list[0])+'_'+str(args.list[1])+'.root'
		else: path_hist = '/nfs/cms/vazqueze/analisisWW/hists/hists_v1v2v3v4_'+s+'_range_'+str(args.list[0])+'_'+str(args.list[1])+'.root'
		myfile = TFile( path_hist, 'RECREATE' )

		hist_nJetGood_M[s].Write()
		hist_nJetGood_E[s].Write()
		hist_nMuoninJet_M[s].Write()
		hist_nMuoninJet_E[s].Write()
		hist_nLooseLepton_M[s].Write()
		hist_nLooseLepton_E[s].Write()
		hist_jet_muon_pt_M[s].Write()
		hist_jet_muon_nmu_M[s].Write()
		hist_jet_muon_eta_M[s].Write()
		hist_jet_muon_pt_E[s].Write()
		hist_jet_muon_eta_E[s].Write()
		hist_jet_muon_nmu_E[s].Write()
		hist_jet_notmuon_pt_M[s].Write()
		hist_jet_notmuon_eta_M[s].Write()
		hist_jet_notmuon_nmu_M[s].Write()
		hist_jet_notmuon_mass_M[s].Write()
		hist_jet_notmuon_qgl_M[s].Write()
		hist_jet_notmuon_pt_E[s].Write()
		hist_jet_notmuon_eta_E[s].Write()
		hist_jet_notmuon_nmu_E[s].Write()
		hist_jet_notmuon_nmu_E[s].Write()
		hist_jet_notmuon_mass_E[s].Write()
		hist_jet_notmuon_qgl_E[s].Write()
		hist_lepton_pt_M[s].Write()
		hist_lepton_eta_M[s].Write()
		hist_lepton_pt_E[s].Write()
		hist_lepton_eta_E[s].Write()
		hist_muon_jet_pt_M[s].Write()
		hist_muon_jet_eta_M[s].Write()
		hist_muon_jet_pt_E[s].Write()
		hist_muon_jet_eta_E[s].Write()
		hist_InvM_2jets_M[s].Write()
		hist_InvM_2jets_E[s].Write()
		hist_InvM_jetM_lepM[s].Write()
		hist_InvM_jetM_lepE[s].Write()
		hist_deltaR_jetM_lepM[s].Write()
		hist_deltaR_jetM_lepE[s].Write()
		hist_deltaR_jetM_jetNM_M[s].Write()
		hist_deltaR_jetM_jetNM_E[s].Write()
		hist_deltaphi_jetM_jetNM_M[s].Write()
		hist_deltaphi_jetM_jetNM_E[s].Write()
		hist_deltaeta_jetM_jetNM_M[s].Write()
		hist_deltaeta_jetM_jetNM_E[s].Write()
		hist_deltapt_jetM_jetNM_M[s].Write()
		hist_deltapt_jetM_jetNM_E[s].Write()
		hist_MET_M[s].Write()
		hist_MET_E[s].Write()
		hist_tranverse_massM[s].Write()
		hist_tranverse_massE[s].Write()
		hist_tracks_jetM_M[s].Write()
		hist_tracks_jetNM_M[s].Write()
		hist_tracks_jetM_E[s].Write()
		hist_tracks_jetNM_E[s].Write()
		hist_EMN_jetM_M[s].Write()
		hist_EMC_jetM_M[s].Write()
		hist_EMN_jetM_E[s].Write()
		hist_EMC_jetM_E[s].Write()
		hist_EMtotal_jetM_M[s].Write()
		hist_EMtotal_jetM_E[s].Write()
		hist_InvM_muon_jet_M[s].Write()
		hist_muon_jet_mva_M[s].Write()
		hist_muon_jet_mva_E[s].Write()
		hist_muon_jet_tight_M[s].Write()
		hist_muon_jet_tight_E[s].Write()
		hist_muon_jet_relpt_M[s].Write()
		hist_muon_jet_relpt_E[s].Write()
		hist_muon_jet_sigxy_M[s].Write()
		hist_muon_jet_sigxy_E[s].Write()
		hist_muon_jet_sigz_M[s].Write()
		hist_muon_jet_sigz_E[s].Write()
		hist_muon_jet_sigr_M[s].Write()
		hist_muon_jet_sigr_E[s].Write()
		hist_SSOS_M[s].Write()
		hist_SSOS_E[s].Write()
		hist_lepton_iso_E[s].Write()
		hist_lepton_mva_E[s].Write()
		hist_pT_sum_M[s].Write()
		hist_pT_sum_E[s].Write()
		hist_pT_product_M[s].Write()
		hist_pT_product_E[s].Write()
		hist_deltaR_lep_2jets_M[s].Write()
		hist_deltaR_lep_2jets_E[s].Write()
		hist_deltaphi_MET_2jets_M[s].Write()
		hist_deltaphi_MET_2jets_E[s].Write()
		hist_deltaphi_lephad_M[s].Write()
		hist_deltaphi_lephad_E[s].Write()
		hist_eta_2jets_M[s].Write()
		hist_eta_2jets_E[s].Write()
		hist_pt_2jets_M[s].Write()
		hist_pt_2jets_E[s].Write()
		hist_pt_Wlep_M[s].Write()
		hist_pt_Wlep_E[s].Write()
		hist_deltaR_lephad_M[s].Write()
		hist_deltaR_lephad_E[s].Write()
		hist_deltaphi_lep_2jets_M[s].Write()
		hist_deltaphi_lep_2jets_E[s].Write()
		hist_deltaeta_lephad_M[s].Write()
		hist_deltaeta_lephad_E[s].Write()
		hist_deltaeta_lep_2jets_M[s].Write()
		hist_deltaeta_lep_2jets_E[s].Write()
		hist_deltapt_lephad_M[s].Write()
		hist_deltapt_lephad_E[s].Write()
		hist_deltapt_lep_2jets_M[s].Write()
		hist_deltapt_lep_2jets_E[s].Write()
		hist_deltaR_lep_jet1_M[s].Write()
		hist_deltaR_lep_jet2_M[s].Write()
		hist_deltaPhi_lep_jet1_M[s].Write()
		hist_deltaPhi_lep_jet2_M[s].Write()
		hist_deltaEta_lep_jet1_M[s].Write()
		hist_deltaEta_lep_jet2_M[s].Write()
		hist_deltaR_lep_jet1_E[s].Write()
		hist_deltaR_lep_jet2_E[s].Write()
		hist_deltaPhi_lep_jet1_E[s].Write()
		hist_deltaPhi_lep_jet2_E[s].Write()
		hist_deltaEta_lep_jet1_E[s].Write()
		hist_deltaEta_lep_jet2_E[s].Write()
		hist_jet_muon_btag_M[s].Write()
		hist_jet_muon_btag_E[s].Write()
		hist_jet_notmuon_btag_M[s].Write()
		hist_jet_notmuon_btag_E[s].Write()
		hist_jet_notmuon_deeptagG_M[s].Write()
		hist_jet_notmuon_deeptagG_E[s].Write()
		hist_jet_notmuon_deeptagC_M[s].Write()
		hist_jet_notmuon_deeptagC_E[s].Write()
		hist_Deta_InvM_M[s].Write()
		hist_Deta_InvM_E[s].Write()
		hist_Deta_qgl_M[s].Write()
		hist_Deta_qgl_E[s].Write()
		hist_DR_InvM_M[s].Write()
		hist_DR_InvM_E[s].Write()
		hist_qgl_InvM_M[s].Write()
		hist_qgl_InvM_E[s].Write()
		hist_pT_proy_M[s].Write()
		hist_pT_proy_E[s].Write()
		hist_pT_sum_2J_M[s].Write()
		hist_pT_sum_2J_E[s].Write()

		myfile.Close()

else:
	for s in samples:
                if args.ssos: path_hist = '/nfs/cms/vazqueze/analisisWW/hists/ssos/hists_v1v2v3v4SSOS_'+s+'.root'
                else: path_hist = '/nfs/cms/vazqueze/analisisWW/hists/hists_v1v2v3v4_'+s+'.root'
                myfile = TFile( path_hist, 'RECREATE' )

                hist_nJetGood_M[s].Write()
                hist_nJetGood_E[s].Write()
                hist_nMuoninJet_M[s].Write()
                hist_nMuoninJet_E[s].Write()
                hist_nLooseLepton_M[s].Write()
                hist_nLooseLepton_E[s].Write()
                hist_jet_muon_pt_M[s].Write()
                hist_jet_muon_nmu_M[s].Write()
                hist_jet_muon_eta_M[s].Write()
                hist_jet_muon_pt_E[s].Write()
                hist_jet_muon_eta_E[s].Write()
                hist_jet_muon_nmu_E[s].Write()
                hist_jet_notmuon_pt_M[s].Write()
                hist_jet_notmuon_eta_M[s].Write()
                hist_jet_notmuon_nmu_M[s].Write()
                hist_jet_notmuon_mass_M[s].Write()
                hist_jet_notmuon_qgl_M[s].Write()
                hist_jet_notmuon_pt_E[s].Write()
                hist_jet_notmuon_eta_E[s].Write()
                hist_jet_notmuon_nmu_E[s].Write()
                hist_jet_notmuon_nmu_E[s].Write()
                hist_jet_notmuon_mass_E[s].Write()
                hist_jet_notmuon_qgl_E[s].Write()
                hist_lepton_pt_M[s].Write()
                hist_lepton_eta_M[s].Write()
                hist_lepton_pt_E[s].Write()
                hist_lepton_eta_E[s].Write()
                hist_muon_jet_pt_M[s].Write()
                hist_muon_jet_eta_M[s].Write()
                hist_muon_jet_pt_E[s].Write()
                hist_muon_jet_eta_E[s].Write()
                hist_InvM_2jets_M[s].Write()
                hist_InvM_2jets_E[s].Write()
                hist_InvM_jetM_lepM[s].Write()
                hist_InvM_jetM_lepE[s].Write()
                hist_deltaR_jetM_lepM[s].Write()
                hist_deltaR_jetM_lepE[s].Write()
                hist_deltaR_jetM_jetNM_M[s].Write()
                hist_deltaR_jetM_jetNM_E[s].Write()
                hist_deltaphi_jetM_jetNM_M[s].Write()
                hist_deltaphi_jetM_jetNM_E[s].Write()
                hist_deltaeta_jetM_jetNM_M[s].Write()
                hist_deltaeta_jetM_jetNM_E[s].Write()
                hist_deltapt_jetM_jetNM_M[s].Write()
                hist_deltapt_jetM_jetNM_E[s].Write()
                hist_MET_M[s].Write()
                hist_MET_E[s].Write()
                hist_tranverse_massM[s].Write()
                hist_tranverse_massE[s].Write()
                hist_tracks_jetM_M[s].Write()
                hist_tracks_jetNM_M[s].Write()
                hist_tracks_jetM_E[s].Write()
                hist_tracks_jetNM_E[s].Write()
                hist_EMN_jetM_M[s].Write()
                hist_EMC_jetM_M[s].Write()
                hist_EMN_jetM_E[s].Write()
                hist_EMC_jetM_E[s].Write()
                hist_EMtotal_jetM_M[s].Write()
                hist_EMtotal_jetM_E[s].Write()
                hist_InvM_muon_jet_M[s].Write()
                hist_muon_jet_mva_M[s].Write()
                hist_muon_jet_mva_E[s].Write()
                hist_muon_jet_tight_M[s].Write()
                hist_muon_jet_tight_E[s].Write()
                hist_muon_jet_relpt_M[s].Write()
                hist_muon_jet_relpt_E[s].Write()
                hist_muon_jet_sigxy_M[s].Write()
                hist_muon_jet_sigxy_E[s].Write()
                hist_muon_jet_sigz_M[s].Write()
                hist_muon_jet_sigz_E[s].Write()
                hist_muon_jet_sigr_M[s].Write()
                hist_muon_jet_sigr_E[s].Write()
                hist_SSOS_M[s].Write()
                hist_SSOS_E[s].Write()
                hist_lepton_iso_E[s].Write()
                hist_lepton_mva_E[s].Write()
                hist_pT_sum_M[s].Write()
                hist_pT_sum_E[s].Write()
                hist_pT_product_M[s].Write()
                hist_pT_product_E[s].Write()
                hist_deltaR_lep_2jets_M[s].Write()                                                                   
                hist_deltaR_lep_2jets_E[s].Write()
                hist_deltaphi_MET_2jets_M[s].Write()
                hist_deltaphi_MET_2jets_E[s].Write()
                hist_deltaphi_lephad_M[s].Write()
                hist_deltaphi_lephad_E[s].Write()
                hist_eta_2jets_M[s].Write()
                hist_eta_2jets_E[s].Write()
                hist_pt_2jets_M[s].Write()
                hist_pt_2jets_E[s].Write()
                hist_pt_Wlep_M[s].Write()
                hist_pt_Wlep_E[s].Write()
                hist_deltaR_lephad_M[s].Write()
                hist_deltaR_lephad_E[s].Write()
                hist_deltaphi_lep_2jets_M[s].Write()
                hist_deltaphi_lep_2jets_E[s].Write()
                hist_deltaeta_lephad_M[s].Write()
                hist_deltaeta_lephad_E[s].Write()
                hist_deltaeta_lep_2jets_M[s].Write()
                hist_deltaeta_lep_2jets_E[s].Write()
                hist_deltapt_lephad_M[s].Write()
                hist_deltapt_lephad_E[s].Write()
                hist_deltapt_lep_2jets_M[s].Write()
                hist_deltapt_lep_2jets_E[s].Write()
                hist_deltaR_lep_jet1_M[s].Write()
                hist_deltaR_lep_jet2_M[s].Write()
                hist_deltaPhi_lep_jet1_M[s].Write()
                hist_deltaPhi_lep_jet2_M[s].Write()
                hist_deltaEta_lep_jet1_M[s].Write()
                hist_deltaEta_lep_jet2_M[s].Write()
                hist_deltaR_lep_jet1_E[s].Write()
                hist_deltaR_lep_jet2_E[s].Write()
                hist_deltaPhi_lep_jet1_E[s].Write()
                hist_deltaPhi_lep_jet2_E[s].Write()
                hist_deltaEta_lep_jet1_E[s].Write()
                hist_deltaEta_lep_jet2_E[s].Write()
                hist_jet_muon_btag_M[s].Write()
                hist_jet_muon_btag_E[s].Write()
                hist_jet_notmuon_btag_M[s].Write()
                hist_jet_notmuon_btag_E[s].Write()
                hist_jet_notmuon_deeptagG_M[s].Write()
                hist_jet_notmuon_deeptagG_E[s].Write()
                hist_jet_notmuon_deeptagC_M[s].Write()
                hist_jet_notmuon_deeptagC_E[s].Write()
                hist_Deta_InvM_M[s].Write()
                hist_Deta_InvM_E[s].Write()
                hist_Deta_qgl_M[s].Write()
                hist_Deta_qgl_E[s].Write()
                hist_DR_InvM_M[s].Write()
                hist_DR_InvM_E[s].Write()
                hist_qgl_InvM_M[s].Write()
                hist_qgl_InvM_E[s].Write()
                hist_pT_proy_M[s].Write()
                hist_pT_proy_E[s].Write()
                hist_pT_sum_2J_M[s].Write()
                hist_pT_sum_2J_E[s].Write()

                myfile.Close()

for s in samples_test:
	print("Expected events for %s sample are %s" %(s,event_test[s]))
	print("Number of events analyzed for %s sample are %s" %(s,df_test[s].Count().GetValue()))
	if mode == "mc": print("Sum of weights of events for %s sample are %s" %(s,df_test[s].Sum("weight_aux").GetValue()))

if args.ssos: print('SSOS version')

print('Ended succesfully')

