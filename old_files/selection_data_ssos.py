###################################                                        
###################################                                        
#######     DATA VERSION    #######
###################################                                        
###################################                                        

## Modified version, meant to be run locally, not sent to condor 

## Selection for data samples, meant to distinguish between SingleMuon or SingeleElectron datasets
## Different histogram files are produced for each situation 

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
parser.add_argument("--process", type=string, default="M",
                    help="SingleMuon or SingleElectron")
parser.add_argument("--year", type=string, default="2016",
                    help="Select data year to run")
parser.add_argument("--notfull", action="store_true", default=False,
                    help="Use just a range of the sample")
parser.add_argument('-l','--list', nargs='+', help='range of sample to use')
# Use like:
# python arg.py -l 1234 2345 3456 4567

args = parser.parse_args()

if args.process == "all": proc = ["M","E"]
elif (args.process == "M" or args.process == "E"): proc = [str(args.process)]
else: raise NameError('Incorrect process name')

if args.year == "all": years = ["2016","2016B","2017","2018"]
elif (args.year == "2016" or args.year == "2016B" or args.year == "2017" or args.year == "2018"): years = [str(args.year)]
else: raise NameError('Incorrect year')

samples = []

for p in proc:
  for y in years:
    samples.append(y+p)

print("the samples treated are",samples)

if args.notfull:
	if len(args.list) != 2: raise NameError('List has to have 2 elements')
	print("the range of files is",args.list)
 

# Create a ROOT dataframe for each dataset
# Note that we load the filenames from the external json file placed in the same folder than this script.
# Example "python analisisWW/selection_v2.py --process="WW" --notfull -l 0 50"

files = json.load(open("/nfs/cms/vazqueze/analisisWW/data_info_v9.json"))
processes = files.keys()

df = {}
xsecs = {}
sumws = {}
archives = {}

for s in samples:
  archives[s]=[]

for p in processes:
    # Construct the dataframes
    folder = files[p]["folder_name"] # Folder name
    filename = files[p]["sample"] # Sample name
    num_events = files[p]["events"] # Number of events
    num_files = files[p]["files"] # Number of files
    list_files = [f for f in listdir(folder) if isfile(join(folder, f))] # Lista de archivos
    #print(len(list_files))
    for s in samples:
      if files[p]["type"]==s:
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
      auto muonInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vbool tID) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
		if (pt[i]>30. && fabs(eta[i])<2.4 && iso[i]<0.20 && tID[i]){
                	vb.push_back(i);
		}
            }
            return vb;
      };
      auto elInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso, Vint cutB, Vbool mva80 , Vbool mva90) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
                if (pt[i]>30. && fabs(eta[i])<2.5 && iso[i]<0.15 && mva80[i]){
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
	    bool cond1 = false;
	    bool cond = false;
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

#############################################################
######################## Definitions ########################
#############################################################

df_muon = {}
df_electron = {}

for s in samples:
	df[s] = df[s].Define('MuonGoodInd','muonInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all, Muon_tightId)')
	df[s] = df[s].Define('ElectronGoodInd','elInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all, Electron_cutBased, Electron_mvaFall17V2noIso_WP80, Electron_mvaFall17V2noIso_WP90)')
	df[s] = df[s].Define('nMuonGood','MuonGoodInd.size()')
	df[s] = df[s].Define('nElectronGood','ElectronGoodInd.size()')
	df[s] = df[s].Define('JetGoodInd','JetInd(nJet, Jet_pt, Jet_eta, Jet_phi, MuonGoodInd, Muon_eta, Muon_phi, ElectronGoodInd, Electron_eta, Electron_phi, Jet_puId, Jet_jetId)')
	df[s] = df[s].Define('nJetGood','JetGoodInd.size()')
	df[s] = df[s].Define('MuonJetInd','JetMuonInd(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx)')
	df[s] = df[s].Define('JetMuonInd','JetMuonIndJet(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, Muon_pfRelIso04_all, ElectronGoodInd, Muon_jetIdx)')
	df[s] = df[s].Define('MuonJetGood','muonCond( MuonJetInd, Muon_pt, Muon_eta, Muon_pfRelIso04_all, Muon_softId)')
	df[s] = df[s].Define('MuonLepSign','SSOS(MuonGoodInd, ElectronGoodInd, MuonJetInd, Electron_charge, Muon_charge)')
	df[s] = df[s].Define('weightSSOS','-1*MuonLepSign')

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
	df[s] = df[s].Define('jet_notmuon_pt','nJetGood>1 ? Jet_pt[JetnotMuonInd] : 0')
	df[s] = df[s].Define('jet_muon_eta','Jet_eta[JetMuonInd[0]]')
	df[s] = df[s].Define('jet_notmuon_eta','nJetGood>1 ? Jet_eta[JetnotMuonInd] : 0')
	df[s] = df[s].Define('muon_jet_pt','Muon_pt[MuonJetInd[0]]')
	df[s] = df[s].Define('muon_jet_eta','Muon_eta[MuonJetInd[0]]')
	df[s] = df[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Jet_pt[JetnotMuonInd],Jet_eta[JetnotMuonInd],Jet_phi[JetnotMuonInd],Jet_mass[JetnotMuonInd]) : 0')
	df[s] = df[s].Define('deltaR_jetM_jetNM','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotMuonInd], Jet_eta[JetMuonInd[0]] , Jet_phi[JetnotMuonInd], Jet_phi[JetMuonInd[0]])  : 10')
	df[s] = df[s].Define('tracks_jetM','Jet_nConstituents[JetMuonInd[0]]')
	df[s] = df[s].Define('tracks_jetNM','nJetGood>1 ? Jet_nConstituents[JetnotMuonInd] : 0')
	df[s] = df[s].Define('EMN_jetM','Jet_neEmEF[JetMuonInd[0]]')
	df[s] = df[s].Define('EMC_jetM','Jet_chEmEF[JetMuonInd[0]]')
	df[s] = df[s].Define('EMtotal_jetM','Jet_chEmEF[JetMuonInd[0]]+Jet_neEmEF[JetMuonInd[0]]')
	df[s] = df[s].Define('muon_jet_mva','Muon_softMva[MuonJetInd[0]]')
	df[s] = df[s].Define('muon_jet_relpt','Muon_pt[MuonJetInd[0]]/Jet_pt[JetMuonInd[0]]')
	df[s] = df[s].Define('jet_muon_btag','Jet_btagDeepFlavB[JetMuonInd[0]]')
	df[s] = df[s].Define('jet_notmuon_btag','Jet_btagDeepFlavB[JetnotMuonInd]')

## Differentiated definitions
for s in samples:
	df_muon[s] = df[s].Filter('nMuonGood>0')
	df_electron[s] = df[s].Filter('nElectronGood >0')
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
	## jet pt
	#df_electron[s] = df_electron[s].Filter('jet_muon_pt>35.')
	#df_muon[s] = df_muon[s].Filter('jet_muon_pt>35.')
	#df_electron[s] = df_electron[s].Filter('jet_notmuon_pt>35.')
	#df_muon[s] = df_muon[s].Filter('jet_notmuon_pt>35.')

############################################################
####################     HISTS    ##########################
############################################################

hist_nJetGood_M = {}
hist_nJetGood_E = {}
hist_jet_muon_pt_M = {}
hist_jet_notmuon_pt_M = {}
hist_jet_muon_eta_M = {}
hist_jet_notmuon_eta_M = {}
hist_jet_muon_pt_E = {}
hist_jet_notmuon_pt_E = {}
hist_jet_muon_eta_E = {}
hist_jet_notmuon_eta_E = {}
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
hist_muon_jet_relpt_M = {}
hist_muon_jet_relpt_E = {}                       
hist_SSOS_M = {}
hist_SSOS_E = {}                          
hist_jet_muon_btag_M = {}
hist_jet_muon_btag_E = {}
hist_jet_notmuon_btag_M = {}
hist_jet_notmuon_btag_E = {}
                          
for s in samples:
	hist_nJetGood_M[s] = df_muon[s].Histo1D(("nJetGood M","",10,0,10),"nJetGood","weightSSOS")
	hist_nJetGood_E[s] = df_electron[s].Histo1D(("nJetGood E","",10,0,10),"nJetGood","weightSSOS")

	hist_jet_muon_pt_M[s] = df_muon[s].Histo1D(("jet muon pt M","",50,20,120),"jet_muon_pt","weightSSOS")
	hist_jet_notmuon_pt_M[s] = df_muon[s].Histo1D(("jet not muon pt M","",50,20,120),"jet_notmuon_pt","weightSSOS")
	hist_jet_muon_eta_M[s] = df_muon[s].Histo1D(("jet muon eta M","",80,-4,4),"jet_muon_eta","weightSSOS")
	hist_jet_notmuon_eta_M[s] = df_muon[s].Histo1D(("jet not muon eta M","",80,-4,4),"jet_notmuon_eta","weightSSOS")

	hist_jet_muon_pt_E[s] = df_electron[s].Histo1D(("jet muon pt E","",50,20,120),"jet_muon_pt","weightSSOS")
	hist_jet_notmuon_pt_E[s] = df_electron[s].Histo1D(("jet not muon pt E","",50,20,120),"jet_notmuon_pt","weightSSOS")
	hist_jet_muon_eta_E[s] = df_electron[s].Histo1D(("jet muon eta E","",80,-4,4),"jet_muon_eta","weightSSOS")
	hist_jet_notmuon_eta_E[s] = df_electron[s].Histo1D(("jet not muon eta E","",80,-4,4),"jet_notmuon_eta","weightSSOS")

	hist_lepton_pt_M[s] = df_muon[s].Histo1D(("lepton pt M","",50,20,120),"lepton_pt","weightSSOS")	
	hist_lepton_eta_M[s] = df_muon[s].Histo1D(("lepton eta M","",80,-4,4),"lepton_eta","weightSSOS")

	hist_lepton_pt_E[s] = df_electron[s].Histo1D(("lepton pt E","",50,20,120),"lepton_pt","weightSSOS")
	hist_lepton_eta_E[s] = df_electron[s].Histo1D(("lepton eta E","",80,-4,4),"lepton_eta","weightSSOS")

	hist_muon_jet_pt_M[s] = df_muon[s].Histo1D(("muon jet pt M","",50,0,25),"muon_jet_pt","weightSSOS")
	hist_muon_jet_eta_M[s] = df_muon[s].Histo1D(("muon jet eta M","",80,-4,4),"muon_jet_eta","weightSSOS")

	hist_muon_jet_pt_E[s] = df_electron[s].Histo1D(("muon jet pt E","",50,0,25),"muon_jet_pt","weightSSOS")
	hist_muon_jet_eta_E[s] = df_electron[s].Histo1D(("muon jet eta E","",80,-4,4),"muon_jet_eta","weightSSOS")

	hist_InvM_2jets_M[s] = df_muon[s].Histo1D(("InvM 2jets M","",100,0,300),"InvM_2jets","weightSSOS")
	hist_InvM_2jets_E[s] = df_electron[s].Histo1D(("InvM 2jets E","",100,0,300),"InvM_2jets","weightSSOS")
	hist_InvM_jetM_lepM[s] = df_muon[s].Histo1D(("InvM jetM lepM","",100,0,300),"InvM_jetM_lep","weightSSOS")
	hist_InvM_jetM_lepE[s] = df_electron[s].Histo1D(("InvM jetM lepE","",100,0,300),"InvM_jetM_lep","weightSSOS")

	hist_MET_M[s] = df_muon[s].Histo1D(("MET pt M","",100,0,150),"MET_pt","weightSSOS")
	hist_MET_E[s] = df_electron[s].Histo1D(("MET pt E","",100,0,150),"MET_pt","weightSSOS")

	hist_deltaR_jetM_lepM[s] = df_muon[s].Histo1D(("deltaR jetM lepM","",100,0,5),"deltaR_jetM_lepM","weightSSOS")
	hist_deltaR_jetM_lepE[s] = df_electron[s].Histo1D(("deltaR jetM lepE","",100,0,5),"deltaR_jetM_lepE","weightSSOS")
	hist_deltaR_jetM_jetNM_M[s] = df_muon[s].Histo1D(("deltaR jetM jetNM M","",100,0,5),"deltaR_jetM_jetNM","weightSSOS")
	hist_deltaR_jetM_jetNM_E[s] = df_electron[s].Histo1D(("deltaR jetM jetNM E","",100,0,5),"deltaR_jetM_jetNM","weightSSOS")

	hist_tranverse_massM[s] = df_muon[s].Histo1D(("transverse massM","",100,0,150),"transverse_mass","weightSSOS")
	hist_tranverse_massE[s] = df_electron[s].Histo1D(("tranverse massE","",100,0,150),"transverse_mass","weightSSOS")

	hist_tracks_jetM_M[s] = df_muon[s].Histo1D(("tracks jetM M","",60,0,60),"tracks_jetM","weightSSOS")
	hist_tracks_jetNM_M[s] = df_muon[s].Histo1D(("tracks jetNM M","",60,0,60),"tracks_jetNM","weightSSOS")

	hist_tracks_jetM_E[s] = df_electron[s].Histo1D(("tracks jetM E","",60,0,60),"tracks_jetM","weightSSOS")
	hist_tracks_jetNM_E[s] = df_electron[s].Histo1D(("tracks jetNM E","",60,0,60),"tracks_jetNM","weightSSOS")

	hist_EMN_jetM_M[s] = df_muon[s].Histo1D(("EMN jetM M","",60,0,1),"EMN_jetM","weightSSOS")
	hist_EMC_jetM_M[s] = df_muon[s].Histo1D(("EMC jetM M","",60,0,1),"EMC_jetM","weightSSOS")
	hist_EMN_jetM_E[s] = df_electron[s].Histo1D(("EMN jetM E","",60,0,1),"EMN_jetM","weightSSOS")
	hist_EMC_jetM_E[s] = df_electron[s].Histo1D(("EMC jetM E","",60,0,1),"EMC_jetM","weightSSOS")

	hist_EMtotal_jetM_M[s] = df_muon[s].Histo1D(("EMtotal jetM M","",60,0,1),"EMtotal_jetM","weightSSOS")
	hist_EMtotal_jetM_E[s] = df_muon[s].Histo1D(("EMtotal jetM E","",60,0,1),"EMtotal_jetM","weightSSOS")

	hist_InvM_muon_jet_M[s] = df_muon[s].Histo1D(("InvM muon jet M","",50,0,200),"InvM_muon_jet","weightSSOS")

	hist_muon_jet_mva_M[s] = df_muon[s].Histo1D(("muon jet mva M","",50,0,1),"muon_jet_mva","weightSSOS")
	hist_muon_jet_mva_E[s] = df_electron[s].Histo1D(("muon jet mva E","",50,0,1),"muon_jet_mva","weightSSOS")

	hist_muon_jet_relpt_M[s] = df_muon[s].Histo1D(("muon jet relpt M","",50,0,1),"muon_jet_relpt","weightSSOS")
	hist_muon_jet_relpt_E[s] = df_electron[s].Histo1D(("muon jet relpt E","",50,0,1),"muon_jet_relpt","weightSSOS")                                                                                             

	hist_SSOS_M[s] = df_muon[s].Histo1D(("SSOS M","",4,-2,2),"MuonLepSign","weightSSOS")
	hist_SSOS_E[s] = df_electron[s].Histo1D(("SSOS E","",4,-2,2),"MuonLepSign","weightSSOS")

	hist_jet_muon_btag_M[s] = df_muon[s].Histo1D(("jet muon btag M","",50,0,1),"jet_muon_btag","weightSSOS")
	hist_jet_muon_btag_E[s] = df_electron[s].Histo1D(("jet muon btag E","",50,0,1),"jet_muon_btag","weightSSOS")
	hist_jet_notmuon_btag_M[s] = df_muon[s].Histo1D(("jet notmuon btag M","",50,0,1),"jet_notmuon_btag","weightSSOS")
	hist_jet_notmuon_btag_E[s] = df_electron[s].Histo1D(("jet notmuon btag E","",50,0,1),"jet_notmuon_btag","weightSSOS")

#############################
####     DATA SAVING     ####
#############################

if args.notfull:
	for s in samples:
		path_hist = '/nfs/cms/vazqueze/analisisWW/hists/ssos/ver9_v1v2v3v4SSOS_'+s+'_range_'+str(args.list[0])+'_'+str(args.list[1])+'.root'
		myfile = TFile( path_hist, 'RECREATE' )

		hist_nJetGood_M[s].Write()
		hist_nJetGood_E[s].Write()
		hist_jet_muon_pt_M[s].Write()
		hist_jet_muon_eta_M[s].Write()
		hist_jet_muon_pt_E[s].Write()
		hist_jet_muon_eta_E[s].Write()
		hist_jet_notmuon_pt_M[s].Write()
		hist_jet_notmuon_eta_M[s].Write()
		hist_jet_notmuon_pt_E[s].Write()
		hist_jet_notmuon_eta_E[s].Write()
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
		hist_muon_jet_relpt_M[s].Write()
		hist_muon_jet_relpt_E[s].Write()                                  
		hist_SSOS_M[s].Write()
		hist_SSOS_E[s].Write()
		hist_jet_muon_btag_M[s].Write()
		hist_jet_muon_btag_E[s].Write()
		hist_jet_notmuon_btag_M[s].Write()
		hist_jet_notmuon_btag_E[s].Write()

		myfile.Close()

else:
	for s in samples:
		path_hist = '/nfs/cms/vazqueze/analisisWW/hists/ssos/ver9_v1v2v3v4SSOS_'+s+'.root'
		myfile = TFile( path_hist, 'RECREATE' )

		hist_nJetGood_M[s].Write()
		hist_nJetGood_E[s].Write()
		hist_jet_muon_pt_M[s].Write()
		hist_jet_muon_eta_M[s].Write()
		hist_jet_muon_pt_E[s].Write()
		hist_jet_muon_eta_E[s].Write()
		hist_jet_notmuon_pt_M[s].Write()
		hist_jet_notmuon_eta_M[s].Write()
		hist_jet_notmuon_pt_E[s].Write()
		hist_jet_notmuon_eta_E[s].Write()
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
		hist_muon_jet_relpt_M[s].Write()
		hist_muon_jet_relpt_E[s].Write()                                  
		hist_SSOS_M[s].Write()
		hist_SSOS_E[s].Write()
		hist_jet_muon_btag_M[s].Write()
		hist_jet_muon_btag_E[s].Write()
		hist_jet_notmuon_btag_M[s].Write()
		hist_jet_notmuon_btag_E[s].Write()

		myfile.Close()

print('Ended succesfully')
