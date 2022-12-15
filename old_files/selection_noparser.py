import ROOT, os, sys
from ROOT import *

import numpy as np
import awkward as ak
import uproot

import json

#if not sys.flags.interactive: ROOT.EnableImplicitMT() 

plotdir = 'analisisWW/plots/' # this is where we'll save your plots                                                                                                                                 
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

ROOT.EnableImplicitMT()

# Create a ROOT dataframe for each dataset
# Note that we load the filenames from the external json file placed in the same folder than this script.

files = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "example.json")))
processes = files.keys()

print(processes)

df = {}
xsecs = {}
sumws = {}

files_WW = set()
files_Wjets = set()
files_ttbar = set()

for p in processes:
    # Construct the dataframes
    folder = files[p]["folder_name"] # Folder name
    filename = files[p]["sample"] # Sample name
    num_events = files[p]["events"] # Number of events
    if files[p]["type"]=="WW":
        files_WW.add(folder+filename)
    if files[p]["type"]=="Wjets":
       	files_Wjets.add(folder+filename)
    if files[p]["type"]=="ttbar":
       	files_ttbar.add(folder+filename)

samples = ["WW","Wjets","ttbar"]

df["WW"] = ROOT.RDataFrame("Events",files_WW)
df["Wjets"] = ROOT.RDataFrame("Events",files_Wjets)
df["ttbar"] = ROOT.RDataFrame("Events",files_ttbar)

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
      auto muonInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
		if (pt[i]>30. && fabs(eta[i])<2.5 && iso[i]<0.20){
                	vb.push_back(i);
		}
            }
            return vb;
      };
      auto elInd(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso) {
            vector<int> vb;
            for (unsigned int i=0; i<nmu; ++i) {
                if (pt[i]>30. && fabs(eta[i])<2.5 && iso[i]<0.15){
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
      auto JetInd(UInt_t njet, Vfloat pt, Vfloat eta, Vfloat phi, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vint el_good, Vfloat el_eta, Vfloat el_phi) {
            vector<int> vb;
	    bool cond = false;
            for (unsigned int i=0; i<njet; ++i) {
		if(mu_good.size()>0){
			cond = ROOT::VecOps::DeltaR(mu_eta[mu_good[0]],eta[i],mu_phi[mu_good[0]],phi[i]) > 0.4;
                        if (pt[i]>30. && fabs(eta[i])<2.5 && cond){
                                vb.push_back(i);
                        }
		}
                if(el_good.size()>0){
                        cond = ROOT::VecOps::DeltaR(el_eta[el_good[0]],eta[i],el_phi[el_good[0]],phi[i]) > 0.4;
                        if (pt[i]>30. && fabs(eta[i])<2.5 && cond){
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
      auto JetMuonIndJet(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_pt, Vint el_good) {
            vector<int> vb;
            bool cond = false;
            bool cond1 = true;
            float ptM{-10.};
            float ptJ{-10.};
            if (good.size() == 1){
                for (unsigned int i=0; i<nmu; ++i){
                        cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[0]],mu_phi[i],phi[good[0]]) < 0.4;
                        if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                        if(cond && mu_pt[i] > ptM && cond1){
                                vb.push_back(0);
                                ptM = mu_pt[i];
                        }
                }
            }
            if (good.size() > 1){
                for (unsigned int j=0; j<good.size(); ++j){
                        for (unsigned int i=0; i<nmu; ++i){
                                cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[j]],mu_phi[i],phi[good[j]]) < 0.4;
                                if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                                if(cond && mu_pt[i] > ptM && pt[good[j]] > ptJ && cond1){
                                        vb.push_back(j);
                                        ptM = mu_pt[i];
                                        ptJ = pt[good[j]];
                                }
                        }
                }
            }
            return vb;
      };
      auto JetMuonInd(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_pt, Vint el_good) {
            vector<int> vb;
            bool cond = false;
	    bool cond1 = true;
            float ptM{-10.};
            float ptJ{-10.};
	    if (good.size() == 1){
		for (unsigned int i=0; i<nmu; ++i){
			cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[0]],mu_phi[i],phi[good[0]]) < 0.4;
			if (mu_good.size() > 0) cond1 = mu_good[0] != i;
			if(cond && mu_pt[i] > ptM && cond1){
				vb.push_back(i);
				ptM = mu_pt[i];
			}
		}
	    }
	    if (good.size() > 1){
		for (unsigned int j=0; j<good.size(); ++j){
                	for (unsigned int i=0; i<nmu; ++i){
                        	cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[j]],mu_phi[i],phi[good[j]]) < 0.4;
				if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                        	if(cond && mu_pt[i] > ptM && pt[good[j]] > ptJ && cond1){
                                	vb.push_back(i);
                                	ptM = mu_pt[i];
                                        ptJ = pt[good[j]];
                        	}
                	}
		}
	    }
	    return vb;
      };
""")

######### pedimos algunas condiciones al muon seleccionado

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto muonCond( Vint mu_jet,Vfloat mu_pt, Vfloat mu_eta, Vfloat mu_iso) {
	    bool cond = false;
	    if (mu_jet.size()>0){ 
	    	if (mu_pt[mu_jet[0]]<25. && fabs(mu_eta[mu_jet[0]])<2.5 && mu_iso[mu_jet[0]] > 0.2) {
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
	df[s] = df[s].Define('MuonGoodInd','muonInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all)')
	df[s] = df[s].Define('ElectronGoodInd','elInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all)')
	df[s] = df[s].Define('nMuonGood','MuonGoodInd.size()')
	df[s] = df[s].Define('nElectronGood','ElectronGoodInd.size()')
	df[s] = df[s].Define('JetGoodInd','JetInd(nJet, Jet_pt, Jet_eta, Jet_phi, MuonGoodInd, Muon_eta, Muon_phi, ElectronGoodInd, Electron_eta, Electron_phi)')
	df[s] = df[s].Define('nJetGood','JetGoodInd.size()')
	df[s] = df[s].Define('MuonJetInd','JetMuonInd(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, ElectronGoodInd)')
	df[s] = df[s].Define('JetMuonInd','JetMuonIndJet(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, ElectronGoodInd)')
	df[s] = df[s].Define('MuonJetGood','muonCond( MuonJetInd, Muon_pt, Muon_eta, Muon_pfRelIso04_all)')
	df[s] = df[s].Define('MuonLepSign','SSOS(MuonGoodInd, ElectronGoodInd, MuonJetInd, Electron_charge, Muon_charge)')
	


#################     FILTERS     #######################

###########    Extra definitions

for s in samples:
	df[s] = df[s].Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>0').Filter('MuonJetInd.size() == 1').Filter('MuonJetGood')
	df[s] = df[s].Define('JetnotMuonInd','nJetGood>1 ? (JetMuonInd[0] == JetGoodInd[0] ? JetGoodInd[1] : JetGoodInd[0]) : -1')
	df_muon[s] = df[s].Filter('nMuonGood>0')
	df_electron[s] = df[s].Filter('nElectronGood >0')

	### hists definitions
	df[s] = df[s].Define('jet_muon_pt','Jet_pt[JetMuonInd[0]]')
	df[s] = df[s].Define('jet_notmuon_pt','nJetGood>1 ? Jet_pt[JetnotMuonInd] : 0')
	df[s] = df[s].Define('jet_muon_eta','Jet_eta[JetMuonInd[0]]')
	df[s] = df[s].Define('jet_notmuon_eta','nJetGood>1 ? Jet_eta[JetnotMuonInd] : 0')
	df_muon[s] = df_muon[s].Define('lepton_pt','Muon_pt[MuonGoodInd[0]]')
	df_muon[s] = df_muon[s].Define('lepton_eta','Muon_eta[MuonGoodInd[0]]')
	df_electron[s] = df_electron[s].Define('lepton_pt','Electron_pt[ElectronGoodInd[0]]')
	df_electron[s] = df_electron[s].Define('lepton_eta','Electron_eta[ElectronGoodInd[0]]')
	df[s] = df[s].Define('muon_jet_pt','Muon_pt[MuonJetInd[0]]')
	df[s] = df[s].Define('muon_jet_eta','Muon_eta[MuonJetInd[0]]')
	df[s] = df[s].Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Jet_pt[JetnotMuonInd],Jet_eta[JetnotMuonInd],Jet_phi[JetnotMuonInd],Jet_mass[JetnotMuonInd]) : 0')
	df_muon[s] = df_muon[s].Define('deltaR_jetM_lep','ROOT::VecOps::DeltaR(Muon_eta[MuonGoodInd[0]],Jet_eta[JetMuonInd[0]] , Muon_phi[MuonGoodInd[0]], Jet_phi[JetMuonInd[0]])')
	df_electron[s] = df_electron[s].Define('deltaR_jetM_lep','ROOT::VecOps::DeltaR(Electron_eta[ElectronGoodInd[0]],Jet_eta[JetMuonInd[0]] , Electron_phi[ElectronGoodInd[0]], Jet_phi[JetMuonInd[0]])')
	df[s] = df[s].Define('deltaR_jetM_jetNM','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotMuonInd], Jet_eta[JetMuonInd[0]] , Jet_phi[JetnotMuonInd], Jet_phi[JetMuonInd[0]])  : 10')
	df_muon[s] = df_muon[s].Define('InvM_jetM_lep', 'InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]])')
	df_electron[s] = df_electron[s].Define('InvM_jetM_lep', 'InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Electron_pt[ElectronGoodInd[0]],Electron_eta[ElectronGoodInd[0]],Electron_phi[ElectronGoodInd[0]], Electron_mass[ElectronGoodInd[0]])')
	df_muon[s] = df_muon[s].Define('transverse_mass','std::sqrt(2*Muon_pt[MuonGoodInd[0]]*MET_pt*(1-std::cos(Muon_phi[MuonGoodInd[0]]-MET_phi)))')
	df_electron[s] = df_electron[s].Define('transverse_mass','std::sqrt(2*Electron_pt[ElectronGoodInd[0]]*MET_pt*(1-std::cos(Electron_phi[ElectronGoodInd[0]]-MET_phi)))')
	df[s] = df[s].Define('tracks_jetM','Jet_nConstituents[JetMuonInd[0]]')
	df[s] = df[s].Define('tracks_jetNM','nJetGood>1 ? Jet_nConstituents[JetnotMuonInd] : 0')

############################################################
####################     HISTS    ##########################
############################################################

hist_nJetGood = {}
hist_jet_muon_pt = {}
hist_jet_notmuon_pt = {}
hist_jet_muon_eta = {}
hist_jet_notmuon_eta = {}
hist_leptonM_pt = {}
hist_leptonM_eta = {}
hist_leptonE_pt = {}
hist_leptonE_eta = {}
hist_muon_jet_pt = {}
hist_muon_jet_eta = {}
hist_InvM_2jets = {}
hist_InvM_jetM_lepM = {}
hist_InvM_jetM_lepE = {}
hist_MET = {}
hist_tranverse_massM = {}
hist_tranverse_massE = {}
hist_tracks_jetM = {}
hist_tracks_jetNM = {}

for s in samples:
	hist_nJetGood[s] = df[s].Histo1D(("nJetGood","",10,0,10),"nJetGood")

	hist_jet_muon_pt[s] = df[s].Histo1D(("jet muon pt","",100,0,500),"jet_muon_pt")
	hist_jet_notmuon_pt[s] = df[s].Histo1D(("jet not muon pt","",100,0,500),"jet_notmuon_pt")
	hist_jet_muon_eta[s] = df[s].Histo1D(("jet muon eta","",80,-4,4),"jet_muon_eta")
	hist_jet_notmuon_eta[s] = df[s].Histo1D(("jet not muon eta","",80,-4,4),"jet_notmuon_eta")

	hist_leptonM_pt[s] = df_muon[s].Histo1D(("leptonM pt","",100,0,500),"lepton_pt")	
	hist_leptonM_eta[s] = df_muon[s].Histo1D(("leptonM eta","",80,-4,4),"lepton_eta")
	hist_leptonE_pt[s] = df_electron[s].Histo1D(("leptonE pt","",100,0,500),"lepton_pt")
	hist_leptonE_eta[s] = df_electron[s].Histo1D(("leptonE eta","",80,-4,4),"lepton_eta")

	hist_muon_jet_pt[s] = df[s].Histo1D(("muon jet pt","",100,0,50),"muon_jet_pt")
	hist_muon_jet_eta[s] = df[s].Histo1D(("muon jet eta","",80,-4,4),"muon_jet_eta")

	hist_InvM_2jets[s] = df[s].Histo1D(("InvM 2jets","",100,0,500),"InvM_2jets")
	hist_InvM_jetM_lepM[s] = df_muon[s].Histo1D(("InvM jetM lepM","",100,0,500),"InvM_jetM_lep")
	hist_InvM_jetM_lepE[s] = df_electron[s].Histo1D(("InvM jetM lepE","",100,0,500),"InvM_jetM_lep")

	hist_MET[s] = df[s].Histo1D(("MET pt","",100,0,500),"MET_pt")

	hist_tranverse_massM[s] = df_muon[s].Histo1D(("transverse massM","",100,0,500),"transverse_mass")
	hist_tranverse_massE[s] = df_electron[s].Histo1D(("tranverse massE","",100,0,500),"transverse_mass")

	hist_tracks_jetM[s] = df[s].Histo1D(("tracks jetM","",60,0,60),"tracks_jetM")
	hist_tracks_jetNM[s] = df[s].Histo1D(("tracks jetNM","",60,0,60),"tracks_jetNM")


#############################
####     DATA SAVING     ####
#############################


for s in samples:
	path_hist = 'analisisWW/hists/hists_'+s+'.root'
	myfile = TFile( path_hist, 'RECREATE' )

	hist_nJetGood[s].Write()
	hist_jet_muon_pt[s].Write()
	hist_jet_muon_eta[s].Write()
	hist_jet_notmuon_pt[s].Write()
	hist_jet_notmuon_eta[s].Write()
	hist_leptonM_pt[s].Write()
	hist_leptonM_eta[s].Write()
	hist_leptonE_pt[s].Write()
	hist_leptonE_eta[s].Write()
	hist_muon_jet_pt[s].Write()
	hist_muon_jet_eta[s].Write()
	hist_InvM_2jets[s].Write()
	hist_InvM_jetM_lepM[s].Write()
	hist_InvM_jetM_lepE[s].Write()
	hist_MET[s].Write()
	hist_tranverse_massM[s].Write()
	hist_tranverse_massE[s].Write()
	hist_tracks_jetM[s].Write()
	hist_tracks_jetNM[s].Write()

	myfile.Close()
