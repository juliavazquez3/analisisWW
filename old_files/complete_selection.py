import ROOT, os, sys
from ROOT import *

import numpy as np
import awkward as ak
import uproot


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

### signal

fileName1 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/*.root"
fileName2 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/410000/*.root"
fileName3 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/40000/*.root"

### background

# W jets

fileNameWJ = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/260000/*.root"

# TTbar

fileNameTtbar1 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/00000/*.root"
fileNameTtbar2 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/10000/*.root"
fileNameTtbar3 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/230000/*.root"
fileNameTtbar4 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/240000/*.root"
fileNameTtbar5 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/250000/*.root"
fileNameTtbar6 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/50000/*.root"

treeName = "Events"

# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.

samples = ['signal','Wjets','ttbar']

d = {}

d['signal'] = ROOT.RDataFrame(treeName, {fileName1,fileName2, fileName3})
d['Wjets'] = ROOT.RDataFrame(treeName, fileNameWJ)
d['ttbar'] = ROOT.RDataFrame(treeName, {fileNameTtbar1, fileNameTtbar2, fileNameTtbar3, fileNameTtbar4, fileNameTtbar5, fileNameTtbar6})

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

# All input variables for the method present in the tree
my_call = "typeWW(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother)"
ww_df1 = d['signal'].Define("typeWW", my_call)


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

# All input variables for the method present in the tree
my_call = "typeC(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother)"
d['signal'] = ww_df1.Define("typeC", my_call)

#entries_semi = d.Filter("typeWW == 2").Count()

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

df5 = {}

for s in samples:
	d[s] = d[s].Define('MuonGoodInd','muonInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all)')
	d[s] = d[s].Define('ElectronGoodInd','elInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all)')
	d[s] = d[s].Define('nMuonGood','MuonGoodInd.size()')
	d[s] = d[s].Define('nElectronGood','ElectronGoodInd.size()')
	d[s] = d[s].Define('JetGoodInd','JetInd(nJet, Jet_pt, Jet_eta, Jet_phi, MuonGoodInd, Muon_eta, Muon_phi, ElectronGoodInd, Electron_eta, Electron_phi)')
	d[s] = d[s].Define('nJetGood','JetGoodInd.size()')
	d[s] = d[s].Define('MuonJetInd','JetMuonInd(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, ElectronGoodInd)')
	d[s] = d[s].Define('JetMuonInd','JetMuonIndJet(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, ElectronGoodInd)')
	d[s] = d[s].Define('MuonJetGood','muonCond( MuonJetInd, Muon_pt, Muon_eta, Muon_pfRelIso04_all)')
	df5[s] = d[s].Define('MuonLepSign','SSOS(MuonGoodInd, ElectronGoodInd, MuonJetInd, Electron_charge, Muon_charge)')



#################     FILTERS     #######################

###########    Extra definitions

df = df5['ttbar'].Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>0').Filter('MuonJetInd.size() == 1').Filter('MuonJetGood')

df_muon = df.Filter('nMuonGood>0')
df_electron = df.Filter('nElectronGood >0')

df = df.Define('JetnotMuonInd','nJetGood>1 ? (JetMuonInd[0] == JetGoodInd[0] ? JetGoodInd[1] : JetGoodInd[0]) : -1')

### hists definitions

df = df.Define('jet_muon_pt','Jet_pt[JetMuonInd[0]]')
df = df.Define('jet_notmuon_pt','nJetGood>1 ? Jet_pt[JetnotMuonInd] : 0')

df = df.Define('jet_muon_eta','Jet_eta[JetMuonInd[0]]')
df = df.Define('jet_notmuon_eta','nJetGood>1 ? Jet_eta[JetnotMuonInd] : 0')

df_muon = df_muon.Define('lepton_pt','Muon_pt[MuonGoodInd[0]]')
df_muon = df_muon.Define('lepton_eta','Muon_eta[MuonGoodInd[0]]')

df_electron = df_electron.Define('lepton_pt','Electron_pt[ElectronGoodInd[0]]')
df_electron = df_electron.Define('lepton_eta','Electron_eta[ElectronGoodInd[0]]')

df = df.Define('muon_jet_pt','Muon_pt[MuonJetInd[0]]')
df = df.Define('muon_jet_eta','Muon_eta[MuonJetInd[0]]')

df = df.Define('InvM_2jets','nJetGood>1 ? InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Jet_pt[JetnotMuonInd],Jet_eta[JetnotMuonInd],Jet_phi[JetnotMuonInd],Jet_mass[JetnotMuonInd]) : 0')

df_muon = df_muon.Define('deltaR_jetM_lep','ROOT::VecOps::DeltaR(Muon_eta[MuonGoodInd[0]],Jet_eta[JetMuonInd[0]] , Muon_phi[MuonGoodInd[0]], Jet_phi[JetMuonInd[0]])')
df_electron = df_electron.Define('deltaR_jetM_lep','ROOT::VecOps::DeltaR(Electron_eta[ElectronGoodInd[0]],Jet_eta[JetMuonInd[0]] , Electron_phi[ElectronGoodInd[0]], Jet_phi[JetMuonInd[0]])')

df = df.Define('deltaR_jetM_jetNM','nJetGood>1 ? ROOT::VecOps::DeltaR(Jet_eta[JetnotMuonInd], Jet_eta[JetMuonInd[0]] , Jet_phi[JetnotMuonInd], Jet_phi[JetMuonInd[0]])  : 10')

df_muon = df_muon.Define('InvM_jetM_lep', 'InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Muon_pt[MuonGoodInd[0]],Muon_eta[MuonGoodInd[0]],Muon_phi[MuonGoodInd[0]],Muon_mass[MuonGoodInd[0]])')
df_electron = df_electron.Define('InvM_jetM_lep', 'InvariantM(Jet_pt[JetMuonInd[0]],Jet_eta[JetMuonInd[0]],Jet_phi[JetMuonInd[0]],Jet_mass[JetMuonInd[0]],Electron_pt[ElectronGoodInd[0]],Electron_eta[ElectronGoodInd[0]],Electron_phi[ElectronGoodInd[0]], Electron_mass[ElectronGoodInd[0]])')

df_muon = df_muon.Define('transverse_mass','std::sqrt(2*Muon_pt[MuonGoodInd[0]]*MET_pt*(1-std::cos(Muon_phi[MuonGoodInd[0]]-MET_phi)))')
df_electron = df_electron.Define('transverse_mass','std::sqrt(2*Electron_pt[ElectronGoodInd[0]]*MET_pt*(1-std::cos(Electron_phi[ElectronGoodInd[0]]-MET_phi)))')

df = df.Define('tracks_jetM','Jet_nConstituents[JetMuonInd[0]]')
df = df.Define('tracks_jetNM','nJetGood>1 ? Jet_nConstituents[JetnotMuonInd] : 0')

############################################################
####################     HISTS    ##########################
############################################################

hist_nJetGood = df.Histo1D(("nJetGood","",10,0,10),"nJetGood")

hist_jet_muon_pt = df.Histo1D(("jet muon pt","",100,0,500),"jet_muon_pt")
hist_jet_notmuon_pt = df.Histo1D(("jet not muon pt","",100,0,500),"jet_notmuon_pt")
hist_jet_muon_eta = df.Histo1D(("jet muon eta","",80,-4,4),"jet_muon_eta")
hist_jet_notmuon_eta = df.Histo1D(("jet not muon eta","",80,-4,4),"jet_notmuon_eta")

hist_leptonM_pt = df_muon.Histo1D(("leptonM pt","",100,0,500),"lepton_pt")
hist_leptonM_eta = df_muon.Histo1D(("leptonM eta","",80,-4,4),"lepton_eta")
hist_leptonE_pt = df_electron.Histo1D(("leptonE pt","",100,0,500),"lepton_pt")
hist_leptonE_eta = df_electron.Histo1D(("leptonE eta","",80,-4,4),"lepton_eta")

hist_muon_jet_pt = df.Histo1D(("muon jet pt","",100,0,50),"muon_jet_pt")
hist_muon_jet_eta = df.Histo1D(("muon jet eta","",80,-4,4),"muon_jet_eta")

hist_InvM_2jets = df.Histo1D(("InvM 2jets","",100,0,500),"InvM_2jets")
hist_InvM_jetM_lepM = df_muon.Histo1D(("InvM jetM lepM","",100,0,500),"InvM_jetM_lep")
hist_InvM_jetM_lepE = df_electron.Histo1D(("InvM jetM lepE","",100,0,500),"InvM_jetM_lep")

hist_MET = df.Histo1D(("MET pt","",100,0,500),"MET_pt")

hist_tranverse_massM = df_muon.Histo1D(("transverse massM","",100,0,500),"transverse_mass")
hist_tranverse_massE = df_electron.Histo1D(("tranverse massE","",100,0,500),"transverse_mass")

hist_tracks_jetM = df.Histo1D(("tracks jetM","",60,0,60),"tracks_jetM")
hist_tracks_jetNM = df.Histo1D(("tracks jetNM","",60,0,60),"tracks_jetNM")


#############################
####     DATA SAVING     ####
#############################


myfile = TFile( 'analisisWW/hists/hists_ttbar.root', 'RECREATE' )

hist_nJetGood.Write()
hist_jet_muon_pt.Write()
hist_jet_muon_eta.Write()
hist_jet_notmuon_pt.Write()
hist_jet_notmuon_eta.Write()
hist_leptonM_pt.Write()
hist_leptonM_eta.Write()
hist_leptonE_pt.Write()
hist_leptonE_eta.Write()
hist_muon_jet_pt.Write()
hist_muon_jet_eta.Write()
hist_InvM_2jets.Write()
hist_InvM_jetM_lepM.Write()
hist_InvM_jetM_lepE.Write()
hist_MET.Write()
hist_tranverse_massM.Write()
hist_tranverse_massE.Write()
hist_tracks_jetM.Write()
hist_tracks_jetNM.Write()



myfile.Close()
