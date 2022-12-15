import ROOT, os, sys
from ROOT import *

import numpy as np
import awkward as ak
import uproot


#if not sys.flags.interactive: ROOT.EnableImplicitMT() 

plotdir = 'plots/' # this is where we'll save your plots                                                                                                                                 
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


# A simple helper function to fill a test tree: this makes the example stand-alone.
#def fill_tree(treeName, fileName):
 #   df = ROOT.RDataFrame(10)
  #  df.Define("b1", "(double) rdfentry_")\
   #   .Define("b2", "(int) rdfentry_ * rdfentry_").Snapshot(treeName, fileName)

# We prepare an input tree to run on
fileName1 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/FB4B351A-BD0F-CF4A-9FAF-100A2A3A86BC.root"
fileName2 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/FAF7F5F5-8750-8347-A890-F6106D178BB6.root"
fileName3 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/F44E06AC-57B3-2E42-9B5D-6E82C4443387.root"
fileName4 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/E4130C64-2E68-A649-80FC-2161759A5174.root"
fileName5 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/E32B21CA-53DA-814F-9F2B-2149A0E76712.root"
fileName6 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/DA100535-8EE3-7749-BBDC-9F75EBA8D7F5.root"
fileName7 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/D86B5C0C-B2FA-9A45-8CFD-F6BDFE6260EB.root"
fileName8 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/D5C526BF-95A2-F248-BBA6-82AAC42CD06B.root"
fileName9 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/C8277FEB-771C-B34D-B638-335566291F6A.root"
fileName10 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/C28F9BCD-5A01-D645-8759-EA075D16E17D.root"
fileName11 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/A85AC2BA-F1DD-FB4B-ADAB-A9C4C8AE35E6.root"
fileName12 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/A4640085-D831-AF4F-8661-70A108B00584.root"
fileName13 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/A19D9628-DA0B-0B44-ACCC-5597E7F8DF9C.root"
fileName14 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/93D0C189-5EE1-0140-9307-9F90661FFF90.root"
fileName15 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/864423BD-4FA5-FE4E-A67E-44B94331731C.root"
fileName16 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/76541735-6144-DC45-8D92-A68A9A8A51FE.root"
fileName17 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/5925A9F5-0C75-B54F-950A-501D77409970.root"
fileName18 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/5093F917-D106-F342-AA0E-9BCA519D87C0.root"
fileName19 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/3E60A247-EEB8-D54C-98DC-1AE489D38F19.root"
fileName20 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/37C01D05-8328-9A4F-AD4B-1AD09B4288AA.root"
fileName21 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/309C4E04-65DC-9B42-995B-5D45916B2D06.root"
fileName22 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/237E819E-2AFE-A445-98E4-1FC0CB5E1D61.root"
fileName23 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/1D6F3E2F-7A4E-6841-926F-800FDA9C16E7.root"
fileName24 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/13A6B559-A234-3446-97C0-1AD9550E912B.root"
fileName25 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/10DA70C0-DD84-2049-82FB-DCC828A74F5D.root"
fileName26 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/410000/85D3A6CC-3F6A-B946-A843-1853EE2B96D0.root"
fileName27 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/40000/D6C803BD-AF23-BB4F-8D52-A248F2186127.root"
fileName28 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/40000/A1BCCAA8-14FA-0845-BAEE-031B0EDA1D19.root"
fileName29 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/40000/9B8578C7-49D7-8844-B082-8172C1FB534D.root"
fileName30 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/40000/914978FA-D96B-E04C-9D19-3849F81B3EE0.root"
fileName31 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/40000/6F56FDE1-5554-4B49-8B64-A841D42BD8C3.root"
fileName32 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/40000/6DD0AF9E-308B-1C40-B1E3-EDD5BF96DE91.root"
fileName33 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/40000/4A13C66D-3984-7148-BAA2-D01EA8AA2685.root"
fileName34 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/40000/21F66F87-19EE-404F-93C5-17526EA1C0ED.root"
fileName35 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/40000/044FC323-464D-464F-A7B3-D2BEA2EC53C7.root"
fileName36 = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/40000/004E4C32-9FC4-A24C-AF08-374E95DF0B4C.root"


treeName = "Events"
#ROOT.fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, {fileName1,fileName2,fileName3,fileName4,fileName5,fileName6,fileName7,fileName8,fileName9,fileName10,fileName11,fileName12,fileName13,fileName14,fileName15,fileName16,fileName17,fileName18,fileName19,fileName20,fileName21,fileName22,fileName23,fileName24,fileName25,fileName26,fileName27,fileName28,fileName29,fileName30,fileName31,fileName32,fileName33,fileName34,fileName35,fileName36})

########GEN_PART###########

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
ww_df1 = d.Define("typeWW", my_call)
#filtered_good_df = good_df.Filter("nJetGood>0")
#filtered_good_df.Report().Print()


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
ww_df = ww_df1.Define("typeC", my_call)
#filtered_good_df = good_df.Filter("nJetGood>0")
#filtered_good_df.Report().Print()

####### aux filter

gInterpreter.Declare("""
      using VboolR = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using namespace std;
      auto MuonAuxPT(UInt_t nmu,VboolR tId, Vfloat pt, Vfloat eta, Vfloat iso) {
	    float pT = -10.;
            for (unsigned int i=0; i<nmu; ++i) {
                if(pt[i]>pT){
		//if(pt[i]>pT && fabs(eta[i])<2.5 && iso[i]<0.15){
			pT = pt[i];
		}
            }
            return pT;
      };
      auto MuonAuxIso(UInt_t nmu,VboolR tId, Vfloat pt, Vfloat eta, Vfloat iso) {
            float pT = -10.;
	    float IS = -10.;
            for (unsigned int i=0; i<nmu; ++i) {
		if(pt[i]>pT){
                //if(pt[i]>30. && fabs(eta[i])<2.5 && pt[i]>pT){
                        pT = pt[i];
			IS = iso[i];
                }
            }
            return IS;
      };
      auto MuonAuxEta(UInt_t nmu,VboolR tId, Vfloat pt, Vfloat eta, Vfloat iso) {
            float pT = -10.;
	    float ET = -10.;
            for (unsigned int i=0; i<nmu; ++i) {
		if(pt[i]>pT){
                //if(pt[i]>pT && pt[i]>30. && iso[i]<0.15){
                        pT = pt[i];
			ET = eta[i];
                }
            }
            return ET;
      };

""")


ww_df = ww_df.Define("muonptaux", 'MuonAuxPT(nMuon,Muon_tightId,Muon_pt,Muon_eta,Muon_pfRelIso04_all)')
ww_df = ww_df.Define("muonisoaux", 'MuonAuxIso(nMuon,Muon_tightId,Muon_pt,Muon_eta,Muon_pfRelIso04_all)')
ww_df = ww_df.Define("muonetaaux", 'MuonAuxEta(nMuon,Muon_tightId,Muon_pt,Muon_eta,Muon_pfRelIso04_all)')

semi = ww_df.Filter("typeWW == 2")

print("semi events count on %s" %semi.Count().GetValue())

brlist1 = ["typeWW","typeC","nJet","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop", 
"FatJet_pt","nSV","nGenJetAK8","GenJetAK8_partonFlavour","GenJetAK8_hadronFlavour", 
"FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta","GenPart_mass","GenPart_phi","GenPart_pt", 
"GenPart_genPartIdxMother","GenPart_pdgId","GenJet_partonFlavour","GenJet_hadronFlavour","GenJetAK8_eta","GenJetAK8_mass","GenJetAK8_phi","GenJetAK8_pt", 
"Muon_eta","Muon_mass","Muon_phi","Muon_pt","Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area", 
"MET_phi","MET_pt","GenPart_statusFlags","Jet_partonFlavour","Jet_puId","Jet_nConstituents","Jet_nMuons","Jet_muonIdx1","Jet_muonIdx2", 
"nElectron","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","Electron_pfRelIso03_all","Muon_genPartIdx"]

ww_df.Snapshot("Events","analisisWW/files_mine/WW_complete.root",brlist1) 
