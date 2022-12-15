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
fileName = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/C28F9BCD-5A01-D645-8759-EA075D16E17D.root"
treeName = "Events"
#ROOT.fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, fileName)

# Create a df by selecting all events with at least 1 good muon using C++ code
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using namespace std;
      int nGoodJets(UInt_t njet, Vfloat pt, Vfloat eta) {
            unsigned int ngood = 0;
            for (unsigned int i=0; i<njet; ++i) {
                  if (pt[i]<30.) continue;
                  if (fabs(eta[i])>2.5) continue;
                  ngood++;
            }
            return ngood;
      };
      vector<bool> isGoodJets(UInt_t njet, Vfloat pt, Vfloat eta) {
            vector<bool> vb;
	    for (unsigned int i=0; i<njet; ++i) {
		vb.push_back(pt[i]>30. && fabs(eta[i])<2.5);
	    }
            return vb;
      };

""")


my_call = "isGoodJets(nJet,Jet_pt,Jet_eta)"
df_vb  = d.Define("isJetGood", my_call)


# All input variables for the method present in the tree
my_call = "nGoodJets(nJet,Jet_pt,Jet_eta)"
df1  = df_vb.Define("nJetGood", my_call)
#filtered_good_df = good_df.Filter("nJetGood>0")
#filtered_good_df.Report().Print()

# Create a df by selecting all events with at least 1 good muon using C++ code
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      int FirstJetPt(UInt_t njet, Vfloat pt, Vfloat eta) {
            float jet_pt = -1.;
	    if (njet>0){
		int i=0;
		while(i<njet && pt[i]>30.){
             	   jet_pt=pt[i];
		   i++;
            	}
            }
            return jet_pt;
      };
""")
# All input variables for the method present in the tree
my_call = "FirstJetPt(nJet,Jet_pt,Jet_eta)"
df2  = df1.Define("FirstJetPt", my_call)

########### FAT JETS #############

my_call = "isGoodJets(nFatJet,FatJet_pt,FatJet_eta)"
df3  = df2.Define("isFatJetGood", my_call)


# All input variables for the method present in the tree
my_call = "nGoodJets(nFatJet,FatJet_pt,FatJet_eta)"
df4  = df3.Define("nFatJetGood", my_call)

########## MUON ###########

# At some point we will need a much more sophisticated selection procedure
# A C++ function is needed for that
# Create a df by selecting all events with at least 1 good muon using C++ code
gInterpreter.Declare("""
      using VboolR = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using namespace std;
      int nGoodMuons(UInt_t nmu, Vfloat pt, Vfloat eta, Vfloat iso) {
            unsigned int ngood = 0;
            for (unsigned int i=0; i<nmu; ++i) {
                  if (pt[i]<30.) continue;
                  if (fabs(eta[i])>2.5) continue;
                  if (iso[i]>0.15) continue;
                  ngood++;
            }
            return ngood;
      };
      vector<bool> isGoodMuons(UInt_t nmuon, VboolR tId, Vfloat pt, Vfloat eta, Vfloat iso) {
            vector<bool> vb;
            for (unsigned int i=0; i<nmuon; ++i) {
                vb.push_back(pt[i]>30. && fabs(eta[i])<2.5 && iso[i]<0.15 && tId[i]);
            }
            return vb;
      };

""")

# All input variables for the method present in the tree
my_call = "nGoodMuons(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all)"
df3 = df4.Define("nMuonGood", my_call)
my_call2 = "isGoodMuons(nMuon,Muon_tightId,Muon_pt,Muon_eta,Muon_pfRelIso04_all)"
dg = df3.Define("isMuonGood", my_call2)


#########ELECTRON###########

gInterpreter.Declare("""
      using VboolR = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;
      int nGoodEl(UInt_t nel, Vfloat pt, Vfloat eta, Vint cutB) {
            unsigned int ngood = 0;
            for (unsigned int i=0; i<nel; ++i) {
                  if (pt[i]<30.) continue;
                  if (fabs(eta[i])>2.5) continue;
                  if (cutB[i]<3) continue;
                  ngood++;
            }
            return ngood;
      };
      vector<bool> isGoodEl(UInt_t nel, Vfloat pt, Vfloat eta, Vint cutB) {
            vector<bool> vb;
            for (unsigned int i=0; i<nel; ++i) {
                vb.push_back(pt[i]>30. && fabs(eta[i])<2.5 && cutB[i]>2);
            }
            return vb;
      };

""")

# All input variables for the method present in the tree
my_call = "nGoodEl(nElectron,Electron_pt,Electron_eta,Electron_cutBased)"
dg2 = dg.Define("nElGood", my_call)
my_call2 = "isGoodEl(nElectron, Electron_pt, Electron_eta, Electron_cutBased)"
good_df = dg2.Define("isElGood", my_call2)

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
ww_df1 = good_df.Define("typeWW", my_call)
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

c = ROOT.TCanvas('c','c',800,700)
c.cd()

# How the analysis is stacked and scheduled to run: twice slower
hmet = semi.Histo1D("muonptaux")

hmet.Draw()
plot_filename = plotdir+'/muonptaux2.pdf'
c.SetLogy()
c.Print(plot_filename)

c1 = ROOT.TCanvas('c','c',800,700)
c1.cd()

# How the analysis is stacked and scheduled to run: twice slower
hmet = semi.Histo1D("muonisoaux")

hmet.Draw()
plot_filename = plotdir+'/muonisoaux2.pdf'
c1.SetLogy()
c1.Print(plot_filename)

c2 = ROOT.TCanvas('c','c',800,700)
c2.cd()

# How the analysis is stacked and scheduled to run: twice slower
hmet = semi.Histo1D("muonetaaux")

hmet.Draw()
plot_filename = plotdir+'/muonetaaux2.pdf'
c2.SetLogy()
c2.Print(plot_filename)


#print(d.GetColumnNames())

c = ROOT.TCanvas('c','c',800,700)
c.cd()

# How the analysis is stacked and scheduled to run: twice slower
hmet = ww_df.Histo1D(("typeWW","typeWW",50,-25,25),"typeWW")

hmet.Draw()
plot_filename = plotdir+'/typeWW.pdf'

c.Print(plot_filename)

c1 = ROOT.TCanvas('c','c',800,700)
c1.cd()

# How the analysis is stacked and scheduled to run: twice slower
hmet = ww_df.Histo1D(("typeC","typeC",10,-5,5),"typeC")

hmet.Draw()
plot_filename = plotdir+'/typeC.pdf'

c1.Print(plot_filename)

c2 = ROOT.TCanvas('c','c',800,700)
c2.cd()

# How the analysis is stacked and scheduled to run: twice slower
hmet = ww_df.Histo1D("nFatJetGood")

hmet.Draw()
plot_filename = plotdir+'/nFatJetGood.pdf'

c2.Print(plot_filename)

brlist1 = ["typeWW","typeC","nJet","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop",
"FatJet_phi","FatJet_pt","nSV","nGenJetAK8","GenJetAK8_partonFlavour","GenJetAK8_hadronFlavour",
"FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta","GenPart_mass","GenPart_phi","GenPart_pt",
"GenPart_genPartIdxMother","GenPart_pdgId","GenJet_partonFlavour","GenJet_hadronFlavour","GenJetAK8_eta","GenJetAK8_mass","GenJetAK8_phi","GenJetAK8_pt",
"Muon_eta","Muon_mass","Muon_phi","Muon_pt","Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area",
"FirstJetPt","MET_phi","MET_pt","GenPart_statusFlags","Jet_partonFlavour","Jet_puId","Jet_nConstituents","Jet_nMuons","Jet_muonIdx1","Jet_muonIdx2",
"nElectron","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","nFatJetGood","isFatJetGood"]

ww_df.Snapshot("Events", "analisisWW/files_mine/WW_complete.root",brlist1)

brlist = ["typeWW","typeC","nJet","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop",
"FatJet_phi","FatJet_pt","nSV","nGenJetAK8","GenJetAK8_partonFlavour","GenJetAK8_hadronFlavour",
"isMuonGood","isJetGood","FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta","GenPart_mass","GenPart_phi","GenPart_pt",
"GenPart_genPartIdxMother","GenPart_pdgId","GenJet_partonFlavour","GenJet_hadronFlavour","GenJetAK8_eta","GenJetAK8_mass","GenJetAK8_phi","GenJetAK8_pt",
"Muon_eta","Muon_mass","Muon_phi","Muon_pt","Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area",
"FirstJetPt","MET_phi","MET_pt","GenPart_statusFlags","Jet_partonFlavour","Jet_puId","Jet_nConstituents","Jet_nMuons","Jet_muonIdx1","Jet_muonIdx2",
"nElectron","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","nElGood","isElGood","nFatJetGood","isFatJetGood"]

#ww_df.Snapshot("Events", "MyWW.root",brlist)

Hadronic = ww_df.Filter("typeWW == 1")
Semi = ww_df.Filter("typeWW == 2")
semicharm = Semi.Filter("typeC == 1")
Scrap = ww_df.Filter("typeWW == -1")

#Hadronic.Snapshot("Events", "MyWWHadronic.root",brlist)
Semi.Snapshot("Events", "MyWWSemi.root",brlist)
semicharm.Snapshot("Events", "MyWWSemiCharm.root",brlist)

#Scrap.Snapshot("Events", "Scrap.root",brlist)

