import ROOT, os, sys 
from ROOT import *

#from ROOT import VecOps 
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
gStyle.SetCanvasDefH(1600)

# We prepare an input tree to run on
fileName = "/nfs/cms/vazqueze/analisisWW/files_mine/WW_complete.root" 
treeName = "Events"

#ROOT.fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that 
#allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, fileName)

d = d.Filter('typeWW==2')

################################################
################   LEPTONES   ##################
################################################


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

d1 = d.Define('MuonGoodInd','muonInd(nMuon,Muon_pt,Muon_eta,Muon_pfRelIso04_all)')
d2 = d1.Define('ElectronGoodInd','elInd(nElectron, Electron_pt, Electron_eta, Electron_pfRelIso03_all)')
d3 = d2.Define('nMuonGood','MuonGoodInd.size()')
d4 = d3.Define('nElectronGood','ElectronGoodInd.size()')

########################################
##############   JETS   ################
########################################


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
""")

df1 = d4.Define('JetGoodInd','JetInd(nJet, Jet_pt, Jet_eta, Jet_phi, MuonGoodInd, Muon_eta, Muon_phi, ElectronGoodInd, Electron_eta, Electron_phi)')
df2 = df1.Define('nJetGood','JetGoodInd.size()')


## Muon dentro del jet

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
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

df3 = df2.Define('MuonJetInd','JetMuonInd(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, ElectronGoodInd)')

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

df4 = df3.Define('MuonJetGood','muonCond( MuonJetInd, Muon_pt, Muon_eta, Muon_pfRelIso04_all)')

##################
####   SSOS   ####
##################


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

df5 = df4.Define('MuonLepSign','SSOS(MuonGoodInd, ElectronGoodInd, MuonJetInd, Electron_charge, Muon_charge)')

##### Filtros

entries1 = df5.Count()

## leptones

entries_lepton = df5.Filter('!(nMuonGood==1) != !(nElectronGood==1)').Count()

# comprobamos que efectivamente se hace bien
entries_muon = df5.Filter('(nMuonGood==1) && (nElectronGood!=1)').Count()
entries_electron = df5.Filter('(nMuonGood!=1) && (nElectronGood==1)').Count()

## jets

entries_jet = df5.Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>0').Count()

## filter of this characteristic

entries_jet_muon = df5.Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>0').Filter('MuonJetInd.size() == 1').Count()

## filter

entries_jet_muon_good = df5.Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>0').Filter('MuonJetInd.size() == 1').Filter('MuonJetGood').Count()

entries_jet_muon_good_charm = df5.Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>0').Filter('MuonJetInd.size() == 1').Filter('MuonJetGood').Filter('typeC==1').Count()
entries_jet_muon_good_nocharm = df5.Filter('!(nMuonGood==1) != !(nElectronGood==1)').Filter('nJetGood>0').Filter('MuonJetInd.size() == 1').Filter('MuonJetGood').Filter('typeC==-1').Count()

########### prints

print("%s semi events" %entries1.GetValue())

print("%s is the quantity with 1 good muon or 1 good electron " %entries_lepton.GetValue())

print("%s is the quantity with 1 good muon" %entries_muon.GetValue())
print("%s is the quantity with 1 good electron" %entries_electron.GetValue())

print("%s is left when requesting at least 1 good jet " %entries_jet.GetValue())

print("%s charm events with 1 good muon inside jet " %entries_jet_muon_good_charm.GetValue())
print("%s no charm events with 1 good muon inside jet " %entries_jet_muon_good_nocharm.GetValue())

#####################################

#print("%s charm events with 0 muons inside jet " %df3.Filter('MuonJetInd.size() == 0 && typeC == 1').Count().GetValue())
#print("%s charm events with 1 muon inside jet " %df3.Filter('MuonJetInd.size() == 1 && typeC == 1').Count().GetValue())
#print("%s charm events with 2 muons inside jet " %df3.Filter('MuonJetInd.size() == 2 && typeC == 1').Count().GetValue())

#print("%s not charm events with 0 muons inside jet " %df3.Filter('MuonJetInd.size() == 0 && typeC == -1').Count().GetValue())
#print("%s not charm events with 1 muon inside jet " %df3.Filter('MuonJetInd.size() == 1 && typeC == -1').Count().GetValue())
#print("%s not charm events with 2 muons inside jet " %df3.Filter('MuonJetInd.size() == 2 && typeC == -1').Count().GetValue())

#print("%s charm events with muon from jet and lepton with same sign" %df5.Filter('MuonLepSign>0 && typeC == 1').Count().GetValue())
#print("%s charm events with muon from jet and lepton with different sign" %df5.Filter('MuonLepSign<0 && typeC == 1').Count().GetValue())

#print("%s not charm events with muon from jet and lepton with same sign" %df5.Filter('MuonLepSign>0 && typeC == -1').Count().GetValue())
#print("%s not charm events with muon from jet and lepton with different sign" %df5.Filter('MuonLepSign<0 && typeC == -1').Count().GetValue())


########################################################################################################################################################################


#############################
####     DATA SAVING     ####
#############################

brlist = ["typeWW","typeC","nJet","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop",
"FatJet_pt","nSV","nGenJetAK8","GenJetAK8_partonFlavour","GenJetAK8_hadronFlavour",
"FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta","GenPart_mass","GenPart_phi","GenPart_pt",
"GenPart_genPartIdxMother","GenPart_pdgId","GenJet_partonFlavour","GenJet_hadronFlavour","GenJetAK8_eta","GenJetAK8_mass","GenJetAK8_phi","GenJetAK8_pt",
"Muon_eta","Muon_mass","Muon_phi","Muon_pt","Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area",
"MET_phi","MET_pt","GenPart_statusFlags","Jet_partonFlavour","Jet_puId","Jet_nConstituents","Jet_nMuons","Jet_muonIdx1","Jet_muonIdx2",
"nElectron","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","Electron_pfRelIso03_all",
"MuonGoodInd","nMuonGood","ElectronGoodInd","nElectronGood","JetGoodInd","nJetGood","MuonJetInd","Muon_genPartIdx","MuonLepSign","MuonJetGood"]


#df5.Snapshot("Events", "analisisWW/files_mine/semi_filtered_reco.root",brlist)

