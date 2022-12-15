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
fileName = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/260000/062085CD-8DF4-1D40-8042-63998B6A3A95.root "
#fileName = "/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL16NanoAODv2/WW_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v15-v1/270000/10DA70C0-DD84-2049-82FB-DCC828A74F5D.root"
treeName = "Events"
#ROOT.fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, fileName)
#d = d.Range(0,3)

########GEN_PART###########

### status flag

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
            //std::cout << "statusflags " << statusflags_string << std::endl;
            for(unsigned int j=statusflags_string.length()-1; j>=statusflags_string.length()-(n+1); j--)
            {
                   //std::cout << "statusflags bit " << statusflags_string.length()-j-1 << " " << statusflags_string[j] <<std::endl;
                   ind = statusflags_string.at(j);
            }
            if(ind=='1') hardP = true;
            return hardP;
      };
      auto charmHard(UInt_t nPart, Vint status, Vint pdg) {
            int indC = 0;
	    int indCH = 0;
            for (unsigned int i=0; i<nPart; ++i) {
		if (fabs(pdg[i])==4 && testSF(i,status[i],13)) indC++;
		if (fabs(pdg[i])==4 && testSF(i,status[i],8) && testSF(i,status[i],13)) indCH++;
            }
	    vector<int> vb;
	    vb.push_back(indC);
	    vb.push_back(indCH);
            return vb;
      };
      auto charmIndex(UInt_t nPart, Vint status, Vint pdg) {
            vector<int> vb;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==4 && testSF(i,status[i],13)) vb.push_back(i);
            }
            return vb;
      };
      auto LeptonIndex(UInt_t nPart, Vint status, Vint pdg, UInt_t n) {
            vector<int> vb;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==n && testSF(i,status[i],13)) vb.push_back(i);
            }
            return vb;
      };
      auto LeptonIndexGood(UInt_t nPart, Vint status, Vint pdg, UInt_t n,Vfloat pt, Vfloat eta) {
            vector<int> vb;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==n && testSF(i,status[i],13) && pt[i]>30 && fabs(eta[i])<2.5) vb.push_back(i);
            }
            return vb;
      };
      auto vectorHP(UInt_t nPart, Vint status, Vint pdg, UInt_t n) {
            vector<int> vb;
            for (unsigned int i=0; i<nPart; ++i) {
            	vb.push_back(testSF(i,status[i],n));
            }
            return vb;
      };
      auto wpluscbool(UInt_t nPart, Vint status, Vint pdg, Vint hardP, Vint firstC) {
            int typeC = 0;
            int indC = 0;
            for (unsigned int i=0; i<nPart; ++i) {
                if (fabs(pdg[i])==5 && hardP[i]==true) {
                          typeC = 3;
                } else {
                     if (fabs(pdg[i])==4 && hardP[i]==true && firstC[i]==true) indC++;
                }
            }
            if (typeC!=3 && (indC % 2 != 0)) typeC = 2;
            if (typeC!=3 && (indC % 2 == 0) && (indC > 0)) typeC = 1;
            return typeC;
      };

""")
d = d.Define('ishard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,7)')
d = d.Define('fromhard','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,8)')
d = d.Define('firstC','vectorHP(nGenPart,GenPart_statusFlags,GenPart_pdgId,12)')
d = d.Define('isWplusc','wpluscbool(nGenPart,GenPart_statusFlags,GenPart_pdgId,ishard,firstC)')

#d = d.Filter('isWplusc == 1')

brlist1 = ["charm_status","charm_status_13","nGenPart","GenPart_genPartIdxMother","GenPart_pdgId","CharmJetDeltaR"]

brlist = ["typeWW","typeC","nJet","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop",
"FatJet_phi","FatJet_pt","nSV","nGenJetAK8","GenJetAK8_partonFlavour","GenJetAK8_hadronFlavour",
"isMuonGood","isJetGood","FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta","GenPart_mass","GenPart_phi","GenPart_pt",
"GenPart_genPartIdxMother","GenPart_pdgId","GenJet_partonFlavour","GenJet_hadronFlavour","GenJetAK8_eta","GenJetAK8_mass","GenJetAK8_phi","GenJetAK8_pt",
"Muon_eta","Muon_mass","Muon_phi","Muon_pt","Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area",
"FirstJetPt","MET_phi","MET_pt",
"nElectron","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","nElGood","isElGood","nFatJetGood","isFatJetGood"]

brlist2 = ["typeWW","typeC","nJet","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop",
"FatJet_phi","FatJet_pt","nSV","isMuonGood","isJetGood","FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta","GenJet_partonFlavour",
"GenPart_mass","GenPart_phi","GenPart_pt","GenPart_genPartIdxMother","GenPart_pdgId","Muon_eta","Muon_mass","Muon_phi","Muon_pt",
"Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area","FirstJetPt","MET_phi","MET_pt",
"nElectron","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","Jet_partonFlavour","charm_pt_1","charm_quantity","charm_status_quantity",
"nCharmGenJets","nGoodGenJets","nRawCharmGenJets","nSemiRawCharmGenJets","RawCharmJet_pt","muon_jetGen","GoodGenJetsID","W_hadronic_charge","ishard","fromhard"]

brlist3 = ["nGenPart","GenPart_genPartIdxMother","GenPart_pdgId","ishard","fromhard","isWplusc","firstC","GenJet_partonFlavour"]

dir = "/nfs/cms/vazqueze/analisisWW/old_files/"

d.Snapshot("Events", dir+"Wjets.root",brlist3)

