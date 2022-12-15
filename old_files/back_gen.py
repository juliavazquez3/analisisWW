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

# We prepare an input tree to run on
fileName = "/nfs/cms/vazqueze/WJets.root"
treeName = "Events"
#ROOT.fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, fileName)

## status flag

gInterpreter.Declare("""
      #include <bitset>
      #include <string>
      #include <iostream>
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
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
""")

df_aux = d.Define('aux_charm_status','charmHard(nGenPart,GenPart_statusFlags,GenPart_pdgId)')
df_aux = df_aux.Define('charm_quantity','aux_charm_status[0]')
df = df_aux.Define('charm_status_quantity','aux_charm_status[1]')

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
      };
      auto partsIDs(UInt_t nGenPart, Vint partId, Vint motherId) {
            int type = -1;
            auto c = (partId == 24)||(partId == -24);
            std::vector<int> v1(size(partId));
            std::iota(std::begin(v1),std::end(v1),0);
            ROOT::VecOps::RVec<int> myRVec(v1.data(), v1.size());
            int v2 = -1;
            auto Windexes = ROOT::VecOps::Where(c,myRVec,v2);
            int idWlast1 = ROOT::VecOps::Max(Windexes);
            int idpart1 = getIndex(nGenPart,motherId, idWlast1);

            std::vector<int> v;
            v.push_back(idWlast1);
            v.push_back(idpart1);
            return idWlast1;
      };
      auto firstW(UInt_t nGenPart, Vint mother, Vint partId, Vfloat part_pt) {
	    float pt{-10.};
	    int ind = 0;
	    if (nGenPart>0){
		    ind = partsIDs(nGenPart,partId,mother);
		    pt = part_pt[ind];
	    }
            return pt;
      };
""")

df_aux = df.Define("W_pt",'firstW(nGenPart,GenPart_genPartIdxMother,GenPart_pdgId,GenPart_pt)')


########################################################################

def plot(hists,labels,filename,logY=False):

        c = ROOT.TCanvas('c','c',800,700)
        c.cd()
        colores = [ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kMagenta,ROOT.kOrange,ROOT.kViolet]

        if len(hists)>len(colores):
                raise ValueError

        color = {}
        label = {}
        for i in range(len(hists)):
                color[str(hists[i])]=colores[i]
                label[str(hists[i])]=labels[i]

        hists[0].Draw()
        hists[0].SetName(str(hists[0]))
        hists[0].SetLineColor(colores[0])

        # Add legend
        legend = ROOT.TLegend(0.62, 0.70, 0.82, 0.88)
        legend.SetFillColor(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.03)
        legend.AddEntry(str(hists[0]),labels[0],"L")

        for hist in hists[1:]:
                hist.Draw("Sames")
                hist.SetName(str(hist))
                hist.SetLineColor(color[str(hist)])
                legend.AddEntry(str(hist),label[str(hist)],"L")

        legend.Draw()

        if logY:
                c.SetLogy()

        plot_filename = plotdir+filename

        c.Print(plot_filename)
        c.Delete()


############################################
###########    PLOTS   #####################
############################################

########## gen

hmet1 = df_aux.Histo1D(("","",10,0,10),"charm_quantity")
hmet2 = df_aux.Histo1D(("","",10,0,10),"charm_status_quantity")
plot([hmet1,hmet2],["Charm part","Hard process charm part"],"/WJets_charm_quantity.pdf")

hmet1 = df_aux.Histo1D(("","",100,0,500),"W_pt")
plot([hmet1],["W pt"],"/WJets_W_pt.pdf",logY=True)

brlist = ["nJetGood","nMuonGood","nElectronGood","nJet","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop",
"FatJet_phi","FatJet_pt","nSV","nGenJetAK8","GenJetAK8_partonFlavour","GenJetAK8_hadronFlavour",
"isMuonGood","isJetGood","FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta","GenPart_mass","GenPart_phi","GenPart_pt",
"GenPart_genPartIdxMother","GenPart_pdgId","GenJet_partonFlavour","GenJet_hadronFlavour","GenJetAK8_eta","GenJetAK8_mass","GenJetAK8_phi","GenJetAK8_pt",
"Muon_eta","Muon_mass","Muon_phi","Muon_pt","Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area",
"MET_phi","MET_pt","GenPart_statusFlags","Jet_partonFlavour","pruebaSF_8",
"nElectron","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","isElGood"]

#d.Snapshot("Events", "WJets_testSF.root",brlist)


