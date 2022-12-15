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
fileName = "/nfs/cms/vazqueze/analisisWW/files_mine/semi_filtered_reco.root" 
treeName = "Events"

#ROOT.fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that 
#allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, fileName)

print("%s semi events" %d.Count().GetValue())

## pequeña comprobacion con Jet_muonid

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto JetMuonInd(UInt_t njet, Vint good, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vint mu_good, Vfloat mu_eta, Vfloat mu_phi, Vfloat mu_pt, Vint el_good, Vint jet_muon, Vint jet_nmu) {
            vector<bool> vb;
            bool cond = false;
            bool cond1 = true;
	    bool cond2 = false;
            float ptM{-10.};
            float ptJ{-10.};
            if (good.size() == 1){
                for (unsigned int i=0; i<nmu; ++i){
                        cond = ROOT::VecOps::DeltaR(mu_eta[i],eta[good[0]],mu_phi[i],phi[good[0]]) < 0.4;
                        if (mu_good.size() > 0) cond1 = mu_good[0] != i;
                        if(cond && mu_pt[i] > ptM && cond1){
                                if (jet_nmu[good[0]]>0) cond2 = (i == jet_muon[good[0]]); 
				vb.push_back(cond2);
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
					if (jet_nmu[good[j]]>0) cond2 = (i == jet_muon[good[j]]);
	                                vb.push_back(cond2);
                                        ptM = mu_pt[i];
                                        ptJ = pt[good[j]];
                                }
                        }
                }
            }
            return vb;
      };
""")

df2 = d.Define('MuonJetCorr','JetMuonInd(nJet, JetGoodInd, Jet_pt, Jet_eta, Jet_phi, nMuon, MuonGoodInd, Muon_eta, Muon_phi, Muon_pt, ElectronGoodInd, Jet_muonIdx1, Jet_nMuons)')

print("%s events with 1 muon corresponding to nano " %df2.Filter('MuonJetInd.size() == 1').Filter('MuonJetCorr[0]').Count().GetValue())
print("%s events with 1 muon not corresponding to nano" %df2.Filter('MuonJetInd.size() == 1').Filter('!MuonJetCorr[0]').Count().GetValue())

## Origen de los muones SS tanto para charm como no charm

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto MuonMother(Vint mu_jet, Vint pdgID, Vint motherID, Vint mu_genId) {

	int mother = pdgID[motherID[mu_genId[mu_jet[0]]]];
	return mother;

      };
""")

df1 = df2.Define('MuonJetMother','MuonMother(MuonJetInd,GenPart_pdgId, GenPart_genPartIdxMother, Muon_genPartIdx)')

###########plots

def plot(hists,labels,filename,logY=False,lim1=True):

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

        if lim1:
                hists[0].GetYaxis().SetRangeUser(0,1)

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
                if lim1:
                        hist.GetYaxis().SetRangeUser(0,1)

        legend.Draw()

        if logY:
                c.SetLogy()

        plot_filename = plotdir+filename

        c.Print(plot_filename)
        c.Delete()

######################################################################################################3

charm_SS = df1.Filter('typeC == 1 && MuonLepSign ==1')
nocharm_SS = df1.Filter('typeC == -1 && MuonLepSign ==1')

charm_OS = df1.Filter('typeC == 1 && MuonLepSign ==-1')
nocharm_OS = df1.Filter('typeC == -1 && MuonLepSign ==-1')

## muon from jet mother ID

hmet1 = charm_SS.Histo1D(("","",200,-500,500),"MuonJetMother")
hmet2 = charm_OS.Histo1D(("","",200,-500,500),"MuonJetMother")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

plot([hmet2,hmet1],["charm OS","charm SS"],"/muon_jet_mother_charm.pdf",lim1=False)

hmet1 = nocharm_SS.Histo1D(("","",200,-500,500),"MuonJetMother")
hmet2 = nocharm_OS.Histo1D(("","",200,-500,500),"MuonJetMother")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

plot([hmet2,hmet1],["no charm OS","no charm SS"],"/muon_jet_mother_nocharm.pdf",lim1=False)


#############################
####     DATA SAVING     ####
#############################

brlist = ["typeWW","typeC","nJet","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop",
"FatJet_phi","FatJet_pt","nSV","nGenJetAK8","GenJetAK8_partonFlavour","GenJetAK8_hadronFlavour",
"FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta","GenPart_mass","GenPart_phi","GenPart_pt",
"GenPart_genPartIdxMother","GenPart_pdgId","GenJet_partonFlavour","GenJet_hadronFlavour","GenJetAK8_eta","GenJetAK8_mass","GenJetAK8_phi","GenJetAK8_pt",
"Muon_eta","Muon_mass","Muon_phi","Muon_pt","Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area",
"MET_phi","MET_pt","GenPart_statusFlags","Jet_partonFlavour","Jet_puId","Jet_nConstituents","Jet_nMuons","Jet_muonIdx1","Jet_muonIdx2",
"nElectron","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","Electron_pfRelIso03_all",
"MuonGoodInd","nMuonGood","ElectronGoodInd","nElectronGood","JetGoodInd","nJetGood","MuonJetInd"]

#df_jet_muon.Snapshot("Events", "analisisWW/files_mine/semi_filtered_reco.root",brlist)

