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
fileName_reco = "/nfs/cms/vazqueze/analisisWW/old_files/SemiNew.root"
fileName_gen = "/nfs/cms/vazqueze/analisisWW/old_files/SemiGen.root"

treeName = "Events"
#ROOT.fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.

d_reco = ROOT.RDataFrame(treeName, fileName_reco)
d_gen = ROOT.RDataFrame(treeName, fileName_gen)

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto charmJets(UInt_t njet,Vint flavour, Vfloat pt, Vfloat eta) {
            int ind = 0;
            for (unsigned int i=0; i<njet; ++i) {
                if(fabs(flavour[i])==4 && pt[i]>30 && fabs(eta[i])<2.5) {
                        ind++;                
		}
            }
            return ind;
      };
      auto rawFlavour(UInt_t njet,Vint flavour, UInt_t n) {
            int ind = 0;
            for (unsigned int i=0; i<njet; ++i) {
                if(fabs(flavour[i])==n) {
                        ind++;
                }
            }
            return ind;
      };
      auto charmJetPT(UInt_t njet,Vint flavour, Vfloat pt, Vfloat eta) {
	    float pT{-10.};
            for (unsigned int i=0; i<njet; ++i) {
                if(fabs(flavour[i])==4 && pt[i]>pT) {
                        pT = pt[i];
                }
            }
            return pT;
      };
      auto strangeJetPT(UInt_t njet,Vint flavour, Vfloat pt, Vfloat eta) {
            float pT{-10.};
            for (unsigned int i=0; i<njet; ++i) {
                if(fabs(flavour[i])==3 && pt[i]>pT) {
                        pT = pt[i];
                }
            }
            return pT;
      };
      auto jetsMuons(UInt_t njet,Vint flavour, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t nmu, Vfloat mu_eta, Vfloat mu_phi) {
	    int indJ = 0;
            int	indC = 0;
	    bool cond = true;
	    int c = 0;
	    int	c1 = 0;
            for (unsigned int i=0; i<njet; ++i) {
		c = 0;
		c1 = 0;
		for (unsigned int j=0; j<nmu; ++j) {
			cond = ROOT::VecOps::DeltaR(eta[i],mu_eta[j],phi[i],mu_phi[j])<0.4;
                	if(fabs(eta[i])<2.5 && pt[i]>30 && cond) {
                        	c = 1;
                        	if(fabs(flavour[i])==4) {
                                	c1 = 1;
                        	}
                	}
            	}
                indJ = indJ + c;
       	       	indC = indC + c1;
	    }
	    vector<int> vb;
	    vb.push_back(indJ);
            vb.push_back(indC);
            return vb;
      };
      auto jetsGenMuons(UInt_t njet,Vint flavour, Vfloat pt, Vfloat eta, Vfloat phi, UInt_t npart, Vfloat part_eta, Vfloat part_phi, Vint pdg, Vint mother) {
            int indJ = 0;
            int indC = 0;
            bool cond = true;
            int c = 0;
            int c1 = 0;
            for (unsigned int i=0; i<njet; ++i) {
                c = 0;
                c1 = 0;
                for (unsigned int j=0; j<npart; ++j) {
                        cond = ROOT::VecOps::DeltaR(eta[i],part_eta[j],phi[i],part_phi[j])<0.4;
                        if(fabs(eta[i])<2.5 && pt[i]>30 && cond && fabs(pdg[j])==13 && (fabs(pdg[mother[j]]) != 24)) {
                                c = 1;
                                if(fabs(flavour[i])==4) {
                                        c1 = 1;
                                }
                        }
                }
                indJ = indJ + c;
                indC = indC + c1;
            }
            vector<int> vb;
            vb.push_back(indJ);
            vb.push_back(indC);
            return vb;
      };

""")

d_reco = d_reco.Define('nCharmJets','charmJets(nJet,Jet_partonFlavour, Jet_pt, Jet_eta)')
d_gen = d_gen.Define('RawCharmRECO_pt','charmJetPT(nJet,Jet_partonFlavour, Jet_pt, Jet_eta)')
d_gen = d_gen.Define('RawStrangeRECO_pt','strangeJetPT(nJet,Jet_partonFlavour, Jet_pt, Jet_eta)')
d_gen = d_gen.Define('RawCharmGEN_pt','charmJetPT(nGenJet,GenJet_partonFlavour, GenJet_pt, GenJet_eta)')
d_gen = d_gen.Define('RawStrangeGEN_pt','strangeJetPT(nGenJet,GenJet_partonFlavour, GenJet_pt, GenJet_eta)')

d1 = d_reco.Define('auxRawMuonRecoJets','jetsMuons(nJet,Jet_partonFlavour, Jet_pt, Jet_eta, Jet_phi, nMuon, Muon_eta, Muon_phi)')
d2 = d1.Define('nJetRawMuonReco','auxRawMuonRecoJets[0]')
d_reco = d2.Define('nCharmJetRawMuonReco','auxRawMuonRecoJets[1]')

d1 = d_gen.Define('auxRawMuonGenJets','jetsGenMuons(nGenJet,GenJet_partonFlavour, GenJet_pt, GenJet_eta, GenJet_phi, nGenPart, GenPart_eta, GenPart_phi,GenPart_pdgId, GenPart_genPartIdxMother)')
d2 = d1.Define('nJetRawMuonGen','auxRawMuonGenJets[0]')
d_gen = d2.Define('nCharmJetRawMuonGen','auxRawMuonGenJets[1]')

d_gen = d_gen.Define('nRawGENCharm','rawFlavour(nGenJet,GenJet_partonFlavour,4)')
d_gen = d_gen.Define('nRawGENStrange','rawFlavour(nGenJet,GenJet_partonFlavour,3)')

d_reco = d_reco.Define('nRawRECOCharm','rawFlavour(nJet,Jet_partonFlavour,4)')
d_reco = d_reco.Define('nRawRECOStrange','rawFlavour(nJet,Jet_partonFlavour,3)')

## pile up test

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto pileupID(UInt_t njet, Vfloat pt, Vint pileup) {
            int ind1 = -1;
	    int ind2 = -1;
            if(njet>0) {
            	if(pt[0]<50) {
			ind1 = pileup[0];
		}
                if(pt[0]>50) {
                        ind2 = pileup[0];
                }
            }
	    vector<int> vb;
	    vb.push_back(ind1);
            vb.push_back(ind2);
            return vb;
      };
""")

d_reco=d_reco.Define('pileupGood','pileupID(nJet,Jet_pt,Jet_puId)[0]')
d_reco=d_reco.Define('pileupTry','pileupID(nJet,Jet_pt,Jet_puId)[1]')

## muon test

d_reco=d_reco.Define("muon_iso_test",'nMuon>0 ? Muon_pfRelIso04_all[0] : -1')

print("%s Reco events" %d_reco.Count().GetValue())
print("%s Gen events" %d_gen.Count().GetValue())

charm_reco = d_reco.Filter('typeC == 1')
nocharm_reco = d_reco.Filter('typeC == -1')

charm_gen = d_gen.Filter('typeC == 1')
nocharm_gen = d_gen.Filter('typeC == -1')

print("%s test" %d_reco.Filter('nJetRawMuonReco>0').Count().GetValue())

### jets pt

####### PT jet study  #############

jet_none = d_reco.Filter('nJet == 0')
jet_one = d_reco.Filter('nJet == 1')
jet_two = d_reco.Filter('nJet == 2')
jet_three = d_reco.Filter('nJet == 3')
jet_four = d_reco.Filter('nJet == 4')

print("%s entries dont have gen jets" %jet_none.Count().GetValue())
print("%s entries have 1 gen jet" %jet_one.Count().GetValue())
print("%s entries have 2 gen jets" %jet_two.Count().GetValue())
print("%s entries have 3 gen jets" %jet_three.Count().GetValue())
print("%s entries have 4 gen jets" %jet_four.Count().GetValue())

jet_one = jet_one.Define('jetG_pt1','Jet_pt[0]')
jet_two = jet_two.Define('jetG_pt1','Jet_pt[0]')
jet_two = jet_two.Define('jetG_pt2','Jet_pt[1]')
jet_three = jet_three.Define('jetG_pt1','Jet_pt[0]')
jet_three = jet_three.Define('jetG_pt2','Jet_pt[1]')
jet_three = jet_three.Define('jetG_pt3','Jet_pt[2]')
jet_four = jet_four.Define('jetG_pt1','Jet_pt[0]')
jet_four = jet_four.Define('jetG_pt2','Jet_pt[1]')
jet_four = jet_four.Define('jetG_pt3','Jet_pt[2]')
jet_four = jet_four.Define('jetG_pt4','Jet_pt[3]')

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


############################################
###########    PLOTS   #####################
############################################

############## RECO ############

### Leptones

#hmet1 = charm_reco.Histo1D(("","",10,0,10),"nElectron")
#hmet2 = nocharm_reco.Histo1D(("","",10,0,10),"nElectron")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nElectron.pdf")

#hmet3 = charm_reco.Histo1D(("","",10,0,10),"nMuonGood")
#hmet4 = nocharm_reco.Histo1D(("","",10,0,10),"nMuonGood")

#hmet3.Scale(1/hmet3.Integral())
#hmet4.Scale(1/hmet4.Integral())

#plot([hmet3,hmet4],["Charm","NoCharm"],"/normed_nMuonGood.pdf")

#hmet1 = charm_reco.Histo1D(("","",10,0,10),"nElectronGood")
#hmet2 = nocharm_reco.Histo1D(("","",10,0,10),"nElectronGood")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nElectronGood.pdf")

## muon iso test

#hmet1 = charm_reco.Histo1D(("","",30,0,15),"muon_iso_test")
#hmet2 = nocharm_reco.Histo1D(("","",30,0,15),"muon_iso_test")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_muon_iso_test.pdf",logY=True,lim1=False)


##### Jets

#hmet1 = charm_reco.Histo1D(("","",10,0,10),"nJet")
#hmet2 = nocharm_reco.Histo1D(("","",10,0,10),"nJet")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nJetRaw.pdf")


#hmet1 = charm_reco.Histo1D(("","",10,0,10),"nJetGood")
#hmet2 = nocharm_reco.Histo1D(("","",10,0,10),"nJetGood")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nJetGood.pdf")

## flavour
#hmet1 = charm_reco.Histo1D(("","",10,0,10),"nCharmJets")
#hmet2 = nocharm_reco.Histo1D(("","",10,0,10),"nCharmJets")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nCharmJetGood.pdf")

#hmet1 = charm_reco.Histo1D(("","",10,0,10),"nRawRECOCharm")
#hmet2 = charm_reco.Histo1D(("","",10,0,10),"nRawRECOStrange")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","Strange"],"/normed_nCharmRECORAW.pdf")

## muon
#hmet1 = charm_reco.Histo1D(("","",10,0,10),"nJetRawMuonReco")
#hmet2 = nocharm_reco.Histo1D(("","",10,0,10),"nJetRawMuonReco")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nJetswithRawMuon_reco.pdf")

#hmet1 = charm_reco.Histo1D(("","",10,0,10),"nCharmJetRawMuonReco")
#hmet2 = nocharm_reco.Histo1D(("","",10,0,10),"nCharmJetRawMuonReco")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nCharmJetswithRawMuon_reco.pdf")

## muon from nanoAOD variable

hmet1 = d_reco.Histo1D(("","",10,0,10),"jetsraw_nMuNano")
hmet2 = d_reco.Histo1D(("","",10,0,10),"charmjetsraw_nMuNano")

hmet1.Scale(1/hmet1.Integral())
hmet2.Scale(1/hmet2.Integral())

plot([hmet1,hmet2],["not charm jet","charm jet"],"/normed_nRawJetwithMuonNano_reco.pdf")

hmet1 = d_reco.Histo1D(("","",10,0,10),"jetsgood_nMuNano")
hmet2 = d_reco.Histo1D(("","",10,0,10),"charmjetsgood_nMuNano")

hmet1.Scale(1/hmet1.Integral())
hmet2.Scale(1/hmet2.Integral())

plot([hmet1,hmet2],["not charm jet","charm jet"],"/normed_nGoodJetwithMuonNano_reco.pdf")

## pileup test

#hmet1 = d_reco.Histo1D(("","",10,-1,9),"pileupGood")
#hmet2 = d_reco.Histo1D(("","",10,-1,9),"pileupTry")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["less than 50 pt","more than 50 pt"],"/normed_pileupJets.pdf")

## flavour pt

#hmet1 = charm_gen.Histo1D(("","",100,-10,400),"RawCharmRECO_pt")
#hmet2 = charm_gen.Histo1D(("","",100,-10,400),"RawStrangeRECO_pt")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","Strange"],"/normed_RecoJetsPT_charmevents.pdf",logY=True,lim1=False)


############## GEN ############

##### Jets

#hmet1 = charm_gen.Histo1D(("","",10,0,10),"nGenJet")
#hmet2 = nocharm_gen.Histo1D(("","",10,0,10),"nGenJet")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nGenJetRaw.pdf")


#hmet1 = charm_gen.Histo1D(("","",10,0,10),"nGoodGenJets")
#hmet2 = nocharm_gen.Histo1D(("","",10,0,10),"nGoodGenJets")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nGenJetGood.pdf")

## flavour
#hmet1 = charm_gen.Histo1D(("","",10,0,10),"charm_status_quantity")
#hmet2 = nocharm_gen.Histo1D(("","",10,0,10),"charm_status_quantity")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nCharmPartsHardP.pdf")

#hmet1 = charm_gen.Histo1D(("","",10,0,10),"nRawGENCharm")
#hmet2 = charm_gen.Histo1D(("","",10,0,10),"nRawGENStrange")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","Strange"],"/normed_nCharmGENRAW.pdf")

## muon
#hmet1 = charm_gen.Histo1D(("","",10,0,10),"nJetRawMuonGen")
#hmet2 = nocharm_gen.Histo1D(("","",10,0,10),"nJetRawMuonGen")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nJetswithRawMuon_gen.pdf")

#hmet1 = charm_gen.Histo1D(("","",10,0,10),"nCharmJetRawMuonGen")
#hmet2 = nocharm_gen.Histo1D(("","",10,0,10),"nCharmJetRawMuonGen")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","NoCharm"],"/normed_nCharmJetswithRawMuon_gen.pdf")


## flavour pt 

#hmet1 = charm_gen.Histo1D(("","",100,-10,400),"RawCharmGEN_pt")
#hmet2 = charm_gen.Histo1D(("","",100,-10,400),"RawStrangeGEN_pt")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())

#plot([hmet1,hmet2],["Charm","Strange"],"/normed_GenJetsPT_charmevents.pdf",logY=True,lim1=False)


# jets pt

#hmet1 = jet_one.Histo1D(("","",100,0,400),"jetG_pt1")
#hmet2 = jet_two.Histo1D(("","",100,0,400),"jetG_pt1")
#hmet3 = jet_two.Histo1D(("","",100,0,400),"jetG_pt2")
#hmet4 = jet_three.Histo1D(("","",100,0,400),"jetG_pt1")
#hmet5 = jet_three.Histo1D(("","",100,0,400),"jetG_pt2")
#hmet6 = jet_three.Histo1D(("","",100,0,400),"jetG_pt3")
#hmet7 = jet_four.Histo1D(("","",100,0,400),"jetG_pt1")
#hmet8 = jet_four.Histo1D(("","",100,0,400),"jetG_pt2")
#hmet9 = jet_four.Histo1D(("","",100,0,400),"jetG_pt3")
#hmet10 = jet_four.Histo1D(("","",100,0,400),"jetG_pt4")

#hmet1.Scale(1/hmet1.Integral())
#hmet2.Scale(1/hmet2.Integral())
#hmet3.Scale(1/hmet3.Integral())
#hmet4.Scale(1/hmet4.Integral())
#hmet5.Scale(1/hmet5.Integral())
#hmet6.Scale(1/hmet6.Integral())
#hmet7.Scale(1/hmet7.Integral())
#hmet8.Scale(1/hmet8.Integral())
#hmet9.Scale(1/hmet9.Integral())
#hmet10.Scale(1/hmet10.Integral())


#plot([hmet1],["first"],"/normed_recojet_pt1.pdf",logY=True,lim1=False)
#plot([hmet2,hmet3],["first","second"],"/normed_recojet_pt2.pdf",logY=True,lim1=False)
#plot([hmet4,hmet5,hmet6],["first","second","third"],"/normed_recojet_pt3.pdf",logY=True,lim1=False)
#plot([hmet7,hmet8,hmet9,hmet10],["first","second","third","fourth"],"/normed_recojet_pt4.pdf",logY=True,lim1=False)

########### scatters

#c3 = ROOT.TCanvas('c1','c1',800,700)

#h2 = d_gen.Histo2D(("h2", "scatter plot", 100,0 , 400, 100, 0, 400),"charm_pt_1","RawCharmRECO_pt")
#h2.Draw("COLZ")
#h2.SetStats(0);
#h2.GetXaxis().SetTitle( 'Charm Particle pT' )
#h2.GetYaxis().SetTitle( 'Charm Jet pT (RECO)' )
#c3.SaveAs(plotdir+"scatterCharmPTpartVSjetRECO.pdf")

#c2 = ROOT.TCanvas('c1','c1',800,700)

#h2 = d_gen.Histo2D(("h2", "scatter plot", 100,0 , 400, 100, 0, 400),"RawCharmRECO_pt","RawCharmJet_pt")
#h2.Draw("COLZ")
#h2.SetStats(0);
#h2.GetXaxis().SetTitle( 'Charm Jet pT (RECO)' )
#h2.GetYaxis().SetTitle( 'Charm Jet pT (GEN)' )
#c2.SaveAs(plotdir+"scatterCharmjetRECOVSjetGEN.pdf")

