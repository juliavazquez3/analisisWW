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
fileName = "/nfs/cms/vazqueze/MyWW.root"
treeName = "Events"
#ROOT.fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, fileName)

df = d.Filter('typeWW == 2')


##df = d_filt.Define('deltaR_jetmuon','ROOT::VecOps::DeltaR(jet_eta1,muon_eta1,jet_phi1,muon_phi1)')

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
            int idWlast2 = idWlast1-1;
            int idpart1 = getIndex(nGenPart,motherId, idWlast1);
            int idpart2 = getIndex(nGenPart,motherId, idWlast2);

	    std::vector<int> v;
	    v.push_back(idWlast1);
	    v.push_back(idWlast2);
	    v.push_back(idpart1);
	    v.push_back(idpart2);
            v.push_back(idpart1+1);
            v.push_back(idpart2+1);
            return v;
      };
      auto tauMother(UInt_t nGenPart, Vint partId, Vint motherId, Vint parts) {
	    int ind = -1;
	    if (fabs(partId[parts[2]])==15) {
 	       ind = parts[2];
	    } else if (fabs(partId[parts[3]])==15) {
  	       ind = parts[3];
            } else if (fabs(partId[parts[4]])==15) {
               ind = parts[4];
	    } else {
  	       ind = parts[5];
	    }

	    auto c = (motherId == ind);
            std::vector<int> v1(size(motherId));
            std::iota(std::begin(v1),std::end(v1),0);
            ROOT::VecOps::RVec<int> myRVec(v1.data(), v1.size());
	    int v2 = -1;
            auto TAUind = ROOT::VecOps::Where(c,myRVec,v2);

	    int idlep = -1; 
	    for (unsigned int i=0; i<nGenPart; ++i) {
		if (TAUind[i]>0 && (fabs(partId[i])==11 || fabs(partId[i])==13)) {
			idlep = i;
		}
	    }
            return idlep;
      };
""")

my_call = "partsIDs(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother)"
df1 = df.Define("partID",my_call)

df2 = df1.Define("W_hadronic_pt",'fabs(GenPart_pdgId[partID[2]])<9 ? GenPart_pt[partID[0]] : GenPart_pt[partID[1]]')
df3 = df2.Define("W_leptonic_pt",'fabs(GenPart_pdgId[partID[2]])<9 ? GenPart_pt[partID[1]] : GenPart_pt[partID[0]]')

df4 = df3.Define("quarks_dr",'fabs(GenPart_pdgId[partID[2]])<9 ? ROOT::VecOps::DeltaR(GenPart_eta[partID[2]],GenPart_eta[partID[2]+1],GenPart_phi[partID[2]],GenPart_phi[partID[2]+1]) : ROOT::VecOps::DeltaR(GenPart_eta[partID[3]],GenPart_eta[partID[3]+1],GenPart_phi[partID[3]],GenPart_phi[partID[3]+1])')

df5 = df4.Define("Genmuon_pt",'fabs(GenPart_pdgId[partID[2]])==13 ?  GenPart_pt[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==13 ?  GenPart_pt[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==13 ?  GenPart_pt[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==13 ?  GenPart_pt[partID[5]] :-20)))')
df6 = df5.Define("Genelectron_pt",'fabs(GenPart_pdgId[partID[2]])==11 ?  GenPart_pt[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==11 ?  GenPart_pt[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==11 ?  GenPart_pt[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==11 ?  GenPart_pt[partID[5]] :-20)))')
df4 = df6.Define("Gentau_pt",'fabs(GenPart_pdgId[partID[2]])==15 ?  GenPart_pt[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==15 ?  GenPart_pt[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==15 ?  GenPart_pt[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==15 ?  GenPart_pt[partID[5]] :-20)))')

df5 = df4.Define("Genmuon_eta",'fabs(GenPart_pdgId[partID[2]])==13 ?  GenPart_eta[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==13 ?  GenPart_eta[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==13 ?  GenPart_eta[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==13 ?  GenPart_eta[partID[5]] :-20)))')
df6 = df5.Define("Genelectron_eta",'fabs(GenPart_pdgId[partID[2]])==11 ?  GenPart_eta[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==11 ?  GenPart_eta[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==11 ?  GenPart_eta[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==11 ?  GenPart_eta[partID[5]] :-20)))')
df4 = df6.Define("Gentau_eta",'fabs(GenPart_pdgId[partID[2]])==15 ?  GenPart_eta[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==15 ?  GenPart_eta[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==15 ?  GenPart_eta[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==15 ?  GenPart_eta[partID[5]] :-20)))')

entries1 = df4.Count()
print("%s entries passed semi filter" %entries1.GetValue())

testTAU = df4.Filter('Gentau_pt>0')
testMUON = df4.Filter('Genmuon_pt>0')
testEL = df4.Filter('Genelectron_pt>0')

entries1 = testTAU.Count()
entries2 = testMUON.Count()
entries3 = testEL.Count()
#print("%s entries are tau" %entries1.GetValue())
#print("%s entries are muon" %entries2.GetValue())
#print("%s entries are electron" %entries3.GetValue())

test1 = df4.Filter('Gentau_pt>0 && Genmuon_pt>0')
#print("%s entries are tau and muon" %(test1.Count()).GetValue())

test1 = df4.Filter('Gentau_pt>0 && Genmuon_pt>0&& Genelectron_pt>0')
#print("%s entries are tau and muon and el" %(test1.Count()).GetValue())


#### TAU decays, we take the set test tau

testTAU = testTAU.Define("leptonTAU_pt",'tauMother(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother, partID) != -1 ? GenPart_pt[tauMother(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother, partID)] : -1')

## Prints

#print("%s entries are from W decaying to a muon with pt>30" %df4.Filter('Genmuon_pt>30').Count().GetValue())
#print("%s entries are from W decaying to a electron with pt>30" %df4.Filter('Genelectron_pt>30').Count().GetValue())
#print("%s entries are from W -> tau then decaying to lepton with pt>30" %testTAU.Filter('leptonTAU_pt>30').Count().GetValue())


############ GEN JETS ##############

gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using namespace std;
      int nGoodJets(UInt_t njet, Vfloat pt, Vfloat eta) {
            unsigned int ngood = 0;
            for (unsigned int i=0; i<njet; ++i) {
                  if (pt[i]<25.) continue;
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
      auto jetFlavour(UInt_t nGenJ,Vint flavour,Vfloat pt, Vfloat eta){
            int ind = 0;
            for (unsigned int i=0; i<nGenJ; ++i) {
                if (fabs(flavour[i])==4 && pt[i]>30 && fabs(eta[i])<2.4) {
                        ind++;
                }
            }
            return ind;
      };
      auto CharmJetInd(UInt_t nGenJ,Vint flavour,Vfloat pt, Vfloat eta){
	    float pT{-10};
	    int ind = -1;
            for (unsigned int i=0; i<nGenJ; ++i) {
                if (pt[i]>pT && fabs(flavour[i])==4 && pt[i]>25. && fabs(eta[i])<2.4) {
                        pT = pt[i];
			ind = i;
                }
            }
            return ind;
      };
      auto StrangeJetInd(UInt_t nGenJ,Vint flavour,Vfloat pt, Vfloat eta){
            float pT{-10};
	    int ind = -1;
	    auto ptC = CharmJetInd(nGenJ,flavour,pt,eta);
            for (unsigned int i=0; i<nGenJ; ++i) {
                if (ptC!=-1 && pt[i]>pT && (fabs(flavour[i])==3 || fabs(flavour[i])==0) && pt[i]>25. && fabs(eta[i])<2.4) {
                        pT = pt[i];
			ind = i;
                }
            }
            return ind;
      };
      auto ExtraInd(UInt_t nGenJ,Vint flavour,Vfloat pt, Vfloat eta){
            float pT{-10};
            int ind = -1;
            auto ptC = CharmJetInd(nGenJ,flavour,pt,eta);
            auto ptS = StrangeJetInd(nGenJ,flavour,pt,eta);
            for (unsigned int i=0; i<nGenJ; ++i) {
                if (ptC!=-1 && ptS!=i && i!=ptC && pt[i]>pT && pt[i]>15. && fabs(flavour[i])<9 && fabs(eta[i])<2.4) {
                        pT = pt[i];
                        ind = i;
                }
            }
            return ind;
      };
      auto CharmJetPt(UInt_t nGenJ,Vint flavour,Vfloat pt, Vfloat eta){
            float pT{-10};
	    auto ind = CharmJetInd(nGenJ,flavour,pt,eta);            
            if (ind != -1) {
                   pT = pt[ind];
            }
            return pT;
      };
      auto StrangeJetPt(UInt_t nGenJ,Vint flavour,Vfloat pt, Vfloat eta){
            float pT{-10};
            auto ind = StrangeJetInd(nGenJ,flavour,pt,eta);
            if (ind != -1) {
                   pT = pt[ind];
            }
            return pT;
      };
      auto CSdeltaR(UInt_t nGenJ,Vint flavour,Vfloat phi, Vfloat eta, Vfloat pt){
            float dR{10.};
            auto indC = CharmJetInd(nGenJ,flavour,pt,eta);
            auto indS = StrangeJetInd(nGenJ,flavour,pt,eta);
            if (indC != -1 && indS != -1) {
                   dR = ROOT::VecOps::DeltaR(eta[indC],eta[indS],phi[indC],phi[indS]);
            }
            return dR;
      };
      auto CextraDeltaR(UInt_t nGenJ,Vint flavour,Vfloat phi, Vfloat eta, Vfloat pt){
            float dR{10.};
            auto indC = CharmJetInd(nGenJ,flavour,pt,eta);
            auto indE = ExtraInd(nGenJ,flavour,pt,eta);
            if (indC != -1 && indE != -1) {
                   dR = ROOT::VecOps::DeltaR(eta[indC],eta[indE],phi[indC],phi[indE]);
            }
            return dR;
      };
""")

my_call = 'jetFlavour(nGenJet,GenJet_partonFlavour, GenJet_pt, GenJet_eta)'
my_call1 = 'nGoodJets(nGenJet, GenJet_pt, GenJet_eta)'
my_call2 = 'isGoodJets(nGenJet,GenJet_pt,GenJet_eta)' 
my_call3 = 'CharmJetPt(nGenJet,GenJet_partonFlavour, GenJet_pt, GenJet_eta)'
my_call4 = 'StrangeJetPt(nGenJet,GenJet_partonFlavour, GenJet_pt, GenJet_eta)'
my_call5 = 'CSdeltaR(nGenJet,GenJet_partonFlavour, GenJet_phi, GenJet_eta, GenJet_pt)'
my_call6 = 'CextraDeltaR(nGenJet,GenJet_partonFlavour, GenJet_phi, GenJet_eta, GenJet_pt)'


df5 = df4.Define("nCharmJets",my_call)
df6 = df5.Define("nGoodGenJets",my_call1)
df4 = df6.Define("isGenJetGood",my_call2)

df5 = df4.Define("charmJet_pt",my_call3)
df6 = df5.Define("strangeJet_pt",my_call4)
df4 = df6.Define("deltaRCS",my_call5)

df5 = df4.Define("deltaRCExtra",my_call6)
df6 = df5.Define("nCharmFatJets",'jetFlavour(nGenJetAK8,GenJetAK8_partonFlavour, GenJetAK8_pt, GenJetAK8_eta)')
df4 = df6.Define("nGoodFatGenJets",'nGoodJets(nGenJetAK8, GenJetAK8_pt, GenJetAK8_eta)')

entries1 = df4.Filter('typeC==1 && nCharmJets>0').Count()
entries2 = df4.Filter('typeC==-1 && nCharmJets>0').Count()

print("%s charm entries have gen jet charm" %entries1.GetValue())
print("%s not charm entries have gen jet charm" %entries2.GetValue())

######### MUON FROM JET #########

gInterpreter.Declare("""
        using Vfloat = const ROOT::RVec<float>&;
	using Vint = const ROOT::RVec<int>&;
        using namespace std;
	bool checkprefix(int A, int B)
	{
    		// Convert numbers into strings
    		string s1 = to_string(A);
    		string s2 = to_string(B);
  
       		if (s1[0]!= s2[0]) {
         		return false;
       		}
  
    		// Return true
    		return true;
	};
	auto indMuonJet(UInt_t nGen, Vint pdg, Vint mother) {
		int ind = -1;
                for (unsigned int i=0; i<nGen; ++i) {

                // If at any index characters
                // are unequals then return false
                        if (fabs(pdg[i])==13 && checkprefix(fabs(pdg[mother[i]]),4) && fabs(pdg[mother[mother[i]]]) == 4) {
                                ind = i;
                        }
                }
		return ind;
	};
""")

df5 = df4.Define("muon_jetGen",'indMuonJet(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother)')
df4 = df5.Define("muon_jetGen_pt",'muon_jetGen != -1 ? GenPart_pt[muon_jetGen] : -10')

entries1 = df4.Filter('typeC==1 && muon_jetGen != -1').Count()
entries2 = df4.Filter('typeC==-1 && muon_jetGen != -1').Count()

print("%s charm entries have muon from gen jet charm" %entries1.GetValue())
print("%s not charm entries have muon from gen jet charm" %entries2.GetValue())

############## PLOTS ################

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


charm = df4.Filter('typeC==1')
notCharm = df4.Filter('typeC==-1')

#hmet1 = df4.Histo1D("W_hadronic_pt")
#hmet2 = df4.Histo1D("W_leptonic_pt")
#hmet3 = df4.Histo1D("quarks_dr")

hmet1 = df4.Histo1D("Lepton_ID")

#plot([hmet1,hmet2],["W_had","W_lep"],"/gen_ptDR.pdf",logY=True)
plot([hmet1],["lepton_id"],"/prueba.pdf")

######## GEN JETS #########

#hmet1 = testMUON.Histo1D(("","",10,0,10),"GenJet_partonFlavour")
#hmet2 = testEL.Histo1D(("","",10,0,10),"GenJet_partonFlavour")
#hmet3 = testTAU.Histo1D(("","",10,0,10),"GenJet_partonFlavour")
#plot([hmet1,hmet2,hmet3],["Muon","Electron","Tau"],"/GenJet_flavour.pdf",logY=True)

#hmet1 = charm.Histo1D(("","",10,0,10),"nCharmJets")
#hmet2 = notCharm.Histo1D(("","",10,0,10),"nCharmJets")
#plot([hmet1,hmet2],["Charm","Not charm"],"/nCharmJets.pdf",logY=True)

#hmet1 = charm.Histo1D(("","",10,0,10),"nGoodGenJets")
#hmet2 = notCharm.Histo1D(("","",10,0,10),"nGoodGenJets")
#plot([hmet1,hmet2],["Charm","Not charm"],"/nGoodGenJets.pdf",logY=True)

#hmet1 = charm.Histo1D(("","",100,-10,500),"charmJet_pt")
#hmet2 = charm.Histo1D(("","",100,-10,500),"strangeJet_pt")
#hmet3 = notCharm.Histo1D(("","",100,-10,500),"charmJet_pt")
#hmet4 = notCharm.Histo1D(("","",100,-10,500),"strangeJet_pt")
#plot([hmet1,hmet2,hmet3,hmet4],["Charm C","Charm S","No Charm C","Not charm S"],"/quarks_pt.pdf",logY=True)

hmet1 = charm.Histo1D(("","",80,-10,150),"muon_jetGen_pt")
hmet2 = notCharm.Histo1D(("","",80,-10,150),"muon_jetGen_pt")
plot([hmet1,hmet2],["Charm","Not charm"],"/muon_jetGen_pt.pdf",logY=True)

######## muon electron pt

#hmet1 = testMUON.Histo1D(("","",200,-20,600),"Genmuon_pt")
#hmet2 = testEL.Histo1D(("","",200,-20,600),"Genelectron_pt")
#hmet3 = testTAU.Histo1D(("","",200,-20,600),"Gentau_pt")
#plot([hmet1,hmet2,hmet3],["Muon_pt","Electron_pt","Tau_pt"],"/gen_leptonPT.pdf",logY=True)

#hmet1 = testMUON.Histo1D(("","",100,-20,20),"Genmuon_eta")
#hmet2 = testEL.Histo1D(("","",100,-20,20),"Genelectron_eta")
#hmet3 = testTAU.Histo1D(("","",100,-20,20),"Gentau_eta")
#plot([hmet1,hmet2,hmet3],["Muon_eta","Electron_eta","Tau_eta"],"/gen_leptonETA.pdf",logY=True)

#hmet1 = testTAU.Histo1D(("","",100,0,300),"leptonTAU_pt")
#plot([hmet1],["Lepton_tau"],"/gen_leptonTAUpt.pdf",logY=True)


######## Scatter plot

c = ROOT.TCanvas('c','c',800,700)

h2 = df4.Histo2D(("h2", "scatter plot", 100, 0.0, 400.0, 50, 0.0, 10.0),"W_hadronic_pt","quarks_dr")
h2.Draw("COLZ")
h2.SetStats(0);
h2.GetXaxis().SetTitle( 'Hadronic W pt' )
h2.GetYaxis().SetTitle( 'Quarks deltaR' )
c.SaveAs(plotdir+"scatterWptQdr.pdf")

######## W PT filters ###########

W_normal = df4.Filter('W_hadronic_pt < 400')
W_boosted = df4.Filter('W_hadronic_pt>300')

entries1 = df4.Filter('typeC==1 && muon_jetGen != -1').Count()
entries2 = df4.Filter('typeC==-1 && muon_jetGen != -1').Count()

print("%s semi entries have hadronic W with pT under 400" %W_normal.Count().GetValue())
print("%s semi entries have hadronic W with pT over 300" %W_boosted.Count().GetValue())

print("%s semi charm entries normal W have at least one charm jet" %W_normal.Filter('typeC==1').Filter('nCharmJets>0').Count().GetValue())
print("%s semi charm entries boosted W have at least one charm jet" %W_boosted.Filter('typeC==1').Filter('nCharmJets>0').Count().GetValue())

print("%s semi charm entries normal W have at least one charm fat jet" %W_normal.Filter('typeC==1').Filter('nCharmFatJets>0').Count().GetValue())
print("%s semi charm entries boosted W have at least one charm fat jet" %W_boosted.Filter('typeC==1').Filter('nCharmFatJets>0').Count().GetValue())

##### plots

hmet1 = W_normal.Filter('typeC==1').Histo1D("nGoodGenJets")
hmet2 = W_boosted.Filter('typeC==1').Histo1D("nGoodGenJets")
plot([hmet1,hmet2],["Charm Normal","Charm Boosted"],"/nGenJets_boosted.pdf",logY=True)

hmet1 = W_normal.Filter('typeC==1').Histo1D("nGoodFatGenJets")
hmet2 = W_boosted.Filter('typeC==1').Histo1D("nGoodFatGenJets")
plot([hmet1,hmet2],["Charm Normal","Charm Boosted"],"/nGenFatJets_boosted.pdf",logY=True)

brlist = ["typeWW","typeC","nJet","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop",
"FatJet_phi","FatJet_pt","nSV","nGenJetAK8","GenJetAK8_partonFlavour","GenJetAK8_hadronFlavour",
"isMuonGood","isJetGood","FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta","GenPart_mass","GenPart_phi","GenPart_pt",
"GenPart_genPartIdxMother","GenPart_pdgId","GenJet_partonFlavour","GenJet_hadronFlavour","GenJetAK8_eta","GenJetAK8_mass","GenJetAK8_phi","GenJetAK8_pt",
"Muon_eta","Muon_mass","Muon_phi","Muon_pt","Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area",
"FirstJetPt","MET_phi","MET_pt",
"nElectron","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","nElGood","isElGood","nFatJetGood","isFatJetGood"]

#W_normal.Snapshot("Events", "W_normal.root",brlist)
#W_boosted.Snapshot("Events", "W_boosted.root",brlist)
