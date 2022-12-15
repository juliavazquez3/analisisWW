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

df5 = df4.Define("Genmuon_phi",'fabs(GenPart_pdgId[partID[2]])==13 ?  GenPart_phi[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==13 ?  GenPart_phi[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==13 ?  GenPart_phi[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==13 ?  GenPart_phi[partID[5]] :-20)))')
df6 = df5.Define("Genelectron_phi",'fabs(GenPart_pdgId[partID[2]])==11 ?  GenPart_phi[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==11 ?  GenPart_phi[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==11 ?  GenPart_phi[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==11 ?  GenPart_eta[partID[5]] :-20)))')
df4 = df6.Define("Gentau_phi",'fabs(GenPart_pdgId[partID[2]])==15 ?  GenPart_phi[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==15 ?  GenPart_phi[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==15 ?  GenPart_phi[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==15 ?  GenPart_phi[partID[5]] :-20)))')

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

entries1 = testTAU.Filter('fabs(Gentau_eta)<2.4 && Gentau_pt>30').Count()
entries2 = testMUON.Filter('fabs(Genmuon_eta)<2.4 && Genmuon_pt>30').Count()
entries3 = testEL.Filter('fabs(Genelectron_eta)<2.4 && Genelectron_pt>30').Count()
#print("%s entries have a proper tau" %entries1.GetValue())
#print("%s entries have a proper muon" %entries2.GetValue())
#print("%s entries have a proper electron" %entries3.GetValue())

## We filter to maintain only events with a proper muon or electron

df_lepton = df4.Filter('(fabs(Genmuon_eta)<2.4 && Genmuon_pt>30) || (fabs(Genelectron_eta)<2.4 && Genelectron_pt>30)')

print("%s is the quantity we are left with, one good muon or electron" %df_lepton.Count().GetValue())

######### charm particle study #########

df_aux = df_lepton.Define("charm_pt",'fabs(GenPart_pdgId[partID[2]])==4 ?  GenPart_pt[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==4 ?  GenPart_pt[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==4 ?  GenPart_pt[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==4 ?  GenPart_pt[partID[5]] :-20)))')
df_aux = df_aux.Define("charm_eta",'fabs(GenPart_pdgId[partID[2]])==4 ?  GenPart_eta[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==4 ?  GenPart_eta[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==4 ?  GenPart_eta[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==4 ?  GenPart_eta[partID[5]] :-20)))')
df_aux = df_aux.Define("charm_phi",'fabs(GenPart_pdgId[partID[2]])==4 ?  GenPart_phi[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==4 ?  GenPart_phi[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==4 ?  GenPart_phi[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==4 ?  GenPart_phi[partID[5]] :-20)))')
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto BinStatus( Vint pdgId, Vint partID, Vint GenPart_statusFlags) {
            string binS;
            if (fabs(pdgId[partID[2]])==4) {
                binS = std::bitset<15>(GenPart_statusFlags[partID[2]]).to_string();
            }
            else if (fabs(pdgId[partID[3]])==4) {
                binS = std::bitset<15>(GenPart_statusFlags[partID[3]]).to_string();
            }
            else if (fabs(pdgId[partID[4]])==4) {
                binS = std::bitset<15>(GenPart_statusFlags[partID[4]]).to_string();
            }
            else if (fabs(pdgId[partID[5]])==4) {
                binS = std::bitset<15>(GenPart_statusFlags[partID[5]]).to_string();
            }
            return binS;
      };
""")

df_aux = df_aux.Define("charm_status",'BinStatus(GenPart_pdgId,partID, GenPart_statusFlags)')
df_aux = df_aux.Define("charm_status_13",'charm_status[0]')
############ GEN JETS ##############

gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using namespace std;
      vector<bool> isGoodJets(UInt_t njet, Vfloat pt, Vfloat eta) {
            vector<bool> vb;
            for (unsigned int i=0; i<njet; ++i) {
                vb.push_back(pt[i]>30. && fabs(eta[i])<2.5);
            }
            return vb;
      };
      auto jetFlavour(UInt_t nGenJ,Vint flavour,Vfloat pt, Vfloat eta, Vfloat jet_phi, float muon_pt, float el_pt, float muon_eta, float muon_phi, float el_eta, float el_phi){
            int ind = 0;
            for (unsigned int i=0; i<nGenJ; ++i) {
                if (fabs(flavour[i])==4 && pt[i]>30 && fabs(eta[i])<2.5) {
                        if (muon_pt < 0 && ROOT::VecOps::DeltaR(el_eta,eta[i] , el_phi, jet_phi[i])>0.5) {
                                ind++;
                        }
                        else if (el_pt < 0 && ROOT::VecOps::DeltaR(muon_eta,eta[i] , muon_phi, jet_phi[i])>0.5) {
                                ind++;
                        }
                }
            }
            return ind;
      };
      auto GoodJetsIDs(UInt_t njet, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, float muon_pt, float el_pt, float muon_eta, float muon_phi, float el_eta, float el_phi) {
            vector<int> vb;
            for (unsigned int i=0; i<njet; ++i) {
                if (jet_pt[i]>30 && fabs(jet_eta[i])<2.5) {
			if (muon_pt < 0 && ROOT::VecOps::DeltaR(el_eta,jet_eta[i] , el_phi, jet_phi[i])>0.5) {
				vb.push_back(i);
	                }
                        else if (el_pt < 0 && ROOT::VecOps::DeltaR(muon_eta,jet_eta[i] , muon_phi, jet_phi[i])>0.5) {
                                vb.push_back(i);
                        }
                }
            }
            return vb;
      };
      auto charmJetMinDelta(UInt_t njet, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, float charm_pt, float charm_eta, float charm_phi, float muon_pt, float el_pt, float muon_eta, float muon_phi, float el_eta, float el_phi) {
            auto dR{10.};
	    int ind = -1;
            for (unsigned int i=0; i<njet; ++i) {
                if (jet_pt[i]>30 && fabs(jet_eta[i])<2.5 && (ROOT::VecOps::DeltaR(charm_eta,jet_eta[i] , charm_phi, jet_phi[i])<dR)) {
                        if (muon_pt < 0 && ROOT::VecOps::DeltaR(el_eta,jet_eta[i] , el_phi, jet_phi[i])>0.5) {
				dR = ROOT::VecOps::DeltaR(charm_eta,jet_eta[i] , charm_phi, jet_phi[i]);
				ind = i;
                        }
                        else if (el_pt < 0 && ROOT::VecOps::DeltaR(muon_eta,jet_eta[i] , muon_phi, jet_phi[i])>0.5) {
				dR = ROOT::VecOps::DeltaR(charm_eta,jet_eta[i] , charm_phi, jet_phi[i]);
				ind = i;
                        }
                }
            }
            return ind;
      };

""")

my_call = 'jetFlavour(nGenJet,GenJet_partonFlavour, GenJet_pt, GenJet_eta, GenJet_phi, Genmuon_pt, Genelectron_pt, Genmuon_eta, Genmuon_phi, Genelectron_eta, Genelectron_phi)' 
my_call2 = 'isGoodJets(nGenJet,GenJet_pt,GenJet_eta)' 


df5 = df_aux.Define("nCharmGenJets",my_call) 
df6 = df5.Define("nGoodGenJets",'GoodJetsIDs(nGenJet, GenJet_pt, GenJet_eta, GenJet_phi, Genmuon_pt, Genelectron_pt, Genmuon_eta, Genmuon_phi, Genelectron_eta, Genelectron_phi).size()') 
df4 = df6.Define("isGenJetGood",my_call2)

df5 = df4.Define("nCharmFatJets",'jetFlavour(nGenJetAK8,GenJetAK8_partonFlavour, GenJetAK8_pt, GenJetAK8_eta, GenJetAK8_phi, Genmuon_pt, Genelectron_pt, Genmuon_eta, Genmuon_phi, Genelectron_eta, Genelectron_phi)')
df6 = df5.Define("nGoodFatGenJets",'GoodJetsIDs(nGenJetAK8, GenJetAK8_pt, GenJetAK8_eta, GenJetAK8_phi, Genmuon_pt, Genelectron_pt, Genmuon_eta, Genmuon_phi, Genelectron_eta, Genelectron_phi).size()')
df4 = df6.Define('GoodGenJetsID','GoodJetsIDs(nGenJet, GenJet_pt, GenJet_eta, GenJet_phi, Genmuon_pt, Genelectron_pt, Genmuon_eta, Genmuon_phi, Genelectron_eta, Genelectron_phi)')

df4 = df4.Filter('nGoodGenJets > 0')

print("%s entries have at least 1 good jet" %df4.Count().GetValue())

#print("%s entries have at least 1 good charm jet" %df4.Filter('nCharmGenJets > 0').Count().GetValue())

df5 = df4.Define('GenJetdeltaRlepton','Genmuon_pt < 0 ? (ROOT::VecOps::DeltaR(Genelectron_eta, GenJet_eta[GoodGenJetsID[0]], Genelectron_phi, GenJet_phi[GoodGenJetsID[0]])) : (ROOT::VecOps::DeltaR(Genmuon_eta, GenJet_eta[GoodGenJetsID[0]], Genmuon_phi, GenJet_phi[GoodGenJetsID[0]]))')
df6 = df5.Define('CharmJetInd', 'charmJetMinDelta(nGenJet, GenJet_pt, GenJet_eta, GenJet_phi,charm_pt, charm_eta, charm_phi, Genmuon_pt, Genelectron_pt, Genmuon_eta, Genmuon_phi, Genelectron_eta, Genelectron_phi)')
df7 = df6.Define('CharmJetDeltaR','CharmJetInd<0 ? 10 : ROOT::VecOps::DeltaR(charm_eta,GenJet_eta[CharmJetInd] , charm_phi, GenJet_phi[CharmJetInd])')
df8 = df7.Define('CharmJetPartFlavour','CharmJetInd<0 ? 10 : (CharmJetDeltaR<0.4 ? GenJet_partonFlavour[CharmJetInd] : 10)')
df4 = df8.Define('CharmJetHadrFlavour','CharmJetInd<0 ? 10 : (CharmJetDeltaR<0.4 ? GenJet_hadronFlavour[CharmJetInd] : 10)')

#print("%s entries have at least 1 good jet associated to a charm gen particle" %df4.Filter('CharmJetPartFlavour < 6').Count().GetValue())
#print("%s entries have at least 1 good charm jet associated to a charm gen particle" %df4.Filter('fabs(CharmJetPartFlavour) == 4').Count().GetValue())
####### Flavour var study  #############
# I personally would not use it.

#jet_none = df4.Filter('nGoodGenJets == 0')
#jet_one = df4.Filter('nGoodGenJets == 1')
#jet_two = df4.Filter('nGoodGenJets == 2')
#jet_three = df4.Filter('nGoodGenJets == 3')

#print("%s entries dont have good jets" %jet_none.Count().GetValue())
#print("%s entries have 1 good jet" %jet_one.Count().GetValue())
#print("%s entries have 2 good jets" %jet_two.Count().GetValue())
#print("%s entries have 3 good jets" %jet_three.Count().GetValue())

#jet_one = jet_one.Define('jetG_flavour1','GenJet_partonFlavour[GoodGenJetsID[0]]')
#jet_two = jet_two.Define('jetG_flavour1','GenJet_partonFlavour[GoodGenJetsID[0]]')
#jet_two = jet_two.Define('jetG_flavour2','GenJet_partonFlavour[GoodGenJetsID[1]]')
#jet_three = jet_three.Define('jetG_flavour1','GenJet_partonFlavour[GoodGenJetsID[0]]')
#jet_three = jet_three.Define('jetG_flavour2','GenJet_partonFlavour[GoodGenJetsID[1]]')
#jet_three = jet_three.Define('jetG_flavour3','GenJet_partonFlavour[GoodGenJetsID[2]]')

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
	auto indMuonJet(UInt_t nGen, Vint pdg, Vint mother, Vfloat part_pt, Vfloat part_eta, Vfloat part_phi, Vfloat jet_eta, Vfloat jet_phi, Vint GoodJets) {
		int ind = -1;
	        float pT{-10.};
                for (unsigned int i=0; i<nGen; ++i) {
			for (unsigned int j=0; j<(GoodJets.size()); ++j) {
			auto cond1 = ROOT::VecOps::DeltaR(part_eta[i],jet_eta[GoodJets[j]] , part_phi[i], jet_phi[GoodJets[j]])<0.25;
                        	if (fabs(pdg[i])==13 && cond1 && pT < part_pt[i]) {
                                	ind = i;
					pT = part_pt[i];
                        	}
			}
                }
		return ind;
	};
""")

df5 = df4.Define("muon_jetGen",'indMuonJet(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, GenJet_eta, GenJet_phi, GoodGenJetsID)')
df4 = df5.Define("muon_jetGen_pt",'muon_jetGen != -1 ? GenPart_pt[muon_jetGen] : -10')

#entries1 = df4.Filter('typeC==1 && muon_jetGen_pt > 0').Count()
#entries2 = df4.Filter('typeC==-1 && muon_jetGen_pt > 0').Count()

#print("%s charm entries have muon from gen jet charm" %entries1.GetValue())
#print("%s not charm entries have muon from gen jet charm" %entries2.GetValue())

############### Hadronic W pt dependent study #################

W_normal = df4.Filter('W_hadronic_pt <200')
W_boosted = df4.Filter('W_hadronic_pt>200')

print("%s semi entries with 1 good jet have hadronic W with pT under 200" %W_normal.Count().GetValue())
print("%s semi entries with 1 good jet have hadronic W with pT over 200" %W_boosted.Count().GetValue())

print("%s semi charm entries normal W have muon inside jet" %W_normal.Filter('typeC==1').Filter('muon_jetGen_pt >0').Count().GetValue())
print("%s semi not charm entries normal W have muon inside jet" %W_normal.Filter('typeC==-1').Filter('muon_jetGen_pt >0').Count().GetValue())

print("%s semi charm entries boosted W have muon inside jet" %W_boosted.Filter('typeC==1').Filter('muon_jetGen_pt >0').Count().GetValue())
print("%s semi not charm entries boosted W have muon inside jet" %W_boosted.Filter('typeC==-1').Filter('muon_jetGen_pt >0').Count().GetValue())

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


charm = df_aux.Filter('typeC==1')
notCharm = df_aux.Filter('typeC==-1')

########### charms ##############3

#hmet1 = charm.Histo1D(("","",50,0,600),"charm_pt")
#hmet2 = notCharm.Histo1D(("","",50,0,600),"charm_pt")
#plot([hmet1,hmet2],["Charm","Not charm"],"/charm_pt_filtered.pdf",logY=True)

#hmet1 = charm.Histo1D("charm_eta")
#hmet2 = notCharm.Histo1D("charm_eta")
#plot([hmet1,hmet2],["Charm","Not charm"],"/charm_eta_filtered.pdf",logY=True)

#hmet1 = charm.Histo1D("charm_status_13")
#hmet2 = notCharm.Histo1D("charm_status_13")
#plot([hmet1,hmet2],["Charm","Not charm"],"/charm_status_filtered.pdf",logY=True)

charm = df4.Filter('typeC==1')
notCharm = df4.Filter('typeC==-1')


######## GEN JETS #########

#hmet1 = charm.Histo1D(("","",10,0,10),"nCharmJets")
#hmet2 = notCharm.Histo1D(("","",10,0,10),"nCharmJets")
#plot([hmet1,hmet2],["Charm","Not charm"],"/nCharmJets.pdf",logY=True)

#hmet1 = charm.Histo1D(("","",10,0,10),"nGoodGenJets")
#hmet2 = notCharm.Histo1D(("","",10,0,10),"nGoodGenJets")
#plot([hmet1,hmet2],["Charm","Not charm"],"/nGoodGenJets_filtered.pdf",logY=True)

#hmet = df4.Histo1D("GenJetdeltaRlepton")
#plot([hmet],["deltaR from lepton"],"/GenJetdeltaRlepton_firstjet_filtered.pdf", logY = True)

#hmet1 = charm.Histo1D(("","",50,0,6),"CharmJetDeltaR")
#hmet2 = notCharm.Histo1D(("","",50,0,6),"CharmJetDeltaR")
#plot([hmet1,hmet2],["Charm","Not charm"],"/Gencharm_jet_deltaR.pdf",logY=True)

#hmet1 = charm.Histo1D(("","",12,-6,6),"CharmJetPartFlavour")
#hmet2 = notCharm.Histo1D(("","",12,-6,6),"CharmJetPartFlavour")
#plot([hmet1,hmet2],["Charm","Not charm"],"/Gencharm_jet_flavourP.pdf",logY=True)


# Flavour var tests

#hmet1 = jet_one.Histo1D(("","",20,-10,10),"jetG_flavour1")
#hmet2 = jet_two.Histo1D(("","",20,-10,10),"jetG_flavour1")
#hmet3 = jet_two.Histo1D(("","",20,-10,10),"jetG_flavour2")
#hmet4 = jet_three.Histo1D(("","",20,-10,10),"jetG_flavour1")
#hmet5 = jet_three.Histo1D(("","",20,-10,10),"jetG_flavour2")
#hmet6 = jet_three.Histo1D(("","",20,-10,10),"jetG_flavour3")

#plot([hmet4,hmet5,hmet6],["3 jets 1","3 jets 2","3 jets 3"],"/GoodJetsFlavours_3jets_filtered_Charm.pdf")
#plot([hmet2,hmet3],["2 jets 1","2 jets 2"],"/GoodJetsFlavours_2jets_filtered_Charm.pdf")

#c = ROOT.TCanvas('c','c',800,700)

#h2 = jet_two.Histo2D(("h2", "scatter plot", 20,-10 , 10, 20, -10, 10),"jetG_flavour1","jetG_flavour2")
#h2.Draw("COLZ")
#h2.SetStats(0);
#h2.GetXaxis().SetTitle( 'Higher pt good jet' )
#h2.GetYaxis().SetTitle( 'Lower pt gen jet' )
#c.SaveAs(plotdir+"scatter2genjets.pdf")

######## muon from jet

hmet1 = charm.Histo1D(("","",50,0,150),"muon_jetGen_pt")
hmet2 = notCharm.Histo1D(("","",50,0,150),"muon_jetGen_pt")
plot([hmet1,hmet2],["Charm","Not charm"],"/jetmuon_pt_filtered.pdf",logY=True)

brlist1 = ["charm_status","charm_status_13","nGenPart","GenPart_genPartIdxMother","GenPart_pdgId","CharmJetDeltaR"]

brlist = ["typeWW","typeC","nJet","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop",
"FatJet_phi","FatJet_pt","nSV","nGenJetAK8","GenJetAK8_partonFlavour","GenJetAK8_hadronFlavour",
"isMuonGood","isJetGood","FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta","GenPart_mass","GenPart_phi","GenPart_pt",
"GenPart_genPartIdxMother","GenPart_pdgId","GenJet_partonFlavour","GenJet_hadronFlavour","GenJetAK8_eta","GenJetAK8_mass","GenJetAK8_phi","GenJetAK8_pt",
"Muon_eta","Muon_mass","Muon_phi","Muon_pt","Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area",
"FirstJetPt","MET_phi","MET_pt",
"nElectron","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","nElGood","isElGood","nFatJetGood","isFatJetGood"]

#W_normal.Snapshot("Events", "W_normal.root",brlist)
#W_boosted.Snapshot("Events", "W_boosted.root",brlist)

#df4.Snapshot("Events", "charm_aux.root",brlist1)
