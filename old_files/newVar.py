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
fileName = "/nfs/cms/vazqueze/analisisWW/old_files/MyWW.root"
treeName = "Events"
#ROOT.fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, fileName)

#df = d.Filter('nMuonGood>0 && nJetGood >0')


##df = d_filt.Define('deltaR_jetmuon','ROOT::VecOps::DeltaR(jet_eta1,muon_eta1,jet_phi1,muon_phi1)') 

## Funcion auxiliar para DeltaR 

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      auto InvariantM(const float pt, const float eta, const float phi, const float mass) {
            auto x = pt*std::cos(phi);
            auto y = pt*std::sin(phi);
            auto z = pt*std::sinh(eta);
            auto e = std::sqrt(x*x+y*y+z*z+mass*mass);
            auto mJet = std::sqrt(e*e-x*x-y*y-z*z);
            return mJet;
      };
      auto invMassJet(UInt_t nmu, UInt_t njet, Vfloat muon_eta,Vfloat muon_phi, Vfloat jet_eta, Vfloat jet_phi,Vfloat muon_pt,Vfloat muon_mass, Vfloat jet_pt, Vfloat jet_mass) {
            auto mJet{-10.};
	    if (njet>0){
                int i = ROOT::VecOps::ArgMax(jet_pt);
                mJet = InvariantM(jet_pt[i],jet_eta[i],jet_phi[i],jet_mass[i]);
            }
            return mJet;
      };
      auto invMass2ptJet(UInt_t njet,Vfloat jet_eta, Vfloat jet_phi,Vfloat jet_pt, Vfloat jet_mass, Vbool jet_good) {
            auto mJet{-10.};
            if (njet>1){
		int i1 = ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(jet_pt))[0];
                int i2 = ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(jet_pt))[1];
		if (jet_good[i1] && jet_good[i2]){
			auto x1 = jet_pt[i1]*std::cos(jet_phi[i1]);
            		auto y1 = jet_pt[i1]*std::sin(jet_phi[i1]);
            		auto z1 = jet_pt[i1]*std::sinh(jet_eta[i1]);
            		auto e1 = std::sqrt(x1*x1+y1*y1+z1*z1+jet_mass[i1]*jet_mass[i1]);

                	auto x2 = jet_pt[i2]*std::cos(jet_phi[i2]);
                        auto y2 = jet_pt[i2]*std::sin(jet_phi[i2]);
                        auto z2 = jet_pt[i2]*std::sinh(jet_eta[i2]);
                        auto e2 = std::sqrt(x2*x2+y2*y2+z2*z2+jet_mass[i2]*jet_mass[i2]);
			
			mJet = std::sqrt((e1+e2)*(e1+e2) - (x1+x2)*(x1+x2) - (y1+y2)*(y1+y2) - (z1+z2)*(z1+z2));
		}
	    }
            return mJet;
      };
""")

######## MUONES ##########

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      auto invMassMuon(UInt_t nmu, UInt_t njet, Vfloat muon_eta,Vfloat muon_phi, Vfloat jet_eta, Vfloat jet_phi,Vfloat muon_pt,Vfloat muon_mass, Vfloat jet_pt, Vfloat jet_mass) {
            auto mMuon{-10.};
            if (nmu>0){
                int i = ROOT::VecOps::ArgMax(muon_pt);
                mMuon = InvariantM(muon_pt[i],muon_eta[i],muon_phi[i],muon_mass[i]);
            }
            return mMuon;
      };
      auto TransverseMass(UInt_t nmu,Vbool muon_good,Vfloat muon_phi, Vfloat muon_pt,Vfloat muon_mass, float met_phi, float met_pt) {
            auto mt{-10.};
	    auto pt{-10.};
	    for (unsigned int i=0; i<nmu; ++i) {
            	if (muon_good[i]&&(pt<muon_pt[i])){
			pt = muon_pt[i];
			mt = std::sqrt(2*muon_pt[i]*met_pt*(1-std::cos(muon_phi[i]-met_phi)));
            	}
	    }
            return mt;
      };
      auto countG(UInt_t n, Vbool good) {
	    int c = 0;
	    for (unsigned int i=0; i<n; ++i) {
                if (good[i]){
		    c++;
                }
            }
            return c;
      };
""")

my_call1 = "countG(nMuon,isMuonGood)"
my_call2 = "countG(nElectron, isElGood)"
my_call3 = "countG(nJet,isJetGood)"

dft2 = d.Define("nMuonGood",my_call1)
dft = dft2.Define("nElectronGood",my_call2)
d = dft.Define("nJetGood",my_call3)


########### ELECTRONES ##########

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      auto TransverseMassEl(UInt_t nel,Vbool el_good,Vfloat el_phi, Vfloat el_pt,Vfloat el_mass, float met_phi, float met_pt) {
            auto mt{-10.};
            auto pt{-10.};
            for (unsigned int i=0; i<nel; ++i) {
                if (el_good[i]&&(pt<el_pt[i])){
                        pt = el_pt[i];
                        mt = std::sqrt(2*el_pt[i]*met_pt*(1-std::cos(el_phi[i]-met_phi)));
                }
            }
            return mt;
      };
""")

# All input variables for the method present in the tree
my_call1 = "invMassMuon(nMuon,nJet,Muon_eta,Muon_phi,Jet_eta,Jet_phi,Muon_pt,Muon_mass,Jet_pt,Jet_mass)"
my_call2 = "invMassJet(nMuon,nJet,Muon_eta,Muon_phi,Jet_eta,Jet_phi,Muon_pt,Muon_mass,Jet_pt,Jet_mass)"
my_call3 = "invMass2ptJet(nJet,Jet_eta,Jet_phi,Jet_pt,Jet_mass,isJetGood)"
my_call4 = "TransverseMass(nMuon,isMuonGood,Muon_phi,Muon_pt,Muon_mass,MET_phi,MET_pt)"
my_call5 = "TransverseMassEl(nElectron,isElGood,Electron_phi,Electron_pt,Electron_mass,MET_phi,MET_pt)"
my_call6 = "invMassJet(nMuon,nFatJet,Muon_eta,Muon_phi,FatJet_eta,FatJet_phi,Muon_pt,Muon_mass,FatJet_pt,FatJet_mass)"
my_call7 = "invMass2ptJet(nFatJet,FatJet_eta,FatJet_phi,FatJet_pt,FatJet_mass,isFatJetGood)"

dft2 = d.Define("FatJet_invM2",my_call7)
dft = dft2.Define("FatJet_invM",my_call6)
det = dft.Define("TransverseMassEl",my_call5)
dmt = det.Define("TransverseMass",my_call4)
d1 = dmt.Define("ImassMuon", my_call1)
d2 = d1.Define("ImassJet", my_call2)
d4 = d2.Define("Inv2massJetFiltered", my_call3)

d3 = d4.Define("FatJet_ptmax",'nFatJet>0 ? ROOT::VecOps::Max(FatJet_pt) : -10')

#good_df.Snapshot("Events","Filtered.root")

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      auto deltaR(UInt_t nmu, UInt_t njet, Vfloat muon_eta,Vfloat muon_phi, Vfloat jet_eta, Vfloat jet_phi, Vfloat jet_pt,Vbool jet_good, Vbool muon_good) {
            std::vector<float> vR;
	    float dR{10.};
	    int ind = -1;
            for (unsigned int i=0; i<nmu; ++i) {
                  dR = 10.;
		  for (unsigned int j=0; j<njet; ++j) {
			if(jet_good[j]){
                        	auto delt = ROOT::VecOps::DeltaR(jet_eta[j],muon_eta[i],jet_phi[j],muon_phi[i]);
                        	if (delt < dR){
                                	dR = delt;
					ind = i;
                        	}
			}
                  }
		  vR.push_back(dR);
            }
	    ROOT::VecOps::RVec<float> myRVec(vR.data(), vR.size());
	    return myRVec;
      };
""")

######### muon pt

# All input variables for the method present in the tree
my_call = "deltaR(nMuon,nJet,Muon_eta,Muon_phi,Jet_eta,Jet_phi,Jet_pt,isJetGood,isMuonGood)"
df1 = d3.Define("muonDR", my_call)
my_call = "deltaR(nElectron,nJet,Electron_eta,Electron_phi,Jet_eta,Jet_phi,Jet_pt,isJetGood,isElGood)"
df = df1.Define("electronDR", my_call)

df = df.Define("muondrSize",'muonDR.size()')

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      auto ElInd(UInt_t nel,Vbool electron_good, Vfloat electron_pt) {
            float pT{-10.};
            int ind = -1;
	    for (unsigned int i=0; i<nel; ++i) {
                if((electron_pt[i]>pT)&&(electron_good[i])) {
                        ind = i;
                        pT = electron_pt[i];
                }
            }
	    return ind;
      };
      auto MuInd(UInt_t nmu,Vbool muon_good, Vfloat muon_pt) { 
            float pT{-10.};
            int ind = -1;
            for (unsigned int i=0; i<nmu; ++i) {
                if((muon_pt[i]>pT)&&(muon_good[i])) {
                        ind = i;
                        pT = muon_pt[i];
                }
            }
	    return ind;
      };
      auto GoodJetsIDs(UInt_t njet, UInt_t nmu, UInt_t nel, Vbool jet_good, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vfloat muon_eta, Vfloat muon_phi, Vfloat electron_eta, Vfloat electron_phi, Vbool electron_good, Vbool muon_good, Vfloat electron_pt, Vfloat muon_pt) {
	    float pT{-10.};
            int ind = -1;
            int ind1 = MuInd(nmu,muon_good,muon_pt);
            int ind2 = ElInd(nel,electron_good,electron_pt);
            bool cond1 = false;
            bool cond2 = false;
            vector<int> vb;
	    for (unsigned int i=0; i<njet; ++i) {
                if (ind1 != -1) {
                        cond1 = ROOT::VecOps::DeltaR(muon_eta[ind1],jet_eta[i],muon_phi[ind1],jet_phi[i]) > 0.5;
                }
                if (ind2 != -1) {
                        cond2 = ROOT::VecOps::DeltaR(electron_eta[ind2],jet_eta[i],electron_phi[ind2],jet_phi[i]) > 0.5;
                }
                if((jet_pt[i]>pT)&&(jet_good[i])&&(cond1 || cond2)) {
                        vb.push_back(i);
                        pT = jet_pt[i];
                }
            }
            return vb;
      };
      auto JetsMuonsNAOD(UInt_t njet, Vfloat jet_pt, Vint jet_flavour, Vfloat jet_eta, Vint jet_cons ) {
            vector<int> vb;
            int indJ = -1;
            int indC = -1;
            int indGJ = -1;
            int indGC = -1;
            float jetPT1 = -10;
            float jetPT2 = -10;
            float jetPT3 = -10;
            float jetPT4 = -10;
            bool cond1 = false;
            bool cond2 = false;
            for (unsigned int i=0; i<njet; ++i) {
                cond1 = jet_pt[i]>30. && fabs(jet_eta[i])<2.5;
                cond2 = fabs(jet_flavour[i])==4;
                if (!cond2 && jet_pt[i]>jetPT1) {
                        indJ = jet_cons[i];
                        jetPT1 = jet_pt[i];
                }
                if (cond2 && jet_pt[i]>jetPT2) {
                        indC = jet_cons[i];
                        jetPT2 = jet_pt[i];
                }
                if (cond1 && !cond2 && jet_pt[i]>jetPT3) {
                        indGJ = jet_cons[i];
                        jetPT3 = jet_pt[i];
                }
                if (cond1 && cond2 && jet_pt[i]>jetPT4) {
                        indGC = jet_cons[i];
                        jetPT4 = jet_pt[i];
                }
            }
            vb.push_back(indJ);
            vb.push_back(indC);
            vb.push_back(indGJ);
            vb.push_back(indGC);
            return vb;
      };

      auto JetsTracks(UInt_t njet, Vfloat jet_pt, Vint jet_flavour, Vfloat jet_eta, Vint jet_cons ) {
            vector<int> vb;
            int indJ = 0;
            int indC = 0;
            int indGJ = 0;
            int indGC = 0;
            float jetPT1 = -10;
            float jetPT2 = -10;
            float jetPT3 = -10;
            float jetPT4 = -10;
	    bool cond1 = false;
	    bool cond2 = false;
            for (unsigned int i=0; i<njet; ++i) {
		cond1 = jet_pt[i]>30. && fabs(jet_eta[i])<2.5;
		cond2 = fabs(jet_flavour[i])==4;
		if (!cond2 && jet_pt[i]>jetPT1) {
                	indJ = jet_cons[i];
			jetPT1 = jet_pt[i];      
                }
                if (cond2 && jet_pt[i]>jetPT2) {
                        indC = jet_cons[i];
                        jetPT2 = jet_pt[i];
                }
                if (cond1 && !cond2 && jet_pt[i]>jetPT3) {
                        indGJ = jet_cons[i];
                        jetPT3 = jet_pt[i];
                }
                if (cond1 && cond2 && jet_pt[i]>jetPT4) {
                        indGC = jet_cons[i];
                        jetPT4 = jet_pt[i];
                }
            }
            vb.push_back(indJ);
            vb.push_back(indC);
            vb.push_back(indGJ);
            vb.push_back(indGC);
            return vb;
      };
      auto FinalInd1(UInt_t nmu, Vfloat muon_eta,Vfloat muon_phi, Vbool muon_iso, Vfloat muon_pt, Vfloat muonDR, UInt_t njet, Vfloat jet_pt, Vfloat jet_eta, Vfloat jet_phi, Vbool jet_good, Vint GoodJets) {
            float dR{10.};
            float pT{-10.};
            bool cond1 = false;
	    vector<int> vb;
            int ind = -1;
	    int indJ = -1;
	    float jetPT = -10;
            for (unsigned int i=0; i<nmu; ++i) {
                        for (unsigned int j=0; j<(GoodJets.size()); ++j) {
                        jetPT = jet_pt[GoodJets[j]];
			cond1 = ROOT::VecOps::DeltaR(muon_eta[i],jet_eta[GoodJets[j]] , muon_phi[i], jet_phi[GoodJets[j]])<0.4;
                                //if (cond1 && muon_pt[i]>pT && fabs(muon_eta[i])<2.1 && muon_pt[i]<25. && muon_pt[i]/jetPT < 0.6) {
                                if (cond1 && muon_pt[i]>pT && muon_iso[i]<5.) {
                                        ind = i;
                                        pT = muon_pt[i];
					indJ = j;
                                }
                        }
            }
	    vb.push_back(ind);
            vb.push_back(indJ);
            return vb;
      };
""")

df = df.Define("GoodJets_vector", 'GoodJetsIDs(nJet,nMuon,nElectron,isJetGood, Jet_pt,Jet_eta,Jet_phi,Muon_eta,Muon_phi,Electron_eta,Electron_phi,isElGood,isMuonGood,Electron_pt,Muon_pt)')

#my_call = 'FinalInd(nMuon, nElectron, Electron_eta, Electron_phi, Muon_eta, Muon_phi, isElGood, Electron_pt, isMuonGood, Muon_pt, muonDR, nJet, Jet_pt, Jet_eta, Jet_phi, isJetGood)'
my_call = 'FinalInd1(nMuon, Muon_eta, Muon_phi, Muon_pfRelIso04_all, Muon_pt, muonDR, nJet, Jet_pt, Jet_eta, Jet_phi, isJetGood, GoodJets_vector)'
my_call1 = 'MuInd(nMuon,isMuonGood,Muon_pt)'
my_call2 = 'ElInd(nElectron,isElGood, Electron_pt)'
df1 = df.Define("leadMuon",my_call1)
df2 = df1.Define("leadElectron",my_call2)
df = df2.Define("muonJetAux",my_call)

df = df.Define('tracks_vector','JetsTracks(nJet, Jet_pt, Jet_partonFlavour, Jet_eta, Jet_nConstituents)')
df2 = df.Define("jetsraw_ntracks", 'tracks_vector[0]')
df3 = df2.Define("charmjetsraw_ntracks", 'tracks_vector[1]')
df2 = df3.Define("jetsgood_ntracks", 'tracks_vector[2]')
df = df2.Define("charmjetsgood_ntracks", 'tracks_vector[3]')

df = df.Define('nMuNAOD_vector','JetsMuonsNAOD(nJet, Jet_pt, Jet_partonFlavour, Jet_eta, Jet_nMuons)')
df2 = df.Define("jetsraw_nMuNano", 'nMuNAOD_vector[0]')
df3 = df2.Define("charmjetsraw_nMuNano", 'nMuNAOD_vector[1]')
df2 = df3.Define("jetsgood_nMuNano", 'nMuNAOD_vector[2]')
df = df2.Define("charmjetsgood_nMuNano", 'nMuNAOD_vector[3]')

df1 = df.Define("jet_muonid_pre", 'muonJetAux[0]')
df2 = df1.Define("jet_muonid", 'jet_muonid_pre != -1 ? (Muon_pfRelIso04_all[jet_muonid_pre]>0.2 ? jet_muonid_pre : -1) : -1')
df3 = df2.Define("jet_electronid",'nElectron>0 ? (ROOT::VecOps::Min(electronDR)<0.2 ? ROOT::VecOps::ArgMin(electronDR) : -1) : -1')
df2 = df3.Define("PTmuon_jet", 'jet_muonid != -1 ? Muon_pt[jet_muonid] : -10')
df4 = df2.Define("PTelectron_jet", 'jet_electronid != -1 ? Electron_pt[jet_electronid] : -10')

df1 = df4.Define("ETA_muon_jet", 'jet_muonid != -1 ? Muon_eta[jet_muonid] : -10')
df2 = df1.Define("ISO_muon_jet", 'jet_muonid != -1 ? Muon_pfRelIso04_all[jet_muonid] : -10')
df4 = df2.Define("RELPT_muon_jet", 'jet_muonid != -1 ? Muon_pt[jet_muonid]/Jet_pt[muonJetAux[1]] : -10')

df4 = df4.Define('J_nConstituents_1','nJet>0 ? Jet_nConstituents[0] : -1')

##########

gInterpreter.Declare("""
      using Vbool = const ROOT::RVec<bool>&;
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      std::vector<int> chargeV(UInt_t nmu, UInt_t nel, Vint muon_charge,Vint electron_charge,Vbool electron_good, Vbool muon_good, int jet_muonid, int jet_elid, Vfloat muon_pt, Vfloat electron_pt) {
            auto muon_sub = -10;
	    auto muon_lead = -10;
	    auto electron_sub = -10;
            auto electron_lead = -10;
	    std::vector<int> v;
	    auto pt{-10};
	    for (unsigned int i=0; i<nmu; ++i) {
                if (muon_good[i]&&(pt<muon_pt[i])){
                        pt = muon_pt[i];
			muon_sub = muon_charge[jet_muonid];
	               	muon_lead = muon_charge[i];
                }
            }
	    pt = -10;
	    for (unsigned int i=0; i<nel; ++i) {
                if (electron_good[i]&&(pt<electron_pt[i])){
                        pt = electron_pt[i];
			electron_sub = electron_charge[jet_elid];
                        electron_lead = electron_charge[i];
                }
            }
	    v.push_back(muon_lead);
	    v.push_back(muon_sub);
	    v.push_back(electron_lead);
	    v.push_back(electron_sub);	    
            return v;
      };
      auto SSOS(Vint chargeme){
	auto fin = 0;
	if ((chargeme[0] == -chargeme[1])||(chargeme[0] == -chargeme[3])||(chargeme[2] == -chargeme[1])){
		fin = 1;
	}
	return fin;
      };
""")
my_call = "chargeV(nMuon,nElectron,Muon_charge,Electron_charge,isElGood,isMuonGood,jet_muonid,jet_electronid,Muon_pt,Electron_pt)"
df = df4.Define("chargeME", my_call)

my_call = "SSOS(chargeME)"
df = df.Define("SSOS", my_call)

df = df.Define("charge_leadMuon",'chargeME[0]')
df = df.Define("charge_subMuon",'chargeME[1]')
df = df.Define("charge_leadElectron",'chargeME[2]')
df = df.Define("charge_subElectron",'chargeME[3]')

df = df.Define("leadptmuoncharge",'nMuon>0 ? Muon_charge[ROOT::VecOps::ArgMax(Muon_pt)]:-10')

leptonic = df.Filter("typeWW==3")
hadronic = df.Filter("typeWW==1")
semi = df.Filter("typeWW==2")

entries1 = hadronic.Count()
entries2 = semi.Count()
entries3 = leptonic.Count()
#print("%s entries passed hadronic filter" %entries1.GetValue())
print("%s entries passed semi filter" %entries2.GetValue())
#print("%s entries passed leptonic filter" %entries3.GetValue())

semiCharm = semi.Filter("typeC==1")
semiNoCharm = semi.Filter("typeC==-1")

entries1 = semiCharm.Count()
entries2 = semiNoCharm.Count()

#print("%s entries passed Charm filter" %entries1.GetValue())
#print("%s entries did not pass Charm filter" %entries2.GetValue())


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


## Pintamos varias comparaciones

#hmet1 = semiCharm.Histo1D(("","",10,-1,9),"jet_muonid")
#hmet2 = semiNoCharm.Histo1D(("","",10,-1,9),"jet_muonid")
#hmet3 = semiCharm.Histo1D(("","",10,-1,9),"jet_electronid")
#hmet4 = semiNoCharm.Histo1D(("","",10,-1,9),"jet_electronid")
#plot([hmet1,hmet2,hmet3,hmet4],["CharmM","NoCharmM","CharmE","NoCharmE"],"/deltaSemi.pdf")

#hmet1 = semiCharm.Histo1D(("","",50,0,50),"PTmuon_jet")
#hmet2 = semiNoCharm.Histo1D(("","",50,0,50),"PTmuon_jet")
#hmet3 = semiCharm.Histo1D(("","",310,-10,300),"PTelectron_jet")
#hmet4 = semiNoCharm.Histo1D(("","",310,-10,300),"PTelectron_jet")
#plot([hmet1,hmet2],["CharmM","NoCharmM"],"/jet_muonpt_reco.pdf",logY=True)

#hmet1 = semiCharm.Histo1D(("","",50,-3,3),"ETA_muon_jet")
#hmet2 = semiNoCharm.Histo1D(("","",50,-3,3),"ETA_muon_jet")
#plot([hmet1,hmet2],["CharmM","NoCharmM"],"/jet_muonETA_reco.pdf",logY=True)

#hmet1 = semiCharm.Histo1D(("","",50,0,25),"ISO_muon_jet")
#hmet2 = semiNoCharm.Histo1D(("","",50,0,25),"ISO_muon_jet")
#plot([hmet1,hmet2],["CharmM","NoCharmM"],"/jet_muonISO_reco.pdf",logY=True)

#hmet1 = semi.Histo1D("charge_leadMuon")
#hmet2 = semi.Histo1D("charge_subMuon")
#hmet3 = semi.Histo1D("charge_leadElectron")
#hmet4 = semi.Histo1D("charge_subElectron")
#plot([hmet1,hmet2,hmet3,hmet4],["LeadMuon","SubMuon","LeadElectron","SubElectron"],"/charge2.pdf",logY=True)

#hmet1 = hadronic.Histo1D("SSOS")
#hmet2 = semi.Histo1D("SSOS")
#hmet3 = leptonic.Histo1D("SSOS")
#hmet4 = semiCharm.Histo1D("SSOS")
#hmet5 = semiNoCharm.Histo1D("SSOS")
#plot([hmet1,hmet2,hmet3,hmet4,hmet5],["Hadronic","Semi","Leptonic","SemiCharm","SemiNoCharm"],"/SSOScharm2.pdf",logY=True)

#hmet1 = hadronic.Histo1D("muondrSize")
#hmet2 = semi.Histo1D("muondrSize")
#hmet3 = leptonic.Histo1D("muondrSize")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/muonDR.pdf",logY=True)

#######MUONES#########

#hmet1 = hadronic.Histo1D("nMuonGood")
#hmet2 = semi.Histo1D("nMuonGood")
#hmet3 = leptonic.Histo1D("nMuonGood")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/nGMuons.pdf")

#hmet1 = hadronic.Histo1D(("imassmuon","imassmuon",100,-0.5,0.7),"ImassMuon")
#hmet2 = semi.Histo1D(("imassmuon","imassmuon",100,-0.5,0.7),"ImassMuon")
#hmet3 = leptonic.Histo1D(("imassmuon","imassmuon",100,-0.5,0.7),"ImassMuon")

#hmet1 = hadronic.Histo1D("ImassMuon")
#hmet2 = semi.Histo1D("ImassMuon")
#hmet3 = leptonic.Histo1D("ImassMuon")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/clas_invMuon.pdf",logY=True)

#hmet1 = hadronic.Histo1D("nMuon")
#hmet2 = semi.Histo1D("nMuon")
#hmet3 = leptonic.Histo1D("nMuon")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/nMuon.pdf")


#hmet1 = hadronic.Histo1D(("Tmass","Tmass",200,-10,500),"TransverseMass")
#hmet2 = semi.Histo1D(("Tmass","Tmass",200,-10,500),"TransverseMass")
#hmet3 = leptonic.Histo1D(("Tmass","Tmass",200,-10,500),"TransverseMass")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/TransverseMass.pdf",logY=True)

#hmet1 = hadronic.Histo1D("leadptmuoncharge")
#hmet2 = semi.Histo1D("leadptmuoncharge")
#hmet3 = leptonic.Histo1D("leadptmuoncharge")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/exCharge.pdf")

## Scatter plot

c = ROOT.TCanvas('c','c',800,700)

h2 = semi.Histo2D(("h2", "scatter plot", 10, 0, 5, 10,0,5),"nMuon","nElectron")
h2.Draw("COLZ")
h2.SetStats(0);
h2.GetXaxis().SetTitle( 'Muons' )
h2.GetYaxis().SetTitle( 'Electrons' )
c.SaveAs(plotdir+"scatterUnfilteredME.pdf")

###### ELECTRONES ########

#hmet1 = hadronic.Histo1D("nElectronGood")
#hmet2 = semi.Histo1D("nElectronGood")
#hmet3 = leptonic.Histo1D("nElectronGood")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/nGElectrons.pdf")

#hmet1 = hadronic.Histo1D(("Tmass","Tmass",200,-10,500),"TransverseMassEl")
#hmet2 = semi.Histo1D(("Tmass","Tmass",200,-10,500),"TransverseMassEl")
#hmet3 = leptonic.Histo1D(("Tmass","Tmass",200,-10,500),"TransverseMassEl")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/TransverseMassEl.pdf",logY=True)

####### FAT JETS ########

#hmet1 = hadronic.Histo1D("nFatJet")
#hmet2 = semi.Histo1D("nFatJet")
#hmet3 = leptonic.Histo1D("nFatJet")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/clas_nFatJet.pdf",logY=True)

#hmet1 = hadronic.Histo1D("FatJet_ptmax")
#hmet2 = semi.Histo1D("FatJet_ptmax")
#hmet3 = leptonic.Histo1D("FatJet_ptmax")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/FatJetpt.pdf",logY=True)

#hmet1 = hadronic.Histo1D(("invJet","invJet",210,-10,200),"FatJet_invM")
#hmet2 = semi.Histo1D(("invJet","invJet",210,-10,200),"FatJet_invM")
#hmet3 = leptonic.Histo1D(("invJet","invJet",210,-10,200),"FatJet_invM")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/FatJetINVM.pdf",logY=True)

#hmet1 = hadronic.Histo1D(("invJet","invJet",201,-10,2000),"FatJet_invM2")
#hmet2 = semi.Histo1D(("invJet","invJet",201,-10,200),"FatJet_invM2")
#hmet3 = leptonic.Histo1D(("invJet","invJet",201,-10,2000),"FatJet_invM2")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/FatJet_invm2.pdf",logY=True)


#hmet1 = semiCharm.Histo1D("nFatJetGood")
#hmet2 = semiNoCharm.Histo1D("nFatJetGood")
#plot([hmet1,hmet2],["Charm","NoCharm"],"/clas_nFatJetGood_semiC.pdf",logY=True)

## Saving Fat Jet inv mass

#semi.Snapshot("Events", "FatJet_invM_Semi.root","FatJet_invM")

############# JETS #############

#hmet1 = hadronic.Histo1D("subJet")
#hmet2 = semi.Histo1D("subJet")
#hmet3 = leptonic.Histo1D("subJet")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/Jetsub.pdf")


#hmet1 = hadronic.Histo1D("leadJet")
#hmet2 = semi.Histo1D("leadJet")
#hmet3 = leptonic.Histo1D("leadJet")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/Jetlead.pdf")

#hmet1 = hadronic.Histo1D(("invJet","invJet",100,0,50),"ImassJet")
#hmet2 = semi.Histo1D(("invJet","invJet",100,0,50),"ImassJet")
#hmet3 = leptonic.Histo1D(("invJet","invJet",100,0,50),"ImassJet")

#hmet1 = hadronic.Histo1D(("invJet","invJet",200,-10,500),"Inv2massJetFiltered")
#hmet2 = semi.Histo1D(("invJet","invJet",200,-10,500),"Inv2massJetFiltered")
#hmet3 = leptonic.Histo1D(("invJet","invJet",200,-10,500),"Inv2massJetFiltered")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/clas_invJet2Filtered.pdf",logY=True)

#hmet1 = semiCharm.Histo1D(("imassmuon","imassmuon",100,-0.5,0.7),"ImassMuon")
#hmet2 = semiNoCharm.Histo1D(("imassmuon","imassmuon",100,-0.5,0.7),"ImassMuon")
#plot([hmet1,hmet2],["Charm","NoCharm"],"/clas_invMuon_semiC.pdf",logY=True)

#hmet1 = semiCharm.Histo1D(("","",10,0,10),"nJetGood")
#hmet2 = semiNoCharm.Histo1D(("","",10,0,10),"nJetGood")
#plot([hmet1,hmet2],["Charm","NoCharm"],"/nJetGood_semi.pdf")


############ Save #################

brlist1 = ["muonsize","muonDR","electronDR"]
#df.Snapshot("Events", "wwDR.root",brlist1)

brlist = ["typeWW","typeC","nJet","nMuonGood","nJetGood","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop",
"FatJet_phi","FatJet_pt","nSV","isMuonGood","isJetGood","FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta",
"GenPart_mass","GenPart_phi","GenPart_pt","GenPart_genPartIdxMother","GenPart_pdgId","Muon_eta","Muon_mass","Muon_phi","Muon_pt","Jet_partonFlavour",
"Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area","FirstJetPt","MET_phi","MET_pt","Jet_puId",
"nElectron","nElectronGood","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","nElGood","isElGood","nFatJetGood","isFatJetGood",
"jet_muonid","leadMuon","leadElectron","chargeME","SSOS","FatJet_invM","TransverseMassEl","TransverseMass","ImassJet","Inv2massJetFiltered", "GoodJets_vector",
"jetsraw_nMuNano","charmjetsraw_nMuNano","jetsgood_nMuNano","charmjetsgood_nMuNano"]

semi.Snapshot("Events", "SemiNew.root",brlist)

