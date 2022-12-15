import ROOT, os, sys
from ROOT import *
#from ROOT import VecOps
import argparse
import json 
from os import listdir
from os.path import isfile, join

#if not sys.flags.interactive: ROOT.EnableImplicitMT()

plotdir = '/nfs/cms/vazqueze/analisisWW/plots/genWW/' # this is where we'll save your plots
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

#ROOT.EnableImplicitMT()

files = json.load(open("/nfs/cms/vazqueze/analisisWW/mc_info.json"))
processes = files.keys()

archives=[]

for p in processes:
    # Construct the dataframes
    folder = files[p]["folder_name"] # Folder name
    filename = files[p]["sample"] # Sample name
    num_events = files[p]["events"] # Number of events
    num_files = files[p]["files"] # Number of files
    list_files = [f for f in listdir(folder) if isfile(join(folder, f))] # Lista de archivos
    #print(len(list_files))
    if files[p]["type"]=="WW2016":
      if (num_files == len(list_files)):
        for f in list_files:
          archives.append(join(folder,f))

df = ROOT.RDataFrame("Events",set(archives))
print("Number of files for WW",len(archives))


########GEN_PART###########

# Attempt to classify with GenPart
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;
 
// Function to print the
// index of an element
// 1: Full hadronic 2: One W hadronic, the other leptonic (SEMI) 3: Full leptonic
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


# Attempt to classify with GenPart
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;

// Function to print the
// index of an element
// 1: Charm -1: No Charm
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

# Attempt to classify with GenPart
gInterpreter.Declare("""
      using Vfloat = const ROOT::RVec<float>&;
      using Vint = const ROOT::RVec<int>&;
      using namespace std;

// Function to print the
// index of an element
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
            v.push_back(idpart1+1); // Con la funcion typeC ya se que esta particula tiene como madre la misma que la anterior
            v.push_back(idpart2+1); // Same here
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

df = df.Define('typeWW','typeWW(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother)')
df = df.Define('typeC','typeC(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother)')
hadronic = df.Filter('typeWW == 1')
leptonic = df.Filter('typeWW == 3')
semi_charm = df.Filter('typeWW == 2 && typeC == 1')
semi_nocharm = df.Filter('typeWW == 2 && typeC != 1')

my_call = "partsIDs(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother)"
df = semi_charm.Define("partID",my_call)

df = df.Define("charmID",'fabs(GenPart_pdgId[partID[2]])<9 ? ((fabs(GenPart_pdgId[partID[2]])==4 && fabs(GenPart_pdgId[partID[2]+1])==3) ? partID[2] : partID[2]+1) : ((fabs(GenPart_pdgId[partID[3]])==4 && fabs(GenPart_pdgId[partID[3]+1])==3) ? partID[3] : partID[3]+1 )')
df = df.Define("strangeID",'fabs(GenPart_pdgId[partID[2]])<9 ? ((fabs(GenPart_pdgId[partID[2]])==3 && fabs(GenPart_pdgId[partID[2]+1])==4) ? partID[2] : partID[2]+1) : ((fabs(GenPart_pdgId[partID[3]])==3 && fabs(GenPart_pdgId[partID[3]+1])==4) ? partID[3] : partID[3]+1 )')

df = df.Define("AUX_charm_pdgID",'GenPart_pdgId[charmID]')
df = df.Define("AUX_strange_pdgID",'GenPart_pdgId[strangeID]')

df = df.Define("charm_pt",'GenPart_pt[charmID]')
df = df.Define("charm_eta",'GenPart_eta[charmID]')
df = df.Define("charm_phi",'GenPart_phi[charmID]')
df = df.Define("strange_pt",'GenPart_pt[strangeID]')
df = df.Define("strange_eta",'GenPart_eta[strangeID]')
df = df.Define("strange_phi",'GenPart_phi[strangeID]')

df = df.Define("CS_deltaR",'ROOT::VecOps::DeltaR(GenPart_eta[charmID],GenPart_eta[strangeID],GenPart_phi[charmID],GenPart_phi[strangeID])')
df = df.Define("CS_deltapt",'(GenPart_pt[charmID]-GenPart_pt[strangeID])')
df = df.Define("CS_deltaeta",'(GenPart_eta[charmID]-GenPart_eta[strangeID])')
df = df.Define("CS_deltaphi",'(GenPart_phi[charmID]-GenPart_phi[strangeID])')
df = df.Define("CS_invM",'ROOT::VecOps::InvariantMasses(GenPart_pt[charmID],GenPart_eta[charmID],GenPart_phi[charmID],GenPart_mass[charmID],GenPart_pt[strangeID],GenPart_eta[strangeID],GenPart_phi[strangeID],GenPart_mass[strangeID])')

df = df.Define("Genmuon_pt",'fabs(GenPart_pdgId[partID[2]])==13 ?  GenPart_pt[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==13 ?  GenPart_pt[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==13 ?  GenPart_pt[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==13 ?  GenPart_pt[partID[5]] :-20)))')
df = df.Define("Genelectron_pt",'fabs(GenPart_pdgId[partID[2]])==11 ?  GenPart_pt[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==11 ?  GenPart_pt[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==11 ?  GenPart_pt[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==11 ?  GenPart_pt[partID[5]] :-20)))')
df = df.Define("Gentau_pt",'fabs(GenPart_pdgId[partID[2]])==15 ?  GenPart_pt[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==15 ?  GenPart_pt[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==15 ?  GenPart_pt[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==15 ?  GenPart_pt[partID[5]] :-20)))')

df = df.Define("Genmuon_eta",'fabs(GenPart_pdgId[partID[2]])==13 ?  GenPart_eta[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==13 ?  GenPart_eta[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==13 ?  GenPart_eta[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==13 ?  GenPart_eta[partID[5]] :-20)))')
df = df.Define("Genelectron_eta",'fabs(GenPart_pdgId[partID[2]])==11 ?  GenPart_eta[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==11 ?  GenPart_eta[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==11 ?  GenPart_eta[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==11 ?  GenPart_eta[partID[5]] :-20)))')
df = df.Define("Gentau_eta",'fabs(GenPart_pdgId[partID[2]])==15 ?  GenPart_eta[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==15 ?  GenPart_eta[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==15 ?  GenPart_eta[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==15 ?  GenPart_eta[partID[5]] :-20)))')

df = df.Define("Genmuon_phi",'fabs(GenPart_pdgId[partID[2]])==13 ?  GenPart_phi[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==13 ?  GenPart_phi[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==13 ?  GenPart_phi[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==13 ?  GenPart_phi[partID[5]] :-20)))')
df = df.Define("Genelectron_phi",'fabs(GenPart_pdgId[partID[2]])==11 ?  GenPart_phi[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==11 ?  GenPart_phi[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==11 ?  GenPart_phi[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==11 ?  GenPart_eta[partID[5]] :-20)))')
df = df.Define("Gentau_phi",'fabs(GenPart_pdgId[partID[2]])==15 ?  GenPart_phi[partID[2]] : (fabs(GenPart_pdgId[partID[3]])==15 ?  GenPart_phi[partID[3]] : (fabs(GenPart_pdgId[partID[4]])==15 ?  GenPart_phi[partID[4]] : (fabs(GenPart_pdgId[partID[5]])==15 ?  GenPart_phi[partID[5]] :-20)))')

#entries1 = df.Count()
#print("%s entries passed semi filter" %entries1.GetValue())

testTAU = df.Filter('Gentau_pt>0')
testMUON = df.Filter('Genmuon_pt>0')
testEL = df.Filter('Genelectron_pt>0')

#entries1 = testTAU.Count()
#entries2 = testMUON.Count()
#entries3 = testEL.Count()
#print("%s entries are tau" %entries1.GetValue())
#print("%s entries are muon" %entries2.GetValue())
#print("%s entries are electron" %entries3.GetValue())

## We filter to maintain only events with a proper muon or electron

testMUON = testMUON.Filter('fabs(Genmuon_eta)<2.4 && Genmuon_pt>30')
testEL = testEL.Filter('fabs(Genelectron_eta)<2.4 && Genelectron_pt>30')

#print("%s is the quantity we are left with, one good muon or electron" %df_lepton.Count().GetValue())

#########################
######### HISTS
#########################

histNames = []

histNames.append("AUX_charm_pdgID")
histNames.append("AUX_strange_pdgID")
histNames.append("charm_pt")
histNames.append("charm_eta")
histNames.append("charm_phi")
histNames.append("strange_pt")
histNames.append("strange_eta")
histNames.append("strange_phi")
histNames.append("CS_deltaR")
histNames.append("CS_deltapt")
histNames.append("CS_deltaphi")
histNames.append("CS_deltaeta")
histNames.append("CS_invM")

hists_muon = {}
hists_el = {}

hists_muon["AUX_charm_pdgID"] = testMUON.Histo1D(("AUX_charm_pdgID","",10,-5,5),"AUX_charm_pdgID")
hists_muon["AUX_strange_pdgID"] = testMUON.Histo1D(("AUX_strange_pdgID","",10,-5,5),"AUX_strange_pdgID")
hists_muon["charm_pt"] = testMUON.Histo1D(("charm_pt","",100,0,150),"charm_pt")
hists_muon["charm_eta"] = testMUON.Histo1D(("charm_eta","",100,-4,4),"charm_eta")
hists_muon["charm_phi"] = testMUON.Histo1D(("charm_phi","",100,-1,7),"charm_phi")
hists_muon["strange_pt"] = testMUON.Histo1D(("strange_pt","",100,0,150),"strange_pt")
hists_muon["strange_eta"] = testMUON.Histo1D(("strange_eta","",100,-4,4),"strange_eta")
hists_muon["strange_phi"] = testMUON.Histo1D(("strange_phi","",100,-1,7),"strange_phi")
hists_muon["CS_deltaR"] = testMUON.Histo1D(("CS_deltaR","",100,-1,7),"CS_deltaR")
hists_muon["CS_deltapt"] = testMUON.Histo1D(("CS_deltapt","",100,-50,50),"CS_deltapt")
hists_muon["CS_deltaphi"] = testMUON.Histo1D(("CS_deltaphi","",100,-7,7),"CS_deltaphi")
hists_muon["CS_deltaeta"] = testMUON.Histo1D(("CS_deltaeta","",100,-5,5),"CS_deltaeta")
hists_muon["CS_invM"] = testMUON.Histo1D(("CS_invM","",100,0,300),"CS_invM")

hists_el["AUX_charm_pdgID"] = testEL.Histo1D(("AUX_charm_pdgID","",10,-5,5),"AUX_charm_pdgID")
hists_el["AUX_strange_pdgID"] = testEL.Histo1D(("AUX_strange_pdgID","",10,-5,5),"AUX_strange_pdgID")
hists_el["charm_pt"] = testEL.Histo1D(("charm_pt","",100,0,150),"charm_pt")
hists_el["charm_eta"] = testEL.Histo1D(("charm_eta","",100,-4,4),"charm_eta")
hists_el["charm_phi"] = testEL.Histo1D(("charm_phi","",100,-1,7),"charm_phi")
hists_el["strange_pt"] = testEL.Histo1D(("strange_pt","",100,0,150),"strange_pt")
hists_el["strange_eta"] = testEL.Histo1D(("strange_eta","",100,-4,4),"strange_eta")
hists_el["strange_phi"] = testEL.Histo1D(("strange_phi","",100,-1,7),"strange_phi")
hists_el["CS_deltaR"] = testEL.Histo1D(("CS_deltaR","",100,-1,7),"CS_deltaR")
hists_el["CS_deltapt"] = testEL.Histo1D(("CS_deltapt","",100,-50,50),"CS_deltapt")
hists_el["CS_deltaphi"] = testEL.Histo1D(("CS_deltaphi","",100,-7,7),"CS_deltaphi")
hists_el["CS_deltaeta"] = testEL.Histo1D(("CS_deltaeta","",100,-5,5),"CS_deltaeta")
hists_el["CS_invM"] = testEL.Histo1D(("CS_invM","",100,0,300),"CS_invM")

#####################################
############## PLOTS ################
#####################################

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


########### charms ##############3

#plot([hmet1,hmet2],["Charm","Not charm"],"/charm_pt_filtered.pdf",logY=True)

for name in histNames:
	plot([hists_muon[name],hists_el[name]],["Muon channel","Electron channel"],name+'.pdf')
