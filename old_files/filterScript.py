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
fileName = "/nfs/cms/vazqueze/SemiNew.root"
treeName = "Events"
#ROOT.fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, fileName)

####### FILTERS ########

### Selection of W events

entries1 = d.Count()
print("%s entries passed semi filter" %entries1.GetValue())

d_oneLepton = d.Filter('nMuonGood>0 || nElectronGood>0')

print("%s entries have one good muon or electron" %(d_oneLepton.Count()).GetValue())

print("%s entries have at least one Unfiltered muon or electron" %((d.Filter('nMuon>0 || nElectron>0')).Count()).GetValue())

d_jet = d_oneLepton.Filter('GoodJets_vector.size() >0')

print("%s entries have at least one associated jet" %(d_jet.Count()).GetValue())

print("%s entries have muon in jet and are charm" %d_jet.Filter('jet_muonid != -1 && typeC == 1').Count().GetValue())
print("%s entries have muon in jet and are not charm" %d_jet.Filter('jet_muonid != -1 && typeC == -1').Count().GetValue())

print("%s entries have OS and are charm" %d_jet.Filter('SSOS == 1 && typeC == 1').Count().GetValue())
print("%s entries have OS and are not charm" %d_jet.Filter('SSOS == 1 && typeC == -1').Count().GetValue())

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

### Some plots

#hmet1 = hadronic.Histo1D("nMuon")
#hmet2 = semi.Histo1D("nMuon")
#hmet3 = leptonic.Histo1D("nMuon")

#plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/nMuon.pdf")

############ Save #################


brlist = ["typeWW","typeC","nJet","nMuonGood","nJetGood","nFatJet","nGenJet","nGenPart","nMuon","FatJet_phi","FatJet_area","FatJet_eta","FatJet_mass","FatJet_msoftdrop",
"FatJet_phi","FatJet_pt","nSV","isMuonGood","isJetGood","FatJet_tau1","GenJet_eta","GenJet_mass","GenJet_phi","GenJet_pt","GenPart_eta",
"GenPart_mass","GenPart_phi","GenPart_pt","GenPart_genPartIdxMother","GenPart_pdgId","Muon_eta","Muon_mass","Muon_phi","Muon_pt",
"Muon_pfRelIso04_all","Muon_charge","Muon_jetIdx","Muon_pdgId","Jet_eta","Jet_mass","Jet_phi","Jet_pt","Jet_area","FirstJetPt","MET_phi","MET_pt",
"nElectron","Electron_phi","Electron_pt","Electron_eta","Electron_charge","Electron_mass","nElGood","isElGood","nFatJetGood","isFatJetGood",
"leadMuon","leadElectron","leadJet","chargeME","SSOS","FatJet_invM","TransverseMassEl","TransverseMass","ImassJet","Inv2massJetFiltered"]

#semi.Snapshot("Events", "SemiNew.root",brlist)

