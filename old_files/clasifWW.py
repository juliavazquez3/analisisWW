import ROOT, os, sys
from ROOT import *

import numpy as np
import awkward as ak
import uproot
import pandas as pd
import pyarrow as pa
import urllib.request

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
fileName = "/nfs/cms/vazqueze/MyWW.root"
treeName = "Events"
#ROOT.fill_tree(treeName, fileName)


# We read the tree from the file and create a RDataFrame, a class that
# allows us to interact with the data contained in the tree.
d = ROOT.RDataFrame(treeName, fileName)

leptonic = d.Filter("typeWW==3")
hadronic = d.Filter("typeWW==1")
semi = d.Filter("typeWW==2")

semiCharm = semi.Filter("typeC==1")
semiNoCharm = semi.Filter("typeC==0")

def plot(hists,labels,filename,logY=False):

	c = ROOT.TCanvas('c','c',800,700)
	c.cd()
	colores = [ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kOrange]

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


	if logY:
                c.SetLogy()

	legend.Draw()

	plot_filename = plotdir+filename

	c.Print(plot_filename)


## Pintamos varias comparaciones

hmet1 = hadronic.Histo1D(("nJet","nJet",10,0,10),"nJetGood")
hmet2 = semi.Histo1D(("nJet","nJet",10,0,10),"nJetGood")
hmet3 = leptonic.Histo1D(("nJet","nJet",10,0,10),"nJetGood")

plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/clas_nJet.pdf")

hmet1 = hadronic.Histo1D(("nMuons","nMuons",10,0,10),"nMuonGood")
hmet2 = semi.Histo1D(("nMuons","nMuons",10,0,10),"nMuonGood")
hmet3 = leptonic.Histo1D(("nMuons","nMuons",10,0,10),"nMuonGood")

plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/clas_nMuon.pdf")

##Define new var

hadronic1 = hadronic.Define('jet_pt1',"nJet > 0 ? Jet_pt[0] : -1")
semi1 = semi.Define('jet_pt1',"nJet > 0 ? Jet_pt[0] : -1")
leptonic1 = leptonic.Define('jet_pt1',"nJet > 0 ? Jet_pt[0] : -1")

hadronic2 = hadronic1.Define('Muon_pt1',"nMuon > 0 ? Muon_pt[0] : -1")
semi2 = semi1.Define('Muon_pt1',"nMuon > 0 ? Muon_pt[0] : -1")
leptonic2 = leptonic1.Define('Muon_pt1',"nMuon > 0 ? Muon_pt[0] : -1")

#hmet1 = hadronic1.Histo1D(("Jet_pt1","Jet_pt1",100,0,100),"FirstJetPt")
#hmet2 = semi1.Histo1D(("Jet_pt1","Jet_pt1",100,0,100),"FirstJetPt")
#hmet3 = leptonic1.Histo1D(("Jet_pt1","Jet_pt1",100,0,100),"FirstJetPt")

hmet1 = hadronic1.Histo1D("jet_pt1")
hmet2 = semi1.Histo1D("jet_pt1")
hmet3 = leptonic1.Histo1D("jet_pt1")


plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/clas_JetPt.pdf",logY=True)

hmet1 = hadronic2.Histo1D(("Jet_muon1","Jet_muon1",100,0,200),"Muon_pt1")
hmet2 = semi2.Histo1D(("Jet_muon1","Jet_muon1",100,0,200),"Muon_pt1")
hmet3 = leptonic2.Histo1D(("Jet_muon1","Jet_muon1",100,0,200),"Muon_pt1")

plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/clas_MuonPt.pdf")

hmet1 = hadronic1.Histo1D(("Jet_pt1","Jet_pt1",1000,0,2000),"jet_pt1")
hmet2 = semi1.Histo1D(("Jet_pt1","Jet_pt1",1000,0,2000),"jet_pt1")
hmet3 = leptonic1.Histo1D(("Jet_pt1","Jet_pt1",1000,0,2000),"jet_pt1")


plot([hmet1,hmet2,hmet3],["Hadronic","Semi","Leptonic"],"/clas_JetPt_limited.pdf",logY=True)

