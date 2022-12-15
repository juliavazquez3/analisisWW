import ROOT
import sys
import os

plotdir = 'plots/' # this is where we'll save your plots
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

if len( sys.argv ) != 3:
	print ("USAGE : %s <input file > <output file >"%(sys.argv[0]))
	sys.exit(1)

inFileName = sys.argv[1]
outFileName = sys.argv[2]
print (" Reading from ", inFileName , "and writing to", outFileName)

inFile = ROOT.TFile.Open(inFileName,"READ")
tree = inFile.Get("Events")

mll = ROOT.TH1D("data","m_{ll} , data ",150 ,10.,150.)
mll.Sumw2()

for entryNum in range (0,tree.GetEntries()):
	tree.GetEntry(entryNum)
	FatJet_invM = getattr(tree,"FatJet_invM")
	mll.Fill(FatJet_invM)

mll.SetDirectory(0)

dataHisto = mll

dataHisto.SetDirectory(0)

canvas = ROOT.TCanvas("canvas")
canvas.cd()
canvas.SetLogy(True)

plotFileName = plotdir+'/fatjet_pt_fit.pdf'

dataHisto.SetStats(0)
dataHisto.SetLineColor(ROOT.kBlack)
dataHisto.SetLineWidth(2)
dataHisto.GetYaxis().SetTitle("Number of events")
dataHisto.GetXaxis().SetTitle("")
dataHisto.Draw("pe")
canvas.Print(plotFileName)

fit = ROOT.TF1("myfit","gaus(0)+expo(3)",10.,150.)
gausfit = ROOT.TF1("myfit","gaus(0)",10.,150.)
fit.SetParameters(200,90,10,8,-0.04)
#gausfit.SetParameters(0.5,80,20)

dataHisto.Fit(fit,"L")

dataHisto.Draw("pe")

fit = ROOT.TF1("myfit","gaus(0)+expo(3)",10.,150.)
gausfit = ROOT.TF1("myfit","gaus(0)",10.,150.)
fit.SetParameters(200,90,10,8,-0.04)

#fit.Draw("Sames")
canvas.Print(plotFileName)

