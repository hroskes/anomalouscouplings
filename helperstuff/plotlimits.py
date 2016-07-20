import array
from combinehelpers import Luminosity
from extendedcounter import ExtendedCounter
import ROOT
import stylefunctions as style
import sys

def plotlimits(filename, branchname, outputfilename, xaxislabel):
    f = ROOT.TFile(filename)
    assert f
    t = f.Get("limit")
    assert t

    NLL = ExtendedCounter() 

    for i in range(1, t.GetEntries()):
        t.GetEntry(i)
        fa3 = getattr(t, branchname)
        deltaNLL = t.deltaNLL
        print fa3
        NLL[fa3] = 2*deltaNLL
    if 1 not in NLL and -1 in NLL:
        NLL[1] = NLL[-1]

    c1 = ROOT.TCanvas("c1", "", 8, 30, 800, 800)
    g = NLL.TGraph()
    g.Draw("AC")
    g.GetXaxis().SetTitle(xaxislabel)
    g.GetXaxis().SetRangeUser(-1, 1)
    g.GetYaxis().SetTitle("-2#Deltaln L")

    style.applycanvasstyle(c1)
    style.applyaxesstyle(g)
    style.CMS("Preliminary", float(Luminosity("forexpectedscan")))

    drawlines()
    for ext in "png eps root pdf".split():
        outputfilename = outputfilename.replace("."+ext, "")
    for ext in "png eps root pdf".split():
        c1.SaveAs("{}.{}".format(outputfilename, ext))

def drawlines():
    line68 = ROOT.TLine()
    line68.SetLineStyle(9)
    line68.DrawLine(-1,1,1,1)
    line95 = ROOT.TLine()
    line95.SetLineStyle(9)
    line95.SetLineColor(4)
    line95.DrawLine(-1,3.84,1,3.84)

    oneSig = ROOT.TPaveText(0.18,0.16,0.28,0.19,"NDC")
    oneSig.SetFillColor(0)
    oneSig.SetFillStyle(0)
    oneSig.SetTextFont(42)
    oneSig.SetBorderSize(0)
    oneSig.AddText("68\% CL")
    oneSig.Draw()

    twoSig = ROOT.TPaveText(0.18,0.24,0.28,0.27,"NDC")
    twoSig.SetFillColor(0)
    twoSig.SetFillStyle(0)
    twoSig.SetTextFont(42)
    twoSig.SetBorderSize(0)
    twoSig.AddText("95\% CL")
    twoSig.Draw()
