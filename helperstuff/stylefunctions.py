import ROOT

from ROOT import gStyle, gPad, gROOT, kWhite

def fixOverlay():
  gPad.RedrawAxis()


def setTDRStyle(force):

  gStyle.SetCanvasBorderMode(0)
  gStyle.SetCanvasColor(kWhite)
  gStyle.SetCanvasDefH(600)
  gStyle.SetCanvasDefW(600)
  gStyle.SetCanvasDefX(0)
  gStyle.SetCanvasDefY(0)

  gStyle.SetPadBorderMode(0)
  gStyle.SetPadColor(kWhite)
  gStyle.SetPadGridX(False)
  gStyle.SetPadGridY(False)
  gStyle.SetGridColor(0)
  gStyle.SetGridStyle(3)
  gStyle.SetGridWidth(1)

  gStyle.SetFrameBorderMode(0)
  gStyle.SetFrameBorderSize(1)
  gStyle.SetFrameFillColor(0)
  gStyle.SetFrameFillStyle(0)
  gStyle.SetFrameLineColor(1)
  gStyle.SetFrameLineStyle(1)
  gStyle.SetFrameLineWidth(1)

  if force:
      gStyle.SetHistLineColor(1)
      gStyle.SetHistLineStyle(0)
      gStyle.SetHistLineWidth(1)


  gStyle.SetEndErrorSize(2)
  gStyle.SetErrorX(0.)

  gStyle.SetMarkerStyle(20)

  gStyle.SetOptFit(1)
  gStyle.SetFitFormat("5.4g")
  gStyle.SetFuncColor(2)
  gStyle.SetFuncStyle(1)
  gStyle.SetFuncWidth(1)

  gStyle.SetOptDate(0)

  gStyle.SetOptFile(0)
  gStyle.SetOptStat(0)
  gStyle.SetStatColor(kWhite)
  gStyle.SetStatFont(42)
  gStyle.SetStatFontSize(0.04)
  gStyle.SetStatTextColor(1)
  gStyle.SetStatFormat("6.4g")
  gStyle.SetStatBorderSize(1)
  gStyle.SetStatH(0.1)
  gStyle.SetStatW(0.2)


  gStyle.SetPadTopMargin(0.05)
  gStyle.SetPadBottomMargin(0.13)
  gStyle.SetPadLeftMargin(0.16)
  gStyle.SetPadRightMargin(0.04)


  gStyle.SetOptTitle(0)
  gStyle.SetTitleFont(62, "t")
  gStyle.SetTitleColor(1)
  gStyle.SetTitleTextColor(1)
  gStyle.SetTitleFillColor(10)
  gStyle.SetTitleFontSize(0.05)


  gStyle.SetTitleColor(1, "XYZ")
  gStyle.SetTitleFont(42, "XYZ")
  gStyle.SetTitleSize(0.06, "XYZ")
  gStyle.SetTitleXOffset(0.9)
  gStyle.SetTitleYOffset(1.25)


  gStyle.SetLabelColor(1, "XYZ")
  gStyle.SetLabelFont(42, "XYZ")
  gStyle.SetLabelOffset(0.007, "XYZ")
  gStyle.SetLabelSize(0.05, "XYZ")


  gStyle.SetAxisColor(1, "XYZ")
  gStyle.SetStripDecimals(True)
  gStyle.SetTickLength(0.03, "XYZ")
  gStyle.SetNdivisions(510, "XYZ")
  gStyle.SetPadTickX(1)
  gStyle.SetPadTickY(1)

  gStyle.SetOptLogx(0)
  gStyle.SetOptLogy(0)
  gStyle.SetOptLogz(0)

  gStyle.SetPaperSize(20.,20.)


  gROOT.ForceStyle()



def tdrstyle(force=True):
    setTDRStyle(force)


def asymmerrorsfromhistogram(h, showemptyerrors=False):
    from array import array
    quant = (1.0-0.6827)/2.0
    nbins = h.GetNbinsX()
    x, exu, exd, y, eyu, eyd = array('d', [0]*nbins), array('d', [0]*nbins), array('d', [0]*nbins), array('d', [0]*nbins), array('d', [0]*nbins), array('d', [0]*nbins)
    n = 0
    for i in range(nbins):
        bin = i+1
        x[n] = h.GetBinCenter(bin)
        y[n] = bincontent = h.GetBinContent(bin)
        #print bin, n, x[n], y[n]
        exu[n] = exd[n] = 0
        eyu[n] = ROOT.Math.chisquared_quantile_c(quant,2*(bincontent+1))/2.-bincontent
        eyd[n] = bincontent-ROOT.Math.chisquared_quantile_c(1-quant,2*bincontent)/2. if bincontent else 0
        if bincontent or showemptyerrors:
            n += 1
    return ROOT.TGraphAsymmErrors(n, x, y, exd, exu, eyd, eyu)

def ymax(*args):
    bkpgpad = ROOT.gPad
    result = float("-inf")
    for obj, drawoption in args:
        c = ROOT.TCanvas()
        obj.Draw(drawoption)
        result = max(result, obj.GetYaxis().GetXmax())
    if bkpgpad:
        bkpgpad.cd()
    return result

def applycanvasstyle(c):
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetBorderSize(2)
    c.SetTickx(1)
    c.SetTicky(1)
    c.SetLeftMargin(0.17)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.07)
    c.SetBottomMargin(0.13)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)

def applylegendstyle(l):
    l2D.SetBorderSize(0)
    l2D.SetTextFont(42)
    l2D.SetTextSize(0.04)
    l2D.SetLineColor(1)
    l2D.SetLineStyle(1)
    l2D.SetLineWidth(1)
    l2D.SetFillColor(0)
    l2D.SetFillStyle(0)

def applyaxesstyle(h):
    h.GetXaxis().SetNdivisions(505)
    h.GetXaxis().SetLabelFont(42)
    h.GetXaxis().SetLabelOffset(0.007)
    h.GetXaxis().SetLabelSize(0.04)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(0.9)
    h.GetXaxis().SetTitleFont(42)
    h.GetYaxis().SetNdivisions(505)
    h.GetYaxis().SetLabelFont(42)
    h.GetYaxis().SetLabelOffset(0.007)
    h.GetYaxis().SetLabelSize(0.04)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(1.1)
    h.GetYaxis().SetTitleFont(42)

tokeep = {}
def cuttext(text):
    if (cuttext, text) not in tokeep:
        pt = tokeep[cuttext,text] = TPaveText(0.71,0.84,0.92,0.92,"brNDC")
        pt.SetBorderSize(0)
        pt.SetTextAlign(12)
        pt.SetTextSize(0.04)
        pt.SetFillStyle(0)
        pt.SetTextFont(42)
        pt.AddText(0.01,0.01,text)
    tokeep[cuttext,text].Draw()

def CMS():
    if CMS not in tokeep:
         pass
