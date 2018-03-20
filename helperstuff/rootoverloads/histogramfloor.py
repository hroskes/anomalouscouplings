from itertools import product
import ROOT

assert not hasattr(ROOT.TH1, "Floor")

def Floor(self, floorvalue=1e-18):
    if isinstance(self, ROOT.TH3):
        bins = product(range(self.GetNbinsX()+2), range(self.GetNbinsY()+2), range(self.GetNbinsZ()+2))
    elif isinstance(self, ROOT.TH2):
        bins = product(range(self.GetNbinsX()+2), range(self.GetNbinsY()+2))
    else:
        bins = product(range(self.GetNbinsX()+2))

    for bin in bins:
        bin = self.GetBin(*bin)
        self.SetBinContent(bin, max(self.GetBinContent(bin), floorvalue))

ROOT.TH1.Floor = Floor
