from CJLSTscripts import CJLSTscriptsfolder, convertTGraphstoTH1Fs
import os
import ROOT

def fakeRate13TeV(*args, **kwargs):
    raise ValueError("Have to call setup before calling fakeRate13TeV")

def setup(production):
    global fakeRate13TeV
    status = ROOT.ZXsetup(int(production), CJLSTscriptsfolder)
    if status == ROOT.ZXsetupsuccess:
        pass
    elif status == ROOT.ZXsetupbadproduction:
        raise ValueError("Bad production: {}".format(production))
    elif status == ROOT.ZXsetupfailed:
        raise ValueError("ZX setup failed! {}".format(production))
    else:
        raise ValueError("???????")
    from ROOT import fakeRate13TeV

for rootfile in "FakeRate_SS_Moriond368.root", "FakeRates_SS_Moriond18.root":
    convertTGraphstoTH1Fs.convertTGraphstoTH1Fs(os.path.join(CJLSTscriptsfolder, rootfile))

from ROOT import CRZLLss, test_bit
