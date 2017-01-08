from CJLSTscripts import CJLSTscriptsfolder
import ROOT

def fakeRate13TeV(*args, **kwargs):
    raise ValueError("Have to call setup before calling fakeRate13TeV")

def setup(production):
    global fakeRate13TeV
    success = ROOT.ZXsetup(int(production), CJLSTscriptsfolder)
    if not success:
        raise ValueError("Bad production: {}".format(production))
    from ROOT import fakeRate13TeV

from ROOT import CRZLLss, test_bit
