import ROOT

def seterrortozero(hsmooth, rawprojections):
    for x in xrange(1, hsmooth.GetNbinsX()+1):
      for y in xrange(1, hsmooth.GetNbinsY()+1):
        for z in xrange(1, hsmooth.GetNbinsZ()+1):
          hsmooth.SetBinError(x, y, z, 0)

    return False  #no need for new control plots because we didn't touch the bin content
