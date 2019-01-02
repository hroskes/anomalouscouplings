import ROOT

def seterrorforfloor(hsmooth, rawprojections):
    mincontent = min(
      hsmooth.GetBinContent(x, y, z)
        for x in xrange(1, hsmooth.GetNbinsX()+1)
        for y in xrange(1, hsmooth.GetNbinsY()+1)
        for z in xrange(1, hsmooth.GetNbinsZ()+1)
      if hsmooth.GetBinContent(x, y, z) != 0
    )
    print mincontent

    maxerrorratio, errortoset = max(
      (hsmooth.GetBinError(x, y, z) / hsmooth.GetBinContent(x, y, z), hsmooth.GetBinError(x, y, z))
        for x in xrange(1, hsmooth.GetNbinsX()+1)
        for y in xrange(1, hsmooth.GetNbinsY()+1)
        for z in xrange(1, hsmooth.GetNbinsZ()+1)
      if hsmooth.GetBinContent(x, y, z) != 0
    )
    #the reasoning being that if there's a bin with just one entry 2.3 +/- 2.3, then the zero bin could also have 2.3
    #but we can't draw that conclusion from a bin 1000 +/- 5.5

    for x in xrange(1, hsmooth.GetNbinsX()+1):
      for y in xrange(1, hsmooth.GetNbinsY()+1):
        for z in xrange(1, hsmooth.GetNbinsZ()+1):
          if hsmooth.GetBinError(x, y, z) == 0 or hsmooth.GetBinContent(x, y, z) <= mincontent*1.1:  #buffer for rounding error
            assert hsmooth.GetBinContent(x, y, z) <= mincontent*1.1
            hsmooth.SetBinError(x, y, z, errortoset)
            

    return False  #no need for new control plots because we didn't touch the bin content
