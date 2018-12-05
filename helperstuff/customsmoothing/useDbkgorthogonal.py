from itertools import product
import ROOT

def Integral(hist3d, binx1, binx2, biny1, biny2, binz1, binz2, option=""):
    #wrapper to allow kwargs
    return hist3d.Integral(binx1, binx2, biny1, biny2, binz1, binz2, option)
def SetBinContent(hist3d, binx, biny, binz, content):
    #wrapper to allow kwargs
    return hist3d.SetBinContent(binx, biny, binz, content)

def useDbkgorthogonal(hsmooth, rawprojections):
    """
    Loop over the xy plane.
    For each xy bin that has a z bin with zero content,
    keep the integral over z, but copy the shape from the overall z shape
    """
    totalintegral = hsmooth.Integral()
    if not totalintegral:
        assert hsmooth.GetMaximum() == hsmooth.GetMinimum() == 0
        return True
    zprojection = rawprojections[2]
    for xbin, ybin in product(xrange(1, hsmooth.GetNbinsX()+1), xrange(1, hsmooth.GetNbinsY()+1)):
        for zbin in xrange(1, hsmooth.GetNbinsZ()+1):
            if hsmooth.GetBinContent(xbin, ybin, zbin) == 0:
                break  #process this slice
        else:
            continue   #don't process this slice

        sliceintegral = Integral(hsmooth, xbin, xbin, ybin, ybin, 1, hsmooth.GetNbinsZ())
        for zbin in xrange(1, hsmooth.GetNbinsZ()+1):
            bincontent = zprojection.GetBinContent(zbin) * sliceintegral / totalintegral
            hsmooth.SetBinContent(xbin, ybin, zbin, bincontent)
            hsmooth.SetBinError(xbin, ybin, zbin, bincontent)

    return True #do make new control plots
