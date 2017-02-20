import ROOT

def justcopy(hsmooth, rawprojections, othertemplate):
    integral = hsmooth.Integral()
    hsmooth.Reset("M")
    if othertemplate.hascustomsmoothing: raise ValueError("Can't copy a template that has custom smoothing!")
    f = ROOT.TFile(othertemplate.templatesfile.templatesfile(firststep=True))
    otherh = getattr(f, othertemplate.templatename())
    hsmooth.Add(otherh)
    hsmooth.Scale(integral / hsmooth.Integral())

    return True #do make new control plots
