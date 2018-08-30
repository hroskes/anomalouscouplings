#!/usr/bin/env python
"""
combine ScaleUp, ResUp --> ScaleResUp
ZX systematics templates from qqZZ
"""
import argparse
if __name__ == "__main__":
    def __Analysis(*args, **kwargs):
        from helperstuff.enums import Analysis
        return Analysis(*args, **kwargs)
    parser = argparse.ArgumentParser()
    parser.add_argument("analysis", type=__Analysis)
    args = parser.parse_args()

from itertools import product
import os

import ROOT

import helperstuff.rootoverloads.histogramfloor

from helperstuff import config

from helperstuff.combinehelpers import gettemplate
from helperstuff.enums import analyses, categories, channels, productions
from helperstuff.templates import IntTemplate, Template, TemplatesFile
from helperstuff.utilities import KeepWhileOpenFiles, tfiles

dom4lshapes = any((config.applym4lshapesystematicsUntagged, config.applym4lshapesystematicsVBFVHtagged, config.applym4lshapesystematicsggH, config.applym4lshapesystematicsggHUntagged, config.applym4lshapesystematicsdiagonal))

def combinesystematics(channel, analysis, production, category):
    thetemplatesfiles = []

    if dom4lshapes:
        ScaleAndRes = {}
        for _ in "ggh", "vbf", "zh", "wh", "tth", "bbh":
            for syst in "ScaleUp", "ScaleDown", "ResUp", "ResDown":
                if config.getm4lsystsfromggHUntagged:
                    if _ == "ggh" and category == "Untagged": continue
                elif config.getm4lsystsfromggH:
                    if _ == "ggh": continue
                else:
                    continue
                ScaleAndRes[_,syst] = TemplatesFile(channel, _, syst, analysis, production, category).actualtemplatesfile
            assert all("Scale" in _.templatesfile() or "Res" in _.templatesfile() for _ in ScaleAndRes.values())
        thetemplatesfiles += ScaleAndRes.values()
        if config.combinem4lshapesystematics:
            ScaleResUp, ScaleResDown = {}, {}
            for _ in "ggh", "vbf", "zh", "wh", "tth", "bbh":
                ScaleResUp[_] = TemplatesFile(channel, _, "ScaleResUp", analysis, production, category).actualtemplatesfile
                ScaleResDown[_] = TemplatesFile(channel, _, "ScaleResDown", analysis, production, category).actualtemplatesfile
            assert all("ScaleRes" in _.templatesfile() for _ in ScaleResUp.values()+ScaleResDown.values())
            thetemplatesfiles += ScaleResUp.values() + ScaleResDown.values()

    if config.applyZXshapesystematicsUntagged and category == "Untagged" or config.applyZXshapesystematicsVBFVHtagged and category != "Untagged":
        ZXUp = TemplatesFile(channel, "bkg", "ZXUp", analysis, production, category).actualtemplatesfile
        ZXDown = TemplatesFile(channel, "bkg", "ZXDown", analysis, production, category).actualtemplatesfile
        assert "ZXUp" in ZXUp.templatesfile() and "ZXDown" in ZXDown.templatesfile()
        thetemplatesfiles += [ZXUp, ZXDown]

    if config.applyMINLOsystematics and category != "Untagged":
        MINLOUp = TemplatesFile(channel, "ggh", "MINLOUp", analysis, production, category).actualtemplatesfile
        MINLODn = TemplatesFile(channel, "ggh", "MINLODn", analysis, production, category).actualtemplatesfile
        assert "MINLOUp" in MINLOUp.templatesfile() and "MINLODn" in MINLODn.templatesfile()
        thetemplatesfiles += [MINLOUp, MINLODn]

    assert not any(tf.copyfromothertemplatesfile for tf in thetemplatesfiles)

    with KeepWhileOpenFiles(*[_.templatesfile()+".tmp" for _ in thetemplatesfiles]) as kwofs:
      outfiles = {
        templatesfile: 
          ROOT.TFile(templatesfile.templatesfile(), "CREATE")
            if kwofs[i]
            and not os.path.exists(templatesfile.templatesfile())
          else
            None
        for i, templatesfile in enumerate(thetemplatesfiles)
      }

      store = []

      if dom4lshapes and config.getm4lsystsfromggH:
          if config.getm4lsystsfromggHUntagged: usecategory = "Untagged"
          elif config.getm4lsystsfromggH: usecategory = category
          else: assert False

          ggHSM = {systematic:
                     gettemplate(channel, "ggH", analysis, production, usecategory, analysis.purehypotheses[0], systematic)
                     .ProjectionZ().Clone("projection_{}".format(systematic))
                   for systematic in ("", "ScaleUp", "ScaleDown", "ResUp", "ResDown")}
          ggHnominal = ggHSM[""]
          for syst in "ScaleUp", "ScaleDown", "ResUp", "ResDown":
              ggHsyst = ggHSM[syst]
              for _ in "ggh", "vbf", "zh", "wh", "tth", "bbh":
                  if _ == "ggh" and category == "Untagged" and config.getm4lsystsfromggHUntagged: continue
                  elif _ == "ggh" and not config.getm4lsystsfromggHUntagged: continue
                  outfile = outfiles[ScaleAndRes[_,syst]]
                  if not outfile: continue
                  tf = TemplatesFile(channel, _, analysis, production, category)
                  for t in tf.templates() + tf.inttemplates():
                      h = t.gettemplate().Clone()
                      integral = h.Integral()
                      h.SetDirectory(outfile)
                      for x, y, z in product(xrange(1, h.GetNbinsX()+1), xrange(1, h.GetNbinsY()+1), xrange(1, h.GetNbinsZ()+1)):
                          h.SetBinContent(
                                          x, y, z,
                                            h.GetBinContent(x, y, z)
                                            * ggHsyst.GetBinContent(z)
                                            / ggHnominal.GetBinContent(z)
                                         )
                      store.append(h)

          for templatesfile in ScaleAndRes.values():
              if outfiles[templatesfile]:
                  outfiles[templatesfile].Write()

      if dom4lshapes and config.combinem4lshapesystematics:
          for _ in "ggh", "vbf", "zh", "wh", "tth", "bbh":
              tf = TemplatesFile(channel, _, analysis, production, category)
              for t in tf.templates() + tf.inttemplates():
                  outfile = outfiles[ScaleResDown[_]], outfiles[ScaleResUp[_]]
                  if not any(outfile): continue
                  assert all(outfile)
                  if isinstance(t, Template): hypothesis = t.hypothesis
                  elif isinstance(t, IntTemplate): hypothesis = t.interferencetype
                  h = gettemplate(t.productionmode, hypothesis, channel, analysis, production, category)
                  #ScaleResDown = nominal + (ScaleDown-nominal) + (ResDown-nominal)
                  #           = -nominal + ScaleDown + ResDown
                  hdn = h.Clone()
                  hdn.SetDirectory(outfiles[ScaleResDown[_]])
                  hdn.Scale(-1)
                  hdn.Add(gettemplate(t.productionmode, hypothesis, channel, analysis, "ResDown", production, category))
                  hdn.Add(gettemplate(t.productionmode, hypothesis, channel, analysis, "ScaleDown", production, category))

                  #ScaleResUp = nominal - (ScaleResDown-nominal)
                  #             = 2*nominal - ScaleResDown
                  hup = h.Clone()
                  hup.SetDirectory(outfiles[ScaleResUp[_]])
                  hup.Scale(2)
                  hup.Add(hdn, -1)

                  store += [hup, hdn]

      if config.applyMINLOsystematics and category != "Untagged" and (outfiles[MINLOUp] or outfiles[MINLODn]):
          assert outfiles[MINLOUp] and outfiles[MINLODn]
          tf = TemplatesFile(channel, "ggh", analysis, production, category)
          POWHEG3DSM = gettemplate("ggH", "0+", channel, analysis, production, category)
          MINLO3DSM = gettemplate("ggH", "0+", "MINLO_SM", "MINLO", channel, analysis, production, category)

          for t in tf.templates():
              POWHEG3D = t.gettemplate()
              MINLO3D = POWHEG3D.Clone(Template("ggH", t.hypothesis, "MINLOUp", channel, analysis, production, category).templatename())
              MINLO3D.SetDirectory(outfiles[MINLOUp])
              MINLO3DDn = POWHEG3D.Clone(Template("ggH", t.hypothesis, "MINLODn", channel, analysis, production, category).templatename())
              MINLO3DDn.SetDirectory(outfiles[MINLODn])
              for x, y, z in product(xrange(1, MINLO3D.GetNbinsX()+1), xrange(1, MINLO3D.GetNbinsY()+1), xrange(1, MINLO3D.GetNbinsZ()+1)):
                  try:
                      upbincontent = POWHEG3D.GetBinContent(x, y, z) * MINLO3DSM.GetBinContent(x, y, z) / POWHEG3DSM.GetBinContent(x, y, z)
                  except ZeroDivisionError:
                      if MINLO3DSM.GetBinContent(x, y, z) == POWHEG3D.GetBinContent(x, y, z) == 0:
                          upbincontent = 0
                      else:
                           raise
                  downbincontent = POWHEG3D.GetBinContent(x, y, z) ** 2 / upbincontent
                  MINLO3D.SetBinContent(
                                        x, y, z,
                                        upbincontent
                                       )
                  MINLO3DDn.SetBinContent(
                                          x, y, z,
                                          downbincontent
                                         )
              if production.year == 2016:
                  MINLO3D.Scale(POWHEG3D.Integral()/MINLO3D.Integral())
                  MINLO3DDn.Scale(POWHEG3D.Integral()/MINLO3DDn.Integral())
              elif production.year == 2017:
                  pass
              else:
                  assert False
              MINLO3DDn.Floor()
              store += [MINLO3D, MINLO3DDn]
          for templatesfile in MINLOUp, MINLODn: outfiles[templatesfile].Write()

          for t in tf.inttemplates():
              MINLOintUp = POWHEG3D.Clone(IntTemplate("ggH", t.interferencetype, "MINLOUp", channel, analysis, production, category).templatename())
              MINLOintUp.SetDirectory(outfiles[MINLOUp])
              MINLOintUp.Reset("M")
              for template, factor in t.templatesandfactors:
                  MINLOintUp.Add(getattr(outfiles[MINLOUp], template.templatename()), factor)

              MINLOintDn = POWHEG3D.Clone(IntTemplate("ggH", t.interferencetype, "MINLODn", channel, analysis, production, category).templatename())
              MINLOintDn.SetDirectory(outfiles[MINLODn])
              MINLOintDn.Reset("M")
              for template, factor in t.templatesandfactors:
                  MINLOintDn.Add(getattr(outfiles[MINLODn], template.templatename()), factor)

              store += [MINLOintUp, MINLOintDn]

      if (config.applyZXshapesystematicsUntagged and category == "Untagged" or config.applyZXshapesystematicsVBFVHtagged and category != "Untagged") and (outfiles[ZXUp] or outfiles[ZXDown]):
          assert outfiles[ZXUp] and outfiles[ZXDown], (outfiles[ZXUp], outfiles[ZXDown])
          ZXtemplate = gettemplate("ZX", channel, analysis, production, category)
          qqZZtemplate = gettemplate("qqZZ", channel, analysis, production, category)
          ZXUptemplate = qqZZtemplate.Clone(Template("ZX", channel, analysis, "ZXUp", production, category).templatename())
          ZXUptemplate.Scale(ZXtemplate.Integral() / ZXUptemplate.Integral())
          ZXUptemplate.SetDirectory(outfiles[ZXUp])

          ZXDowntemplate = ZXtemplate.Clone(Template("ZX", channel, analysis, "ZXDown", production, category).templatename())
          ZXDowntemplate.SetDirectory(outfiles[ZXDown])
          ZXDowntemplate.Scale(2)
          ZXDowntemplate.Add(ZXUptemplate, -1)
          ZXDowntemplate.Floor()
          ZXDowntemplate.Scale(ZXtemplate.Integral() / ZXDowntemplate.Integral())

      for f in outfiles.values():
          if not f: continue
          f.Write()
          f.Close()

if __name__ == "__main__":
    for production in productions:
        for analysis in analyses:
            if analysis != args.analysis: continue
            for channel in channels:
                for category in categories:
                    print production, analysis, channel, category
                    combinesystematics(channel, analysis, production, category)
