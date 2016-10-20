#!/usr/bin/env python
from Alignment.OfflineValidation.TkAlAllInOneTool.helperFunctions import replaceByMap
from glob import glob
from helperstuff import config
from helperstuff import style
from helperstuff.enums import Analysis
from helperstuff.filemanager import tfiles
import os
import ROOT
import sys

class Folder(object):
    def __init__(self, folder, title, color, analysis, subdir, plotname):
        self.__folder, self.__title, self.color, self.analysis, self.subdir, self.plotname = folder, title, color, Analysis(analysis), subdir, plotname
        self.repmap = {
                       "analysis": str(self.analysis),
                       "gi": self.analysis.couplingname.replace("g", "").replace("1prime2", "1}^{#prime2"),
                       "intname": "CP" if self.analysis == "fa3" else "int",
                      }
    @property
    def folder(self):
        print config.plotsbasedir, "limits", self.subdir, replaceByMap(self.__folder, self.repmap)
        result = os.path.join(config.plotsbasedir, "limits", self.subdir, replaceByMap(self.__folder, self.repmap))
        gl = glob(result)
        if not gl:
            raise ValueError("{} does not exist!".format(result))
        if len(gl) > 1:
            raise ValueError("{} returns more than one match!".format(result))
        return gl[0]
    @property
    def title(self):
        return replaceByMap(self.__title, self.repmap)
    @property
    def graph(self):
        f = tfiles[os.path.join(self.folder, self.plotname)]
        c = f.c1
        mg = c.GetListOfPrimitives()[1]
        graphs = mg.GetListOfGraphs()
        assert len(graphs) == 1
        graphs[0].SetLineColor(self.color)
        return graphs[0]
    def addtolegend(self, legend):
        legend.AddEntry(self.graph, self.title, "l")

def mergeplots(analysis, subdir="", plotname="limit_nosystematics.root"):
    analysis = Analysis(analysis)
    folders = [
               Folder(".oO[analysis]Oo._discriminants_D_int_decay",    "D_{.oO[intname]Oo.}^{dec}",        1, analysis, subdir, plotname),
               Folder(".oO[analysis]Oo._discriminants_D_int_prod",     "D_{.oO[intname]Oo.}^{prod}",      16, analysis, subdir, plotname),
               Folder(".oO[analysis]Oo._discriminants_D_g11gi3",       "D_{g_{1}^{1}g_{.oO[gi]Oo.}^{3}}",  2, analysis, subdir, plotname),
               Folder(".oO[analysis]Oo._discriminants_D_g11gi3_prime", "D'_{g_{1}^{1}g_{.oO[gi]Oo.}^{3}}", 6, analysis, subdir, plotname),
               Folder(".oO[analysis]Oo._discriminants_D_g12gi2",       "D_{g_{1}^{2}g_{.oO[gi]Oo.}^{2}}",  4, analysis, subdir, plotname),
               Folder(".oO[analysis]Oo._discriminants_D_g12gi2_prime", "D'_{g_{1}^{2}g_{.oO[gi]Oo.}^{2}}", 7, analysis, subdir, plotname),
               Folder(".oO[analysis]Oo._discriminants_D_g13gi1",       "D_{g_{1}^{3}g_{.oO[gi]Oo.}^{1}}",  3, analysis, subdir, plotname),
               Folder(".oO[analysis]Oo._discriminants_D_g13gi1_prime", "D'_{g_{1}^{3}g_{.oO[gi]Oo.}^{1}}", ROOT.kGreen+3, analysis, subdir, plotname),
              ]
    outdir = "{}_discriminants".format(analysis)

    mg = ROOT.TMultiGraph("limit", "limit")
    if analysis == "fa3":
        l = ROOT.TLegend(.6, .2, .9, .5)
    else:
        l = ROOT.TLegend(.6, .6, .9, .9)
    l.SetBorderSize(0)
    l.SetFillStyle(0)
    l.SetNColumns(2)
    for folder in folders:
        mg.Add(folder.graph)
        folder.addtolegend(l)

    c = ROOT.TCanvas()
    mg.Draw("ac")
    l.Draw()
    saveasdir = os.path.join(config.plotsbasedir, "limits", subdir, outdir)
    try:
        os.makedirs(saveasdir)
    except OSError:
        pass
    for ext in "png eps root pdf".split():
        c.SaveAs(os.path.join(saveasdir, plotname.replace("root", ext)))

if __name__ == "__main__":
    args = []
    kwargs = {}
    for arg in sys.argv[1:]:
        if "=" in arg:
            kwargs[arg.split("=")[0]] = arg.split("=", 1)[1]
        else:
            args.append(arg)
    if args or kwargs:
        mergeplots(*args, **kwargs)
    else:
        print """
python mergeplots.py fa3 plotname=limit_VBF_nobkg_nosystematics.root
python mergeplots.py fa2 plotname=limit_VBF_nobkg_nosystematics.root
python mergeplots.py fa3 plotname=limit_WH,ZH_nobkg_nosystematics.root
python mergeplots.py fa2 plotname=limit_WH,ZH_nobkg_nosystematics.root
python mergeplots.py fa3 plotname=limit_VBF_nosystematics.root
python mergeplots.py fa2 plotname=limit_VBF_nosystematics.root
python mergeplots.py fa3 plotname=limit_WH,ZH_nosystematics.root
python mergeplots.py fa2 plotname=limit_WH,ZH_nosystematics.root
mkdir -p {basedir}/limits/discriminantcomparison/
for a in fa2 fa3; do
    for c in png eps root pdf; do
        for file in $(ls {basedir}/limits/${{a}}_discriminants/ | grep ".${{c}}$"); do
            ln -s {basedir}/limits/${{a}}_discriminants/$file {basedir}/limits/discriminantcomparison/${{a}}_${{file}}
        done
    done
done
""".format(basedir=config.plotsbasedir)
