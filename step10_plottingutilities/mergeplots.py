from helperstuff import config
from helperstuff.filemanager import tfiles
from helperstuff import style
import ROOT
import os

class Folder(object):
    def __init__(self, folder, title, color):
        self.folder, self.title, self.color = folder, title, color
    @property
    def graph(self):
        f = tfiles[os.path.join(config.plotsbasedir, "limits", self.folder, plotname)]
        c = f.c1
        mg = c.GetListOfPrimitives()[1]
        graphs = mg.GetListOfGraphs()
        assert len(graphs) == 1
        graphs[0].SetLineColor(self.color)
        return graphs[0]
    def addtolegend(self, legend):
        legend.AddEntry(self.graph, self.title, "l")

plotname = "limit_nosystematics.root"
folders = [
           Folder("fa3_discriminants_D_int_decay_noggH",    "D_{CP}^{dec}",           1),
           Folder("fa3_discriminants_D_int_VBF_noggH",      "D_{CP}^{VBF}",           16),
           Folder("fa3_discriminants_D_g11gi3_noggH",       "D_{g_{1}^{1}g_{4}^{3}}",  2),
           Folder("fa3_discriminants_D_g11gi3_prime_noggH", "D'_{g_{1}^{1}g_{4}^{3}}", 6),
           Folder("fa3_discriminants_D_g12gi2_noggH",       "D_{g_{1}^{2}g_{4}^{2}}",  4),
           Folder("fa3_discriminants_D_g12gi2_prime_noggH", "D'_{g_{1}^{2}g_{4}^{2}}", 7),
           Folder("fa3_discriminants_D_g13gi1_noggH",       "D_{g_{1}^{3}g_{4}^{1}}",  3),
           Folder("fa3_discriminants_D_g13gi1_prime_noggH", "D'_{g_{1}^{3}g_{4}^{1}}", ROOT.kGreen+3),
          ]
outdir = "fa3_discriminants"

mg = ROOT.TMultiGraph("limit", "limit")
#l = ROOT.TLegend(.6, .6, .9, .9)
l = ROOT.TLegend(.6, .2, .9, .5)
l.SetBorderSize(0)
l.SetFillStyle(0)
l.SetNColumns(2)
for folder in folders:
    mg.Add(folder.graph)
    folder.addtolegend(l)

c = ROOT.TCanvas()
mg.Draw("ac")
l.Draw()
saveasdir = os.path.join(config.plotsbasedir, "limits", outdir)
try:
    os.makedirs(saveasdir)
except OSError:
    pass
for ext in "png eps root pdf".split():
    c.SaveAs(os.path.join(saveasdir, plotname.replace("root", ext)))
