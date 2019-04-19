import ROOT

def draw(w, pdfname, discname, fai, saveas="~/www/TEST/test.png"):
    w.obj("CMS_zz4l_fai1").setVal(fai)
    frame = w.obj(discname).frame()
    w.obj(pdfname).plotOn(frame)
    c1 = ROOT.TCanvas()
    frame.Draw()
    c1.SaveAs("~/www/TEST/test.png")

