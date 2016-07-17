from helperstuff import config
from helperstuff.enums import templatesfiles, TemplatesFile
from helperstuff.samples import Sample
import os
import shutil
import subprocess
from time import sleep

cmssw = [int(i) for i in os.environ["CMSSW_VERSION"].split("_")[1:]]
if cmssw[0] == 8:
    raise ValueError("TemplateBuilder does not seem to work in CMSSW_8_X; the templates end up filled with NaNs.  Try CMSSW_7_4_X or CMSSW_7_6_X, I have tested that it works there.")


def buildtemplates(*args):
    templatesfile = TemplatesFile(*args)
    print templatesfile.analysis, templatesfile.channel, templatesfile.signalorbkg, templatesfile.systematic, templatesfile.production
    if os.path.exists(templatesfile.templatesfile()):
        return
    subprocess.call([os.path.join(config.repositorydir, "TemplateBuilder/buildTemplate.exe"), templatesfile.jsonfile()])
    if not os.path.exists(templatesfile.templatesfile()):
        raise RuntimeError("Something is wrong!  {} was not created.".format(templatesfile.templatesfile()))

if __name__ == "__main__":
    for templatesfile in templatesfiles:
        buildtemplates(templatesfile)
        #and copy data
        shutil.copy(Sample("data", templatesfile.production, templatesfile.blindstatus).withdiscriminantsfile(), os.path.join(config.repositorydir, "step7_templates/"))
