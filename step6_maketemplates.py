from helperstuff import config
from helperstuff.enums import channels, treesystematics, TemplatesFile
from helperstuff.samples import Sample
import os
import shutil
import subprocess
from time import sleep

cmssw = [int(i) for i in os.environ["CMSSW_VERSION"].split("_")[1:]]
if cmssw[0] == 8:
    raise ValueError("TemplateBuilder does not seem to work in CMSSW_8_X; the templates end up filled with NaNs.  Try CMSSW_7_4_X, I have tested that it works there.")
if cmssw[0] == 7 and cmssw[1] == 6:
    print "I haven't tested TemplateBuilder in CMSSW_7_6_X.  It definitely works in 7_4_X, but in 8_X the templates end up filled with NaNs.  Proceed at your own risk."
    print "Continuing in 10 seconds..."
    sleep(10)


def buildtemplates(*args):
    templatesfile = TemplatesFile(*args)
    print templatesfile.channel
    if os.path.exists(templatesfile.templatesfile()):
        return
    subprocess.call([os.path.join(config.repositorydir, "TemplateBuilder/buildTemplate.exe"), templatesfile.jsonfile()])
    if not os.path.exists(templatesfile.templatesfile()):
        raise RuntimeError("Something is wrong!  {} was not created.".format(templatesfile.templatesfile()))

if __name__ == "__main__":
    for channel in channels:
        for systematic in treesystematics:
            for analysis in analyses:
                buildtemplates(channel, systematic, "signal", analysis)
            buildtemplates(channel, "", "bkg", analysis)
    #and copy data
    shutil.copy(Sample("data").withdiscriminantsfile(), os.path.join(config.repositorydir, "step7_templates/"))
