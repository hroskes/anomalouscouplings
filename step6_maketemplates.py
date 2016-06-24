from helperstuff import config
from helperstuff.enums import Channel, channels, Systematic, systematics
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


def buildtemplates(flavor, systematic, isbkg):
    flavor = Channel(flavor)
    systematic = Systematic(systematic)
    print flavor
    if os.path.exists(flavor.templatesfile(systematic, isbkg)):
        return
    subprocess.call([os.path.join(config.repositorydir, "TemplateBuilder/buildTemplate.exe"), flavor.jsonfile(systematic, isbkg)])
    if not os.path.exists(flavor.templatesfile(systematic, isbkg)):
        raise RuntimeError("Something is wrong!  {} was not created.".format(flavor.templatesfile(systematic, isbkg)))

if __name__ == "__main__":
    for channel in channels:
        for systematic in systematics:
            buildtemplates(channel, systematic, False)
        buildtemplates(channel, "", True)
    #and copy data
    shutil.copy(Sample("data").withdiscriminantsfile(), os.path.join(config.repositorydir, "step7_templates/"))
