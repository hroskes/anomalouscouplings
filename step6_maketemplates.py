from helperstuff import config
from helperstuff.enums import Channel, channels, Systematic, systematics
from helperstuff.samples import Sample
import os
import shutil
import subprocess

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
