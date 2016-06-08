from helperstuff import config
from helperstuff.enums import Channel
from helperstuff.samples import Sample
import os
import shutil
import subprocess

def buildtemplates(flavor, isbkg):
    flavor = Channel(flavor)
    print flavor
    if os.path.exists(flavor.templatesfile(isbkg)):
        return
    subprocess.call([os.path.join(config.repositorydir, "TemplateBuilder/buildTemplate.exe"), flavor.jsonfile(isbkg)])
    if not os.path.exists(flavor.templatesfile(isbkg)):
        raise RuntimeError("Something is wrong!  {} was not created.".format(flavor.templatesfile(isbkg)))

if __name__ == "__main__":
    buildtemplates("2e2mu", True)
    buildtemplates("4e", True)
    buildtemplates("4mu", True)
    buildtemplates("2e2mu", False)
    buildtemplates("4e", False)
    buildtemplates("4mu", False)
    #and copy data
    shutil.copy(Sample("data").withdiscriminantsfile(), os.path.join(config.repositorydir, "step7_templates/"))
