from helperstuff import config
from helperstuff.enums import Channel
import os
import subprocess

def buildtemplates(flavor):
    flavor = Channel(flavor)
    print flavor
    if os.path.exists(flavor.templatesfile()):
        return
    subprocess.call([os.path.join(config.repositorydir, "TemplateBuilder/buildTemplate.exe"), "step5_json/templates_{}.json".format(flavor)])
    if not os.path.exists(flavor.templatesfile()):
        raise RuntimeError("Something is wrong!  {} was not created.".format(flavor.templatesfile()))

if __name__ == "__main__":
    buildtemplates("2e2mu")
    buildtemplates("4e")
    buildtemplates("4mu")
