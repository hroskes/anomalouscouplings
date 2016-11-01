from copy import deepcopy
from helperstuff.optimizesmoothing import TemplateIterate
from helperstuff.templates import templatesfiles
from helperstuff.submitjob import submitjob
import shutil

def runiteration(iternumber):
    bkptemplatesdir = os.path.join(config.repositorydir, "step7_templates", "bkp_iter{}".format(iternumber-1))
    bkpjsondir = os.path.join(config.repositorydir, "step5_json", "bkp_iter{}".format(iternumber-1))
    try:
        os.makedirs(bkptemplatesdir)
    except OSError:
        pass
    try:
        os.makedirs(bkpjsondir)
    except OSError:
        pass

    smoothingparametersdict = deepcopy(template.getsmoothingparametersdict())

    nchangedfiles = 0
    for templatesfile in templatesfiles:
        anychange = False
        for template in templatesfile.templates():
            template = TemplateIterate(template)
            nextiteration = template.getnextiteration()
            if nextiteration is None:
                continue
            previousparameters = template.smoothingparameters
            if nextiteration == previousparameters:
                continue
            anychange = True

            (smoothingparametersdict
                                    [str(self.productionmode)]
                                    [str(self.category)]
                                    [str(self.analysis)]
                                    [str(self.whichproddiscriminants)]
                                    [str(self.hypothesis)]
            ) = nextiteration

        if anychange and os.path.exists(templatesfile.templatesfile()):
            nchangedfiles += 1
            shutil.move(templatesfile.templatesfile(), bkptemplatesdir)
        else:
            shutil.copy(templatesfile.templatesfile(), bkptemplatesdir)
        shutil.move(templatesfile.jsonfile(), bkpjsondir)

    with cd(config.repositorydir):
        subprocess.check_call(["python", "step4_makejson.py"])
        from helperstuff import ZX, CJLSTscripts #make sure the C++ stuff is compiled
                                                 #otherwise the jobs will try to compile them simultaneously
                                                 #and fail
        jobids = []
        for i in range(nchangedfiles):
            jobids.append(submitjob("python step6_maketemplates.py", jobname="iter{}_job{}".format(iternumber, i), jobtime="2:0:0"))

