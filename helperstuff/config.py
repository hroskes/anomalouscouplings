import getpass
import socket

if "lxplus" in socket.gethostname() and getpass.getuser() == "hroskes":
    CJLSTmaindir = "/afs/cern.ch/work/h/hroskes/reweighting_CJLST/CMSSW_7_6_3_patch2/src/ZZAnalysis/AnalysisStep/test/prod/AnomalousCouplingsReweighting/PT13TeV"
    repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings"
    CJLSTdirnames_ggH = {
                         "0+": "0PM_v2",
                         "a2": "0PH_v2",
                         "0-": "0M_v1",
                         "L1": "0L1_v2",
                         "fa20.5": "0PHf05ph0_v2",
                         "fa30.5": "0Mf05ph0_v2",
                         "fL10.5": "0L1f05ph0_v2",
                        }
    plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings/"

else:
    raise ValueError("Who/where are you?")

usedata = False
luminosity = 10
