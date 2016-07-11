import getpass
import socket

if "lxplus" in socket.gethostname() and getpass.getuser() == "hroskes":
    repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings"
    plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings/"
    CJLSTmaindir_76Xanomalous = "/afs/cern.ch/work/h/hroskes/reweighting_CJLST/CMSSW_7_6_3_patch2/src/ZZAnalysis/AnalysisStep/test/prod/AnomalousCouplingsReweighting/PT13TeV"
    releaseforsignalrates = "76X"
else:
    raise ValueError("Who/where are you?")

usedata = False
luminosity = 10
