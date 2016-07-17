import getpass
import socket

if "lxplus" in socket.gethostname() and getpass.getuser() == "hroskes":
    repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings"
    plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings/"
    productionforsignalrates = "160225"
else:
    raise ValueError("Who/where are you?")

usedata = False
luminosity = 10
blindcut = "D_bkg_0plus < 0.5"
