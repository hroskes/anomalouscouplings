import getpass
import socket

if "lxplus" in socket.gethostname() and getpass.getuser() == "hroskes":
    repositorydir = "/afs/cern.ch/work/h/hroskes/anomalouscouplings"
    plotsbasedir = "/afs/cern.ch/user/h/hroskes/www/anomalouscouplings/"
else:
    raise ValueError("Who/where are you?")

unblinddata = True
expectedscanluminosity = 10
blindcut = lambda self: self.D_bkg_0plus() < 0.5
productionforsignalrates = "160714"
productionforcombine = "160714"
