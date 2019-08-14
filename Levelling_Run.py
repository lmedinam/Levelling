import sys
import itertools
import numpy as np
from Levelling_Beam import Beam
from Levelling_Config import Config
from Levelling_Others import *
from optparse import OptionParser

### General settings:
parser = OptionParser()
parser.add_option("-n", "--name",    dest = "name",    help = "Name of the simulation settings configuration in the Config module to be run", default = "Jun17Paper_Base_nom")
parser.add_option("-v", "--version", dest = "version", help = "Version of the Luminosity module to be used for the simulation",               default = "062")
parser.add_option("-t", "--table",   dest = "table",   help = "Stop simulation after simulation settings and virtual luminosities table",     default = False)

(options, args) = parser.parse_args()

### Create auxiliar classes: bannersclass from Banners module and auxfuncclass from AuxFunc module
bannersclass = Banners()
auxfuncclass = AuxFunc()

### Process the general settings given as input
name    = options.name
version = auxfuncclass.CheckLuminosityVersion(options.version)
table   = auxfuncclass.StringBinaryToLogical(options.table)

### Start of program
bannersclass.header0("\"LUMINOSITY LEVELLING\"")

### Create configclass from Config module
print "\n> Loading CONFIG class...\n"
configclass    = Config(name)

listvers1 = ["011", "012"]
listvers2 = ["011", "012", "013"]
listvers3 = ["013", "021", "031", "032"]
listvers4 = [str(i).zfill(3) for i in range(50,55)]
for i in range(60,63):
    listvers4.append(str(i).zfill(3))

### Check the Luminosity module version and load it
if version == False:
    print "\n", "[!]", "Wrong version of Luminosity. Try one of the following:\n\n", 3*" ", 
    for i in auxfuncclass.ListLuminosityVersions(): print i,
    bannersclass.header0("EXIT")
    sys.exit()
else:
    print "*", "{0:36} {1:10}".format("Version of Luminosity class", "-"), version
    if version in listvers4:  modulename = "Levelling_Luminosity_v" + str(version)
    else:                     modulename = "Luminosity_v" + str(version)
    LumiModule = map(__import__, [modulename])

### Create beamclass class from Beam module, and luminosityclass from the Luminosity module
print "\n> Loading BEAM class..."
beamclass = Beam(configclass)
print "\n> Loading LUMINOSITY class..."
if   version in listvers1: luminosityclass = LumiModule[0].Luminosity(configclass, name)
elif version in listvers3: luminosityclass = LumiModule[0].Luminosity(configclass, name, table)
else:                      luminosityclass = LumiModule[0].Luminosity(configclass, name, table, version)

### Print the simulation settings
bannersclass.header1("SIMULATION SETTINGS")
beamclass.printBeamParam(configclass._longdens)
if   version in listvers2: luminosityclass.printLumiParam()
elif version in listvers4: luminosityclass.PrintLumiParam(beamclass)
else:                      luminosityclass.printLumiParam(beamclass)

### Start luminosityclassnosity levelling simulation
if version in listvers4: luminosityclass.DoFill(beamclass)
else:                    luminosityclass.doFill(beamclass)

### End of program
bannersclass = Banners()
bannersclass.header0("END OF \"LUMINOSITY LEVELLING\"")