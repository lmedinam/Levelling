from Levelling_Densities import Densities
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n", "--name",    dest = "name",    help = "Name of file with the Levelling results from which distributions will be plot (including results_v0...out)", default = None)
parser.add_option("-i", "--ip",      dest = "ip",      help = "IP (0 for IP1 or 2 for IP8)", default = 0)
parser.add_option("-s", "--step",    dest = "step",    help = "Step of levelling: use FIRST or LAST or integer", default = None)

(options, args) = parser.parse_args() 

name = options.name
ip   = options.ip
step = options.step

classdist = Densities()

classdist.GetDistribution(name, ip, step)