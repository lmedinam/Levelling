import numpy as np
from metaclass import *

class Banners:
  
  def __init__(self):
    pass
  
  def header0(self, text):
    print "\n"
    print 112*"*"
    print "*", 108*" ", "*"
    print "*", text.center(108), "*", "\n",
    print "*", 108*" ", "*"
    print 112*"*", "\n"
  
  def header1(self, text):
    print "\n"
    print 112*"*"
    print "*", text.center(108), "*", "\n",
    print 112*"*", "\n"
      
  def header2(self, text):
    print "\n"
    print "+"+108*"-"+"+"
    print "|", text.center(106), "|", "\n",
    print "+"+108*"-"+"+", "\n"
    
  def header3(self, text):
    print "\n"
    print "+"+60*"-"+"+"
    print "|", text.center(58), "|", "\n",
    print "+"+60*"-"+"+", "\n"
  
  def header4(self, text):
    print "\n"
    print "*"+60*"*"+"*"
    print "*", text.center(58), "*", "\n",
    print "*"+60*"*"+"*"
    
  def header5(self, text):
    print "\n"
    print text.center(112), "\n",
    print 112*"-", "\n"
    
  def header6(self, text):
    print "\n"
    print text.center(112), "\n"
  
class AuxFunc:

    def __init__(self):
        pass
    
    def ListLuminosityVersions(self):
        subversions = []
        subversions.append( ["0"+str(x) for x in range(11,14)] )  # 11 to 13
        subversions.append( ["0"+str(x) for x in range(21,22)] )  # 21
        subversions.append( ["0"+str(x) for x in range(31,36)] )  # 31 to 35
        subversions.append( ["0"+str(x) for x in range(40,44)] )  # 40 to 43
        subversions.append( ["0"+str(x) for x in range(50,55)] )  # 50 to 54
        subversions.append( ["0"+str(x) for x in range(60,63)] )  # 60, 62
        allversions = []
        for item in subversions: allversions.extend(item)
        return allversions
    
    def CheckLuminosityVersion(self, version):
        """ Returns version if it matches one of the existing versions of Luminosity module in the list, otherwise it returns logical False. """
        if version in self.ListLuminosityVersions(): return version
        else:                      return False

    def StringBinaryToLogical(self,  value):
        """ Return logical True in the listed cases, otherwise returns logical False. """
        validtrue = []
        validtrue.append(1)
        validtrue.append("1")
        validtrue.append("t")
        validtrue.append("T")
        validtrue.append("y")
        validtrue.append("Y")
        validtrue.append("true")
        validtrue.append("True")
        validtrue.append("yes")
        validtrue.append("Yes")
        if value in validtrue: return True
        else:                  return False
        
    def CountStringInList(self, string, array):
      result = 0
      for i in range(len(array)):
        if array[i] == string: result = result + 1
      return result
    
    def DotToDash(self, string):
      result = ""
      for i in range(len(string)):
        if string[i] == ".": result = result + "_"
        else:                result = result + string[i]
      return result
    
    def OldTable(self, resultsname):
      table = twiss(resultsname)
      NAMEtableold = resultsname + "_old"
      FILEtableold = open(NAMEtableold, "w")
      print >> FILEtableold, "E",                                              "%1.1f" %(table.momeV/1e12), "TeV"
      print >> FILEtableold, "Nb",
      try:    print >> FILEtableold, "%1.1f" %table.ppb
      except: print >> FILEtableold, "%1.1f" %table.intensity
      print >> FILEtableold, "nb",                                             int(table.Nbunch)
      print >> FILEtableold, "#collsIP1&5",                                    int(table.nbunch0), int(table.nbunch1)
      print >> FILEtableold, "Ntot",                                           table.totalpart
      print >> FILEtableold, "beam current",                                   table.current
      print >> FILEtableold, "x-sing angle",                                   table.initphi0, table.initphi1, table.initphi2
      print >> FILEtableold, "beam separation",                                table.sepLR0,   table.sepLR1,   table.sepLR2
      print >> FILEtableold, "beta*",                                          "[", str(table.betx0), str(table.bety0) + "]", "[", str(table.betx1), str(table.bety1) + "]", "[", str(table.betx2), str(table.bety2) + "]"
      print >> FILEtableold, "emit",                                           table.epsxn, table.epsyn
      print >> FILEtableold, "emitL",                                          "0"
      print >> FILEtableold, "Espread",                                        table.initdpp
      print >> FILEtableold, "sigma_s",                                        table.initsigs
      print >> FILEtableold, "IBS [h]",                                        table.virtualtauxIBS
      print >> FILEtableold, "IBS L[h]",                                       table.virtualtauzIBS
      print >> FILEtableold, "Piwinski parameter",                             table.virtpiw0, table.virtpiw1, "planes:", int(table.xplane0), int(table.xplane1)
      print >> FILEtableold, "Total loss factor without crab cavity ",          "["  + str(table.virtredfactornoCC0) + ",", str(table.virtredfactornoCC1) + ",", str(table.virtredfactornoCC2) + "]"
      print >> FILEtableold, "Total loss factor with crab cavity",             "["  + str(table.virtredfactor0)     + ",", str(table.virtredfactor1)     + ",", str(table.virtredfactor2)     + "]"
      print >> FILEtableold, "beam-beam without CC",                           "[[" + str(table.virtualxixnoCC0) + ",", str(table.virtualxiynoCC0) + "],", "[" + str(table.virtualxixnoCC1) + str(table.virtualxiynoCC1) + "]]"
      print >> FILEtableold, "beam-beam with CC",                               "[" + str(table.virtualxix0)     + ",", str(table.virtualxiy0)     + "]"
      print >> FILEtableold, "Peak Lumi without CC",                           "["  + str(table.virtluminoCC0) + ",", str(table.virtluminoCC1) + ",", str(table.virtluminoCC2) + "]"
      print >> FILEtableold, "Virtual lumi with CC",                           "["  + str(table.virtlumi0)     + ",", str(table.virtlumi1)     + ",", str(table.virtlumi2)     + "]"
      print >> FILEtableold, "Events per crossing without lev and without CC",
      try:    print >> FILEtableold, table.virtpileupnoCC0
      except: print >> FILEtableold, table.virtpunoCC0
      print >> FILEtableold, "Leveled lumi",
      try:    print >> FILEtableold, "["  + str(table.level_Lumi0) + ",", str(table.level_Lumi1) + ",", str(table.level_Lumi2) + "]"
      except: print >> FILEtableold, "["  + str(table.levlumi0)    + ",", str(table.levlumi1)    + ",", str(table.levlumi2)    + "]"
      print >> FILEtableold, "Events per crossing with lev and with CC",
      try:    print >> FILEtableold, table.firstpileup0
      except: print >> FILEtableold, table.firstpu0
      print >> FILEtableold, "Peak pile-up density",
      try:    print >> FILEtableold, table.maxpeakpileups0
      except: print >> FILEtableold, table.maxppus0
      print >> FILEtableold, "Leveling time",                                  table.lasttime
      print >> FILEtableold, "Number of collisions IP2/IP2",                   "0", int(table.nbunch2)
      print >> FILEtableold, "nb at injection",                                "0"
      print >> FILEtableold, "nb per injection",                               "0"
      print >> FILEtableold, "Ntot per injection",                             "0"
      print >> FILEtableold, "Emittance at injection",                         "0"
      print "> Writting '" + NAMEtableold + "'..."
      FILEtableold.close()

class Constants:
  
  clight = 2.99792458e8                  # Speed of light [m s-1]
  e      = 1.60217657e-19                # Elementary charge [C]
  eps0   = 8.854187817e-12               # Vacuum permittivity [F m-1 = A2 s4 kg-1 m-3 = C2 N-1 m-2 = C V-1 m-1]
  h      = 6.62606957e-34                # Planck constant [J s-1] 
  hbar   = 6.62606957e-34/(2.0*np.pi)    # Reduced Planck constant [J s-1] 
  m      = 1.67262158e-27                # Proton rest mass [kg]
  pmass  = 0.93827231e9                  # Proton rest mass [eV]
  r0     = 1.535e-18                     # Proton radius [m]
  gamma34 = 1.22541670246517764513       # Gamma(3/4)
  gamma14 = 3.62560990822190831193       # Gamma(1/4)
  gamma54 = 0.90640247705547798267       # Gamma(5/4)
  NormRMSGauss  = 2.**(5./4.) * gamma54                         # Recyprocal of the of the normalization factor of the flat distribution
  FactRMSGauss  = np.sqrt( 1./np.sqrt(2.) * gamma14/gamma34 )   # Factor to take rms sigma as input of the flat distribution (supergaussian of order 4) so that their rms values are equal
  NormFWHMGauss = 5.*np.pi/32.                                  # Recyprocal of normalization factor of the flat distribution
  FactFWHMGauss = 2.*np.sqrt( 2.*np.log(2.)/(1.-2.**(-2/5.)) )  # Factor to take rms sigma as input of the lambda distribution so that their FWHM are equal
  
class qGaussianAux:
  
  def __init__(self, sigs, phi, t1, t2, x=None):
    
    cst = Constants()
    
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    
    self.ctmax  =  (cst.FactFWHMGauss*sigs - (ct2+ct1))/2.
    self.ctmin  = -(cst.FactFWHMGauss*sigs - (ct2+ct1))/2.
    self.ctzero = -(ct2+ct1)/2.
    self.smax   =  (cst.FactFWHMGauss*sigs - (ct2-ct1))/2./np.cos(phi/2.)
    self.smin   = -(cst.FactFWHMGauss*sigs - (ct2-ct1))/2./np.cos(phi/2.)
    self.szero  = -(ct2-ct1)/2./np.cos(phi/2.)
    
    # u  will be s
    self.Act = lambda u:  cst.FactFWHMGauss*sigs/2. + u*np.cos(phi/2.) - ct1
    self.Bct = lambda u: -cst.FactFWHMGauss*sigs/2. + u*np.cos(phi/2.) - ct1
    self.Cct = lambda u:  cst.FactFWHMGauss*sigs/2. - u*np.cos(phi/2.) - ct2
    self.Dct = lambda u: -cst.FactFWHMGauss*sigs/2. - u*np.cos(phi/2.) - ct2
    self.Zct = lambda u: self.ctzero
    
    # u will be ct
    self.As  = lambda u: (-cst.FactFWHMGauss*sigs/2. + u + ct1)/np.cos(phi/2.)
    self.Bs  = lambda u: ( cst.FactFWHMGauss*sigs/2. + u + ct1)/np.cos(phi/2.)
    self.Cs  = lambda u: ( cst.FactFWHMGauss*sigs/2. - u - ct2)/np.cos(phi/2.)
    self.Ds  = lambda u: (-cst.FactFWHMGauss*sigs/2. - u - ct2)/np.cos(phi/2.)
    self.Zs  = lambda u: self.szero

    # u will be ct
    self.Aswithx  = lambda u: (-cst.FactFWHMGauss*sigs/2. - x*np.sin(phi/2.) + u + ct1)/np.cos(phi/2.)
    self.Bswithx  = lambda u: ( cst.FactFWHMGauss*sigs/2. - x*np.sin(phi/2.) + u + ct1)/np.cos(phi/2.)
    self.Cswithx  = lambda u: ( cst.FactFWHMGauss*sigs/2. + x*np.sin(phi/2.) - u - ct2)/np.cos(phi/2.)
    self.Dswithx  = lambda u: (-cst.FactFWHMGauss*sigs/2. + x*np.sin(phi/2.) - u - ct2)/np.cos(phi/2.)
    self.Zswithx  = lambda u: self.szero

