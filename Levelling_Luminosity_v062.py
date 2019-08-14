### Edited by Luis Medina. Jul 17, 2017 ###
### epsxn, epsyn -> epsn
### epsxfact, epsyfact -> tauobs
### epsxfact, _epsyfact -> tauCCs

import os
import sys
import fileinput                                       # To overwrite a line in an existent file
import subprocess                                      # To run external program
import numpy as np
from shutil import copyfile                           # To make copy of input file
from scipy.optimize import fminbound                   # To be used with scipy version "0.15.1" (lxplus, July 2016). Use "from scipy.optimize import minimize_scalar" instead for scipy version "0.17.0" (my personal computer)
from scipy.integrate import quad, dblquad, fixed_quad
from Levelling_Others import Constants as cst
from Levelling_Others import Banners
from Levelling_Others import AuxFunc
from Levelling_Others import qGaussianAux
from Levelling_Densities import Densities

class Luminosity:
  
  #============================================================================================================#
  
  ### Beam parameters inherited from Config:

  # (None)
  
  # Beam parameters NOT inherited:
  # _momeV
  # _ppb
  # _dpp
  # _sigs
  # _epsn
  # _epsn_0
  # _epsn_ibs
  # _epsn_sr
  # _epsn_obs
  # _epsn_cc_ip0
  # _epsn_cc_ip1
  # _beta
  # _betamindict
  # _tau_ibs
  # _tau_sr
  # _tau_obs
  # _tau_cc
  # _kappa
  # _kappac
  # _t1
  # _t2
  
  #------------------------------------------------------------------------------------------------------------#
  
  ### Machine parameters inherited from Config:
  
  _incroncc      = []
  _incronccrate  = 0.0
  _niterlev      = 0
  _IBSRF         = True
    
  _adaptivexsing = []
  
  # Machine parameters NOT inherited:
  # _circ
  # _rho
  # _alfmom
  # _tuneb
  # _tunes
  # _Nbunch
  # _ipnames
  # _nip
  # _nbunch
  # _xplane
  # _sepLR
  # _sepLRsteptime
  # _sepLRstep
  # _wcc
  # _oncc
  # _tcc1
  # _tcc2
  
  #------------------------------------------------------------------------------------------------------------#
  
  ### Levelling parameters inherited from Config:
  
  _levtech         = []
  _levvar          = []
  _levlumi         = []
  _levppus         = []
  _ck              = []
  _constbetar      = []
  _sepLRconst      = []
  _longdens        = "Gaussian"
  _p               = 0.0
  _step            = 0.0
  _optimumfill     = True
  _maxfill         = 0.0
  _penstp          = True
  _timepenstp      = 0.0
  _updatepenstp    = False
  
  # Levelling parameters NOT inherited:
  # (None)
  
  #------------------------------------------------------------------------------------------------------------#
  
  ### Bunch length gymnastics parameters:
  
  # (None)
  
  # Bunch length gymnastics parameters NOT inherited:
  # _constlong
  # _ppblong
  # _redlong
  # _minlong
  
  #------------------------------------------------------------------------------------------------------------#
  
  ### Variable separation parameters inherited from Config:
  
  # (None)
  
  ### Variable separation parameters NOT inherited:
  # _redsepLR  = []
  # _minsepLR  = []
  # _ratesepLR = []

  #------------------------------------------------------------------------------------------------------------#
  
  # Integrated luminosity parameters inherited from Config:
  
  _xsec     = 0.0
  _xsecburn = 0.0
  _days     = 0.0
  _eff      = 0.0
  _turnar   = 0.0
  
  # Integrated luminosity parameters NOT inherited:
  # (None)
  
  #============================================================================================================#
  
  ### New levelling parameters:
  
  _flaglevpu       = []
  _timeppuslev     = None
  _flagtimeppuslev = False
  _steplev2        = None
  _timesteplev2    = 0.0
  _flagsteplev2    = False
  
  ### New independent parameters:
  
  _name    = ""  # Name of simulation settings configuration
  _version = ""  # Number of Luminosity module version
  _ID      = ""  # ID tag
  _table   = ""  # Only parameters table requested
  
  #============================================================================================================#
  
  def __init__(self, config, name, table, version):
    
    # Inherited:
    
    self._levtech         = config._levtech
    self._levvar          = config._levvar
    self._levlumi         = config._levlumi
    self._levppus         = config._levppus
    self._constbetar      = config._constbetar
    self._sepLRconst      = config._sepLRconst
    self._longdens        = config._longdens
    self._ck              = config._ck
    self._p               = config._p
    self._step            = config._step
    self._optimumfill     = config._optimumfill
    self._maxfill         = config._maxfill
    self._penstp          = config._penstp
    if self._penstp == True:
      self._timepenstp   = config._timepenstp
      self._updatepenstp = config._updatepenstp
    else:
      self._timepenstp   = 0.0
      self._updatepenstp = False
    
    self._xsec     = config._xsec
    self._xsecburn = config._xsecburn
    self._days     = config._days
    self._eff      = config._eff
    self._turnar   = config._turnar
    
    self._incroncc = config._incroncc
    self._incronccrate = config._incronccrate
    self._niterlev = config._niterlev
    self._IBSRF    = config._IBSRF

    self._adaptivexsing = config._adaptivexsing
    
    self._timeppuslev  = None 
    self._flaglevpu    = [False, False, False]
    self._steplev2     = None
    self._timesteplev2 = 0.0
    self._flagsteplev2 = False
    
    self._name    = name
    self._version = version
    self._ID      = "v" + self._version + "_" + self._name # + "_" + "p" + str("%1.3f" %(self._p)) + "_" + "oncc" + str("%1.2f" %(config._oncc[0]))   # It should be beam._initoncc, but beam is not loaded
    self._table   = table

  ### Functions to update the values of phi, oncc and phiCR
  
  def GetPhi(self, beam, ip):
    """ Updates in the beam class the value of the crossing angle, computed as function of the beam separation and divergence in the corresponding plane, for a given IP """
    diver = np.sqrt( beam._epsn[beam._xplane[ip]] / (beam._gamma*beam._betarel) / beam._beta[ip][beam._xplane[ip]] )
    beam._phi[ip] = diver*beam._sepLR[ip]  # Maybe this? min(diver*beam._sepLR[ip], beam._initphi[ip])
  
  def GetOnccDueToPhi(self, beam, ip):
    """ Updates in the beam class the value of the ratio between the crab angle and the initial crossing angle, for a given IP. """
    beam._oncc[ip] = min(beam._initoncc[ip]*beam._initphi[ip]/beam._phi[ip], 1.0)
    
  def GetOnccDueToPhiCR(self, beam, ip):
    """ Updates in the beam class the value of the ratio between the crab angle and the crossing angle, for a given IP. """
    beam._oncc[ip] = beam._phiCR[ip]/beam._phi[ip]

  def GetPhiCR(self, beam, ip):
    """" Updates in the beam class the value of the crab angle, computed as function of the crossing angle and the ratio to it, for a given IP"""
    beam._phiCR[ip] = beam._oncc[ip]*beam._phi[ip]

  def GetNewPhi(self, beam, ip, output = True):
    """ Updates in the beam class the value of the crossing angle, crab angle and their ratio, for a given IP. To be used in case beta* for that IP has changed and the beam-beam LR separation has to remains constant. """
    self.GetPhi(beam,ip)
    self.GetOnccDueToPhi(beam,ip)
    self.GetPhiCR(beam,ip)
    if output == True:
      print ""
      print "*", "{0:36} {1:10}".format("Crossing angle", "rad"),                 "%1.4e" %beam._phi[ip]
      print "*", "{0:36} {1:10}".format("BB-LR separation", "sigma"),             "%1.8f" %beam._sepLR[ip] #self.GetNewSepLR(beam, ip)
      print "*", "{0:36} {1:10}".format("Crab cavity angle", "rad"),              "%1.4e" %beam._phiCR[ip]
      print "*", "{0:36} {1:10}".format("Crab cavity ON (oncc)", "1"),            "%1.8f" %beam._oncc[ip]
      print ""

  def GetPhiCRDueToAdaptive(self, beam, ip):
    """"  """
    if beam._phi[ip] <= 380e-6: beam._phiCR[ip] = beam._phi[ip]
    else:                       beam._phiCR[ip] = 380e-6
    
  def GetNewSepLR(self, beam, ip):
    """ Returns the beam-beam LR separation as function of the crossing angle and the divergence in the corresponding plane, for a given IP. """
    diver = np.sqrt( beam._epsn[beam._xplane[ip]] / (beam._gamma*beam._betarel) / beam._beta[ip][beam._xplane[ip]] )
    sepLR = beam._phi[ip]/diver
    #return sepLR
    beam._sepLR[ip] = sepLR
    print "\n> phi is the independent variable. Update sepLR in sigma...\n"
  
  def GetOnckDueToPhiCK(self, beam, ip):
    """ Updates in the beam class the value of the ratio between the crab kissing and the crossing angle, for a given IP. """
    beam._onck[ip] = beam._phiCK[ip]/beam._phi[ip]
    
  #============================================================================================================#
  
  ### Functions to compute luminosity

  def KinFact(self, beam, ip):
    """ Returns the 'kinematic factor' but with cos(phi/2) instead of cos^2(phi/2). """
    phi = beam._phi[ip]
    k   = 2.0*np.cos(phi/2.0)
    return k

  def GetL0(self, beam, ip):
    """ Returns the value of the contribution to the luminosity from the transversal coordinates for a collision with a crossing angle ('L0'), for a given IP. To get the luminosity, it has to be multiplied by a missing factor function of pi, sigs (and values of Gamma, in the case of flat beams), and a reduction factor from GetReduc. """
    sig    = beam.getsigma(ip) 
    ppb    = beam._ppb
    frev   = beam._frev
    nbunch = beam._nbunch[ip]
    k      = self.KinFact(beam,ip)
    L0     = k * ppb**2 * frev * nbunch / (8 * np.pi * sig[0] * sig[1])*1e-4
    return L0

  def GetReduc(self, beam, ip):
    """ Returns the appropiate reduction factor depending of the type of collision with/without crab cavities, and/or crab kissing, for a given IP """
    t1       = beam._t1
    t2       = beam._t2
    sigs     = beam._sigs
    bet      = beam._beta[ip]
    tcc1     = beam._tcc1[ip]
    tcc2     = beam._tcc2[ip]
    xplane   = beam._xplane[ip]
    ck       = self._ck[ip]
    longdens = self._longdens
    betc     = bet[xplane]     # Beta in crossing plane
    betp     = bet[1-xplane]   # Beta in paralell separation plane
    t1c      = tcc1[xplane]
    t1p      = tcc1[1-xplane]
    t2c      = tcc2[xplane]
    t2p      = tcc2[1-xplane]
    sig      = beam.getsigma(ip)
    sigc     = sig[xplane]
    sigp     = sig[1-xplane]
    oncc     = beam._oncc[ip]
    wcc      = beam._wcc
    phi      = beam._phi[ip]
    phiCR    = beam._phiCR[ip]
    phiCK    = beam._phiCK[ip]
    parsep   = beam._parsep
    densities = Densities()
    l0 = self.GetL0(beam,ip)
    printcomp = None
    # if printcomp != None:
    #   if ip == 0:
    #     print ""
    #     result = densities.Evaluate_integrand( 5e-2, -5e-2, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, printcomp)
    #     print result
    #     result = densities.Evaluate_integrand( 5e-2,   0.0, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, printcomp)
    #     print result
    #     result = densities.Evaluate_integrand( 5e-2,  5e-2, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, printcomp)
    #     print result
    #     print ""
    #     result = densities.Evaluate_integrand(  0.0, -5e-2, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, printcomp)
    #     print result
    #     result = densities.Evaluate_integrand(  0.0,   0.0, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, printcomp)
    #     print result
    #     result = densities.Evaluate_integrand(  0.0,  5e-2, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, printcomp)
    #     print result
    #     print ""
    #     result = densities.Evaluate_integrand(-5e-2, -5e-2, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, printcomp)
    #     print result
    #     result = densities.Evaluate_integrand(-5e-2,   0.0, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, printcomp)
    #     print result
    #     result = densities.Evaluate_integrand(-5e-2,  5e-2, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, printcomp)
    #     print result
    #     print ""
    reduc = densities.IntAll(phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) # [1]
    # if printcomp != None:
    #   print reduc
    #   print reduc*l0
    #   print ""
    return reduc

  def getlumi(self, beam):
    """ Returns three arrays: the reduced luminosity in all IPs, L0 in all IPs, and the reduction factors in all IPs. """
    lumi  = []
    lumi0 = []
    reduc = []
    for i in range(beam._nip):
      # Trick to avoid recalculating for ip  = 1 if its settigns are equal to the ones for ip = 0. (That is IP1 = IP5).
      if (i != 1 or not beam._identicalIP01):
        r  = self.GetReduc(beam,i)
        l0 = self.GetL0(beam,i)
      lumi.append(l0*r)
      lumi0.append(l0)
      reduc.append(r)
    return lumi, lumi0, reduc

  def getlumiip(self, beam, ip):
    """ Returns the reduced luminosity, L0, and the reduction factor for a given IP. """
    #print "\nCalling GetReduc..."
    r  = self.GetReduc(beam,ip)
    #print "\nCalling GetL0...", 
    l0 = self.GetL0(beam,ip)
    #print "l0 =", l0
    return l0*r, l0, r
  
  #============================================================================================================#

  ### Pile up

  def GetPileUp(self, frev, levlumi, nbunch):
    """ Returns the pile up, that is, the number of events per bunch crossing. """
    return levlumi * self._xsec * 1.0e-27 / nbunch / frev
  
  #============================================================================================================#
  
  ### Luminous region and time, and peak luminous region and time
  
  def GetLumRegion(self, beam, ip):
    """ Returns the RMS luminous region, for a given IP. """
    t1       = beam._t1
    t2       = beam._t2
    sigs     = beam._sigs
    bet      = beam._beta[ip]
    tcc1     = beam._tcc1[ip]
    tcc2     = beam._tcc2[ip]
    xplane   = beam._xplane[ip]
    ck       = self._ck[ip]
    longdens = self._longdens
    betc     = bet[xplane]     # Beta in crossing plane
    betp     = bet[1-xplane]   # Beta in paralell separation plane
    t1c      = tcc1[xplane]
    t1p      = tcc1[1-xplane]
    t2c      = tcc2[xplane]
    t2p      = tcc2[1-xplane]
    sig      = beam.getsigma(ip)
    sigc     = sig[xplane]
    sigp     = sig[1-xplane]
    oncc     = beam._oncc[ip]
    wcc      = beam._wcc
    phi      = beam._phi[ip]
    phiCR    = beam._phiCR[ip]
    phiCK    = beam._phiCK[ip]
    parsep   = beam._parsep
    densities = Densities() 
    norm = densities.IntAll(phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) # [1]
    if ip == 2:
      if oncc == 0.0:
        if longdens == "Gaussian":
          i    =    quad(lambda s    :     s**2*densities.integrand_noCC_parsep_Gaussian(     s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       -np.inf, np.inf)
        elif longdens == "RF800":
          i    = dblquad(lambda s, ct:     s**2*densities.integrand_noCC_parsep_RF800(    ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
        elif oncc == 0.0 and longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          i1   = dblquad(lambda s, ct:     s**2*densities.integrand_noCC_parsep_qGaussian(ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs )
          i2   = dblquad(lambda s, ct:     s**2*densities.integrand_noCC_parsep_qGaussian(ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs )
          i3   = dblquad(lambda s, ct:     s**2*densities.integrand_noCC_parsep_qGaussian(ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs )
          i4   = dblquad(lambda s, ct:     s**2*densities.integrand_noCC_parsep_qGaussian(ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs )
          i = (i1[0]+i2[0]+i3[0]+i4[0],)
        else:
          sys.exit("\n[!] ERROR.\n")
      else:
        sys.exit("\n[!] SETTINGS HAVE TO BE ONCC = 0 FOR IP8.\n")
    else:
      if ck == False:
        if oncc == 0.0:
          if longdens == "Gaussian":
            i    =    quad(lambda s    :     s**2*densities.integrand_noCC_Gaussian(     s, phi, sigs, betc, sigc, betp, t1, t2),                       -np.inf, np.inf)
          elif longdens == "RF800":
            i    = dblquad(lambda s, ct:     s**2*densities.integrand_noCC_RF800(    ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            i1   = dblquad(lambda s, ct:     s**2*densities.integrand_noCC_qGaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs )
            i2   = dblquad(lambda s, ct:     s**2*densities.integrand_noCC_qGaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs )
            i3   = dblquad(lambda s, ct:     s**2*densities.integrand_noCC_qGaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs )
            i4   = dblquad(lambda s, ct:     s**2*densities.integrand_noCC_qGaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs )
            i = (i1[0]+i2[0]+i3[0]+i4[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
        else:
          if longdens == "Gaussian":
            i    = dblquad(lambda s, ct:     s**2*densities.integrand_CC_Gaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
          elif longdens == "RF800":
            i    = dblquad(lambda s, ct:     s**2*densities.integrand_CC_RF800(      ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            i1   = dblquad(lambda s, ct:     s**2*densities.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs )
            i2   = dblquad(lambda s, ct:     s**2*densities.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs )
            i3   = dblquad(lambda s, ct:     s**2*densities.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs )
            i4   = dblquad(lambda s, ct:     s**2*densities.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs )
            i = (i1[0]+i2[0]+i3[0]+i4[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          i    = dblquad(lambda s, ct:     s**2*densities.integrand_CC_CK_Gaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
          
        elif longdens == "RF800":
          i    = dblquad(lambda s, ct:     s**2*densities.integrand_CC_CK_RF800(     ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          i1   = dblquad(lambda s, ct:     s**2*densities.integrand_CC_CK_qGaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs )
          i2   = dblquad(lambda s, ct:     s**2*densities.integrand_CC_CK_qGaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs )
          i3   = dblquad(lambda s, ct:     s**2*densities.integrand_CC_CK_qGaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs )
          i4   = dblquad(lambda s, ct:     s**2*densities.integrand_CC_CK_qGaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs )
          i = (i1[0]+i2[0]+i3[0]+i4[0],)
        else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
    result = np.sqrt(i[0]/norm)
    return result

  def GetLumTime(self, beam, ip):
    """ Returns the RMS luminous time, for a given IP. """
    t1       = beam._t1
    t2       = beam._t2
    sigs     = beam._sigs
    bet      = beam._beta[ip]
    tcc1     = beam._tcc1[ip]
    tcc2     = beam._tcc2[ip]
    xplane   = beam._xplane[ip]
    ck       = self._ck[ip]
    longdens = self._longdens
    betc     = bet[xplane]     # Beta in crossing plane
    betp     = bet[1-xplane]   # Beta in paralell separation plane
    t1c      = tcc1[xplane]
    t1p      = tcc1[1-xplane]
    t2c      = tcc2[xplane]
    t2p      = tcc2[1-xplane]
    sig      = beam.getsigma(ip)
    sigc     = sig[xplane]
    sigp     = sig[1-xplane]
    oncc     = beam._oncc[ip]
    wcc      = beam._wcc
    phi      = beam._phi[ip]
    phiCR    = beam._phiCR[ip]
    phiCK    = beam._phiCK[ip]
    parsep   = beam._parsep
    densities = Densities()
    norm = densities.IntAll(phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) # [1]
    if ip == 2:
      if oncc == 0.0:
        if longdens == "Gaussian":
          i    = [np.NaN]
        elif longdens == "RF800":
          i    = dblquad(lambda ct, s:    ct**2*densities.integrand_noCC_parsep_RF800(   ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          i1   = dblquad(lambda ct, s:    ct**2*densities.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.szero,  qGAux.smax,   qGAux.Zct, qGAux.Cct)
          i2   = dblquad(lambda ct, s:    ct**2*densities.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.smin,   qGAux.szero,  qGAux.Zct, qGAux.Act)
          i3   = dblquad(lambda ct, s:    ct**2*densities.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.smin,   qGAux.szero,  qGAux.Dct, qGAux.Zct)
          i4   = dblquad(lambda ct, s:    ct**2*densities.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.szero,  qGAux.smax,   qGAux.Bct, qGAux.Zct)
          i = (i1[0]+i2[0]+i3[0]+i4[0],)
        else:
          sys.exit("\n[!] ERROR.\n")
      else:
        sys.exit("\n[!] SETTINGS HAVE TO BE ONCC = 0 FOR IP8.\n")
    else:
      if ck == False:
        if oncc == 0.0:
          if longdens == "Gaussian":
            i    = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_Gaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc,   0.0, 0.0, 0.0), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
          elif longdens == "RF800":
            i    = dblquad(lambda ct, s:    ct**2*densities.integrand_noCC_RF800(   ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            i1   = dblquad(lambda ct, s:    ct**2*densities.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.szero,  qGAux.smax,   qGAux.Zct, qGAux.Cct)
            i2   = dblquad(lambda ct, s:    ct**2*densities.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.smin,   qGAux.szero,  qGAux.Zct, qGAux.Act)
            i3   = dblquad(lambda ct, s:    ct**2*densities.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.smin,   qGAux.szero,  qGAux.Dct, qGAux.Zct)
            i4   = dblquad(lambda ct, s:    ct**2*densities.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.szero,  qGAux.smax,   qGAux.Bct, qGAux.Zct)
            i = (i1[0]+i2[0]+i3[0]+i4[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
        else:
          if longdens == "Gaussian":
            i    = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_Gaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
          elif longdens == "RF800":
            i    = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_RF800(     ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            i1   = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_qGaussian(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.szero,  qGAux.smax,  qGAux.Zct, qGAux.Cct)
            i2   = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_qGaussian(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.smin,   qGAux.szero, qGAux.Zct, qGAux.Act)
            i3   = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_qGaussian(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.smin,   qGAux.szero, qGAux.Dct, qGAux.Zct)
            i4   = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_qGaussian(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.szero,  qGAux.smax,  qGAux.Bct, qGAux.Zct)
            i = (i1[0]+i2[0]+i3[0]+i4[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          i    = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_CK_Gaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
        elif longdens == "RF800":
          i    = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_CK_RF800(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          i1   = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_CK_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.szero,  qGAux.smax,  qGAux.Zct, qGAux.Cct)
          i2   = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_CK_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.smin,   qGAux.szero, qGAux.Zct, qGAux.Act)
          i3   = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_CK_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.smin,   qGAux.szero, qGAux.Dct, qGAux.Zct)
          i4   = dblquad(lambda ct, s:    ct**2*densities.integrand_CC_CK_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.szero,  qGAux.smax,  qGAux.Bct, qGAux.Zct)
          i = (i1[0]+i2[0]+i3[0]+i4[0],)
        else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
    result = np.sqrt(i[0]/norm)/cst.clight
    return result
  
  def GetPeakLumRegion(self, beam, ip):
    """ Returns the peak luminous region (integral on ct), that is, the luminous region assuming s = 0, for a given IP. """
    t1       = beam._t1
    t2       = beam._t2
    sigs     = beam._sigs
    bet      = beam._beta[ip]
    tcc1     = beam._tcc1[ip]
    tcc2     = beam._tcc2[ip]
    xplane   = beam._xplane[ip]
    ck       = self._ck[ip]
    longdens = self._longdens
    betc     = bet[xplane]     # Beta in crossing plane
    betp     = bet[1-xplane]   # Beta in paralell separation plane
    t1c      = tcc1[xplane]
    t1p      = tcc1[1-xplane]
    t2c      = tcc2[xplane]
    t2p      = tcc2[1-xplane]
    sig      = beam.getsigma(ip)
    sigc     = sig[xplane]
    sigp     = sig[1-xplane]
    oncc     = beam._oncc[ip]
    wcc      = beam._wcc
    phi      = beam._phi[ip]
    phiCR    = beam._phiCR[ip]
    phiCK    = beam._phiCK[ip]
    parsep   = beam._parsep
    densities = Densities()
    norm = densities.IntAll(phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) # [1]
    if ip == 2:
      if oncc == 0.0:
        if longdens == "Gaussian":
          i    =                         [densities.integrand_noCC_parsep_Gaussian(    0, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                          None]
        elif longdens == "RF800":
          i    =    quad(lambda ct:    densities.integrand_noCC_parsep_RF800(      ct, 0, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                          -np.inf, np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          i12  =    quad(lambda ct:    densities.integrand_noCC_parsep_qGaussian(     ct, 0, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                          qGAux.ctzero, qGAux.ctmax )
          i34  =    quad(lambda ct:    densities.integrand_noCC_parsep_qGaussian(     ct, 0, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                          qGAux.ctmin,  qGAux.ctzero)
          i = (i12[0]+i34[0],)
        else:
          sys.exit("\n[!] ERROR.\n")
      else:
        sys.exit("\n[!] SETTINGS HAVE TO BE ONCC = 0 FOR IP8.\n")
    else:
      if ck == False:
        if oncc == 0.0:
          if longdens == "Gaussian":
            i    =                         [densities.integrand_noCC_Gaussian(    0, phi, sigs, betc, sigc, betp, t1, t2),                          None]
          elif longdens == "RF800":
            i    =    quad(lambda ct:    densities.integrand_noCC_RF800(      ct, 0, phi, sigs, betc, sigc, betp, t1, t2),                          -np.inf, np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            i12  =    quad(lambda ct:    densities.integrand_noCC_qGaussian(     ct, 0, phi, sigs, betc, sigc, betp, t1, t2),                          qGAux.ctzero, qGAux.ctmax )
            i34  =    quad(lambda ct:    densities.integrand_noCC_qGaussian(     ct, 0, phi, sigs, betc, sigc, betp, t1, t2),                          qGAux.ctmin,  qGAux.ctzero)
            i = (i12[0]+i34[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
        else:
          if longdens == "Gaussian":
            i    =    quad(lambda ct:    densities.integrand_CC_Gaussian(     ct, 0, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c),              -np.inf, np.inf)
          elif longdens == "RF800":
            i    =    quad(lambda ct:    densities.integrand_CC_RF800(        ct, 0, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c),              -np.inf, np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            i12  =    quad(lambda ct:    densities.integrand_CC_qGaussian(        ct, 0, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c),              qGAux.ctzero, qGAux.ctmax )
            i34  =    quad(lambda ct:    densities.integrand_CC_qGaussian(        ct, 0, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c),              qGAux.ctmin,  qGAux.ctzero)
            i = (i12[0]+i34[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          i    =    quad(lambda ct:    densities.integrand_CC_CK_Gaussian(    ct, 0, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf)
        elif longdens == "RF800":
          i    =    quad(lambda ct:    densities.integrand_CC_CK_RF800(       ct, 0, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          i12  =    quad(lambda ct:    densities.integrand_CC_CK_qGaussian(      ct, 0, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.ctzero, qGAux.ctmax )
          i34  =    quad(lambda ct:    densities.integrand_CC_CK_qGaussian(      ct, 0, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.ctmin,  qGAux.ctzero)
          i = (i12[0]+i34[0],)
        else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
    result = i[0]/norm
    return result
  
  def GetPeakLumTime(self, beam, ip):
    """ Returns the peak luminous time (integral on s), that is, the luminous region assuming ct = 0, for a given IP. """
    t1       = beam._t1
    t2       = beam._t2
    sigs     = beam._sigs
    bet      = beam._beta[ip]
    tcc1     = beam._tcc1[ip]
    tcc2     = beam._tcc2[ip]
    xplane   = beam._xplane[ip]
    ck       = self._ck[ip]
    longdens = self._longdens
    betc     = bet[xplane]     # Beta in crossing plane
    betp     = bet[1-xplane]   # Beta in paralell separation plane
    t1c      = tcc1[xplane]
    t1p      = tcc1[1-xplane]
    t2c      = tcc2[xplane]
    t2p      = tcc2[1-xplane]
    sig      = beam.getsigma(ip)
    sigc     = sig[xplane]
    sigp     = sig[1-xplane]
    oncc     = beam._oncc[ip]
    wcc      = beam._wcc
    phi      = beam._phi[ip]
    phiCR    = beam._phiCR[ip]
    phiCK    = beam._phiCK[ip]
    parsep   = beam._parsep
    densities = Densities()
    norm = densities.IntAll(phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) # [1]
    if ip == 2:
      if oncc == 0.0:
        if longdens == "Gaussian":
          i    = [np.NaN]
        elif longdens == "RF800":
          i    =    quad(lambda  s:    densities.integrand_noCC_parsep_RF800(       0, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                        -np.inf, np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          i12  =    quad(lambda  s:    densities.integrand_noCC_parsep_qGaussian(      0, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                        qGAux.szero,  qGAux.smax )
          i34  =    quad(lambda  s:    densities.integrand_noCC_parsep_qGaussian(      0, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                        qGAux.smin,   qGAux.szero)
          i = (i12[0]+i34[0],)
        else:
          sys.exit("\n[!] ERROR.\n")
      else:
        sys.exit("\n[!] SETTINGS HAVE TO BE ONCC = 0 FOR IP8.\n")
    else:
      if ck == False:
        if oncc == 0.0:
          if longdens == "Gaussian":
            i    =    quad(lambda  s:    densities.integrand_CC_Gaussian(      0, s, phi, sigs, betc, sigc, betp, t1, t2, wcc,   0.0, 0.0, 0.0),  -np.inf, np.inf)
          elif longdens == "RF800":
            i    =    quad(lambda  s:    densities.integrand_noCC_RF800(       0, s, phi, sigs, betc, sigc, betp, t1, t2),                        -np.inf, np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            i12  =    quad(lambda  s:    densities.integrand_noCC_qGaussian(      0, s, phi, sigs, betc, sigc, betp, t1, t2),                        qGAux.szero,  qGAux.smax )
            i34  =    quad(lambda  s:    densities.integrand_noCC_qGaussian(      0, s, phi, sigs, betc, sigc, betp, t1, t2),                        qGAux.smin,   qGAux.szero)
            i = (i12[0]+i34[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
        else:
          if longdens == "Gaussian":
            i    =    quad(lambda  s:    densities.integrand_CC_Gaussian(      0, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf)
          elif longdens == "RF800":
            i    =    quad(lambda  s:    densities.integrand_CC_RF800(         0, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            i12  =    quad(lambda  s:    densities.integrand_CC_qGaussian(        0, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.szero,  qGAux.smax )
            i34  =    quad(lambda  s:    densities.integrand_CC_qGaussian(        0, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.smin,   qGAux.szero)
            i = (i12[0]+i34[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          i    =    quad(lambda  s:    densities.integrand_CC_CK_Gaussian(     0, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf)
        elif longdens == "RF800":
          i    =    quad(lambda  s:    densities.integrand_CC_CK_RF800(        0, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          i12  =    quad(lambda  s:    densities.integrand_CC_CK_qGaussian(       0, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.szero,  qGAux.smax )
          i34  =    quad(lambda  s:    densities.integrand_CC_CK_qGaussian(       0, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.smin,   qGAux.szero)
          i = (i12[0]+i34[0],)
        else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
    result = i[0]/norm*cst.clight
    return result
  
  def GetIntSPUS2(self, beam, ip, pu):
    """ Returns the integral on s of the squared of the PU line density for a given IP. """
    t1       = beam._t1
    t2       = beam._t2
    sigs     = beam._sigs
    bet      = beam._beta[ip]
    tcc1     = beam._tcc1[ip]
    tcc2     = beam._tcc2[ip]
    xplane   = beam._xplane[ip]
    ck       = self._ck[ip]
    longdens = self._longdens
    betc     = bet[xplane]
    betp     = bet[1-xplane]
    t1c      = tcc1[xplane]
    t1p      = tcc1[1-xplane]
    t2c      = tcc2[xplane]
    t2p      = tcc2[1-xplane]
    sig      = beam.getsigma(ip)
    sigc     = sig[xplane]
    sigp     = sig[1-xplane]
    oncc     = beam._oncc[ip]
    wcc      = beam._wcc
    phi      = beam._phi[ip]
    phiCR    = beam._phiCR[ip]
    phiCK    = beam._phiCK[ip]
    parsep   = beam._parsep
    densities = Densities()
    norm = densities.IntAll(phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) # [1]
    pus2 = lambda  s: pow(pu * densities.IntCT(s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) / norm, 2) # [event^2/m^2]
    if ck == False:
      if oncc == 0.0:
        if longdens == "Gaussian":
          result = quad(pus2, -np.inf, np.inf)
        elif longdens == "RF800":
          result = quad(pus2, -np.inf, np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          result1 = quad(lambda  s:    pus2(s),  qGAux.szero,  qGAux.smax )
          result2 = quad(lambda  s:    pus2(s),  qGAux.smin,   qGAux.szero)
          result = (result1[0]+result2[0],)
        else:
          sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          result = quad(pus2, -np.inf, np.inf)
        elif longdens == "RF800":
          result = quad(pus2, -np.inf, np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          result1 = quad(lambda  s:    pus2(s),  qGAux.szero,  qGAux.smax )
          result2 = quad(lambda  s:    pus2(s),  qGAux.smin,   qGAux.szero)
          result = (result1[0]+result2[0],)
        else:
          sys.exit("\n[!] Longitudinal density unknown.\n")
    else:
      if longdens == "Gaussian":
        result = quad(pus2, -np.inf, np.inf)
      elif longdens == "RF800":
        result = quad(pus2, -np.inf, np.inf)
      elif longdens == "qGaussian":
        qGAux = qGaussianAux(sigs, phi, t1, t2)
        result1 = quad(lambda  s:    pus2(s),  qGAux.szero,  qGAux.smax )
        result2 = quad(lambda  s:    pus2(s),  qGAux.smin,   qGAux.szero)
        result = (result1[0]+result2[0],)
      else:
        sys.exit("\n[!] Longitudinal density unknown.\n")
    return result[0] # [event^2/m]
  
  def GetIntCTPUT2(self, beam, ip, pu):
    """ Returns the integral on s of the squared of the PU time density for a given IP. """
    t1       = beam._t1
    t2       = beam._t2
    sigs     = beam._sigs
    bet      = beam._beta[ip]
    tcc1     = beam._tcc1[ip]
    tcc2     = beam._tcc2[ip]
    xplane   = beam._xplane[ip]
    ck       = self._ck[ip]
    longdens = self._longdens
    betc     = bet[xplane]
    betp     = bet[1-xplane]
    t1c      = tcc1[xplane]
    t1p      = tcc1[1-xplane]
    t2c      = tcc2[xplane]
    t2p      = tcc2[1-xplane]
    sig      = beam.getsigma(ip)
    sigc     = sig[xplane]
    sigp     = sig[1-xplane]
    oncc     = beam._oncc[ip]
    wcc      = beam._wcc
    phi      = beam._phi[ip]
    phiCR    = beam._phiCR[ip]
    phiCK    = beam._phiCK[ip]
    parsep   = beam._parsep
    densities = Densities()
    norm = densities.IntAll(phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) # [1]
    put2 = lambda ct: pow(pu * densities.IntS(ct, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) / norm * cst.clight, 2) # [event^2/s^2]
    if ck == False:
      if oncc == 0.0:
        if longdens == "Gaussian":
          result = quad(put2, -np.inf, np.inf)
        elif longdens == "RF800":
          result = quad(put2, -np.inf, np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          result1 = quad(lambda ct:    put2(ct), qGAux.ctzero, qGAux.ctmax )
          result2 = quad(lambda ct:    put2(ct), qGAux.ctmin,  qGAux.ctzero)
          result = (result1[0]+result2[0],)
        else:
          sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          result = quad(put2, -np.inf, np.inf)
        elif longdens == "RF800":
          result = quad(put2, -np.inf, np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          result1 = quad(lambda ct:    put2(ct), qGAux.ctzero, qGAux.ctmax )
          result2 = quad(lambda ct:    put2(ct), qGAux.ctmin,  qGAux.ctzero)
          result = (result1[0]+result2[0],)
        else:
          sys.exit("\n[!] Longitudinal density unknown.\n")
    else:
      if longdens == "Gaussian":
        result = quad(put2, -np.inf, np.inf)
      elif longdens == "RF800":
        result = quad(put2, -np.inf, np.inf)
      elif longdens == "qGaussian":
        qGAux = qGaussianAux(sigs, phi, t1, t2)
        result1 = quad(lambda ct:    put2(ct), qGAux.ctzero, qGAux.ctmax )
        result2 = quad(lambda ct:    put2(ct), qGAux.ctmin,  qGAux.ctzero)
        result = (result1[0]+result2[0],)
      else:
        sys.exit("\n[!] Longitudinal density unknown.\n")
    return result[0]/cst.clight # [event^2/m]
    
  #============================================================================================================#
  
  ### Beam-beam parameter

  def GetBBParam(self, beam, ip):
    """ Returns an array with the horizontal, vertical and mean beam-beam parameters (approximation neglecting hour glass effect), for a given IP."""
    sigs   = beam._sigs
    sigs2  = beam._sigs**2
    bet    = beam._beta[ip]
    xplane = beam._xplane[ip]
    phi    = beam._phi[ip]
    betax  = bet[0]
    betay  = bet[1]
    sig    = beam.getsigma(ip)
    sigx   = sig[0]
    sigy   = sig[1]
    sigx2  = sigx**2
    sigy2  = sigy**2
    oncc   = beam._oncc[ip]
    if oncc < 1.0:
      if (xplane == 0):
        geomx = np.sqrt( 1. + (sigs2)/(sigx2) * np.tan( (1.0-oncc)*phi/2. )**2 )
        geomy = 1.0
      else:
        geomx = 1.0
        geomy = np.sqrt( 1. + (sigs2)/(sigy2) * np.tan( (1.0-oncc)*phi/2. )**2 )
    else:
      geomx = 1.0
      geomy = 1.0
    xi = []
    xi.append( beam._ppb*cst.r0*betax/ (2*np.pi * beam._gamma * sigx * geomx * (sigx*geomx + sigy*geomy) ) )
    xi.append( beam._ppb*cst.r0*betay/ (2*np.pi * beam._gamma * sigy * geomy * (sigx*geomx + sigy*geomy) ) )
    xi.append(0.5*(xi[0]+xi[1]))
    return xi

  #============================================================================================================#
  
  ### Levelling functions
  
  def GetBetaMin(self, beam, ip):
    """ Returns the minimum betas at a given levelling step, for a given IP. They might be a general minimum (for that IP), or an specific dictionary of betas as function of the ppb. """
    result = beam._betamin[ip]
    if (beam._betamindict != {} and ip < 2):
      currentint = beam._ppb
      print "beam._betamindict['intensities'] =", beam._betamindict["intensities"]
      print "beam._betamindict['betas'] =",       beam._betamindict["betas"]
      for ints, betas in zip(beam._betamindict["intensities"], beam._betamindict["betas"]):
        if (currentint < ints):
          print "currentint =", currentint, "is less than ints =", ints, "then finalminbetas =", betas
          result = betas
        else:
          print "currentint =", currentint, "is greater than  or equal to ints =", ints
    return result

  def FuncVar(self, lambdavar, beam, ip, levtech, levvar, roundorflat = None, kbeta = None):
    """ Assigns, for a given IP, the previously created lambda variable to the beta in the corresponding plane(s) or phiCk (depending on the request). In the case of beta, it re-assign a new value to the beta in the opposite plane so that the requested relation between them remains fixed. Then, if requested, it updates the angles. In both cases (betas or phiCK) it returns the luminosity or the peak pile up density (depending on the request) computed for these values of the variable (betas or phiCK) at the given IP. """
    if levvar == "Beta":
      xplane = beam._xplane[ip]
      if roundorflat < 1e-3:
        beam._beta[ip] = np.array([lambdavar, lambdavar])
      else:
        if self._constbetar[ip] == True:
          if xplane == 0:
            beam._beta[ip][1] = lambdavar
            beam._beta[ip][0] = lambdavar/kbeta
          else:   # xplane == 1
            beam._beta[ip][0] = lambdavar
            beam._beta[ip][1] = lambdavar*kbeta
        else:
          beam._beta[ip][1-xplane] = lambdavar
      if self._sepLRconst[ip] == True or beam._redsepLR[ip] == True:
        self.GetNewPhi(beam, ip, output = False)
    elif levvar == "PhiCR":
      beam._phiCR[ip] = lambdavar
      self.GetOnccDueToPhiCR(beam, ip)
    elif levvar == "PhiCK":
      beam._phiCK[ip] = lambdavar
      self.GetOnckDueToPhiCK(beam, ip)
    elif levvar == "ParSep":
      beam._parsep = lambdavar
    lumitot = self.getlumiip(beam, ip)
    if levtech == "levlumi":
      result  = lumitot[0]
    elif levtech == "levppus":
      pu     = self.GetPileUp(beam._frev, lumitot[0], beam._nbunch[ip])
      ppus   = pu*self.GetPeakLumRegion(beam, ip)/1e3
      result = ppus
    return result

  def GetConditionsLevel(self, ip, levtech, levvar, lumi, levlumi, beta, bmin, initphi, phiCR, initphiCR, phiCK, ppus, levppus, parsep):
    """ Returns the conditions that either the betas, phiCK, or parsep (depending on request) for a given IP have to fulfill in order to perform or not a levelling step. """
    if ip == 2:
      levlumitol = 0.995 ########################## CHECK 0.9975 before, and for all the Mar18IP8 cases - 0.995 was used for the case mar18IP8 2*sqrt(160**2+135**2) V, 2e34
    else:
      levlumitol = 0.999
      levppustol = 0.999
    condition1 = lumi/levlumi < levlumitol
    
    if levtech == "levlumi":
      condition2 = True
    elif levtech == "levppus":
      condition2 = ppus/levppus < levppustol  # It is expected to have a value != None
      
    if levvar == "Beta":
      condition3 = beta == bmin
      condition4 = beta >  bmin
      condition5 = True
    elif levvar == "PhiCR":
      condition3 = phiCR == initphi
      condition4 = phiCR >= 0.0
    elif levvar == "PhiCK":
      condition3 = phiCK == 0.0
      condition4 = phiCK >= 0.0
    elif levvar == "ParSep":
      condition3 = parsep == 0.0
      condition4 = parsep > 0.0
    else:
      bannersclass = Banners()
      bannersclass.header0("ERROR: VARIABLE NAMES\"")
      sys.exit()
      
    return condition1, condition2, condition3, condition4

  def GetLevel(self, beam, ip, levtech, levvar, t):
    """ Finds by an optimization process (and updates in the beam class), the betas or phiCK (depending on the request) for a given IP such that the desired levelled luminosity is reached or the levelled peak pile-up density is respected (depending on the request). If the minimum betas or a zero phiCK are/is reached, and the levelled value requested has not been reached, these minimum betas or zero phiCK are/is used for the updates in the beam class. """
    
    beta        = beam._beta[ip]
    bmin        = self.GetBetaMin(beam, ip)
    kbeta       = bmin[1]/bmin[0]
    roundorflat = abs(bmin[0] - bmin[1])   # Used to check for round beam
    xplane      = beam._xplane[ip]
    initphi     = beam._initphi[ip]
    initsepLR   = beam._initsepLR[ip]
    if beam._redsepLR[ip] == True:
      #print "> sepLR is the independent variable (to be reduced).\n Angles in radians are updated during the process..."
      self.GetNewPhi(beam, ip, output = False)
    else:
      self.GetNewSepLR(beam, ip)
    phi         = beam._phi[ip]
    sepLR       = beam._sepLR[ip]
    phiCR       = beam._phiCR[ip]
    oncc        = beam._oncc[ip]
    initphiCR   = beam._initphiCR[ip]
    phiCK       = beam._phiCK[ip]
    onck        = beam._onck[ip]
    if levvar == "ParSep":
      parsep    = beam._parsep
    else:
      parsep    = None
    levlumi     = self._levlumi[ip]
    lumitot     = self.getlumiip(beam, ip)
    plregion    = self.GetPeakLumRegion(beam, ip)
    levppus     = self._levppus[ip]
    pu          = self.GetPileUp(beam._frev, lumitot[0], beam._nbunch[ip])
    ppus        = pu*plregion/1e3
    
    print "*", "{0:36} {1:10}".format("Levelling technique ", "-"),             levtech
    print "*", "{0:36} {1:10}".format("Levelling variable ", "-"),              levvar
    print ""
    # 
    print "*", "{0:36} {1:10}".format("Beta horizontal", "m"),                  "%1.8f" %beta[0]
    print "*", "{0:36} {1:10}".format("Beta vertical",   "m"),                  "%1.8f" %beta[1]
    print "*", "{0:36} {1:10}".format("Beta vertical/horizontal ratio", "1"),   "%1.8f" %(beta[1]/beta[0])
    print "*", "{0:36} {1:10}".format("Beta minimum horizontal", "m"),          "%1.8f" %bmin[0]
    print "*", "{0:36} {1:10}".format("Beta minimum vertical",   "m"),          "%1.8f" %bmin[1]
    print "*", "{0:36} {1:10}".format("Constant for beta ratio", "1"),          "%1.8f" %kbeta
    print "*", "{0:36} {1:10}".format("Round beam parameter", "m"),             "%1.8f" %roundorflat
    print ""
    #
    print "*", "{0:36} {1:10}".format("Crossing angle", "rad"),                 "%1.6e" %phi
    print "*", "{0:36} {1:10}".format("BB-LR separation", "sigma"),             "%1.8f" %sepLR
    print "*", "{0:36} {1:10}".format("Crab cavity angle", "rad"),              "%1.6e" %phiCR
    print "*", "{0:36} {1:10}".format("Crab cavity ON (oncc)", "1"),            "%1.8f" %oncc
    print "*", "{0:36} {1:10}".format("Crab kissing angle", "rad"),             "%1.6e" %phiCK
    print "*", "{0:36} {1:10}".format("Crab kissing ON (onck)", "1"),           "%1.8f" %onck
    print "*", "{0:36} {1:10}".format("Parallel separation", "m"),
    if levvar == "ParSep":
      print "%1.6e" %parsep
    else:
      print parsep
    print ""
    #
    print "*", "{0:36} {1:10}".format("Leveled luminosity", "cm-2 s-1"),        "%1.6e" %levlumi
    print "*", "{0:36} {1:10}".format("Luminosity", "cm-2 s-1"),                "%1.6e" %lumitot[0]
    print "*", "{0:36} {1:10}".format("Ratio of luminosities", "1"),            "%1.8f" %(lumitot[0]/levlumi)
    #
    print "*", "{0:36} {1:10}".format("Leveled peak pile-up s-density", "mm-1"),
    if levppus != None: print "%1.8f" %levppus
    else:               print levppus
    print "*", "{0:36} {1:10}".format("Peak pile-up s-density", "mm-1"),         "%1.8f" %ppus
    print "*", "{0:36} {1:10}".format("Ratio of peak pile-up s-densities", "1"),
    if levppus != None: print "%1.8f" %(ppus/levppus)
    else:               print None
    print "*", "{0:36} {1:10}".format("Pile-up", "1"),                          "%1.8f" %pu
    
    # Set the conditions:
    condition1, condition2, condition3, condition4 = self.GetConditionsLevel(ip, levtech, levvar, lumitot[0], levlumi, beta[1-xplane], bmin[1-xplane], initphi, phiCR, initphiCR, phiCK, ppus, levppus, parsep)
    
    # Check the conditions. If not satisfied, no levelling can do anymore...
    if (condition1 and condition2 and condition3):
      print ""
      print 1*" ", "[!!] Current lumi (" + str("%1.6e" %lumitot[0]) + ") is lower than leveled luminosity,"
      if levtech == "levppus":
        print 6*" ", "and peak pile-up s-density has reached the limit (" + str("%1.6f" %ppus) + "),"
      if levvar == "Beta":
        print 6*" ", "but beta (in the opposite plane) has gone lower than beta minimum (in the opposite plane)!"
      elif levvar == "PhiCR":
        print 6*" ", "but phiCR has reached phi!"
        self._niterlev = 1
      elif levvar == "PhiCK":
        print 6*" ", "but phiCK has reached zero!"
      elif levvar == "ParSep":
        print 6*" ", "but parsep has reached zero!"
      print 6*" ", "I can't do anything about it anymore...\n"
    
    else:
      # Set boudaries for the optimization
      if levvar == "Beta":
        boundmin = beam._betamin[ip][1-xplane]
        if t == 0.0: boundmax = 1e2*boundmin
        else:        boundmax = 1.5*beta[1-xplane]  # boundmax = beta[1-xplane]  should work if the betas are always decreasing (it might help converge faster). A 50% margin has been added.
      elif levvar == "PhiCR":
        boundmin = 0.0
        boundmax = phi
      elif levvar == "PhiCK":
        boundmin = 0.0
        if t == 0.0: boundmax = 2.5*initphi# 2.5*phi  # Take 250% of the crossing angle for the optimization at the first step
        else:        boundmax = 1.5*initphi # 1.5*beam._phiCK[ip]  # boundmax = beam._phiCK[ip] should work if the phiCK is always decreasing (it might help converge faster). A 50% margin has been added.
      elif levvar == "ParSep":
        boundmin = 0.0
        boundmax = 1.0
      #print "\nppus levelling in this step?"
      #print ppus
      #print levppus
      #print self._flaglevpu, "(", self._flaglevpu[ip], ")"
      #print t == 0.0
      #print self._flaglevpu[ip] == False or t == 0.0
      #print ""
      
      # Reduce betas/phiCR/phiCK while they/it are/is larger than the minimum values/zero. t = 0 is meant to be for the first levelling step
      if (condition4 or t == 0.0):
        
        if levvar == "Beta":
          print "*", "{0:36} {1:10}".format("Bound for optimization minimum", "m"), "%1.8f" %boundmin
          print "*", "{0:36} {1:10}".format("Bound for optimization maximum", "m"), "%1.8f" %boundmax
        elif levvar == "PhiCR" or levvar == "PhiCK":
          print "*", "{0:36} {1:10}".format("Bound for optimization minimum", "rad"), "%1.6e" %boundmin
          print "*", "{0:36} {1:10}".format("Bound for optimization maximum", "rad"), "%1.6e" %boundmax
        # no parsep
          
        if levtech == "levlumi":
          
          if levvar == "Beta":
            ### Luminosity levelling with betas
            print "\n> Luminosity levelling with betas..."
            F = lambda betavar: (self.FuncVar(betavar, beam, ip, "levlumi", levvar, roundorflat, kbeta) - levlumi)**2
          elif levvar == "PhiCR":
            ### Luminosity levelling with crab crossing
            print "\n> Luminosity levelling with crab crossing..."
            F = lambda phiCRvar: (self.FuncVar(phiCRvar, beam, ip, "levlumi", levvar) - levlumi)**2
          elif levvar == "PhiCK":
            ### Luminosity levelling with crab kissing
            print "\n> Luminosity levelling with crab kissing..."
            F = lambda phiCKvar: (self.FuncVar(phiCKvar, beam, ip, "levlumi", levvar) - levlumi)**2
          elif levvar == "ParSep":
            print "\n> Luminosity levelling with parallel separation..."
            if ip == 2:
              sigp = beam.getsigma(ip)[1]
              #print sigp
              sepfact = np.exp(-beam._parseptmp**2/4./sigp**2)
              lumitottmp = self.getlumiip(beam, ip)[0]/sepfact
              #print beam._parseptmp, sepfact, lumitottmp
              if lumitottmp >= levlumi:
                beam._parsep = 2.*sigp*np.sqrt( np.log(lumitottmp/levlumi) )
              else:
                beam._parsep = 0.
                if self._flagsteplev2 == False:
                  self._flagsteplev2 = True
              sepfact = np.exp(-beam._parsep**2/4./sigp**2)
              lumitot = self.getlumiip(beam, ip)
              #print beam._parsep, sepfact, lumitot
              beam._parseptmp = beam._parsep
            else:
              sys.error("\nTHIS IS ONLY IMPELMENTED FOR IP8\n")
            
        elif levtech == "levppus":
        
          if levvar == "Beta":
            if (self._flaglevpu[ip] == False or t == 0.0):
              ### Pile-up levelling with betas, but we use luminosity levelling since we are below the levelled peak pile-up
              print "\n> Peak pile-up s-density levelling with betas,\n  but we use luminosity levelling since the levelled peak pile-up has not been reached..."
              F = lambda betavar: (self.FuncVar(betavar, beam, ip, "levlumi", levvar, roundorflat, kbeta) - levlumi)**2
            else:
              ### Pile-up levelling with crab kissing, but we use luminosity levelling since we are below the levelled peak pile-up
              print "\n> Peak pile-up s-density levelling with betas..."
              F = lambda betavar: (self.FuncVar(betavar, beam, ip, "levppus", levvar, roundorflat, kbeta) - levppus)**2
            
          elif levvar == "PhiCR":
            if (self._flaglevpu[ip] == False or t == 0.0):
              ### Pile-up levelling with crab crossing
              print "\n> Peak pile-up s-density levelling with crab,\n  but we use luminosity levelling since the levelled peak pile-up has not been reached..."
              F = lambda phiCRvar: (self.FuncVar(phiCKvar, beam, ip, "levlumi", levvar) - levlumi)**2
            else:
              ### Pile-up levelling with crab crossing
              print "\n> Peak pile-up s-density levelling with crab..."
              F = lambda phiCRvar: (self.FuncVar(phiCKvar, beam, ip, "levppus", levvar) - levppus)**2
            
          elif levvar == "PhiCK":
            if (self._flaglevpu[ip] == False or t == 0.0):
              ### Pile-up levelling with crab kissing
              print "\n> Peak pile-up s-density levelling with crab kissing,\n  but we use luminosity levelling since the levelled peak pile-up has not been reached..."
              F = lambda phiCKvar: (self.FuncVar(phiCKvar, beam, ip, "levlumi", levvar) - levlumi)**2
            else:
              ### Pile-up levelling with crab kissing
              print "\n> Peak pile-up s-density levelling with crab kissing..."
              F = lambda phiCKvar: (self.FuncVar(phiCKvar, beam, ip, "levppus", levvar) - levppus)**2
          
          elif levvar == "ParSep":
            sys.exit("\n[!] PEAK PILE-UP S-DENSITY WITH PARALLEL SEPARATION IS NOT SUPPORTED\n")
              
        #optimization = minimize_scalar(F, bounds = (boundmin, boundmax), method = 'bounded')  # For the scipy version in my desktop PC. Use tol instead of xtol
        if levvar != "ParSep": optimization = fminbound(F, boundmin, boundmax, full_output=1, disp=3, xtol=1e-8)      # For the scipy version in LXPLUS. xtol = 1e-5 is default
        
        # Update conditions
        beta     = beam._beta[ip]
        bmin     = self.GetBetaMin(beam, ip)
        xplane   = beam._xplane[ip]
        if beam._redsepLR[ip] == True: self.GetNewPhi(beam, ip, output = False)
        else:                          self.GetNewSepLR(beam, ip)
        phi      = beam._phi[ip]
        sepLR    = beam._sepLR[ip]
        phiCR    = beam._phiCR[ip]
        oncc     = self.GetOnccDueToPhiCR(beam, ip)
        phiCK    = beam._phiCK[ip]
        onck     = self.GetOnckDueToPhiCK(beam, ip)
        if levvar == "ParSep":
          parsep    = beam._parsep
        else:
          parsep    = None
        lumitot  = self.getlumiip(beam, ip)
        plregion = self.GetPeakLumRegion(beam, ip)
        pu       = self.GetPileUp(beam._frev,lumitot[0], beam._nbunch[ip])
        ppus     = pu*plregion/1e3
        
        if levtech == "levppus":
          if (ppus > levppus and self._flaglevpu[ip] == False):
            self._flaglevpu[ip] = True
            if self._timeppuslev == None:
              self._timeppuslev = t/3600.
          
        condition1, condition2, condition3, condition4 = self.GetConditionsLevel(ip, levtech, levvar, lumitot[0], levlumi, beta[1-xplane], bmin[1-xplane], initphi, phiCR, initphiCR, phiCK, ppus, levppus, parsep)
        
        #print condition1, 'and', condition2, 'and not', condition3, '->',  condition1 and condition2 and not(condition3)
        
        # Check the conditions. If not satisfied, the betas or phiCK after optimization, then set them/it to the minimum values / zero. 
        if (condition1 and condition2 and not(condition3)):
          if levvar == "Beta":
            beam._beta[ip] = bmin
          elif levvar == "PhiCR":
            beam._phiCR[ip] = phi
          elif levvar == "PhiCK":
            beam._phiCK[ip] = 0.0
          elif levvar == "ParSep":
            beam._parsep = 0.0
          print ""
          print 3*" ", "[!] Current luminosity (" + str("%1.6e" %lumitot[0]) + ") is lower than the leveled luminosity!"
          if levtech == "levppus":
            print 7*" ", "Current peak pile-up s-density (" + str("%1.6f" %ppus) + ") is lower than the leveled pile-up s-density!"
          if levvar == "Beta": 
            print 7*" ", "Betas are reseted to beta minimum: [" + str("%1.8f" %beam._beta[ip][0]) + ",", str("%1.8f" %beam._beta[ip][1]) + "]"
          elif levvar == "PhiCR":
            print 7*" ", "New phiCR is reseted phi:", str("%1.4e" %beam._phi[ip]) + "rad"
          elif levvar == "PhiCK":
            print 7*" ", "New phiCK is reseted to zero:", str("%1.4e" %beam._phiCK[ip]) + "rad"
          elif levvar== "ParSep":
            print 7*" ", "New parsep is reseted to zero:", str("%1.4e" %beam._parsep) + "rad"
          print ""
    
    return lumitot
  
  #============================================================================================================#

  ### Update parameters due to burn-off, IBS and radiation damping and bunch length gymnastics. 
  
  def GetPPBNew(self, t, lumiall, beam):
    """ Updates the ppb in the beam class after a time t due to burn off. """
    rate = 0.0
    for i in range(beam._nip):
      rate = rate + lumiall[i]*(self._xsecburn*1.0e-27) / beam._nbunch[i]
    beam._ppb = beam._ppb / (1.0 + t * rate/beam._ppb)  # New form, vs. old: beam._ppb = beam._ppb - t*rate
    print  "*", "{0:36} {1:10}".format("New ppb (ppb)", "1"), "%1.6e" %beam._ppb
  
  def GetPhiFromPPBFor6sigmaDA(self, beam, ip):
    """ Updates phiCR with the interpo;ates value between the corresponding pair (ppb, phiCR) to ensure 6sigma DA (Nikos). To be used for IP1 and IP5 """
    print 'IP =', ip
    phiold   = beam._phi[ip]
    sepLRold = beam._sepLR[ip]
    onccold  = beam._oncc[ip]
    phiCRold = beam._phiCR[ip]
    if   beam._ppb > 2.20e11 : # (2.20,  inf)
      phi = 308e-6
    elif beam._ppb > 1.90e11 : # (1.90, 2.20]
      phi = (366e-6 - 308e-6)/(1.90e11 - 2.20e11) * (beam._ppb - 2.20e11)  + 308e-6
    elif beam._ppb > 1.19e11 : # (1.19, 2.20]
      phi = (470e-6 - 366e-6)/(1.19e11 - 1.90e11) * (beam._ppb - 1.90e11)  + 366e-6
    elif beam._ppb > 1.10e11 : # (1.10, 2.20]
      phi = (464e-6 - 470e-6)/(1.10e11 - 1.19e11) * (beam._ppb - 1.19e11)  + 470e-6
    elif beam._ppb > 0.90e11 : # (0.90, 2.20]
      phi = (432e-6 - 464e-6)/(0.90e11 - 1.10e11) * (beam._ppb - 1.10e11)  + 464e-6
    elif beam._ppb > 0.80e11 : # (0.80, 2.20]
      phi = (418e-6 - 418e-6)/(0.80e11 - 0.90e11) * (beam._ppb - 0.90e11)  + 418e-6
    else:                      # (0.0,  0.80]
      phi = 418e-6
    beam._phi[ip] = phi
    self.GetNewSepLR(beam, ip)
    self.GetPhiCRDueToAdaptive(beam, ip)
    self.GetOnccDueToPhiCR(beam,ip)
    print "*", "{0:36} {1:10}".format("Old crossing angle", "rad"),            "%1.4e" %phiold
    print "*", "{0:36} {1:10}".format("New crossing angle", "rad"),            "%1.4e" %beam._phi[ip]
    print "*", "{0:36} {1:10}".format("Old BB-LR separation", "sigma"),        "%1.8f" %sepLRold
    print "*", "{0:36} {1:10}".format("New BB-LR separation", "sigma"),        "%1.8f" %beam._sepLR[ip]
    print "*", "{0:36} {1:10}".format("Old crab cavity angle", "rad"),         "%1.4e" %phiCRold
    print "*", "{0:36} {1:10}".format("New crab cavity angle", "rad"),         "%1.4e" %beam._phiCR[ip]
    print "*", "{0:36} {1:10}".format("Old crab cavity ON (oncc)", "1"),       "%1.8f" %onccold
    print "*", "{0:36} {1:10}".format("New crab cavity ON (oncc)", "1"),       "%1.8f" %beam._oncc[ip]
    print "\n---"
    
  def IBStauxz(self, t, beam):
    """ Returns an array of the three IBS damping times (all in hours). The taux and tauz are obtained 'from 'ibs' module of MAD-X, run on the specified lattice; tauy is obtained from (1/kappa_x)*taux. The old approximation formula
    'tauxtmp = ( (16 * beam._tuneb * (beam._epsxn*beam._epsyn) * np.sqrt(beam._kappa) * np.sqrt(beam._kappa+1.)* beam._gamma * beam._sigs*beam._dpp ) / ( 2.*cst.clight * (cst.r0**2) * beam._ppb * 23 ) ) / 3600.' is deprecated. """
    NAMEupdateIBS = "updateIBS_" + self._ID + ".madx"
    NAMEtauxz     = "tauxz_"     + self._ID + ".out"
    NAMEjobfile   = "jobtau_"  + self._ID + ".madx"
    FILEupdateIBS = open(NAMEupdateIBS, "a")
    print >> FILEupdateIBS, "!", t
    print >> FILEupdateIBS, "tuneb     =", beam._tuneb, ";"
    print >> FILEupdateIBS, "epsxn     =", beam._epsn[0], ";"
    print >> FILEupdateIBS, "epsyn     =", beam._epsn[1], ";"
    print >> FILEupdateIBS, "kappa     =", beam._kappa, ";"
    print >> FILEupdateIBS, "gamma     =", beam._gamma, ";"
    print >> FILEupdateIBS, "sigs      =", beam._rmssigs, ";" # beam._sigs, ";"
    #if self._longdens == "qGaussian": print >> FILEupdateIBS, "*", cst.FactFWHMGauss, "/4./sqrt(2.)", # before, now it has changed in the previous line for sigs -> sigsrms
    #print >> FILEupdateIBS, ";"
    print >> FILEupdateIBS, "dpp       =", beam._dpp,   ";"
    print >> FILEupdateIBS, "ppb       =", beam._ppb, ";"
    print >> FILEupdateIBS, "NRJ       =", beam._momeV/1e9, ";"
    FILEupdateIBS.close()
    # Make a copy of the template script to run the ibs module, and edit it handle input/output files with a unique name given by an ID
    
    if beam._opticsfile == "HELHC-inj":
      copyfile("./Levelling_jobtau_helhc.madx", NAMEjobfile)

    else:
       copyfile("./Levelling_jobtau.madx", NAMEjobfile) # newest hl-lhc sequence by Stefania
      #copyfile("./Levelling_jobtau-old.madx", NAMEjobfile) # to use old hl-lhc sequence
      
    for LINE in fileinput.input(NAMEjobfile, inplace = 1):
      print LINE,
      if LINE.startswith("//call, file = \"slhc/opt.madx\";"):  # not needed for HE-LHC
        print "call, file = \"slhc/" + beam._opticsfile + "\";"
      if LINE.startswith("//call, file = \"updateIBS.madx\";"):
        print "call, file = \"" + NAMEupdateIBS + "\";"
      if LINE.startswith("//assign, echo = \"tauxz.out\";"):
        print "assign, echo = \"" + NAMEtauxz + "\";"
    # Run the MAD-X script and wait until it finishes
    cmd = ["/afs/cern.ch/user/m/mad/bin/madx", NAMEjobfile]
    process = subprocess.Popen(cmd, stdout = subprocess.PIPE)
    process.wait()
    # File where the results for taux and tauz are saved by MAD-X
    FILEtaux = open(NAMEtauxz, "r")  
    line = FILEtaux.readline()
    data   = []
    while line:
      values = line.split()
      data.append(values)
      line = FILEtaux.readline()
    FILEtaux.close()
    try:            os.remove(NAMEtauxz)
    except OSError: pass
    beam._tau_ibs[0] = float(data[0][2])
    beam._tau_ibs[2] = float(data[1][2])
    beam._tau_ibs[1] = (1.0/beam._kappac)*beam._tau_ibs[0]

  def EmitGrowth(self, t, beam):
    """ Updates, after a time t, the values of the emittances, bunch length and energy spread due to IBS, synchrotron radiation damping, an "observed growth" constant factor, and inducded by CC (for crossing and kissing). Contributions from IBS and SR are optional"""
    ## Crossing schemes in IP1&5 (Formally)
    # if   beam._xplane[0] == 0 and beam._xplane[1] == 0: # HH
    #   note_0 = ["CR at IP1 hor", "CK at IP1 ver"]
    #   note_1 = ["CR at IP5 hor", "CK at IP5 ver"]
    #   if beam._phiCR[0] > 0.0: taux0 = beam._tau_cc / (beam._phiCR[0]**2/beam._beta[0][0]) # CR at IP1 hor *
    #   else:                    taux0 = 1e8
    #   if beam._phiCK[0] > 0.0: tauy0 = beam._tau_cc / (beam._phiCK[0]**2/beam._beta[0][1]) # CK at IP1 ver
    #   else:                    tauy0 = 1e8
    #   if beam._phiCR[1] > 0.0: taux1 = beam._tau_cc / (beam._phiCR[1]**2/beam._beta[1][0]) # CR at IP5 hor *
    #   else:                    taux1 = 1e8
    #   if beam._phiCK[1] > 0.0: tauy1 = beam._tau_cc / (beam._phiCK[1]**2/beam._beta[1][1]) # CK at IP5 ver
    #   else:                    tauy1 = 1e8
    # elif beam._xplane[0] == 0 and beam._xplane[1] == 1: # HV
    #   note_0 = ["CR at IP1 hor", "CK at IP1 ver"]
    #   note_1 = ["CK at IP5 hor", "CR at IP5 ver"]
    #   if beam._phiCR[0] > 0.0: taux0 = beam._tau_cc / (beam._phiCR[0]**2/beam._beta[0][0]) # CR at IP1 hor *
    #   else:                    taux0 = 1e8
    #   if beam._phiCK[0] > 0.0: tauy0 = beam._tau_cc / (beam._phiCK[0]**2/beam._beta[0][1]) # CK at IP1 ver
    #   else:                    tauy0 = 1e8
    #   if beam._phiCK[1] > 0.0: taux1 = beam._tau_cc / (beam._phiCK[1]**2/beam._beta[1][0]) # CK at IP5 hor
    #   else:                    taux1 = 1e8
    #   if beam._phiCR[1] > 0.0: tauy1 = beam._tau_cc / (beam._phiCR[1]**2/beam._beta[1][1]) # CR at IP5 ver *
    #   else:                    tauy1 = 1e8
    # elif beam._xplane[0] == 1 and beam._xplane[1] == 0: # VH
      # note_0 = ["CK at IP1 hor", "CR at IP1 ver"]
      # note_1 = ["CR at IP5 hor", "CK at IP5 ver"]
      # if beam._phiCK[0] > 0.0: taux0 = beam._tau_cc / (beam._phiCK[0]**2/beam._beta[0][0]) # CK at IP1 hor
      # else:                    taux0 = 1e8
      # if beam._phiCR[0] > 0.0: tauy0 = beam._tau_cc / (beam._phiCR[0]**2/beam._beta[0][1]) # CR at IP1 ver *
      # else:                    tauy0 = 1e8
      # if beam._phiCR[1] > 0.0: taux1 = beam._tau_cc / (beam._phiCR[1]**2/beam._beta[1][0]) # CR at IP5 hor *
      # else:                    taux1 = 1e8
      # if beam._phiCK[1] > 0.0: tauy1 = beam._tau_cc / (beam._phiCK[1]**2/beam._beta[1][1]) # CK at IP5 ver
      # else:                    tauy1 = 1e8
    # elif beam._xplane[0] == 1 and beam._xplane[1] == 1: # VV
    #   note_0 = ["CK at IP1 hor", "CR at IP1 ver"]
    #   note_1 = ["CK at IP5 hor", "CR at IP5 ver"]
    #   if beam._phiCK[0] > 0.0: taux0 = beam._tau_cc / (beam._phiCK[0]**2/beam._beta[0][0]) # CK at IP1 hor
    #   else:                    taux0 = 1e8
    #   if beam._phiCR[0] > 0.0: tauy0 = beam._tau_cc / (beam._phiCR[0]**2/beam._beta[0][1]) # CR at IP1 ver *
    #   else:                    tauy0 = 1e8
    #   if beam._phiCK[1] > 0.0: taux1 = beam._tau_cc / (beam._phiCK[1]**2/beam._beta[1][0]) # CK at IP5 hor
    #   else:                    taux1 = 1e8
    #   if beam._phiCR[1] > 0.0: tauy1 = beam._tau_cc / (beam._phiCR[1]**2/beam._beta[1][1]) # CR at IP5 ver *
    #   else:                    tauy1 = 1e8
    # else:
    #   sys.exit("\n[!] Crossing scheme not valid! Use HH, HH, HV or VH.\n")
    
    ### Approximation for VH
    note_0 = ["CK at IP1 hor", "CR at IP1 ver"]
    note_1 = ["CR at IP5 hor", "CK at IP5 ver"]
    if beam._phiCK[0] > 0.0: taux0 = beam._tau_cc / (beam._phiCK[0]**2/beam._beta[0][1]) # CK at IP1 hor - small beta
    else:                    taux0 = 1e8
    if beam._phiCR[0] > 0.0: tauy0 = beam._tau_cc / (beam._phiCR[0]**2/beam._beta[0][0]) # CR at IP1 ver * - big beta
    else:                    tauy0 = 1e8
    if beam._phiCR[1] > 0.0: taux1 = beam._tau_cc / (beam._phiCR[1]**2/beam._beta[1][0]) # CR at IP5 hor * - big beta
    else:                    taux1 = 1e8
    if beam._phiCK[1] > 0.0: tauy1 = beam._tau_cc / (beam._phiCK[1]**2/beam._beta[1][1]) # CK at IP5 ver - small beta
    else:                    tauy1 = 1e8
    
    # ### Approximation for HV
    # note_0 = ["CR at IP1 hor", "CK at IP1 ver"]
    # note_1 = ["CK at IP5 hor", "CR at IP5 ver"]
    # if beam._phiCR[0] > 0.0: taux0 = beam._tau_cc / (beam._phiCR[0]**2/beam._beta[0][0]) # CR at IP1 hor *
    # else:                    taux0 = 1e8
    # if beam._phiCK[0] > 0.0: tauy0 = beam._tau_cc / (beam._phiCK[0]**2/beam._beta[0][1]) # CK at IP1 ver
    # else:                    tauy0 = 1e8
    # if beam._phiCK[1] > 0.0: taux1 = beam._tau_cc / (beam._phiCK[1]**2/beam._beta[1][0]) # CK at IP5 hor
    # else:                    taux1 = 1e8
    # if beam._phiCR[1] > 0.0: tauy1 = beam._tau_cc / (beam._phiCR[1]**2/beam._beta[1][1]) # CR at IP5 ver *
    # else:                    tauy1 = 1e8
    
    beam._tau_cc_0 = [taux0, tauy0]
    beam._tau_cc_1 = [taux1, tauy1]
    steph = t/3600.0
    beam._epsn_0      = [beam._epsn[0], 
                         beam._epsn[1]]
    self.IBStauxz(t, beam)
    if self._IBSRF == True:
      beam._epsn_ibs  = [beam._epsn[0]*(steph/beam._tau_ibs[0]),
                         beam._epsn[1]*(steph/beam._tau_ibs[1])]
      beam._epsn_sr   = [-2*steph*(beam._epsn[0] - beam._eps0[0])/beam._tau_sr[0],
                         -2*steph*(beam._epsn[1] - beam._eps0[1])/beam._tau_sr[1]]
    else:
      beam._epsn_ibs  = [0.0,
                         0.0]
      beam._epsn_sr   = [0.0,
                         0.0]
    beam._epsn_obs    = [beam._epsn[0]*(steph/beam._tau_obs[0]),
                         beam._epsn[1]*(steph/beam._tau_obs[1])]
    beam._epsn_cc_ip0 = [beam._initepsn[0]*(steph/beam._tau_cc_0[0]),  # beam._epsn[0] for multiplicative, beam._initepsn[0] for additive
                         beam._initepsn[1]*(steph/beam._tau_cc_0[1])]  # beam._epsn[1] for multiplicative, beam._initepsn[1] for additive
    beam._epsn_cc_ip1 = [beam._initepsn[0]*(steph/beam._tau_cc_1[0]),  # beam._epsn[0] for multiplicative, beam._initepsn[0] for additive
                         beam._initepsn[1]*(steph/beam._tau_cc_1[1])]  # beam._epsn[1] for multiplicative, beam._initepsn[1] for additive
    print ""
    for i in range(2):
      print "EMITTANCE GROWTH TIME -",
      if i == 0:  print "Horizontal\n"
      else:       print "Vertical\n"
      print 3*" ", "{0:20} {1:>14} h or {2:>14}%/h".format("From IBS", "%1.4f" %beam._tau_ibs[i], "%1.4f" %(100.*1./beam._tau_ibs[i]))
      print 3*" ", "{0:20} {1:>14} h or {2:>14}%/h".format("From SR",  "%1.4f" %beam._tau_sr[i],  "%1.4f" %(100.*1./beam._tau_sr[i] ))
      print 3*" ", "{0:20}".format("Extra factor"),
      if beam._tau_obs[i] < 1e6: print "{0:>14} h or {1:>14}%/h".format("%1.4f" %beam._tau_obs[i], "%1.4f" %(100.*1./beam._tau_obs[i]))
      else:                      print "{0:>14} h or {1:>14}%/h".format("INF or None",             "0")
      print 3*" ", "{0:20}".format("From CR/CK at ip0"),
      if beam._tau_cc_0[i] < 1e6: print "{0:>14} h or {1:>14}%/h".format("%1.4f" %beam._tau_cc_0[i], "%1.4f" %(100.*1./beam._tau_cc_0[i])),
      else:                       print "{0:>14} h or {1:>14}%/h".format("INF or None",                 "0"),
      print note_0[i],
      if   note_0[i][0:2] == "CR":
        if beam._initoncc[i] == 0.0: print "--> No CR"
        else:                        print ""
      elif note_0[i][0:2] == "CK":
        if self._ck[i] == False:     print "--> No CK"
        else:                        print ""
      print 3*" ", "{0:20}".format("From CR/CK at ip1"),
      if beam._tau_cc_1[i] < 1e6: print "{0:>14} h or {1:>14}%/h".format("%1.4f" %beam._tau_cc_1[i], "%1.4f" %(100.*1./beam._tau_cc_1[i])),
      else:                       print "{0:>14} h or {1:>14}%/h".format("INF or None",                 "0"),
      print note_1[i],
      if   note_1[i][0:2] == "CR":
        if beam._initoncc[i] == 0.0: print "--> No CR"
        else:                        print ""
      elif note_1[i][0:2] == "CK":
        if self._ck[i] == False:     print "--> No CK"
        else:                        print ""
      print "\nEMITTANCE GROWTH CONTRIBUTIONS -",
      if i == 0:  print "Horizontal\n"
      else:       print "Vertical\n"
      print 3*" ", "{0:20} {1:>14}".format("Initial",           "%1.4e" %beam._initepsn[i]   )
      print 3*" ", "{0:20} {1:>14}".format("Previous step",     "%1.4e" %beam._epsn_0[i]     )
      print 3*" ", "{0:20} {1:>14}".format("From IBS",          "%1.4e" %beam._epsn_ibs[i]   )
      print 3*" ", "{0:20} {1:>14}".format("From SR",           "%1.4e" %beam._epsn_sr[i]    )
      print 3*" ", "{0:20} {1:>14}".format("Extra factor",      "%1.4e" %beam._epsn_obs[i]   )
      print 3*" ", "{0:20} {1:>14}".format("From CR/CK at ip0", "%1.4e" %beam._epsn_cc_ip0[i])
      print 3*" ", "{0:20} {1:>14}".format("From CR/CK at ip1", "%1.4e" %beam._epsn_cc_ip1[i])
      epsni = beam._epsn_0[i] + beam._epsn_ibs[i] + beam._epsn_sr[i] + beam._epsn_obs[i] + beam._epsn_cc_ip0[i] + beam._epsn_cc_ip1[i]
      beam._epsn[i] = epsni
      print ""
    #epsx1 = beam._epsxn * (1 + steph/taux + steph/beam._epsxfact + steph/(beam._epsxfactCC*beam._beta[0])) - 2*steph / beam._taux_sr * (beam._epsxn       - beam._epsx0)
    #epsy1 = beam._epsyn * (1 + steph/tauy + steph/beam._epsyfact + steph/(beam._epsyfactCC*beam._beta[1])) - 2*steph / beam._tauy_sr * (beam._epsyn       - beam._epsy0)
    if self._longdens == "qGaussian":
      sigsrms0 = beam._sigs0 * cst.FactFWHMGauss/4./np.sqrt(2.)
      sigsrms  = beam._sigs  * cst.FactFWHMGauss/4./np.sqrt(2.)
      epss     = sigsrms * beam._dpp  * (1 + steph/beam._tau_ibs[2]) -  2*steph / beam._tau_sr[2] * (sigsrms*beam._dpp - sigsrms0*beam._dpp0)
      sigsrms1 = np.sqrt( epss*sigsrms   /beam._dpp )
      dpp1     = np.sqrt( epss*beam._dpp/sigsrms    )
      beam._sigs  = sigsrms1 / cst.FactFWHMGauss*4.*np.sqrt(2.) 
      beam._dpp   = dpp1
    else:
      epss     = beam._sigs * beam._dpp * (1 + steph/beam._tau_ibs[2]) -  2*steph / beam._tau_sr[2] * (beam._sigs*beam._dpp - beam._sigs0*beam._dpp0)
      sigs1    = np.sqrt( epss*beam._sigs/beam._dpp )
      dpp1     = np.sqrt( epss*beam._dpp/beam._sigs )
      beam._sigs  = sigs1
      beam._dpp   = dpp1
    print "*", "{0:36} {1:10}".format("Normalized emittance horizontal", "m"), "%1.4e" %beam._epsn[0]
    print "*", "{0:36} {1:10}".format("Normalized emittance vertical",   "m"), "%1.4e" %beam._epsn[1]
    print "*", "{0:36} {1:10}".format("Bunch length (initial)",          "m"), "%1.4e" %beam._sigs
    print "*", "{0:36} {1:10}".format("RMS bunch length (initial)",      "m"), "%1.4e" %beam._rmssigs
    print "*", "{0:36} {1:10}".format("FWHM (initial)",                  "m"), "%1.4e" %beam._fwhm
    print "*", "{0:36} {1:10}".format("Energy spread (initial)",         "1"), "%1.4e" %beam._dpp
  
  def BunchLengthGymnastics(self, beam, typestep):
    """ Adjusts the beam length according to the specified settings. """
    if (typestep != "PENAL"):
      print "\nBUNCH LENGTH GYMNASTICS\n"
      if (beam._sigs > beam._minlong):
        if beam._ppb > beam._ppblong:
          if beam._constlong == True:
            print "> Bunch length is larger than the minimum, but we haven't reached the trigger ppb."
            print " ", "Constant bunch length was requested...\n"
            beam._sigs = beam._initsigs
            beam._dpp  = beam._initdpp
          else:
            print "> Bunch length is larger than the minimum, but we haven't reached the trigger ppb."
            print " ", "Constant bunch length was not requested...\n"
            beam._sigs = beam._sigs
            beam._dpp  = beam._dpp
        else:
          if beam._flagminlong == False:
            print "> Bunch length is larger than the minimum and we have passed the trigger ppb."
            print " ", "Shortening was requested...\n"
            print(beam._sigs)
            print(beam._dpp )
            beam._sigs = beam._sigs * (1. + beam._redlong)
            beam._dpp  = beam._dpp  * (1. + beam._redlong)
            print(beam._sigs)
            print(beam._dpp )
            if beam._sigs < beam._minlong:
              print 2*" ", "[!]", "Minimum bunch length has been reached.",
              reductsigs   = beam._minlong/beam._sigs
              beam._sigs   = beam._sigs * reductsigs
              beam._dpp    = beam._dpp  * reductsigs
              beam._mindpp = beam._dpp
              beam._flagminlong = True
              beam._constlong  = beam._constlong0
              print "We redefine mindpp =", "%1.2e" %beam._mindpp, "\n"
          else:
            if beam._constlong == True:
              print "> Bunch length is larger than the minimum and we have passed the trigger ppb and reached the minimum."
              print " ", "Constant bunch length was requested...\n"
              beam._sigs = beam._minlong
              beam._dpp  = beam._mindpp
            else:
              print "> Bunch length is larger than the minimum and we have passed the trigger ppb and reached the minimum."
              print " ", "Constant bunch length was not requested...\n"
              beam._sigs = beam._sigs
              beam._dpp  = beam._dpp
      else:
        if beam._flagminlong == False:
          print 2*" ", "[!]", "Minimum bunch length has been reached.",
          reductsigs   = beam._minlong/beam._sigs
          beam._sigs   = beam._sigs * reductsigs
          beam._dpp    = beam._dpp  * reductsigs
          beam._mindpp = beam._dpp
          beam._flagminlong = True
          beam._constlong  = beam._constlong0
          print "We redefine mindpp =", "%1.2e" %beam._mindpp, "\n"
        if beam._constlong == True:
          print "> Bunch length is equal or smaller than the minimum."
          print " ", "Constant bunch length was requested...\n"
          beam._sigs = beam._minlong
          beam._dpp  = beam._mindpp
        else:
          print "> Bunch length is equal or smaller than the minimum."
          print " ", "Constant bunch length was not requested...\n"
          beam._sigs = beam._sigs
          beam._dpp  = beam._dpp
      print "*", "{0:36} {1:10}".format("Bunch length",       "m"), "%1.4e" %beam._sigs
      print "*", "{0:36} {1:10}".format("RMS bunch length",   "m"), "%1.4e" %beam._rmssigs
      print "*", "{0:36} {1:10}".format("FWHM",               "m"), "%1.4e" %beam._fwhm
      print "*", "{0:36} {1:10}".format("Energy spread dp/p", "1"), "%1.4e" %beam._dpp
  
  def StepReductionSeparation(self, beam, time, ip, flag):
    """ Reduces sepLR by a given amount at a given time. """
    if (beam._sepLRsteptime[ip] != False):
      sepLRold = beam._sepLR[ip]
      phiold   = beam._phi[ip]
      onccold  = beam._oncc[ip]
      phiCRold = beam._phiCR[ip]
      if (time >= beam._sepLRsteptime[ip] and flag == False):
        beam._sepLR[ip] = beam._sepLR[ip] + beam._sepLRstep[ip]
        self.GetNewPhi(beam, ip)
        print "\n> Reduction of BB-LR separation for ip" + str(ip), "at t =", time, "(Time requested:", str(beam._sepLRsteptime[ip]) + ")..."
        print "*", "{0:36} {1:10}".format("Old crossing angle", "rad"),            "%1.4e" %phiold
        print "*", "{0:36} {1:10}".format("New crossing angle", "rad"),            "%1.4e" %beam._phi[ip]
        print "*", "{0:36} {1:10}".format("Old BB-LR separation", "sigma"),        "%1.8f" %sepLRold
        print "*", "{0:36} {1:10}".format("New BB-LR separation", "sigma"),        "%1.8f" %beam._sepLR[ip]
        print "*", "{0:36} {1:10}".format("Final crab cavity angle", "rad"),       "%1.4e" %phiCR
        print "*", "{0:36} {1:10}".format("Final crab cavity angle", "rad"),       "%1.4e" %beam._phiCR[ip]
        print "*", "{0:36} {1:10}".format("Old crab cavity ON (oncc)", "1"),       "%1.8f" %onccold
        print "*", "{0:36} {1:10}".format("New crab cavity ON (oncc)", "1"),       "%1.8f" %beam._oncc[ip]
        print ""
        flag = True
    return flag

  def ReductionSeparation(self, beam, time, ip):
    """ Reduces sepLR at a given rate down to a minimum value. """
    if beam._redsepLR[ip] == True and time != 0.0:
      if beam._sepLR[ip] > beam._minsepLR[ip]:
        print "> BB-LR separation is reduced to", 
        beam._sepLR[ip] = beam._sepLR[ip] - beam._ratesepLR[ip]*time/3600.0
        if beam._sepLR[ip] < beam._minsepLR[ip]:
          print "%1.4f" %beam._sepLR[ip], "sigma\n"
          print 2*" ", "[!]", "Minimum BB-LR separation has been reached."
          print 6*" ", "Reset to",
          beam._sepLR[ip] = beam._minsepLR[ip]
        print "%1.4f" %beam._sepLR[ip], "sigma after", "%1.2f" %time, "s.\n  Angles in radians adjusted accordingly.\n"
      else:
        print 2*" ", "[!]", "Minimum BB-LR separation has been reached.\n"
        beam._sepLR[ip] = beam._minsepLR[ip]
        #beam._redsepLR[ip] = False
      self.GetNewPhi(beam, ip, output = False)
  
  def IncrOncc(self, beam, j):
    """ Increases oncc linearly until reaching 1.0. """
    if beam._oncc[j] < 1.0:
      print "> oncc is increased from", "%1.6f" %beam._oncc[j], "to",
      beam._oncc[j] = beam._oncc[j] * (1 + self._incronccrate)
      print "%1.6f" %beam._oncc[j]
      print "> phiCR changes from", "%1.4e" %beam._phiCR[j], "to",
    else:
      beam._oncc[j] = 1.0
      print 2*" ", "[!] oncc has reached maximum. It is reset to", "%1.6f" %beam._oncc[j], "and phiCR =",
    self.GetPhiCR(beam, j)
    print "%1.4e" %beam._phiCR[j]
    print ""
  
  def ParametersEvolution(self, beam, typestep, timedecay, listLUMITOTi):
    """ Launches the update of the ppb, and the emittances/bunch length/energy spread, with their respective functions, after the appropiate time which depends of the type of step. """
    print "\nBURN-OFF AND INTRABEAM SCATTERING\n"
    if (typestep == "PENAL"):
      ### Penalty steps ALWAYS follow level or decay steps, therefore parameters always have to be updated after t = timedecay
      print "> ppb and emittance taus are updated after a previous decay of", "%1.2f" %timedecay, "s...\n"
      self.GetPPBNew(timedecay, listLUMITOTi[-1], beam)
      self.EmitGrowth(timedecay, beam)
    elif (typestep == "LEVEL"):
      if (self._penstp == True):
        if (self._updatepenstp == True):
          ### A level step following a penalty step (where parameters were updated), with a request to update the parameters again
          print "> ppb and emittance taus are updated after a penalty step of", "%1.2f" %self._timepenstp, "s...\n"
          self.GetPPBNew(self._timepenstp, listLUMITOTi[-2], beam)
          self.EmitGrowth(self._timepenstp, beam)
        else:
          ### A level step following a penalty step (where parameters were updated), where the parameters do not need to be updated again
          print "> ppb and emittance taus remain constant after a penalty step of", "%1.2f" %self._timepenstp, "s...\n"
          print  "*", "{0:36} {1:10}".format("Same ppb (ppb)", "1"), "%1.6e" %beam._ppb
      else:
        ### A level step with no previous penalty step, therefore parameters have to be updated always after t = timedecay
        print "> ppb and emittance taus are updated after a previous decay of", "%1.2f" %timedecay, "s...\n"
        self.GetPPBNew(timedecay, listLUMITOTi[-1], beam)
        self.EmitGrowth(timedecay, beam)
    elif (typestep == "DECAY"):
      if (self._penstp == True):
        ### A decay step following a penalty step (where parameters were updated). It is not needed to update, since the penalty step after the last step of levelling have all zero length, and they always update parameters
        print  "*", "{0:36} {1:10}".format("Same ppb (ppb)", "1"), "%1.6e" %beam._ppb
      else:
        ### A decay step with no previous step, therefore parameters have to be updated always after t = timedecay
        print "> ppb and emittance taus are updated after a previous decay of", "%1.2f" %timedecay, "s...\n"
        self.GetPPBNew(timedecay, listLUMITOTi[-1], beam)
        self.EmitGrowth(timedecay, beam)
    #else: This option does not exist,since typestep == "FIRST" is not possible
  
  #============================================================================================================#
  
  ### Print functions
  
  def PrintLumiParam(self, beam, tofile = False, foutname = None):
    """ Prints to the standard output or to file 'foutname' the Luminosity class parameters. """
    
    if tofile == False:
      
      print "LEVELLING PARAMETERS"
      print ""
      print "*", "{0:36} {1:10}".format("Levelling technique ", "-"),
      for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:16}".format(self._levtech[i]),
      print ""
      print "*", "{0:36} {1:10}".format("Levelling variable ", "-"),
      for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:16}".format(self._levvar[i]),
      print ""
      print "*", "{0:36} {1:10}".format("Leveled luminosity", "cm-2 s-1"),
      for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:16}".format("%1.6e" %self._levlumi[i]),
      print ""
      print "*", "{0:36} {1:10}".format("Leveled peak pile-up s-density", "mm-1"),
      for i in range(beam._nip):
        print beam._ipnames[i] + ":", 
        if self._levppus[i] != None: print "{0:16}".format("%1.4f" %self._levppus[i]),
        else:                        print self._levppus[i]
      print ""
      print "*", "{0:36} {1:10}".format("CK in the separation plane", "T/F"),
      for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:16}".format(str(self._ck[i])),
      print ""
      print "*", "{0:36} {1:10}".format("Constant beta H/V ratio", "T/F"),
      for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:16}".format(str(self._constbetar[i])),
      print ""
      print "*", "{0:36} {1:10}".format("Constant crossing angle in sigmas", "T/F"),
      for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:16}".format(str(self._sepLRconst[i])),
      print ""
      print "*", "{0:36} {1:10}".format("Flat longitudinal density", "G/RF800/q"),     self._longdens
      print ""
      #
      print "*", "{0:36} {1:10}".format("Luminosity step for levelling", "1"),       "%1.4f" %self._p
      print "*", "{0:36} {1:10}".format("Time step (in case p = 1.0)", "s"),         "%1.4f" %self._step
      print "*", "{0:36} {1:10}".format("Simulation for optimum fill", "T/F"),       self._optimumfill
      print "*", "{0:36} {1:10}".format("Maximum fill length (time)", "s"),          "%1.4f" %self._maxfill
      print "*", "{0:36} {1:10}".format("Penalty step", "T/F"),                       self._penstp
      print "*", "{0:36} {1:10}".format("Length of penalty steps", "s"),             "%1.4f" %self._timepenstp
      print "*", "{0:36} {1:10}".format("Update parameters in penalty step ", "T/F"), self._updatepenstp
      print ""
      #
      print "INTEGRATED LUMINOSITY PARAMETERS"
      print ""
      print "*", "{0:36} {1:10}".format("Cross section for pile-up", "mb"),          "%1.4f" %self._xsec
      print "*", "{0:36} {1:10}".format("Cross section for burn-off", "mb"),         "%1.4f" %self._xsecburn
      print "*", "{0:36} {1:10}".format("Dedicated time to physics", "days yr-1"),   "%1.4f" %self._days
      print "*", "{0:36} {1:10}".format("Efficiency", "1"),                          "%1.4f" %self._eff
      print "*", "{0:36} {1:10}".format("Turn-around time", "h"),                    "%1.4f" %self._turnar
      print ""
      #
      print "OTHER PARAMETERS"
      print ""
      print "*", "{0:36} {1:10}".format("Name of simulation settings config.", "-"), self._name
      print "*", "{0:36} {1:10}".format("Number of Luminosity module version", "-"), self._version
      print "*", "{0:36} {1:10}".format("ID tag", "-"),                              self._ID
      print "*", "{0:36} {1:10}".format("Only parameters table requested", "-"),     self._table
      print ""
    
    else:
      
      for i in range(beam._nip):
        for j in range(len(self._levtech[i])):
          if j == 0: print >> foutname, "@", "{0:22} {1:4}".format("initlevtech" + str(i),                "%s"),    self._levtech[i][j]
          else:      print >> foutname, "@", "{0:22} {1:4}".format(    "levtech" + str(i) + "_" + str(j), "%s"),    self._levtech[i][j]
      for i in range(beam._nip):
        for j in range(len(self._levvar[i])):
          if j == 0: print >> foutname, "@", "{0:22} {1:4}".format("initlevvar" + str(i),                "%s"),     self._levvar[i][j]
          else:      print >> foutname, "@", "{0:22} {1:4}".format(    "levvar" + str(i) + "_" + str(j), "%s"),     self._levvar[i][j]
      
      for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("levlumi" + str(i), "%le"),       self._levlumi[i]
      for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("levppus" + str(i), "%le"),       self._levppus[i]
      for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("ck" + str(i), "%b"),             self._ck[i]
      for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("constbetar" + str(i), "%b"),     self._constbetar[i]
      for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("sepLRconst" + str(i), "%b"),     self._sepLRconst[i]
      print >> foutname, "@", "{0:22} {1:4}".format("longdens", "%s"),                  self._longdens
      #
      print >> foutname, "@", "{0:22} {1:4}".format("p", "%le"),                        self._p
      print >> foutname, "@", "{0:22} {1:4}".format("step", "%le"),                     self._step
      print >> foutname, "@", "{0:22} {1:4}".format("optimumfill", "%s"),               self._optimumfill
      print >> foutname, "@", "{0:22} {1:4}".format("maxfill", "%le"),                  self._maxfill
      print >> foutname, "@", "{0:22} {1:4}".format("penstp", "%b"),                    self._penstp
      print >> foutname, "@", "{0:22} {1:4}".format("timepenstp", "%le"),               self._timepenstp
      print >> foutname, "@", "{0:22} {1:4}".format("updatepenstp", "%b"),              self._updatepenstp
      #
      print >> foutname, "@", "{0:22} {1:4}".format("xsec", "%le"),                     self._xsec
      print >> foutname, "@", "{0:22} {1:4}".format("xsecburn", "%le"),                 self._xsecburn
      print >> foutname, "@", "{0:22} {1:4}".format("days", "%le"),                     self._days
      print >> foutname, "@", "{0:22} {1:4}".format("eff", "%le"),                      self._eff
      print >> foutname, "@", "{0:22} {1:4}".format("turnar", "%le"),                   self._turnar
      #
      print >> foutname, "@", "{0:22} {1:4}".format("name", "%s"),                      self._name
      print >> foutname, "@", "{0:22} {1:4}".format("version", "%s"),                   self._version
      print >> foutname, "@", "{0:22} {1:4}".format("ID", "%s"),                        self._ID
      print >> foutname, "@", "{0:22} {1:4}".format("table", "%b"),                     self._table
  
  def PrintVirtualLumiParam(self, beam, foutname):
    """ Prints to the standard output and to file 'foutname' the virtual values of luminosities, pile-up and beam-beam parameters for all IPs, wih and without CCs, as well as beam sizes and Piwinski angle for all IPs. It also computes (and returns) the IBS tau in the three coordinates. """
    
    # Beam sized and Piwinski paramters for all IPs
    
    sig    = []
    xplane = []
    piw    = []
    for i in range(beam._nip):
      sig.append(beam.getsigma(i))
      xplane.append(beam._xplane[i])
      if self._longdens == "qGaussian":
        sigsrms = beam._sigs * cst.FactFWHMGauss/4./np.sqrt(2.)
        piw.append(sigsrms*beam._phi[i]/sig[i][xplane[i]]/2.)
      else:
        piw.append(beam._sigs*beam._phi[i]/sig[i][xplane[i]]/2.)

    print "VIRTUAL PARAMETERS"
    print ""
    print "*", "{0:36} {1:10}".format("Virtual sigma horizontal", "m"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:12}".format("%1.4e" %sig[i][0]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual sigma vertical", "m"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:12}".format("%1.4e" %sig[i][1]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual Piwinski parameter", "1"),    
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:12}".format("%1.4f" %(piw[i])),
    print ""
    
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtsigx" + str(i), "%le"),           sig[i][0]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtsigy" + str(i), "%le"),           sig[i][1]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtpiw" + str(i), "%le"),            piw[i]
    
    # IBS tau
    
    self.IBStauxz("N/A", beam)  # [h]
    
    print "*", "{0:36} {1:10}".format("Virtual taux IBS", "h"),                            "%1.4f" %beam._tau_ibs[0]
    print "*", "{0:36} {1:10}".format("Virtual tauy IBS", "h"),                            "%1.4f" %beam._tau_ibs[1]
    print "*", "{0:36} {1:10}".format("Virtual tauz IBS", "h"),                            "%1.4f" %beam._tau_ibs[2]
    
    print >> foutname, "@", "{0:22} {1:4}".format("virttauxIBS", "%le"),                "%1.4f" %beam._tau_ibs[0]
    print >> foutname, "@", "{0:22} {1:4}".format("virttauyIBS", "%le"),                "%1.4f" %beam._tau_ibs[1]
    print >> foutname, "@", "{0:22} {1:4}".format("virttauzIBS", "%le"),                "%1.4f" %beam._tau_ibs[2]
 
    # Virtual tau from CC missing
    
    lumiall     = []
    lregionall  = []
    plregionall = []
    ltimeall    = []
    pltimeall   = []
    puall       = []
    ppusall     = []
    pputall     = []
    xiall       = []
    for i in range(beam._nip):
      lumiall.append(self.getlumiip(beam, i))
      lregionall.append(self.GetLumRegion(beam, i))
      plregionall.append(self.GetPeakLumRegion(beam, i))
      ltimeall.append(self.GetLumTime(beam, i))
      pltimeall.append(self.GetPeakLumTime(beam, i))
      puall.append(self.GetPileUp(beam._frev, lumiall[i][0], beam._nbunch[i]))
      ppusall.append(puall[i]*plregionall[i]/1e3)
      pputall.append(puall[i]*pltimeall[i]/1e9)
      xiall.append(self.GetBBParam(beam, i))
    
    print ""
    print "*", "{0:36} {1:10}".format("Virtual luminosity", "cm-2 s-1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4e" %lumiall[i][0]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual luminosity no reduct.", "cm-2 s-1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4e" %lumiall[i][1]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual reduction factor", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.6f" %lumiall[i][2]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual RMS luminous region", "mm", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %(1e3*lregionall[i])),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual peak luminous region", "mm-1", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %(1e-3*plregionall[i])),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual RMS luminous time", "ns", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %(1e9*ltimeall[i])),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual peak luminous time", "ns-1", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %(1e-9*pltimeall[i])),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual pile-up", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %puall[i]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual peak pile-up s-density", "mm-1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %ppusall[i]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual peak pile-up t-density", "ns-1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %pputall[i]),
    print "\n"
    
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtlumi" + str(i), "%le"),           lumiall[i][0]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtluminored" + str(i), "%le"),      lumiall[i][1]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtredfactor" + str(i), "%le"),      lumiall[i][2]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtlregion" + str(i), "%le"),        lregionall[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtplregion" + str(i), "%le"),       plregionall[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtltime" + str(i), "%le"),          ltimeall[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtpltime" + str(i), "%le"),         pltimeall[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtpu" + str(i), "%le"),             puall[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtppus" + str(i), "%le"),           ppusall[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtpput" + str(i), "%le"),           pputall[i]
    
    print "*", "{0:36} {1:10}".format("Virtual BB horizontal", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4e" %xiall[i][0]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual BB vertical", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4e" %xiall[i][1]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual BB mean", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4e" %xiall[i][2]),
    print ""
   
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtxix" + str(i), "%le"),         xiall[i][0]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtxiy" + str(i), "%le"),         xiall[i][1]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtxim" + str(i), "%le"),         xiall[i][2]
    
    # Virtual luminosities, pile-up and beam-beam without CC: 
    
    lumiallNoCC     = []
    lregionallNoCC  = []
    plregionallNoCC = []
    ltimeallNoCC    = []
    pltimeallNoCC   = []
    puallNoCC       = []
    ppusallNoCC     = []
    pputallNoCC     = []
    xiallNoCC       = []
    print ""
    for i in range(beam._nip):
      if (beam._initoncc[i] == 0.0):
        print "> CCs for ip" + str(i), "are already off..."
        lumiallNoCC.append(lumiall[i])
        lregionallNoCC.append(lregionall[i])
        plregionallNoCC.append(plregionall[i])
        ltimeallNoCC.append(ltimeall[i])
        pltimeallNoCC.append(pltimeall[i])
        puallNoCC.append(puall[i])
        ppusallNoCC.append(ppusall[i])
        pputallNoCC.append(pputall[i])
        xiallNoCC.append(xiall[i])
      else:
        print "> Computation of parameters for ip" + str(i), "with CCs off..."
        # CC state will be saved for reference, then turned off (it is need to recompute oncc)
        save_oncc   = beam._oncc[i]
        beam._oncc[i] = 0.0 # Switch off CCs
        self.GetPhiCR(beam, i)
        # Compute paramters now with the CCs off
        lumiallNoCC.append(self.getlumiip(beam, i))
        lregionallNoCC.append(self.GetLumRegion(beam, i))
        plregionallNoCC.append(self.GetPeakLumRegion(beam, i))
        ltimeallNoCC.append(self.GetLumTime(beam, i))
        pltimeallNoCC.append(self.GetPeakLumTime(beam, i))
        puallNoCC.append(self.GetPileUp(beam._frev, lumiallNoCC[i][0], beam._nbunch[i]))
        ppusallNoCC.append(puallNoCC[i]*plregionallNoCC[i]/1e3)
        pputallNoCC.append(puallNoCC[i]*pltimeallNoCC[i]/1e9)
        xiallNoCC.append(self.GetBBParam(beam, i))
        # Getting back to original CC state...
        beam._oncc[i]  = save_oncc
        self.GetPhiCR(beam, i)
    print ""

    print "*", "{0:36} {1:10}".format("Virtual luminosity w/o CC", "cm-2 s-1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4e" %lumiallNoCC[i][0]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual luminosity no reduct. w/o CC", "cm-2 s-1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4e" %lumiallNoCC[i][1]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual reduction factor w/o CC", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.6f" %lumiallNoCC[i][2]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual RMS luminous region w/o CC", "mm", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %(1e3*lregionallNoCC[i])),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual peak lum region w/o CC", "mm-1", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %(1e-3*plregionallNoCC[i])),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual RMS luminous time w/o CC", "ns", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %(1e9*ltimeallNoCC[i])),
    print "" 
    print "*", "{0:36} {1:10}".format("Virtual peak lum. time w/o CC", "ns-1", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %(1e-9*pltimeallNoCC[i])),
    print "" 
    print "*", "{0:36} {1:10}".format("Virtual pile-up w/o CC", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %puallNoCC[i]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual peak pile-up s-dens. w/o CC", "mm-1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %ppusallNoCC[i]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual peak pile-up t-dens. w/o CC", "ns-1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4f" %pputallNoCC[i]),
    print "\n"
    
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtluminoCC" + str(i), "%le"),        lumiallNoCC[i][0]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtluminorednoCC" + str(i), "%le"),   lumiallNoCC[i][1]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtredfactornoCC" + str(i), "%le"),   lumiallNoCC[i][2]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtlregionnoCC" + str(i), "%le"),     lregionallNoCC[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtplregionnoCC" + str(i), "%le"),    plregionallNoCC[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtltimenoCC" + str(i), "%le"),       ltimeallNoCC[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtpltimenoCC" + str(i), "%le"),      pltimeallNoCC[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtpunoCC" + str(i), "%le"),          puallNoCC[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtppusnoCC" + str(i), "%le"),        ppusallNoCC[i]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtpputnoCC" + str(i), "%le"),        pputallNoCC[i]
    
    print "*", "{0:36} {1:10}".format("Virtual BB horizontal w/o CC", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4e" %xiallNoCC[i][0]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual BB vertical w/o CC", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4e" %xiallNoCC[i][1]),
    print ""
    print "*", "{0:36} {1:10}".format("Virtual BB mean w/o CC", "1"),
    for i in range(beam._nip): print beam._ipnames[i] + ":", "{0:14}".format("%1.4e" %xiallNoCC[i][2]),
    print "\n"
    
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtxixnoCC" + str(i), "%le"),       xiallNoCC[i][0]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtxiynoCC" + str(i), "%le"),       xiallNoCC[i][1]
    for i in range(beam._nip): print >> foutname, "@", "{0:22} {1:4}".format("virtximnoCC" + str(i), "%le"),       xiallNoCC[i][2]
  
  def WriteFILElevelHeaders(self, FILElevel, beam):
    """ Write the column names in the output file of the full levelling simulation. The names are, in order (for 3 IPs):
        #
        1.  TIME         [h]           Time
        2.  PPB          [1]           Bunch population
        3.  SIGS         [m]           Gaussian RMS bunch length
        4.  RMSSIGS      [m]           RMS bunch length (different to SIGS for q-Gaussian)
        5.  FWHM         [m]           Full width at half maximum
        6.  DPP          [1]           Energy spread
        7.  EPSXN        [m]           Horizontal normalized emittance
        8.  EPSYN        [m]           Vertical   normalized emittance
        9.  EPSXNIBS     [m]           Contribution to horizontal normalized emittance due to IBS
        10. EPSYNIBS     [m]           Contribution to vertical normalized emittance due to IBS
        11. EPSXNSR      [m]           Contribution to horizontal normalized emittance due to SR
        12. EPSYNSR      [m]           Contribution to vertical   normalized emittance due to SR
        13. EPSXNOBS     [m]
        14. EPSYNOBS     [m]
        15. EPSXNCC0     [m]           ### Note: In results, EPSXN[i-1] / ((EPSXNCC0[i]-EPSXNCC0[i-1])/(TIME[i]-TIME[i-1]) * 0.15/BETA0[i]) = emittance growth from CC in um/h
        16. EPSYNCC0     [m]
        17. EPSXNCC1     [m]
        18. EPSYNCC1     [m]
        #
        19. TAUXIBS      [h]
        20. TAUYIBS      [h]
        21. TAUZIBS      [h]
        22. TAUXCC0      [h]
        23. TAUYCC0      [h]
        24. TAUXCC1      [h]
        25. TAUYCC1      [h]
        #
        26. BETX0        [m]           Horizontal beta* at IP1
        27. BETY0        [m]           Vertical   beta* at IP1
        28. LUMITOT0     [cm-2 s-1]    Instantaneous luminosity at IP1
        29. REDUC0       [1]           Luminosity reduction factor at IP1
        30. BETX1        [m]           Horizontal beta* at IP5
        31. BETY1        [m]           Vertical   beta* at IP5
        32. LUMITOT1     [cm-2 s-1]    Instantaneous luminosity at IP5
        33. REDUC1       [1]           Luminosity reduction factor at IP5
        34. BETX2        [m]           Horizontal beta* at IP8
        35. BETY2        [m]           Vertical   beta* at IP8
        36. LUMITOT2     [cm-2 s-1]    Instantaneous luminosity at IP8
        37. REDUC2       [1]           Luminosity reduction factor at IP8
        #
        38. PHI0         [rad]         Full crossing angle at IP1
        39. SEPLR0       [sigma]       Normalized beam-beam separation at IP1
        40. PHICR0       [rad]         Crab cavity angle at IP1
        41. ONCC0        [1]
        42. PHICK0       [rad]
        43. ONCK0        [1]
        44. PHI1         [rad]         Full crossing angle at IP5
        45. SEPLR1       [sigma]       Normalized beam-beam separation at IP5
        46. PHICR1       [rad]         Crab cavity angle at IP5
        47. ONCC1        [1]
        48. PHICK1       [rad]
        49. ONCK1        [1]
        50. PHI2         [rad]         Full crossing angle at IP8
        51. SEPLR2       [sigma]       Normalized beam-beam separation at IP8
        52. PHICR2       [rad]         Crab cavity angle at IP8
        53. ONCC2        [1]
        54. PHICK2       [rad]
        55. ONCK2        [1]
        56. PARSEP2      [1]
        #
        57. LREGION0     [m]
        58. PLREGION0    [1/m]
        59. LTIME0       [s]
        60. PLTIMEN0     [1/s]
        61. PU0          [1]
        62. PPUS0        [mm-1]
        63. PPUT0        [ns-1]
        64. LREGION1     [m]
        65. PLREGION1    [1/m]
        66. LTIME1       [s]
        67. PLTIMEN1     [1/s]
        68. PU1          [1]
        69. PPUS1        [mm-1]
        70. PPUT1        [ns-1]
        71. LREGION2     [m]
        72. PLREGION2    [1/m]
        73. LTIME2       [s]
        74. PLTIMEN2     [1/s]
        75. PU2          [1]
        76. PPUS2        [mm-1]
        77. PPUT2        [ns-1]
        #
        78. XIX0         [1]
        79. XIY0         [1]
        80. XIM0         [1]
        81. XIX1         [1]
        82. XIY1         [1]
        83. XIM1         [1]
        84. XIX2         [1]
        85. XIY2         [1]
        86. XIM2         [1]
        #
        87. XIXTOT       [1]
        88. XIYTOT       [1]
        89. XIMTOT       [1]
        #
        90.  LINT0        [1e-39 cm-2 = fb-1]
        91.  LINT1        [1e-39 cm-2 = fb-1]
        92.  LINT2        [1e-39 cm-2 = fb-1]
        #
        93.  LINT         [1e-39 cm-2 = fb-1]
        94.  RATE         [s-1]
        95.  LLIFETIME    [h]
        96.  STEPLENGTH   [s]
        97.  LUMIINTSTEP  [cm-2]
        #
        98.  PUINT0       [s]
        99.  EPPUS0       [mm-1]
        100. EPPUT0       [ns-1]
        101. PUINT1       [s]
        102. EPPUS1       [mm-1]
        103. EPPUT1       [ns-1]
        104. PUINT2       [s]
        105. EPPUS2       [mm-1]
        106. EPPUT2       [ns-1]
        #
        107. STEP         -
        108. TYPESTEP     -
    """
    # print >> FILElevel, "@ name     %s", self._name
    # print >> FILElevel, "@ version  %s", self._version
    # print >> FILElevel, "@ ID       %s", self._ID
    print >> FILElevel, "*",
    print >> FILElevel, "{0:14}".format("TIME"),
    print >> FILElevel, "{0:16}".format("PPB"),
    print >> FILElevel, "{0:16}".format("SIGS"),
    print >> FILElevel, "{0:16}".format("RMSSIGS"),
    print >> FILElevel, "{0:16}".format("FWHM"),
    print >> FILElevel, "{0:16}".format("DPP"),
    print >> FILElevel, "{0:16}".format("EPSXN"),
    print >> FILElevel, "{0:16}".format("EPSYN"),
    print >> FILElevel, "{0:16}".format("EPSXNIBS"),
    print >> FILElevel, "{0:16}".format("EPSYNIBS"),
    print >> FILElevel, "{0:16}".format("EPSXNSR"),
    print >> FILElevel, "{0:16}".format("EPSYNSR"),
    print >> FILElevel, "{0:16}".format("EPSXNOBS"),
    print >> FILElevel, "{0:16}".format("EPSYNOBS"),
    print >> FILElevel, "{0:16}".format("EPSXNCC0"),
    print >> FILElevel, "{0:16}".format("EPSYNCC0"),
    print >> FILElevel, "{0:16}".format("EPSXNCC1"),
    print >> FILElevel, "{0:16}".format("EPSYNCC1"),
    #
    print >> FILElevel, "{0:16}".format("TAUXIBS"),
    print >> FILElevel, "{0:16}".format("TAUYIBS"),
    print >> FILElevel, "{0:16}".format("TAUZIBS"),
    print >> FILElevel, "{0:16}".format("TAUXCC0"),
    print >> FILElevel, "{0:16}".format("TAUYCC0"),
    print >> FILElevel, "{0:16}".format("TAUXCC1"),
    print >> FILElevel, "{0:16}".format("TAUYCC1"),
    #
    for i in range(beam._nip):
      print >> FILElevel, "{0:16}".format("BETX"    + str(i)),
      print >> FILElevel, "{0:16}".format("BETY"    + str(i)),
      print >> FILElevel, "{0:16}".format("LUMITOT" + str(i)),
      print >> FILElevel, "{0:16}".format("REDUC"   + str(i)),
    #
    for i in range(beam._nip):
      print >> FILElevel, "{0:16}".format("PHI"   + str(i)),
      print >> FILElevel, "{0:16}".format("SEPLR" + str(i)),
      print >> FILElevel, "{0:16}".format("PHICR" + str(i)),
      print >> FILElevel, "{0:16}".format("ONCC"  + str(i)),
      print >> FILElevel, "{0:16}".format("PHICK" + str(i)),
      print >> FILElevel, "{0:16}".format("ONCK" + str(i)),
    print >> FILElevel, "{0:16}".format("PARSEP2"),
    #
    for i in range(beam._nip):
      print >> FILElevel, "{0:16}".format("LREGION"  + str(i)),
      print >> FILElevel, "{0:16}".format("PLREGION" + str(i)),
      print >> FILElevel, "{0:16}".format("LTIME"    + str(i)),
      print >> FILElevel, "{0:16}".format("PLTIME"   + str(i)),
      print >> FILElevel, "{0:16}".format("PU"   + str(i)),
      print >> FILElevel, "{0:16}".format("PPUS"     + str(i)),
      print >> FILElevel, "{0:16}".format("PPUT"     + str(i)),
    #
    for i in range(beam._nip):
      print >> FILElevel, "{0:16}".format("XIX" + str(i)),
      print >> FILElevel, "{0:16}".format("XIY" + str(i)),
      print >> FILElevel, "{0:16}".format("XIM" + str(i)),
    #
    print >> FILElevel, "{0:16}".format("XIXTOT"),
    print >> FILElevel, "{0:16}".format("XIYTOT"),
    print >> FILElevel, "{0:16}".format("XIMTOT"),
    #
    for i in range(beam._nip):
      print >> FILElevel, "{0:16}".format("LINT" + str(i)),
    #
    print >> FILElevel, "{0:16}".format("LINT"),
    print >> FILElevel, "{0:16}".format("RATE"),
    print >> FILElevel, "{0:16}".format("LLIFETIME"),
    print >> FILElevel, "{0:16}".format("STEPLENGTH"),
    print >> FILElevel, "{0:16}".format("LUMIINTSTEP"),
    #
    for i in range(beam._nip):
      print >> FILElevel, "{0:16}".format("PUINT" + str(i)),
      print >> FILElevel, "{0:16}".format("EPPUS" + str(i)),
      print >> FILElevel, "{0:16}".format("EPPUT" + str(i)),
    #
    print >> FILElevel, "{0:16}".format("STEP"),
    print >> FILElevel, "{0:16}".format("TYPESTEP")
    print >> FILElevel, "$",
    print >> FILElevel, "{0:13} ".format("%le"),
    print >> FILElevel, (34+24*beam._nip)*"{0:16} ".format("%le"),   # 33 parameters + 24 parameters for each IP
    print >> FILElevel, "{0:16}".format("%s")

  #============================================================================================================#
  
  # Do Fill:

  def DoFill(self, beam):
    """ Step-based simulation of the fill. """
    
    bannersclass = Banners()
    
    # Write table to file:
    NAMEtable = "table_" + self._ID + ".out"
    FILEtable = open(NAMEtable, "w")
    beam.printBeamParam(self._longdens, tofile = True, foutname = FILEtable)
    # Large numbers to ignore their effect when None
    for i in range(2):
      if beam._tau_obs[i] == None: beam._tau_obs[i] = 1e8
      if beam._tau_cc  == None:    beam._tau_cc     = 1e8
    self.PrintLumiParam(beam, tofile = True,  foutname = FILEtable)
    self.PrintVirtualLumiParam(beam, FILEtable)
    print "> Writting '" + NAMEtable + "'..."
    FILEtable.close()
    
    # Exit if only the table was requested
    if self._table == True:
      bannersclass.header0("END OF \"LUMINOSITY LEVELLING\"")
      sys.exit()
    
    ##############
    #    EXIT    #
    ##############
    
    bannersclass.header1("START OF FILL")
    
    # Empty file to store IBS results per step:
    NAMEupdateIBS = "updateIBS_" + self._ID + ".madx"
    FILEupdateIBS = open(NAMEupdateIBS, "w")
    FILEupdateIBS.close()
   
    # File to store levelling results per step:
    NAMElevel = "level_" + self._ID + ".out"
    print "\n> Writting '" + NAMElevel + "'..."
    FILElevel = open(NAMElevel, "w")
    self.WriteFILElevelHeaders(FILElevel, beam)
    FILElevel.close()
    
    # File to store the luminosity as a piece-wise function in gnuplot-format ("lumiplot")
    NAMElumiplot = "lumiplot_" + self._ID + ".gnu"
    FILElumiplot = open(NAMElumiplot, "w")
    auxfuncclass = AuxFunc()
    print >> FILElumiplot, "L_" + auxfuncclass.DotToDash(self._ID) + "(x,SECtoHR) = ",
    print >> FILElumiplot, "x < 0 ? NaN : ",
    FILElumiplot.close()
    
    # File to save the beam properties to be load in MAD-X each step to recompute IBS:
    FILEupdateIBS = open(NAMEupdateIBS, "a")
    print >> FILEupdateIBS, "! 0.0"
    print >> FILEupdateIBS, "tuneb     =", beam._tuneb, ";"
    print >> FILEupdateIBS, "epsxn     =", beam._epsn[0], ";"
    print >> FILEupdateIBS, "epsyn     =", beam._epsn[1], ";"
    print >> FILEupdateIBS, "kappa     =", beam._kappa, ";"
    print >> FILEupdateIBS, "gamma     =", beam._gamma, ";"
    print >> FILEupdateIBS, "sigs      =", beam._rmssigs, ";" # beam._sigs, ";"
    #if self._longdens == "qGaussian": print >> FILEupdateIBS, "*", cst.FactFWHMGauss, "/4./sqrt(2.)", # before, now it has changed in the previous line for sigs -> sigsrms
    #print >> FILEupdateIBS, ";"
    print >> FILEupdateIBS, "dpp       =", beam._dpp,   ";"
    print >> FILEupdateIBS, "ppb       =", beam._ppb, ";"
    print >> FILEupdateIBS, "NRJ       =", beam._momeV/1e9, ";"
    FILEupdateIBS.close()
    
    # Initial parameters and create empy list to store their evolution:
    
    i          = 0.    
    countsteps = 0
    typestep   = "FIRST"
    timeinleveling   = -1
    flagendlevel     = False
    flagpointaft6h   = False
    timepenstptmp    = self._timepenstp
    flagsepLRstep    = [False]*beam._nip
    
    #TAUIBS  = [tauIBS[0], tauIBS[1], tauIBS[2]]
    
    lintiall    = []
    lintialltmp = [0.0]*beam._nip
    
    lint       = 0.0
    linttmp    = 1e-9 # small number to pass the first while evaluation
    ltot       = self._days *24*3600. * self._eff  # Time doing physics in a year [s yr-1]
    alllint    = []
    linttmpall = [0.0]*beam._nip
    
    puint = [0.0]*beam._nip
    eppus = [0.0]*beam._nip
    epput = [0.0]*beam._nip
    
    listTIME        = []
    listPPB         = []
    listSIGS        = []
    listRMSSIGS     = []
    listFWHM        = []
    listDPP         = []
    listEPSN        = []
    listEPSNIBS     = []
    listEPSNSR      = []
    listEPSNOBS     = []
    listEPSNCC0     = []
    listEPSNCC1     = []
    #
    listTAUIBS      = []
    listTAUCC0      = []
    listTAUCC1      = []
    #
    listBETXi       = []
    listBETYi       = []
    listLUMITOTi    = []
    listREDUCi      = []
    #
    listPHIi        = []
    listSEPLRi      = []
    listPHICRi      = []
    listONCCi       = []
    listPHICKi      = []
    listONCKi       = []
    listPARSEP2     = []
    #
    listLREGIONi    = []
    listPLREGIONi   = []
    listLTIMEi      = []
    listPLTIMEi     = []
    listPUi         = []
    listPPUSi       = []
    listPPUTi       = []
    #
    listXIXi        = []
    listXIYi        = []
    listXIMi        = []
    #
    listXIXTOT      = []
    listXIYTOT      = []
    listXIMTOT      = []
    #
    listLINTi       = []
    #
    listLINT        = []
    listRATE        = []
    listLLIFETIME   = []
    listSTEPLENGTH  = []
    listLUMIINTSTEP = []
    #
    listPUINTi      = []
    listEPPUSi      = []
    listEPPUTi      = []
    #
    listSTEP        = []
    listTYPESTEP    = []
    
    if self._optimumfill == True: cond1 = linttmp > linttmpall[-3] or linttmp == 0.0
    else:                         cond1 = not self._optimumfill
    
    while (cond1 and i <= int(self._maxfill)):
      
      # Keep levelling and lumiplot files open to append data
      FILElevel    = open(NAMElevel,    "a")
      FILElumiplot = open(NAMElumiplot, "a")
      
      # Print the number of step and the time
      textbanner = "STEP: " + str(countsteps) + " (" + "%1.2f" %i + " s)"
      bannersclass.header2(textbanner)
      print "*", "{0:36} {1:10}".format("Type of step", "-"), typestep
      print "*", "{0:36} {1:10}".format("Time", "h"),         "%1.4f" %(i/3600.)
      
      # Adaptive xsing:
      print "\nADAPTIVE CROSSING\n"
      for j in range(beam._nip):
        if self._adaptivexsing[j] == True:
          self.GetPhiFromPPBFor6sigmaDA(beam, j)
      print ""
      
      # Print the ppb at the start of the fill and, for subsequent steps, evolve the parameters (emittances, cbunch length, energy sperand and IBS tau). Also, the integrated luminosit y at this step comes from the previous step 
      if (i == 0.0):
        print  "\n", "*", "{0:36} {1:10}".format("Initial ppb (ppb)", "1"), "%1.6e" %beam._ppb
      else:
        self.BunchLengthGymnastics(beam, typestep)
        self.ParametersEvolution(beam, typestep, timedecay, listLUMITOTi)
        
        lint  = linttmp

      beam.getrmsfwhm(self._longdens)
      
      # Save values into their corresponding lists
      listTIME.append(i/3600.0)
      listPPB.append(beam._ppb)
      listSIGS.append(beam._sigs)
      listRMSSIGS.append(beam._rmssigs)
      listFWHM.append(beam._fwhm)
      listDPP.append(beam._dpp)
      listEPSN.append([   beam._epsn[0],     beam._epsn[1]])
      listEPSNIBS.append([beam._epsn_ibs[0], beam._epsn_ibs[1]])
      listEPSNSR.append([ beam._epsn_sr[0],  beam._epsn_sr[1]])
      listEPSNOBS.append([beam._epsn_obs[0], beam._epsn_obs[1]])
      listEPSNCC0.append([ beam._epsn_cc_ip0[0],  beam._epsn_cc_ip0[1]])
      listEPSNCC1.append([ beam._epsn_cc_ip1[0],  beam._epsn_cc_ip1[1]])
      #
      listTAUIBS.append(beam._tau_ibs)
      listTAUCC0.append([beam._tau_cc_0[0], beam._tau_cc_0[1]])
      listTAUCC1.append([beam._tau_cc_0[1], beam._tau_cc_1[1]])
      #
      
      # List to contain the value of a parameter at all the IPs for a given step
      betxall     = []
      betyall     = []
      lumiall     = []
      reducall    = []
      #
      phiall      = []
      sepLRall    = []
      phiCRall    = []
      onccall     = []
      phiCKall    = []
      onckall     = []
      #
      lregionall  = []
      plregionall = []
      ltimeall    = []
      pltimeall   = []
      puall       = []
      ppusall     = []
      pputall     = []
      #
      xixall      = []
      xiyall      = []
      ximall      = []
      #
      lintiall    = []
      
      # Do the computations for all the IPs:
      for j in range(beam._nip):
        
        bannersclass.header3(beam._ipnames[j] + " (nip = " + str(j) + ")")
        
        flagsepLRstep[j] = self.StepReductionSeparation(beam, i, j, flagsepLRstep[j])
        
        # Run the appropiate levelling routine: with beta or with phiCK, for luminosity or peak pile-up levelling
        if (typestep != "PENAL"):
          
          if i > 0.0:
            if self._penstp == True: self.ReductionSeparation(beam, listSTEPLENGTH[-1] + listSTEPLENGTH[-2], j)
            else:                    self.ReductionSeparation(beam, listSTEPLENGTH[-1], j)
            
            if self._incroncc[j] == True: self.IncrOncc(beam, j)
          
          # Start of peak pile-up s-density levelling (do not take into accout penal steps)
          if (j == 0):
            # Check if peak pile-up s-density levelling hast started
            if self._timeppuslev != None and self._flagtimeppuslev == False:
              #self._timeppuslev = i/3600.
              bannersclass.header4("Start peak pile-up s-density lev.: " + "%1.4f" %self._timeppuslev + " h")
              self._flagtimeppuslev = True
              typestep = "PPUSL"
              print ""
              
          # Copy results to ip1 from ip0 if they are identical
          if (j != 1 or not beam._identicalIP01):
            
            if len(self._levtech[j]) > 1:
              for m in range(self._niterlev):
                print "Iteration:", m
                print 12*"-"
                for k in range(len(self._levtech[j])):
                  lumitot = self.GetLevel(beam, j, self._levtech[j][k], self._levvar[j][k], i)
                print ""
            else:
              lumitot = self.GetLevel(beam, j, self._levtech[j][0], self._levvar[j][0], i)
              #print 'ip (j) =', j, self._flagsteplev2, self._steplev2
              if self._flagsteplev2 == True and self._steplev2 == None: # ip8 only works with one variable (1 iteration)  --- there seems to be a bug here
                self._steplev2 = countsteps
                self._timesteplev2 = i/3600.
                bannersclass.header4("Leveling time at IP8: " + "%1.4f" %self._timesteplev2 + " h")
            
            lregion  = self.GetLumRegion(beam, j)     # [m]
            plregion = self.GetPeakLumRegion(beam, j) # [m-1]
            ltime    = self.GetLumTime(beam, j)       # [s]
            pltime   = self.GetPeakLumTime(beam, j)   # [s-1]
            pu       = self.GetPileUp(beam._frev, lumitot[0], beam._nbunch[j])
            ppus     = pu*plregion/1e3  # [event mm-1]
            pput     = pu*pltime/1e9    # [event ns-1]
            
          else:
            beam._beta[j][0] = beam._beta[0][0]
            beam._beta[j][1] = beam._beta[0][1]
            beam._phi[j]   = beam._phi[0]
            beam._sepLR[j] = beam._sepLR[0]
            beam._phiCR[j] = beam._phiCR[0]
            beam._oncc[j]  = beam._oncc[0]
            beam._phiCK[j] = beam._phiCK[0]
            beam._onck[j]  = beam._onck[0]
            if self._sepLRconst[j] == True:
              self.GetNewPhi(beam, j, output = False)
            self._flaglevpu[j] = self._flaglevpu[0]
            # lumitot, lregionj, ltimej, puj, puj currentyly hold the values at ip0, so we can re-use them 
            print 2*" ", "[!]", beam._ipnames[j], "is identical to IP1."
        
        else:
          
          # Emittances have changed, and thus the crossing angle (either in sigma or radians, keeping the other constant as requested)
          if beam._redsepLR[j] == True: self.GetNewPhi(beam, j, output = False)
          else:                         self.GetNewSepLR(beam, j)
          
          if self._incroncc[j] == True: self.GetPhiCR(beam, j)
          
          lumitot  = [0.0, 0.0, 0.0]  # l0, l0*r, r (per IP)
          lregion  = 0.0
          plregion = 0.0 # Formally np.inf
          ltime    = 0.0
          pltime   = 0.0 # Formally np.inf
          pu       = 0.0
          ppus     = 0.0
          pput     = 0.0
        
        # Results from this levelling step: betas, angles
        # print "\nRESULTS FOR THIS IP AT THIS STEP"
        print ""
        print "Results:"
        print ""
        print "*", "{0:36} {1:10}".format("Beta horizontal", "m"),              "%1.8f" %beam._beta[j][0]
        print "*", "{0:36} {1:10}".format("Beta vertical",   "m"),              "%1.8f" %beam._beta[j][1]
        print ""
        print "*", "{0:36} {1:10}".format("Crossing angle", "rad"),             "%1.6e" %beam._phi[j]
        print "*", "{0:36} {1:10}".format("BB-LR separation", "sigma"),         "%1.8f" %beam._sepLR[j] #self.GetNewSepLR(beam, j)
        print "*", "{0:36} {1:10}".format("Crab cavity angle", "rad"),          "%1.6e" %beam._phiCR[j]
        print "*", "{0:36} {1:10}".format("Crab cavity ON (oncc)", "1"),        "%1.8f" %beam._oncc[j]
        print "*", "{0:36} {1:10}".format("Crab kissing angle", "rad"),         "%1.6e" %beam._phiCK[j]
        print "*", "{0:36} {1:10}".format("Crab kissing ON (onck)", "1"),       "%1.8f" %beam._onck[j]
        #print ''
        #print 'self._levvar =', self._levvar
        #print 'j = ', j
        #print 'self._levvar[j][0]', self._levvar[j][0]
        #print 'self._steplev2 =', self._steplev2
        #print 'self._flagsteplev2 =', self._flagsteplev2
        #print ''
        if self._levvar[j][0] == "ParSep":
          print "*", "{0:36} {1:10}".format("Parallel separation", "m"),
          print "%1.6e" %beam._parsep
        print ""
        
        # Save values into their corresponding lists
        betxall.append(beam._beta[j][0])
        betyall.append(beam._beta[j][1])
        phiall.append(beam._phi[j])
        sepLRall.append(beam._sepLR[j]) #self.GetNewSepLR(beam, j))
        phiCRall.append(beam._phiCR[j])
        onccall.append(beam._oncc[j])
        phiCKall.append(beam._phiCK[j])
        onckall.append(beam._onck[j])
        onckall.append(beam._onck[j])
        # no parsep
        
        # Results from this levelling step: luminosity, pile-up, luminous region and time
        print "*", "{0:36} {1:10}".format("Luminosity", "cm-2 s-1"),            "%1.6e" %lumitot[0]
        print "*", "{0:36} {1:10}".format("Reduction factor", "1"),             "%1.6f" %lumitot[2]
        print "*", "{0:36} {1:10}".format("RMS luminous region", "mm"),         "%1.8f" %(1e3*lregion)
        print "*", "{0:36} {1:10}".format("Peak luminous region", "mm-1"),      "%1.8f" %(1e-3*plregion)
        print "*", "{0:36} {1:10}".format("RMS luminous time", "ns"),           "%1.8f" %(1e9*ltime)
        print "*", "{0:36} {1:10}".format("Peak luminous time", "ns-1"),        "%1.8f" %(1e-9*pltime)
        print "*", "{0:36} {1:10}".format("Pile-up ", "1"),                     "%1.8f" %pu
        print "*", "{0:36} {1:10}".format("Peak pile-up s-density", "mm-1"),    "%1.8f" %ppus
        print "*", "{0:36} {1:10}".format("Peak pile-up t-density", "ns-1"),    "%1.8f" %pput
        print ""
        
        # Save values into their corresponding lists
        lumiall.append(lumitot[0])
        reducall.append(lumitot[2])
        lregionall.append(lregion)
        plregionall.append(plregion)
        ltimeall.append(ltime)
        pltimeall.append(pltime)
        puall.append(pu)
        ppusall.append(ppus)
        pputall.append(pput)
        
        # Beam-beam tune-shift parameters
        xi = self.GetBBParam(beam, j)
        print "*", "{0:36} {1:10}".format("BB tune-shift horizontal", "1"),    "%1.4e" %xi[0]
        print "*", "{0:36} {1:10}".format("BB tune-shift vertical",   "1"),    "%1.4e" %xi[1]
        print "*", "{0:36} {1:10}".format("BB tune-shift mean",       "1"),    "%1.4e" %xi[2]
        print ""
        
        # Save values into their corresponding lists
        xixall.append(xi[0])
        xiyall.append(xi[1])
        ximall.append(xi[2])
        
        #print "CHECK CONDTION"
        #print self._levvar[j]
        #print self._levvar[j][-1]
        
        # ip0 drives the levelling routine
        if (j == 0):
          
          # Check if still leveling:
          condition = timeinleveling < 0
          #print "timeinleveling < 0 :", condition
          
          for k in range(len(self._levtech[j])):
            
            if self._levtech[j][k] == "levlumi":
              condition = condition and lumitot[0] < self._levlumi[j]
            elif self._levtech[j][k] == "levppus":
              condition = condition and ppusall[j] < self._levppus[j]
              
            if self._levvar[j][k] == "Beta":
              condition = condition and beam._beta[j][1-beam._xplane[j]] <= beam._betamin[j][1-beam._xplane[j]] 
            elif self._levvar[j][k] == "PhiCR":
              condition = condition and beam._phiCR[j] == beam._phi[j]
            elif self._levvar[j][k] == "PhiCK":
              condition = condition and beam._phiCK[j] == 0.0
            # parsep does not apply for ip0
          
          # The time where the reduction in beta or phiCK is not longer possible
          if condition:
            timeinleveling = i/3600.0  # In [h]
            bannersclass.header4("Leveling time: " + "%1.4f" %timeinleveling + " h")
            flagendlevel = True
            self._timepenstp = 0.0
            typestep = "LAST"
        
        lintiall.append(lintialltmp[j])
      
      # Save values into their corresponding lists
      listBETXi.append([])
      listBETYi.append([])
      listLUMITOTi.append([])
      listREDUCi.append([])
      for j in range(beam._nip):
        listBETXi[-1].append(betxall[j])
        listBETYi[-1].append(betyall[j])
        listLUMITOTi[-1].append(lumiall[j])
        listREDUCi[-1].append(reducall[j])
      #
      listPHIi.append([])
      listSEPLRi.append([])
      listPHICRi.append([])
      listONCCi.append([])
      listPHICKi.append([])
      listONCKi.append([])
      for j in range(beam._nip):
        listPHIi[-1].append(phiall[j])
        listSEPLRi[-1].append(sepLRall[j])
        listPHICRi[-1].append(phiCRall[j])
        listONCCi[-1].append(onccall[j])
        listPHICKi[-1].append(phiCKall[j])
        listONCKi[-1].append(onckall[j])
      #
      listPARSEP2.append(beam._parsep)
      #
      listLREGIONi.append([])
      listPLREGIONi.append([])
      listLTIMEi.append([])
      listPLTIMEi.append([])
      listPUi.append([])
      listPPUSi.append([])
      listPPUTi.append([])
      for j in range(beam._nip):
        listLREGIONi[-1].append(lregionall[j])
        listPLREGIONi[-1].append(plregionall[j])
        listLTIMEi[-1].append(ltimeall[j])
        listPLTIMEi[-1].append(pltimeall[j])
        listPUi[-1].append(puall[j])
        listPPUSi[-1].append(ppusall[j])
        listPPUTi[-1].append(pputall[j])
      #
      listXIXi.append([])
      listXIYi.append([])
      listXIMi.append([])
      for j in range(beam._nip):
        listXIXi[-1].append(xixall[j])
        listXIYi[-1].append(xiyall[j])
        listXIMi[-1].append(ximall[j])
      
      # The total beam-beam tune-shifts:
      xixtot = sum(listXIXi[-1])
      xiytot = sum(listXIYi[-1])
      ximtot = sum(listXIMi[-1])
      print "\n", "* * *".center(58), "\n"
      print "*", "{0:36} {1:10}".format("Total BB tune-shift horizontal", "1"),        "%1.4e" %xixtot
      print "*", "{0:36} {1:10}".format("Total BB tune-shift vertical",   "1"),        "%1.4e" %xiytot
      print "*", "{0:36} {1:10}".format("Total BB tune-shift mean",       "1"),        "%1.4e" %ximtot
      
      # Save values into their corresponding lists
      listXIXTOT.append(xixtot)
      listXIYTOT.append(xiytot)
      listXIMTOT.append(ximtot)
      
      # Save values into their corresponding lists
      listLINTi.append([])
      for j in range(beam._nip):
        listLINTi[-1].append(lintiall[j])
      #
      listLINT.append(lint)
      
      # Write level file:
      print >> FILElevel, "{0:16}".format("%1.8f" %listTIME[-1]),
      print >> FILElevel, "{0:16}".format("%1.8e" %listPPB[-1]),
      print >> FILElevel, "{0:16}".format("%1.8e" %listSIGS[-1]),
      print >> FILElevel, "{0:16}".format("%1.8e" %listRMSSIGS[-1]),
      print >> FILElevel, "{0:16}".format("%1.8e" %listFWHM[-1]),
      print >> FILElevel, "{0:16}".format("%1.8e" %listDPP[-1]),
      for k in range(2): print >> FILElevel, "{0:16}".format("%1.8e" %listEPSN[-1][k]),
      for k in range(2): print >> FILElevel, "{0:16}".format("%1.8e" %listEPSNIBS[-1][k]),
      for k in range(2): print >> FILElevel, "{0:16}".format("%1.8e" %listEPSNSR[-1][k]),
      for k in range(2): print >> FILElevel, "{0:16}".format("%1.8e" %listEPSNOBS[-1][k]),
      for k in range(2): print >> FILElevel, "{0:16}".format("%1.8e" %listEPSNCC0[-1][k]),
      for k in range(2): print >> FILElevel, "{0:16}".format("%1.8e" %listEPSNCC1[-1][k]),
      #
      for k in range(3):  print >> FILElevel, "{0:16}".format("%1.8f" %(listTAUIBS[-1][k])),
      for k in range(2):  print >> FILElevel, "{0:16}".format("%1.8e" %(listTAUCC0[-1][k])),
      for k in range(2):  print >> FILElevel, "{0:16}".format("%1.8e" %(listTAUCC1[-1][k])),
      #
      for j in range(beam._nip):
        print >> FILElevel, "{0:16}".format("%1.8f" %(listBETXi[-1][j])),
        print >> FILElevel, "{0:16}".format("%1.8f" %(listBETYi[-1][j])),
        print >> FILElevel, "{0:16}".format("%1.8e" %listLUMITOTi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8f" %listREDUCi[-1][j]),
      #
      for j in range(beam._nip):
        print >> FILElevel, "{0:16}".format("%1.8e" %listPHIi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8f" %listSEPLRi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8e" %listPHICRi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8e" %listONCCi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8e" %listPHICKi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8e" %listONCKi[-1][j]),
      print >> FILElevel, "{0:16}".format("%1.8e" %listPARSEP2[-1]),
      #
      for j in range(beam._nip):
        print >> FILElevel, "{0:16}".format("%1.8f" %listLREGIONi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8f" %listPLREGIONi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8e" %listLTIMEi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8e" %listPLTIMEi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8f" %listPUi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8f" %listPPUSi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8f" %listPPUTi[-1][j]),
      #
      for j in range(beam._nip):
        print >> FILElevel, "{0:16}".format("%1.8e" %listXIXi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8e" %listXIYi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8e" %listXIMi[-1][j]),
      #
      print >> FILElevel, "{0:16}".format("%1.8e" %listXIXTOT[-1]),
      print >> FILElevel, "{0:16}".format("%1.8e" %listXIYTOT[-1]),
      print >> FILElevel, "{0:16}".format("%1.8e" %listXIMTOT[-1]),
      #
      for j in range(beam._nip):
        print >> FILElevel, "{0:16}".format("%1.8f" %(listLINTi[-1][j]*1e-39)),
      #
      print >> FILElevel, "{0:16}".format("%1.8f" %(listLINT[-1]*1e-39)),
      
      # Special interest in the integrated luminosity at 6 h:
      if (i < 6*3600):
        pointbef6h = [i, lint*1.0e-39]
      elif (i == 6*3600):
        lumifill6h = lint*1.0e-39  # In [fb-1 yr-1]
      else:
        if (flagpointaft6h == False):
          pointaft6h = [i, lint*1.0e-39]
          flagpointaft6h = True
      
      #%%%%%%%%%%%%%
      
      rate = 0.0
      if (typestep != "PENAL"): # LEVEL steps (including FIRST and LAST)
        
        print "\n", "* * *".center(58)
        print "\n> Compute the time it takes the luminosity at ip0 to decay and move there...\n"
        
        # Loss rate and time decaying, which depends on the requested luminosity step. Two extreme cases:
        # 1) Keep 100% of luminosity at all times - this cannot be acomplished (number of steps tends to infinity), therefore fixed steps are performed instead.
        # 2) Let the luminosity decay util it reaches 0% - this cannot be acomplished (luminosity approaces asymptotically to zero),therefore the maximum fill length is taken. 
        for j in range(beam._nip):
          rate = rate + lumiall[j]*(self._xsecburn*1.0e-27) / beam._nbunch[j]
        print "*", "{0:36} {1:10}".format("Loss rate", "s-1"),                         "%1.6e" %rate
        print "*", "{0:36} {1:10}".format("Luminosity step for levelling", "1"),       "%1.4f" %self._p
        if (self._p != 0.0):
          if (self._p != 1.0):
            timedecay  = (np.sqrt(1.0/self._p) - 1.0) * beam._ppb / rate
            print "*", "{0:36} {1:10}".format("Time decaying", "s"),                   "%1.2f" %timedecay
          else:
            timedecay = self._step
            print ""
            print 2*" ", "[!] A 100% of the luminosity corresponds to not doing leveling"
            print 6*" ", "(otherwise it would stay in an infinite loop of zero-length step)."
            print 6*" ", "Fixed steps of", timedecay, "s will be performed instead."
        else:
          timedecay = self._maxfill
          print ""
          print 2*" ", "[!] A 0% of the luminosity can not be reached analytically,"
          print 6*" ", "since the ppb decay tends to it, but it does not reach it."
          print 6*" ", "The maximum time of", timedecay, "s will be used instead."
        steplength = timedecay
        
        # Move forward virtually 
        print "\n> We move from current time =", "%1.2f" %i, "s by a decay time =", "%1.2f" %timedecay, "s to", "%1.2f" %(i+timedecay),"s (virtually)...\n"
        textbanner = "VIRTUAL STEP: " + str(countsteps+1) + "* (" + "%1.2f" %(i+timedecay) + " s)"
        bannersclass.header5(textbanner)
        print ""
        print "*", "{0:36} {1:10}".format("Time", "h"),                                "%1.4f" %((i+timedecay)/3600.)
        print ""
        
        # Save a new branch for the step-wise luminosity function
        print >> FILElumiplot, "x <=", str(i+timedecay)+"/SECtoHR", "?", lumiall[0], "/ ( 1 + ( x -", str(i)+"/SECtoHR", ")*", str(rate)+"*SECtoHR", "/", beam._ppb, ")**2", ": "  ,
        
        # Integrated luminosity for each IP
        #print "\n\nLINT\n"
        for j in range(beam._nip):
          #print "lumiall[", j, "] =", lumiall[j]
          integrand = lambda t: lumiall[j] / (1 + (t-i)*rate/beam._ppb)**2
          linti = quad(integrand, i, i+timedecay)
          lintiall[j] = lintiall[j] + linti[0]
        # print "\nwith integrand =  lumiall[j] / (1 + (t-i)*rate/beam._ppb)**2 = lumiall[j] / (1 + (t-", i, ") *", "1.6e" %rate, "/", beam._ppb, ")**2,\nfrom i =", i, "to i+timedecay =", i, "+", timedecay, "=", i+timedecay, "\nThen\n"
        # for j in range(beam._nip):
        #   print "lintiall[", j, "] =", lintiall[j]
        # print "\n"
        
        # We compute the contribution to integrated luminosity from this step
        lumiint = lumiall[0]*beam._ppb/rate * ( 1 - 1/(1 + timedecay*rate/beam._ppb) )
        
        # We add up the value to the (cumulative) integrated luminosity from previous steps:
        factnow = ltot/ (i + timedecay + self._turnar*3600.)  # [yr-1]
        factbef = ltot/ (i             + self._turnar*3600.)  # [yr-1]
      
      else: # PENAL steps
        
        print "\n", "* * *".center(58)
        print "\n> Simulate a penalty step...\n"
        
        # Penalty step length
        print "*", "{0:36} {1:10}".format("Time in penalty step", "s"),                "%1.2f" %self._timepenstp
        steplength = self._timepenstp
        
        # Move forward virtually
        print "\n> We move from current time =", "%1.2f" %i, "s by timepenstp =", "%1.2f" %self._timepenstp, "s to", "%1.2f" %(i+self._timepenstp),"s (virtually)...\n"
        textbanner = "VIRTUAL STEP: " + str(countsteps+1) + "* (" + "%1.2f" %(i+self._timepenstp) + " s)"
        bannersclass.header5(textbanner)
        print "*", "{0:36} {1:10}".format("Time", "h"),                                "%1.4f" %((i+self._timepenstp)/3600.)
        
        # Save a new branch for the step-wise luminosity function
        print >> FILElumiplot, "x <", str(i+self._timepenstp)+"/SECtoHR", "?", "0.0", ": ",
        
        # Integrated luminosity for each IP is the same to previous step since there's no new contribution
        # lintiall = lintiall 
        
        # Integrated luminosity for each IP
        #print "\n\nLINT\n"
        for j in range(beam._nip):
          #print "lumiall[", j, "] =", lumiall[j]
          #integrand = lambda t: 0.0
          #linti = quad(integrand, i, i+self._timepenstp)
          #lintiall[j] = linti[0]
          lintiall[j] = lintiall[j] + 0.0
        # print "\nwith integrand = 0,\nfrom i =", i, "to i+timepenstp =", i, "+", self._timepenstp, "=", i+self._timepenstp, "\nThen\n"
        # for j in range(beam._nip):
        #   print "lintiall[", j, "] =", lintiall[j]
        # print "\n"
        
        # Null contribution to integrated luminosity from this step
        lumiint = 0.0   
        
        # We add up the value to the (cummulative) integrated luminosity from previous steps:
        if (i <= 0.0):
          factnow = 1.0
          factbef = 1.0
        else:
          factnow = ltot/ (i + self._timepenstp)  # [yr-1]
          factbef = ltot/ (i)                     # [yr-1]
      
      # The computation of integrated luminosity goes one step ahead:
      # print "\nBefore reasignment"
      # print "listLINTi   =", listLINTi
      # print "lintiall    =", lintiall
      # print "lintialltmp =", lintialltmp
      print "*", "{0:36} {1:10}".format("Integrated luminosity per IP (tmp)", "cm-2"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.6e" %lintialltmp[j]),
      print ""
      print "*", "{0:36} {1:10}".format("Integrated luminosity per IP (tmp)", "fb-1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %(lintialltmp[j]*1e-39)),
      print ""
      
      lintialltmp = []
      for j in range(len(lintiall)):
        lintialltmp.append(lintiall[j])
      # print "\nAfter reasignment"
      # print "listLINTi   =", listLINTi
      # print "lintiall    =", lintiall
      # print "lintialltmp =", lintialltmp
      # print "\n"
      print "*", "{0:36} {1:10}".format("Integrated luminosity per IP", "cm-2"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.6e" %lintiall[j]),
      print ""
      print "*", "{0:36} {1:10}".format("Integrated luminosity per IP", "fb-1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %(lintiall[j]*1e-39)),
      print ""
      
      print "\n"
      print "*", "{0:36} {1:10}".format("Integrated luminosity (this step)", "cm-2"), "%1.6e" %lumiint
      
      if rate != 0.0: llifetime = beam._ppb/2./rate/3600.  # [h]
      else:           llifetime = 0.0
      
      # Save values into their corresponding lists, and write level file:
      listRATE.append(rate)
      listLLIFETIME.append(llifetime)
      listSTEPLENGTH.append(steplength)
      listLUMIINTSTEP.append(lumiint)
      print >> FILElevel, "{0:16}".format("%1.8e" %listRATE[-1]),
      print >> FILElevel, "{0:16}".format("%1.8f" %listLLIFETIME[-1]),
      print >> FILElevel, "{0:16}".format("%1.8f" %listSTEPLENGTH[-1]),
      print >> FILElevel, "{0:16}".format("%1.8e" %listLUMIINTSTEP[-1]),
      
      # The computation of integrated luminosity goes one step ahead:
      linttmp = factnow*(linttmp/factbef) + factnow*lumiint
      linttmpall.append(linttmp)
      print "*", "{0:36} {1:10}".format("Integrated luminosity (tmp)", "cm-2 yr-1"), "%1.6e" %linttmp, 4*" ", "{0:10}".format("fb-1 yr-1"), "%1.4f" %(linttmp*1e-39)
      print "*", "{0:36} {1:10}".format("Integrated luminosity",       "cm-2 yr-1"), "%1.6e" %lint,    4*" ", "{0:10}".format("fb-1 yr-1"), "%1.4f" %(lint*1e-39)
      
      # A simple value to see how far we are from the optimum integrated luminosity
      print "*", "{0:36} {1:10}".format("Fraction of integrated luminosity", "%"),   "%1.2f" %(100.*linttmpall[-3]/linttmp)
      print ""
      
      # Integrated pile-up and effective pile-up densities:
      puintall = []
      eppusall = []
      epputall = []
      for j in range(beam._nip):
        intpus2 = self.GetIntSPUS2( beam, j, listPUi[-1][j]) /1e3 # [event^2/mm]
        intput2 = self.GetIntCTPUT2(beam, j, listPUi[-1][j]) /1e9 # [event^2/ns]
        oldeppu = False
        if oldeppu == True:
          # Old: (only for Gaussian, it starts at second step)
          if i == 0.0:
            puint = 0.0
            eppus = 0.0
            epput = 0.0
          else:
            puint = listPUINTi[-1][j] + listSTEPLENGTH[-2]*listPUi[-2][j]  # [event s]
            eppus = ( listEPPUSi[-1][j]*listPUINTi[-1][j] + listSTEPLENGTH[-2]* listPUi[-2][j]*listPPUSi[-2][j]/np.sqrt(2.) )/puint # [event/mm]
            epput = ( listEPPUTi[-1][j]*listPUINTi[-1][j] + listSTEPLENGTH[-2]* listPUi[-2][j]*listPPUTi[-2][j]/np.sqrt(2.) )/puint # [event/ns]
          print 2*"*\n"
          print "puint =", puint
          print "eppus =", eppus, " ( from intpus2 =", intpus2, ")"
          print "epput =", epput, " ( from intput2 =", intput2, ")"
          print 2*"*\n"
        else:
          if i == 0.0:
            puint =   listSTEPLENGTH[-1]*listPUi[-1][j]    # [event s]
            eppus = ( listSTEPLENGTH[-1]*intpus2 ) / puint # [event/mm]
            epput = ( listSTEPLENGTH[-1]*intput2 ) / puint # [event/ns]
          else:
            puint =                     listPUINTi[-1][j] + listSTEPLENGTH[-1]*listPUi[-1][j]
            eppus = ( listEPPUSi[-1][j]*listPUINTi[-1][j] + listSTEPLENGTH[-1]*intpus2 ) / puint
            epput = ( listEPPUTi[-1][j]*listPUINTi[-1][j] + listSTEPLENGTH[-1]*intput2 ) / puint
        puintall.append(puint)
        eppusall.append(eppus)
        epputall.append(epput)
      
      listPUINTi.append([])
      listEPPUSi.append([])
      listEPPUTi.append([])
      for j in range(beam._nip):
        listPUINTi[-1].append(puintall[j])
        listEPPUSi[-1].append(eppusall[j])
        listEPPUTi[-1].append(epputall[j])
      
      for j in range(beam._nip):
        print >> FILElevel, "{0:16}".format("%1.8e" %listPUINTi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8f" %listEPPUSi[-1][j]),
        print >> FILElevel, "{0:16}".format("%1.8f" %listEPPUTi[-1][j]),
      
      print "*", "{0:36} {1:10}".format("Integrated pile-up", "s"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4e" %listPUINTi[-1][j]),
      print ""
      print "*", "{0:36} {1:10}".format("Eff. peak pile-up s-density", "mm-1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %listEPPUSi[-1][j]),
      print ""
      print "*", "{0:36} {1:10}".format("Eff. peak pile-up t-density", "ns-1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %listEPPUTi[-1][j]),
      print ""
      
      # Save number and type of step
      listSTEP.append(countsteps)
      listTYPESTEP.append(typestep)      
      print >> FILElevel, "{0:16}".format("%1.0f" %listSTEP[-1]),
      print >> FILElevel, listTYPESTEP[-1]
      
      # We now prepare formally for the new step: increase time and assign the corresponding type of step label
      if (typestep != "PENAL"):
        i = i + timedecay
        if (self._penstp == True):
          typestep = "PENAL"
        else:
          if flagendlevel == False:
            typestep = "LEVEL"
          else:
            typestep = "DECAY"
      else:
        i = i + self._timepenstp
        if (flagendlevel == False):
          typestep = "LEVEL"
        else:
          typestep = "DECAY"
      countsteps = countsteps + 1
      
      # Update condition
      if self._optimumfill == True: cond1 = linttmp > linttmpall[-3] or linttmp == 0.0
      else:                         cond1 = not self._optimumfill
      
      # Close the levelling and lumiplot files every step so the appending is done line by line
      FILElevel.close()
      FILElumiplot.close()
      
    bannersclass.header1("END OF FILL")
    
    listTYPESTEP[-2] = "END"
    
    # Restore the value of the time in penalty step for printting purpoises
    self._timepenstp = timepenstptmp
    
    # Finish the step-wise luminosity function and close the file
    FILElumiplot = open(NAMElumiplot, "a")
    print >> FILElumiplot, "NaN"
    FILElumiplot.close()
    
    # Remove file that stores IBS results per step:
    try:            os.remove(NAMEupdateIBS)
    except OSError: pass
    
    # Count each type of steps:
    auxfuncclass = AuxFunc()
    nls = auxfuncclass.CountStringInList("LEVEL", listTYPESTEP) + auxfuncclass.CountStringInList("FIRST", listTYPESTEP) + auxfuncclass.CountStringInList("LAST", listTYPESTEP) + auxfuncclass.CountStringInList("PPUSL", listTYPESTEP)
    nps = auxfuncclass.CountStringInList("PENAL", listTYPESTEP)
    nds = auxfuncclass.CountStringInList("DECAY", listTYPESTEP)
    nts = nls + nps + nds
    
    # Levelling step at IP8:
    
    if self._steplev2 == None:
      self._steplev2 = nts
      self._timesteplev2 = listTIME[-1]
      
    # Some list having of parameters for each IP separately
    listLREGIONsep  = []
    listPLREGIONsep = []
    listLTIMEsep    = []
    listPLTIMEsep   = []
    listPUsep       = []
    listPPUSsep     = []
    listPPUTsep     = []
    listPUINTsep    = []
    listEPPUSsep    = []
    listEPPUTsep    = []
    listLINTsep     = []
    for j in range(beam._nip):
      listLREGIONsep.append([])
      listPLREGIONsep.append([])
      listLTIMEsep.append([])
      listPLTIMEsep.append([])
      listPUsep.append([])
      listPPUSsep.append([])
      listPPUTsep.append([])
      listPUINTsep.append([])
      listEPPUSsep.append([])
      listEPPUTsep.append([])
      listLINTsep.append([])
      for k in range(len(listPPUSi)):
        listLREGIONsep[j].append(listLREGIONi[k][j])
        listPLREGIONsep[j].append(listPLREGIONi[k][j])
        listLTIMEsep[j].append(listLTIMEi[k][j])
        listPLTIMEsep[j].append(listPLTIMEi[k][j])
        listPUsep[j].append(listPUi[k][j])
        listPPUSsep[j].append(listPPUSi[k][j])
        listPPUTsep[j].append(listPPUTi[k][j])
        listPUINTsep[j].append(listPUINTi[k][j])
        listEPPUSsep[j].append(listEPPUSi[k][j])
        listEPPUTsep[j].append(listEPPUTi[k][j])
        listLINTsep[j].append(listLINTi[k][j])
    
    # Maximum value of each of the parameters from the lists above
    maxlregion  = []
    maxplregion = []
    maxltime    = []
    maxpltime   = []
    maxpu       = []
    maxppus     = []
    maxpput     = []
    for j in range(beam._nip):
      maxlregion.append(max(listLREGIONsep[j]))
      maxplregion.append(max(listPLREGIONsep[j]))
      maxltime.append(max(listLTIMEsep[j]))
      maxpltime.append(max(listPLTIMEsep[j]))
      maxpu.append(max(listPUsep[j]))
      maxppus.append(max(listPPUSsep[j]))
      maxpput.append(max(listPPUTsep[j]))
    
    # Maximum total beam-beam tune-shifts:
    maxxixtot = max(listXIXTOT)
    maxxiytot = max(listXIYTOT)
    maxximtot = max(listXIMTOT)
    
    FILEtable = open(NAMEtable, "a")
    
    # We are interested in the step before the last levelling step (where the betas/crab kissing angle have/has been reseted to their minimum/zero. This is the last step where the lumininosity or peak pile-up density reached the requested levelled value. We re-define this as the "last" step
    
    if nls < 2:
        print "\n[!] LEVELLING NOT POSSIBLE. DECAY STARTS FROM FIRST STEP.\n" 
        nls             = 0
        lastbutoneindex = 0
        lastindex       = 0
    else:
      if self._penstp == True:
          lastbutoneindex = 2*(nls-2)
          lastindex       = 2*(nls-1)
      else:
          lastbutoneindex = nls-2
          lastindex       = nls-1
    endindex        = nts
    
    for k in [0, lastbutoneindex, lastindex, endindex]:
      
      if k == 0:
        print "PARAMETERS AT FIRST STEP OF FILL/LEVELLING"
        firstlast = "first"
      elif k == lastbutoneindex:
        print "PARAMETERS AT LAST-BUT-ONE STEP OF LEVELLING"
        firstlast = "lastbutone"
      elif k == lastindex:
        print "PARAMETERS AT LAST STEP OF LEVELLING"
        firstlast = "last"
      else:
        print "PARAMETERS AT END OF FILL"
        firstlast = "end"
        
      print ""
      print "*", "{0:36} {1:10}".format("Type of step", "-"),                          listTYPESTEP[k]
      print "*", "{0:36} {1:10}".format("Time", "h"),                                  "%1.4f" %listTIME[k]
      print "*", "{0:36} {1:10}".format("Intensity (ppb)", "1"),                       "%1.4e" %listPPB[k]
      print "*", "{0:36} {1:10}".format("Step length", "s"),                           "%1.2f" %listSTEPLENGTH[k]
      print "*", "{0:36} {1:10}".format("Energy spread (initial)", "1"),               "%1.4e" %listDPP[k]
      print "*", "{0:36} {1:10}".format("Bunch length (initial)", "m"),                "%1.4e" %listSIGS[k]
      print "*", "{0:36} {1:10}".format("Normalized emittance horizontal", "m"),       "%1.4e" %listEPSN[k][0]
      print "*", "{0:36} {1:10}".format("Normalized emittance vertical", "m"),         "%1.4e" %listEPSN[k][1]
      
      print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "time",       "%le"), listTIME[k]
      print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "ppb",  "%le"),       listPPB[k]
      print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "steplength", "%le"), listSTEPLENGTH[k]
      print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "dpp",        "%le"), listDPP[k]
      print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "sigs",       "%le"), listSIGS[k]
      print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "epsxn",      "%le"), listEPSN[k][0]
      print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "epsyn",      "%le"), listEPSN[k][1]
      
      print "*", "{0:36} {1:10}".format("Beta* horizontal", "m"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.8f" %listBETXi[k][j]),
      print ""
      print "*", "{0:36} {1:10}".format("Beta* vertical", "m"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.8f" %listBETYi[k][j]),
      print ""
      print "*", "{0:36} {1:10}".format("Crossing angle", "rad"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4e" %listPHIi[k][j]),
      print ""
      print "*", "{0:36} {1:10}".format("BB-LR separation", "sigma"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %listSEPLRi[k][j]),
      print ""
      print "*", "{0:36} {1:10}".format("Crab cavity angle", "rad"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4e" %listPHICRi[k][j]),
      print ""
      print "*", "{0:36} {1:10}".format("Crab cavity ON", "1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %listONCCi[k][j]),
      print ""
      print "*", "{0:36} {1:10}".format("Crab kissing angle", "rad"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4e" %listPHICKi[k][j]),
      print ""
      print "*", "{0:36} {1:10}".format("Crab kissing ON", "1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %listONCKi[k][j]),
      print ""
      #print 'self._levvar =', self._levvar
      #print 'j = ', j
      #print 'self._levvar[j][0]', self._levvar[j][0]
      #print 'self._steplev2 =', self._steplev2
      #print 'self._flagsteplev2 =', self._flagsteplev2
      #print ''
      if self._levvar[j][0] == "ParSep":
        print "*", "{0:36} {1:10}".format("Parallel separation", "m"),
        print "%1.6e" %beam._parsep
      #else:
      #  print beam._parsep
      print ""
      print "*", "{0:36} {1:10}".format("Luminosity", "cm-2 s-1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.6e" %listLUMITOTi[k][j]),
      print "\n"
      #
      print "*", "{0:36} {1:10}".format("RMS luminous region", "mm"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %(1e3*listLREGIONi[k][j])),
      print ""
      print "*", "{0:36} {1:10}".format("Peak luminous region", "mm-1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %(1e-3*listPLREGIONi[k][j])),
      print ""
      print "*", "{0:36} {1:10}".format("RMS luminous time", "ns"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %(1e9*listLTIMEi[k][j])),
      print ""
      print "*", "{0:36} {1:10}".format("Peak luminous time", "ns-1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %(1e-9*listPLTIMEi[k][j])),
      print ""
      print "*", "{0:36} {1:10}".format("Pile-up ", "1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %listPUi[k][j]),
      print ""
      print "*", "{0:36} {1:10}".format("Peak pile-up s-density", "mm-1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %listPPUSi[k][j]),
      print ""
      print "*", "{0:36} {1:10}".format("Peak pile-up t-density", "ns-1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %listPPUTi[k][j]),
      print "\n"
      #
      print "*", "{0:36} {1:10}".format("Integrated luminosity per IP", "cm-2"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.6e" %listLINTi[k][j]),
      print ""
      print "*", "{0:36} {1:10}".format("Integrated luminosity per IP", "fb-1"),
      for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %(listLINTi[k][j]*1e-39)),
      print "\n"
      
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "betx"     + str(j), "%le"), listBETXi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "bety"     + str(j), "%le"), listBETYi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "phi"      + str(j), "%le"), listPHIi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "sepLR"    + str(j), "%le"), listSEPLRi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "phiCR"    + str(j), "%le"), listPHICRi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "oncc"     + str(j), "%le"), listONCCi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "phiCK"    + str(j), "%le"), listPHICKi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "onck"     + str(j), "%le"), listONCKi[k][j]
      print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "parsep2", "%le"), listPARSEP2[k]
      #
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "lumi"     + str(j), "%le"), listLUMITOTi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "lregion"  + str(j), "%le"), listLREGIONi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "plregion" + str(j), "%le"), listPLREGIONi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "ltime"    + str(j), "%le"), listLTIMEi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "pltime"   + str(j), "%le"), listPLTIMEi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "pu"       + str(j), "%le"), listPUi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "ppus"     + str(j), "%le"), listPPUSi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "pput"     + str(j), "%le"), listPPUTi[k][j]
      #
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "xix"     + str(j), "%le"), listXIXi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "xiy"     + str(j), "%le"), listXIYi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "xim"     + str(j), "%le"), listXIMi[k][j]
      #
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "lintcms"  + str(j), "%le"), listLINTi[k][j]
      for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format(firstlast + "lintfb"   + str(j), "%le"), listLINTi[k][j]*1.0e-39
      
    print ""
    
    #
    
    print "FINAL RESULTS"
    print ""
    print "*", "{0:36} {1:10}".format("Fill time", "h"),                               "%1.4f" %listTIME[-1] # (listTIME[-1]+listSTEPLENGTH[-1])
    print "*", "{0:36} {1:10}".format("Leveling time", "h"),                           "%1.4f" %listTIME[lastbutoneindex]
    print "*", "{0:36} {1:10}".format("Start peak pile-up s-density lev.", "h"),
    if self._timeppuslev != None:
      print "%1.4f" %self._timeppuslev
    else:
      print self._timeppuslev
    print "*", "{0:36} {1:10}".format("Integrated luminosity", "cm-2 s-1"),            "%1.6e" %listLINT[-1]
    print "*", "{0:36} {1:10}".format("Integrated luminosity", "fb-1 yr-1"),           "%1.4f" %(listLINT[-1]*1.0e-39)
    print "*", "{0:36} {1:10}".format("Number of LEVEL steps", "1"),                   nls
    print "*", "{0:36} {1:10}".format("Number of PENALTY steps", "1"),                 nps
    print "*", "{0:36} {1:10}".format("Number of DECAY steps", "1"),                   nds
    print "*", "{0:36} {1:10}".format("Total number of steps", "1"),                   nts
    print "*", "{0:36} {1:10}".format("Last-but-one step index", "1"),                 lastbutoneindex
    print "*", "{0:36} {1:10}".format("Last step index", "1"),                         lastindex
    print ""
    print "*", "{0:36} {1:10}".format("Last step index (at IP8)", "1"),                self._steplev2
    print "*", "{0:36} {1:10}".format("Leveling time (at IP8) ", "h"),                 "%1.4f" %self._timesteplev2
    print ""
    
    print >> FILEtable, "@", "{0:22} {1:4}".format("filltime",        "%le"),           listTIME[-1] # +listSTEPLENGTH[-1]
    print >> FILEtable, "@", "{0:22} {1:4}".format("timeinleveling",  "%le"),           listTIME[lastbutoneindex]
    print >> FILEtable, "@", "{0:22} {1:4}".format("timeppuslev",     "%le"),           self._timeppuslev
    print >> FILEtable, "@", "{0:22} {1:4}".format("lintcms",         "%le"),           listLINT[-1]
    print >> FILEtable, "@", "{0:22} {1:4}".format("lintfb",          "%le"),           listLINT[-1]*1.0e-39
    print >> FILEtable, "@", "{0:22} {1:4}".format("nls",             "%d"),            nls
    print >> FILEtable, "@", "{0:22} {1:4}".format("nps",             "%d"),            nps
    print >> FILEtable, "@", "{0:22} {1:4}".format("nds",             "%d"),            nds
    print >> FILEtable, "@", "{0:22} {1:4}".format("nts",             "%d"),            nts
    print >> FILEtable, "@", "{0:22} {1:4}".format("lastbutoneindex", "%d"),            lastbutoneindex
    print >> FILEtable, "@", "{0:22} {1:4}".format("lastindex",       "%d"),            lastindex
    print >> FILEtable, "@", "{0:22} {1:4}".format("steplev2",        "%d"),            self._steplev2
    print >> FILEtable, "@", "{0:22} {1:4}".format("timesteplev2",    "%le"),           self._timesteplev2
    
    #
    
    print "*", "{0:36} {1:10}".format("Integrated pile-up", "s"),
    for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4e" %listPUINTi[-1][j]),
    print ""
    print "*", "{0:36} {1:10}".format("Eff. peak pile-up s-density", "mm-1"),
    for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %listEPPUSi[-1][j]),
    print ""
    print "*", "{0:36} {1:10}".format("Eff. peak pile-up t-density", "ns-1"),
    for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %listEPPUTi[-1][j]),
    print ""
    
    for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format("puint" + str(j),  "%le"), listPUINTi[-1][j]
    for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format("eppus" + str(j),  "%le"), listEPPUSi[-1][j]
    for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format("epput" + str(j),  "%le"), listEPPUTi[-1][j]
    
    #
    
    print "*", "{0:36} {1:10}".format("Maximum RMS luminous region", "mm"),
    for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %(1e3*maxlregion[j])),
    print ""
    print "*", "{0:36} {1:10}".format("Maximum peak luminous region", "mm-1"),
    for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %(1e-3*maxplregion[j])),
    print ""
    print "*", "{0:36} {1:10}".format("Maximum RMS luminous time", "ns"),
    for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %(1e9*maxltime[j])),
    print ""
    print "*", "{0:36} {1:10}".format("Maximum peak luminous time", "ns-1"),
    for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %(1e-9*maxpltime[j])),
    print ""
    print "*", "{0:36} {1:10}".format("Maximum pile-up", "1"),
    for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %maxpu[j]),
    print ""
    print "*", "{0:36} {1:10}".format("Maximum peak pile-up s-density", "mm-1"),
    for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %maxppus[j]),
    print ""
    print "*", "{0:36} {1:10}".format("Maximum peak pile-up t-density", "ns-1"),
    for j in range(beam._nip): print beam._ipnames[j] + ":", "{0:16}".format("%1.4f" %maxpput[j]),
    print "\n"
    
    for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format("maxlregion" + str(j),  "%le"), maxlregion[j]
    for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format("maxplregion" + str(j), "%le"), maxplregion[j]
    for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format("maxltime" + str(j),    "%le"), maxltime[j]
    for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format("maxpltime" + str(j),   "%le"), maxpltime[j]
    for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format("maxpu" + str(j),       "%le"), maxpu[j]
    for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format("maxppus" + str(j),     "%le"), maxppus[j]
    for j in range(beam._nip): print >> FILEtable, "@", "{0:22} {1:4}".format("maxpput" + str(j),     "%le"), maxpput[j]
    
    print "*", "{0:36} {1:10}".format("Max. total BB tune-shift hor.", "1"),           "%1.4e" %maxxixtot
    print "*", "{0:36} {1:10}".format("Max. total BB tune-shift ver.", "1"),           "%1.4e" %maxxiytot
    print "*", "{0:36} {1:10}".format("Max. total BB tune-shift mean", "1"),           "%1.4e" %maxximtot
    print ""
    print "*", "{0:36} {1:10}".format("Energy spread (minimum)", "1"),                 "%1.4e" %beam._mindpp
    print "*", "{0:36} {1:10}".format("Minimum bunch length reached", "T/F"),          beam._flagminlong
    print ""
    
    print >> FILEtable, "@", "{0:22} {1:4}".format("maxxixtot", "%le"),                maxxixtot
    print >> FILEtable, "@", "{0:22} {1:4}".format("maxxiytot", "%le"),                maxximtot
    #
    print >> FILEtable, "@", "{0:22} {1:4}".format("mindpp", "%le"),                   beam._mindpp
    print >> FILEtable, "@", "{0:22} {1:4}".format("flagminlong", "%b"),               beam._flagminlong
    
    FILEtable.close()
    
    # Results file:
    NAMEresults = "results_" + self._ID + ".out"
    print "> Writting '" + NAMEresults + "'..."
    FILEresults = open(NAMEresults, "w")
    FILEtable   = open(NAMEtable,   "r")
    FILElevel   = open(NAMElevel,   "r")
    for FILEtmp in [FILEtable, FILElevel]:
      for line in FILEtmp:
        FILEresults.write(line)
    FILElevel.close()
    FILEtable.close()
    FILEresults.close()

