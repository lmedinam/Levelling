import numpy as np
from Levelling_Others import Constants as cst

class Beam:
  
  #============================================================================================================#
  
  ### Beam parameters inherited from Config:
  
  _momeV         = 0.0
  _ppb           = 0.0
  _dpp           = 0.0
  _sigs          = 0.0
  _epsn          = []
  _beta          = []
  _betamindict   = {}
  _tau_obs       = []
  _tau_cc        = 0.0
  _kappa         = 0.0
  _kappac        = 0.0
  _t1            = 0.0
  _t2            = 0.0
  
  # Beam parameters NOT inherited:
  
  # _incroncc
  # _incronccrate
  # _niterlev
  # _IBSRF
  
  # _adaptivexsing
  
  #------------------------------------------------------------------------------------------------------------#
  
  ### Machine parameters inherited from Config:
  
  _opticsfile    = ""
  _circ          = 0.0
  _rho           = 0.0
  _alfmom        = 0.0
  _tuneb         = 0.0
  _tunes         = 0.0
  _Nbunch        = 0
  _ipnames       = []
  _nip           = 0
  _nbunch        = []
  _xplane        = []
  _sepLR         = []
  _sepLRsteptime = []
  _sepLRstep     = []
  _wcc           = 0.0
  _oncc          = []
  _tcc1          = []
  _tcc2          = []
  
  # Machine parameters NOT inherited:
  # (None)
  
  #------------------------------------------------------------------------------------------------------------#
  
  ### Levelling parameters inherited from Config:
  
  # (None)
  
  # Levelling parameters NOT inherited:
  # _levtech
  # _levvar
  # _levlumi
  # _levppus
  # _ck
  # _constbetar
  # _sepLRconst
  # _longdens  
  # _p
  # _step
  # _optimumfill
  # _maxfill
  # _penstp
  # _timepenstp
  # _updatepenstp
  
  #------------------------------------------------------------------------------------------------------------#
  
  ### Bunch length gymnastics parameters inherited from Config:
  
  _constlong     = True
  _ppblong       = 0.0
  _redlong       = 0.0
  _minlong       = 0.0
  
  # Bunch length gymnastics parameters NOT inherited:
  # (None)
  
  #------------------------------------------------------------------------------------------------------------#
  
  ### BB-LR Separation reduction parameters inherited from Config:
  
  _redsepLR      = []  # Keep crossing angle constant [T/F]
  _minsepLR      = []  # Minimum crossing angle [sigma]
  _ratesepLR     = []  # Rate for crossing angle reduction [sigma/h]
  
  ### Variable separation parameters NOT inherited:
  
  # (None)
  
  #------------------------------------------------------------------------------------------------------------#
  
  # Integrated luminosity parameters inherited from Config:
  
  # (None)
  
  # Integrated luminosity parameters NOT inherited:
  # _xsec
  # _xsecburn
  # _days
  # _eff
  # _turnar
  
  #============================================================================================================#
  
  ### New Machine parameters:

  _frev          = 0.0    # Revolution frequency [Hz]
  _betaver       = 0.0    # Average horizontal beta [m]
  _dx            = 0.0    # Average dispersion [m]
  _initphi       = []     # Initial crossing angle per IP [rad]
  _phi           = []     # Crossing angle per IP [rad]
  _initoncc      = []     # Initial crab cavity ON per IP [1]
  _initphiCR     = []     # initial crab cavity angle per IP [rad]
  _phiCR         = []     # Crab cavity  angle per IP [rad]
  _parsep        = 0.0    # Parallel separation for IP8 [m]
  _parseptmp     = 0.0    # Parallel separation for IP8 [m]
  
  #------------------------------------------------------------------------------------------------------------#
  
  ### New Levelling parameters:
  
  _phiCK         = []  # Crab kissing angle per IP [rad]
  _onck          = []
  
  #------------------------------------------------------------------------------------------------------------#
  
  ### New beam parameters:
  
  _gamma         = 0.0  # Relativistic gamma [1]
  _betarel       = 0.0  # Relativistic beta [1]
  _totalpart     = 0.0  # Total number of particles in the beam [1]
  _current       = 0.0  # Total current of the beam [A]
  _rmssigs       = 0.0
  _fwhm          = 0.0
  _tau_ibs       = []  # [Horizontal, vertical, longitudinal] emittance growth from IBS [h]
  _tau_sr        = []  # [Horizontal, vertical, longitudinal] radiation damping time [h]
  _dpp0          = 0.0  # ?
  _sigs0         = 0.0  # ?
  _eps0          = []   # ? [Horizontal, vertical]
  _initepsn      = []
  _epsn_0        = []
  _epsn_ibs      = []
  _epsn_sr       = []
  _epsn_obs      = []
  _epsn_cc_ip0   = []
  _epsn_cc_ip1   = []
  _tau_cc_mh     = []   # Emittance growth from CCs for beta* = 15 cm and phiCR = 380 urad (exact), in m/h (not um/h)
  _betamin       = []   # Array containing the minimum betas per IP.
  
  #------------------------------------------------------------------------------------------------------------#
  
  ### New bunch length gymnastics parameters:
  
  _initsigs      = 0.0
  _initdpp       = 0.0
  _mindpp        = 0.0
  _flagminlong   = False
  _constlong0    = True

  ### New BB-LR Separation reduction parameters:
  
  _initsepLR            = []

  #------------------------------------------------------------------------------------------------------------#

  ### New independent parameters:
  
  _identicalIP01 = False  # Identical ip = 0 and ip = 1? [T/F]
  
  #============================================================================================================#
  
  def __init__(self,config):
    
    # Inherited:
    
    self._momeV        = config._momeV
    self._ppb          = config._ppb
    self._dpp          = config._dpp
    self._sigs         = config._sigs
    self._epsn         = config._epsn
    self._beta         = np.array(config._beta)
    self._tau_obs      = config._tau_obs
    self._tau_cc       = config._tau_cc
    self._kappa        = config._kappa
    self._kappac       = config._kappac
    self._t1           = config._t1
    self._t2           = config._t2
    
    self._opticsfile    = config._opticsfile
    self._circ          = config._circ
    self._rho           = config._rho
    self._alfmom        = config._alfmom
    self._tuneb         = config._tuneb
    self._tunes         = config._tunes
    self._Nbunch        = config._Nbunch
    self._ipnames       = config._ipnames
    self._nip           = config._nip
    self._nbunch        = config._nbunch
    self._xplane        = config._xplane
    self._sepLR         = np.array(config._sepLR)
    self._sepLRsteptime = np.array(config._sepLRsteptime)
    self._sepLRstep     = np.array(config._sepLRstep)
    self._wcc           = config._wcc
    self._oncc          = config._oncc
    self._tcc1          = np.array(config._tcc1)
    self._tcc2          = np.array(config._tcc2)
    
    self._constlong     = config._constlong
    self._ppblong       = config._ppblong
    self._redlong       = config._redlong
    self._minlong       = config._minlong
    
    self._redsepLR      = config._redsepLR
    self._minsepLR      = config._minsepLR
    self._ratesepLR     = config._ratesepLR
  
    ### New:
    
    self._frev      = cst.clight / config._circ
    self._betaver   = config._circ/(2*np.pi) / config._tuneb
    self._dx        = config._circ/(2*np.pi) * config._alfmom
    
    self._gamma     = config._momeV/cst.pmass + 1.0
    self._betarel   = np.sqrt(1.0 - 1.0/(self._gamma**2.0))
    self._totalpart = self._Nbunch*self._ppb
    self._current   = self._totalpart * 1.815e-15
    
    self._initsepLR = []
    for i in range(self._nip):
      self._initsepLR.append(self._sepLR[i])
      
    self._initphi = []
    self._phi     = []
    for i in range(self._nip):
      diver = np.sqrt( self._epsn[self._xplane[i]] / (self._gamma*self._betarel) / self._beta[i][self._xplane[i]] )
      self._initphi.append(diver*self._sepLR[i])
      self._phi.append(self._initphi[i])
    
    self._initoncc  = []
    for i in range(self._nip):
      self._initoncc.append(self._oncc[i])
    
    self._initphiCR = []
    self._phiCR     = []
    for i in range(self._nip):
      self._initphiCR.append(self._initoncc[i]*self._initphi[i])
      self._phiCR.append(self._initphiCR[i])
  
    self._phiCK     = [0.0 for i in range(self._nip)]
    self._onck      = [0.0 for i in range(self._nip)]
    self._initphiCK = [0.0 for i in range(self._nip)]
    self._initonck  = [0.0 for i in range(self._nip)]
    
    self._identicalIP01 = (self._beta[0][0] == self._beta[1][0]) * (self._beta[0][1] == self._beta[1][1]) * (self._nbunch[0] == self._nbunch[1]) * (self._xplane[0] == self._xplane[1]) * (self._sepLR[0] == self._sepLR[1]) * (self._sepLRsteptime[0] == self._sepLRsteptime[1]) * (self._sepLRstep[0] == self._sepLRstep[1]) * (self._initoncc[0] == self._initoncc[1]) * (config._levtech[0] == config._levtech[1]) * (config._levlumi[0] == config._levlumi[1]) * (config._levppus[0] == config._levppus[1]) * (config._ck[0] == config._ck[1]) * (config._constbetar[0] == config._constbetar[1]) * (config._sepLRconst[0] == config._sepLRconst[1]) * (self._redsepLR[0] == self._redsepLR[1]) * (self._minsepLR[0] == self._minsepLR[1]) * (self._ratesepLR[0] == self._ratesepLR[1]) 
    
    self._parsep = 0.0
    self._parsep = self._parseptmp
    
    self._initsigs      = self._sigs
    self._initdpp       = self._dpp
    self._mindpp        = self._initdpp  # It needs to be redefined with the correct value once sigs reaches minlong
    self.flagminlong    = False
    self._constlong0    = config._constlong
    
    # hor, ver
    self._tau_ibs  = [0.0, 0.0, 0.0]
    if self._tau_cc == None: self._tau_cc_mh = None
    else:                    self._tau_cc_mh = self._epsn[0]/(self._tau_cc * 0.15/380e-6**2)
    self._tau_cc_0 = [0.0, 0.0] # ip0
    self._tau_cc_1 = [0.0, 0.0] # ip1
  
    dEsr = cst.e**2 * self._betarel**3 * self._gamma**4 / (3*cst.eps0*self._rho) / cst.e / self._momeV
    
    taux_sr = 2.0/(dEsr*self._frev*3600.0)
    tauy_sr = 2.0/(dEsr*self._frev*3600.0)
    tauz_sr = 1.0/(dEsr*self._frev*3600.0)
    
    print dEsr, taux_sr, tauz_sr
    quit()
    
    self._tau_sr = [taux_sr, tauy_sr, tauz_sr]

    I1 = self._dx*2*np.pi
    I2 = 2*np.pi/self._rho
    I3 = 2*np.pi/self._rho**2
    I4 = self._circ*self._alfmom/self._rho**2
    I5 = 1/self._betaver*self._dx**2/self._rho**2*2*np.pi
    
    cq = 55/(32*np.sqrt(3.0))*cst.hbar/(cst.m*cst.clight)
    
    self._dpp0  = np.sqrt(cq*self._gamma**2*I3/(2*I2))
    self._sigs0 = self._alfmom*cst.clight/(2*np.pi*self._tunes*cst.clight/self._circ)*self._dpp0
    epsx0 = cq*self._gamma**2*I5/I2*self._betarel*self._gamma
    epsy0 = 13.0/55.0*cq/I2*self._betaver/self._rho**2*2*np.pi*self._betarel*self._gamma
    self._eps0 = [epsx0, epsy0]
    
    # hor, ver
    self._initepsn = [self._epsn[0], self._epsn[1]]
    self._epsn_0   = [self._epsn[0], self._epsn[1]]
    self._epsn_ibs = [0.0, 0.0]
    self._epsn_sr  = [0.0, 0.0]
    self._epsn_obs = [0.0, 0.0]
    
    self._epsn_cc_ip0  = [0.0, 0.0]
    self._epsn_cc_ip1  = [0.0, 0.0]
  
    self._betamin   = np.array(config._beta)
  
  def getsigma(self, ip):
    """ Returns the horizontal and vertical beam sizes. """
    betx = self._beta[ip][0]
    bety = self._beta[ip][1]
    sigx = np.sqrt( self._epsn[0]*betx/(self._gamma*self._betarel) )
    sigy = np.sqrt( self._epsn[1]*bety/(self._gamma*self._betarel) )
    return sigx, sigy

  def getrmsfwhm(self, longdens):
    if longdens == "Gaussian":
      self._rmssigs = self._sigs
      self._fwhm    = 2.*self._sigs*np.sqrt(2.*np.log(2.))
    elif longdens == "qGaussian":
      self._rmssigs = self._sigs*cst.FactFWHMGauss / 4./np.sqrt(2.)
      self._fwhm    = 2.*self._sigs*np.sqrt(2.*np.log(2.))
    elif longdens == "RF800":
      self._rmssigs = self._sigs
      self._fwhm    = self._sigs*cst.FactRMSGauss * 2.*np.sqrt(np.sqrt(2.*np.log(2.)))
  
  def printBeamParam(self, longdens, tofile = False, foutname = None):
    """ Prints to the standard output or to file 'foutname' the Beam class parameters. """
    
    if tofile == False:
      
      print "BEAM PARAMETERS"
      print ""
      print "*", "{0:36} {1:10}".format("Energy", "eV"),                               "%1.4e" %self._momeV
      print "*", "{0:36} {1:10}".format("Intensity (ppb)", "1"),                       "%1.4e" %self._ppb
      
      self.getrmsfwhm(longdens)
      self._initrmssigs = self._rmssigs
      self._initfwhm    = self._fwhm
      
      print "*", "{0:36} {1:10}".format("Energy spread (initial)", "1"),               "%1.4e" %self._initdpp
      print "*", "{0:36} {1:10}".format("Gaussian RMS bunch length (initial)", "m"),   "%1.4e" %self._initsigs
      print "*", "{0:36} {1:10}".format("RMS bunch length (initial)", "m"),            "%1.4e" %self._initrmssigs
      print "*", "{0:36} {1:10}".format("FWHM (initial)", "m"),                        "%1.4e" %self._initfwhm
      
      print "*", "{0:36} {1:10}".format("Normalized emittance horizontal", "m"),       "%1.4e" %self._epsn[0]
      print "*", "{0:36} {1:10}".format("Normalized emittance vertical", "m"),         "%1.4e" %self._epsn[1]
      print "*", "{0:36} {1:10}".format("Beta* horizontal", "m"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.8f" %self._beta[i][0]),
      print ""
      print "*", "{0:36} {1:10}".format("Beta* vertical", "m"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.8f" %self._beta[i][1]),
      print ""

      sig = []
      print "*", "{0:36} {1:10}".format("Beam size horizontal", "m"),
      for i in range(self._nip):
        sig.append(self.getsigma(i))
        print self._ipnames[i] + ":", "{0:16}".format("%1.4e" %sig[i][0]),
      print ""
      print "*", "{0:36} {1:10}".format("Beam size vertical", "m"),
      for i in range(self._nip):
        print self._ipnames[i] + ":", "{0:16}".format("%1.4e" %sig[i][1]),
      print ""
      
      print "*", "{0:36} {1:10}".format("IBS kappa", "1"),                             "%1.2f" %self._kappa
      print "*", "{0:36} {1:10}".format("IBS kappac", "1"),                            "%1.2f" %self._kappac
      print "*", "{0:36} {1:10}".format("Bunch time delay beam 1", "s"),               "%1.2e" %self._t1
      print "*", "{0:36} {1:10}".format("Bunch time delay beam 2", "s"),               "%1.2e" %self._t2
      #
      print "*", "{0:36} {1:10}".format("Relativistic gamma", "1"),                    "%1.2f" %self._gamma
      print "*", "{0:36} {1:10}".format("Relativistic beta", "1"),                     "%1.12f" %self._betarel
      print "*", "{0:36} {1:10}".format("Total no. of particles in the beam", "1"),    "%1.4e" %self._totalpart
      print "*", "{0:36} {1:10}".format("Total current of the beam", "A"),             "%1.4e" %self._current
      print "*", "{0:36} {1:10}".format("Radiation damp. time horizontal", "h"),       "%1.4f" %self._tau_sr[0]
      print "*", "{0:36} {1:10}".format("Radiation damp. time vertical", "h"),         "%1.4f" %self._tau_sr[1]
      print "*", "{0:36} {1:10}".format("Radiation damp. time longitud.", "h"),        "%1.4f" %self._tau_sr[2]
      print "*", "{0:36} {1:10}".format("Norm. emit. growth time obs. hor.", "h"),
      if self._tau_obs[0] != None: print "%1.4e" %self._tau_obs[0]
      else:                        print self._tau_obs[0]
      print "*", "{0:36} {1:10}".format("Norm. emit. growth time obs. ver.", "h"),
      if self._tau_obs[1] != None: print "%1.4e" %self._tau_obs[1]
      else:                        print self._tau_obs[1]
      print "*", "{0:36} {1:10}".format("Norm. emit. growth time CC/CK", "h*m/1^2"),
      if self._tau_cc != None:     print "%1.4e" %self._tau_cc
      else:                        print self._tau_cc
      print "*", "{0:36} {1:10}".format("... at beta=15cm, phiCR=380 urad", "m/h"),
      if self._tau_cc != None:     print "%1.4e" %self._tau_cc_mh
      else:                        print self._tau_cc_mh
      print ""
      print ""
      
      #
      print "MACHINE PARAMETERS"
      print ""
      print "*", "{0:36} {1:10}".format("Circumference", "m"),                          "%1.4f" %self._circ
      print "*", "{0:36} {1:10}".format("Curvature radius", "m"),                       "%1.4f" %self._rho
      print "*", "{0:36} {1:10}".format("Momentum compaction", "1"),                    "%1.4e" %self._alfmom
      print "*", "{0:36} {1:10}".format("Horizontal betatron tune", "1"),               "%1.4e" %self._tuneb
      print "*", "{0:36} {1:10}".format("Synchrotron tune", "1"),                       "%1.4e" %self._tunes
      print "*", "{0:36} {1:10}".format("Total number of bunches", "1"),                self._Nbunch
      print "*", "{0:36} {1:10}".format("Number of IPs", "1"),                          self._nip
      print "*", "{0:36} {1:10}".format("IP names", "-"),
      for i in range(self._nip): print "ip" + str(i) + ":", repr(self._ipnames[i]).ljust(16),
      print ""
      print "*", "{0:36} {1:10}".format("Colliding pairs", "1"),
      for i in range(self._nip): print self._ipnames[i] + ":", repr(self._nbunch[i]).ljust(16),
      print ""
      print "*", "{0:36} {1:10}".format("Crossing plane", "H=0/V=1"),
      for i in range(self._nip): print self._ipnames[i] + ":", repr(self._xplane[i]).ljust(16),
      print ""
      print "*", "{0:36} {1:10}".format("BB-LR separation (initial)", "sigma"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.4f" %self._initsepLR[i]),
      print ""
      print "*", "{0:36} {1:10}".format("Trigger time for sepLR reduction", "s"), 
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.4f" %self._sepLRsteptime[i]),
      print ""
      print "*", "{0:36} {1:10}".format("Reduction step for sepLR", "sigma"), 
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.4f" %self._sepLRstep[i]),
      print ""
      print "*", "{0:36} {1:10}".format("Crab cavity frequency", "Hz"),                 "%1.4e" %self._wcc
      print "*", "{0:36} {1:10}".format("Time delay of horizontal CCs beam 1", "s"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.8f" %self._tcc1[i][0]),
      print ""
      print "*", "{0:36} {1:10}".format("Time delay of horizontal CCs beam 1", "s"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.8f" %self._tcc1[i][1]),
      print ""
      print "*", "{0:36} {1:10}".format("Time delay of horizontal CCs beam 2", "s"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.8f" %self._tcc2[i][0]),
      print ""
      print "*", "{0:36} {1:10}".format("Time delay of horizontal CCs beam 2", "s"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.8f" %self._tcc2[i][1]),
      print ""
      #
      print "*", "{0:36} {1:10}".format("Revolution frequency", "Hz"),                  "%1.4e" %self._frev    
      print "*", "{0:36} {1:10}".format("Average horizontal beta", "m"),                "%1.4f" %self._betaver
      print "*", "{0:36} {1:10}".format("Average dispersion", "m"),                     "%1.4f" %self._dx
      print "*", "{0:36} {1:10}".format("Crossing angle (initial)", "rad"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.4e" %self._initphi[i]),
      print ""
      print "*", "{0:36} {1:10}".format("Crab cavity angle (initial)", "rad"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.4e" %self._initphiCR[i]),
      print ""
      print "*", "{0:36} {1:10}".format("Crab cavity ON (initial)", "1"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.4f" %self._initoncc[i]),
      print ""
      print "*", "{0:36} {1:10}".format("Crab kissing angle (inital)", "rad"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.4e" %self._initphiCK[i]),
      print ""
      print "*", "{0:36} {1:10}".format("Crab kissing ON (initial)", "1"),
      for i in range(self._nip): print self._ipnames[i] + ":", "{0:16}".format("%1.4f" %self._initonck[i]),
      print "\n"
      #
      print "BUNCH LENGTH GYMNASTICS"
      print ""
      print "*", "{0:36} {1:10}".format("Bunch length const. throughout fill", "T/F"),  self._constlong
      print "*", "{0:36} {1:10}".format("Intensity that triggers shortening", "1"),     "%1.4e" %self._ppblong
      print "*", "{0:36} {1:10}".format("Fractional reduction for sigs, dpp", "1"),     "%1.4e" %self._redlong
      print "*", "{0:36} {1:10}".format("Minimum Gaussian bunch length", "m"),          "%1.4e" %self._minlong
      print ""
      #
      print "SEPARATION REDUCTION PARAMETERS"
      print ""
      print "*", "{0:36} {1:10}".format("sepLR reduction", "T/F"),
      for i in range(self._nip): print self._ipnames[i] + ":", repr(self._redsepLR[i]).ljust(16),
      print ""
      print "*", "{0:36} {1:10}".format("Minimum sepLR", "sigma"),
      for i in range(self._nip):
        print self._ipnames[i] + ":",
        if self._minsepLR[i] != None: print "{0:16}".format("%1.4f" %self._minsepLR[i]),
        else:                         print repr(self._minsepLR[i]).ljust(16),
      print ""
      print "*", "{0:36} {1:10}".format("sepLR reduction rate", "sigma/h"),
      for i in range(self._nip):
        print self._ipnames[i] + ":",
        if self._ratesepLR[i] != None: print "{0:16}".format("%1.4f" %self._ratesepLR[i]),
        else:                          print repr(self._ratesepLR[i]).ljust(16),
      print ""
      print ""
      #
      print "*", "{0:36} {1:10}".format("Identical ip = 0 and ip = 1", "T/F"),          self._identicalIP01
      print ""

    else:
      
      print >> foutname, "@", "{0:22} {1:4}".format("momeV", "%le"),                    self._momeV
      print >> foutname, "@", "{0:22} {1:4}".format("ppb", "%le"),                      self._ppb
      print >> foutname, "@", "{0:22} {1:4}".format("initdpp", "%le"),                  self._initdpp
      print >> foutname, "@", "{0:22} {1:4}".format("initsigs", "%le"),                 self._initsigs
      print >> foutname, "@", "{0:22} {1:4}".format("initrmssigs", "%le"),              self._initrmssigs
      print >> foutname, "@", "{0:22} {1:4}".format("initfwhm", "%le"),                 self._initfwhm
      
      print >> foutname, "@", "{0:22} {1:4}".format("epsxn", "%le"),                    self._epsn[0]
      print >> foutname, "@", "{0:22} {1:4}".format("epsyn", "%le"),                    self._epsn[1]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("betx" + str(i), "%le"),          self._beta[i][0]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("bety" + str(i), "%le"),          self._beta[i][1]
      
      sig = []
      for i in range(self._nip):
        sig.append(self.getsigma(i))
      
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("sigx" + str(i), "%le"),          sig[i][0]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("sigy" + str(i), "%le"),          sig[i][1]
      
      print >> foutname, "@", "{0:22} {1:4}".format("taux_sr", "%le"),                  self._tau_sr[0]
      print >> foutname, "@", "{0:22} {1:4}".format("tauy_sr", "%le"),                  self._tau_sr[1]
      print >> foutname, "@", "{0:22} {1:4}".format("tauz_sr", "%le"),                  self._tau_sr[2]
      print >> foutname, "@", "{0:22} {1:4}".format("taux_obs", "%le"),                 self._tau_obs[0]
      print >> foutname, "@", "{0:22} {1:4}".format("tauy_obs", "%le"),                 self._tau_obs[1]
      print >> foutname, "@", "{0:22} {1:4}".format("tau_cc", "%le"),                   self._tau_cc
      print >> foutname, "@", "{0:22} {1:4}".format("tau_cc_mh", "%le"),                self._tau_cc_mh
      print >> foutname, "@", "{0:22} {1:4}".format("kappa", "%le"),                    self._kappa
      print >> foutname, "@", "{0:22} {1:4}".format("kappac", "%le"),                   self._kappac
      print >> foutname, "@", "{0:22} {1:4}".format("t1", "%le"),                       self._t1
      print >> foutname, "@", "{0:22} {1:4}".format("t2", "%le"),                       self._t2
      #
      print >> foutname, "@", "{0:22} {1:4}".format("gamma", "%le"),                    self._gamma
      print >> foutname, "@", "{0:22} {1:4}".format("betarel", "%le"),                  self._betarel
      print >> foutname, "@", "{0:22} {1:4}".format("totalpart", "%le"),                self._totalpart
      print >> foutname, "@", "{0:22} {1:4}".format("current", "%le"),                  self._current

      #
      print >> foutname, "@", "{0:22} {1:4}".format("circ", "%le"),                     self._circ
      print >> foutname, "@", "{0:22} {1:4}".format("rho", "%le"),                      self._rho
      print >> foutname, "@", "{0:22} {1:4}".format("alfmom", "%le"),                   self._alfmom
      print >> foutname, "@", "{0:22} {1:4}".format("tuneb", "%le"),                    self._tuneb
      print >> foutname, "@", "{0:22} {1:4}".format("tunes", "%le"),                    self._tunes
      print >> foutname, "@", "{0:22} {1:4}".format("Nbunch", "%le"),                   self._Nbunch
      print >> foutname, "@", "{0:22} {1:4}".format("nip", "%d"),                       self._nip
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("ipnames" + str(i), "%s"),        self._ipnames[i]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("nbunch" + str(i), "%le"),        self._nbunch[i]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("xplane" + str(i), "%d"),         self._xplane[i]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("initsepLR" + str(i), "%le"),     self._sepLR[i]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("sepLRsteptime" + str(i), "%b"),  self._sepLRsteptime[i]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("sepLRstep" + str(i), "%le"),     self._sepLRstep[i]
      print >> foutname, "@", "{0:22} {1:4}".format("wcc", "%le"),                      self._wcc
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("tcc1x" + str(i), "%le"),         self._tcc1[i][0]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("tcc1y" + str(i), "%le"),         self._tcc1[i][1]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("tcc2x" + str(i), "%le"),         self._tcc2[i][0]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("tcc2y" + str(i), "%le"),         self._tcc2[i][1]
      #
      print >> foutname, "@", "{0:22} {1:4}".format("frev", "%le"),                     self._frev    
      print >> foutname, "@", "{0:22} {1:4}".format("betaver", "%le"),                  self._betaver
      print >> foutname, "@", "{0:22} {1:4}".format("dx", "%le"),                       self._dx
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("initphi" + str(i), "%le"),       self._initphi[i]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("initphiCR" + str(i), "%le"),     self._initphiCR[i]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("initoncc" + str(i), "%le"),      self._initoncc[i]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("initphiCK" + str(i), "%le"),     self._initphiCK[i]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("initonck" + str(i), "%le"),      self._initonck[i]
      print >> foutname, "@", "{0:22} {1:4}".format("constlong", "%b"),                 self._constlong
      print >> foutname, "@", "{0:22} {1:4}".format("ppblong", "%le"),                  self._ppblong
      print >> foutname, "@", "{0:22} {1:4}".format("redlong", "%le"),                  self._redlong
      print >> foutname, "@", "{0:22} {1:4}".format("minlong", "%le"),                  self._minlong
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("redsepLR" + str(i), "%s"),       self._redsepLR[i]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("minsepLR" + str(i), "%le"),      self._minsepLR[i]
      for i in range(self._nip): print >> foutname, "@", "{0:22} {1:4}".format("ratesepLR" + str(i), "%le"),     self._ratesepLR[i]
      #
      print >> foutname, "@", "{0:22} {1:4}".format("identicalIP01", "%b"),             self._identicalIP01


