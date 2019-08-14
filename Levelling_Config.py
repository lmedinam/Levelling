import numpy as np
from Levelling_Others import AuxFunc
from Levelling_Others import Constants as cst

auxfuncclass = AuxFunc()

class Config:
  
  ### Beam parameters:
  
  _momeV         = 0.0   # Momentum [eV]
  _ppb           = 0.0   # bunch intensity or popualtion (ppb) [1]
  _dpp           = 0.0   # Energy spread dp/p [1]
  _sigs          = 0.0   # Bunch length [m]
  _epsn          = []    # Normalized [horizontal, vertical] emittance [m]
  _beta          = []    # [betax, betay] per IP [m]
  _betamindict   = {}    # Dictionary storing a list of specific minimum betas for ip = 0 and ip = 1, as a function of a list of bunch intensities. Structure and units: { "intensities": [1, 1, ...], "betas": [[m, m], [m,m], ...] }
  _tau_obs       = []    # [Horizontal, vertical] normalized emittance growth time (observed)   [None or h]
  _tau_cc        = 0.0   # Inverse of the normalized emittance growth time due to CCs (for crab crossing and crab kissing) in the plane with CCs [None or %/h*(m/1^2)]
  _kappa         = 0.0   # ? [1] - No longer used since IBS is done by MAD-X and not approximate formula
  _kappac        = 0.0   # ? [1] - Still used to get tauy from taux
  _t1            = 0.0   # Time delay of bunches of beam 1 [s]
  _t2            = 0.0   # Time delay of bunches of beam 2 [s]
  
  _incroncc      = []    # Increase crabbing angle from an initial value up to a value equal to the crossing angle (oncc = 1) [rad]
  _incronccrate  = 0.0   # Fractional increment of oncc [1], e.g. = 0.005 = 0.12/24 will result in an increment of 0.12 (as in 0.88 to 1.00) in 24 steps (not counting penalty steps)
  _niterlev      = 0     # Number of iterations for dobule levelling [1]
  _IBSRF         = True  # Emittance growth due to IBS and SR [T/F]
  
  _adaptivexsing = []    # Adaptive crossing angle per IP, as function of ppb (to ensure 6sigma DA, by Nikos) [T/F]
  
  ### Machine parameters:

  _opticsfile    = ""    # Optics file ["opt_round.madx" or "opt_flat.madx"]
  _circ          = 0.0   # Circunference [m]
  _rho           = 0.0   # Curvature radius [m]
  _alfmom        = 0.0   # Momentum compaction [1]
  _tuneb         = 0.0   # Horizontal betatron tune [1]
  _tunes         = 0.0   # Synchrotron tune [1]
  _Nbunch        = 0     # Total number of bunches [1]
  _ipnames       = []    # IP names ["IPname1", ...]
  _nip           = 0     # Number of IPs, automatically determined from the length of ipnames [1]
  _nbunch        = []    # Colliding pairs per IP [1]
  _xplane        = []    # Crossing plane per IP [H=0,V=1]
  _sepLR         = []    # BB-LR separation per IP [sigma]
  _sepLRsteptime = []    # Trigger time for sepLR reduction per IP [s]
  _sepLRstep     = []    # Reduction step for sepLR per IP [sigma]
  _wcc           = 0.0   # Crab cavity frequency [Hz]
  _oncc          = []    # Crab cavity ON per IP [1]
  _tcc1          = []    # Time delay of the CCs for beam 1 in the horizontal and vertical planes per IP [s]
  _tcc2          = []    # Time delay of the CCs for beam 2 in the horizontal and vertical planes per IP [s]
  
  ### Levelling parameters: 
  
  _levtech       = []    # Levelling technique per IP ["levlumi" or "levppus"]
  _levvar        = []    # Levelling variable per IP ["Beta", "PhiCR", "PhiCK", or "ParSep"]
  _levlumi       = []    # Leveled luminosity per IP [cm-2 s-1]
  _levppus       = []    # Leveled peak pile-up per IP [1]
  _ck            = []    # CK in the separation plane per IP [T/F]
  _constbetar    = []    # Constant beta H/V ratio throughout levelling per IP [T/F]
  _sepLRconst    = []    # Crossing angle (in rad) constant and sepLR (in sigma) variable if False, Crossing angle (in rad) variable and sepLR (in sigma) constant if True. per IP [T/F]
  _longdens      = ""    # Longitudinal density ["Gaussian", "RF800", "qGaussian"]
  _p             = 0.0   # Luminosity step for levelling (Reduction in percentage of the luminosity that triggers a new levelling step), same for all IPs [1]
  _step          = 0.0   # Time step (in case p = 1.0) [s]
  _optimumfill   = True  # Run for optimum fill or a given maxfill time [T/F]
  _maxfill       = 0.0   # Maximum fill length (time) [s]
  _penstp        = True  # Add penalty steps (for example, time to align beams after squeeze) [T/F]
  _timepenstp    = 0.0   # Length of penalty steps [s]
  _updatepenstp  = False # Update parameters (emittances and sigz) in penalty step [T/F]
  
  ### Bunch lenght gymnastics parameters:
  
  _constlong     = True  # Keep bunch length or FWHM constant throughout the fill [T/F]
  _ppblong       = 0.0   # Intensity (ppb) that triggers bunch length sortening [1]
  _redlong       = 0.0   # Fractional reduction for sigs and dpp [1]
  _minlong       = 0.0   # Minimum bunch length or FWHM [m]

  ### BB-LR Separation reduction parameters:
  
  _redsepLR      = []    # sepLR reduction [T/F]
  _minsepLR      = []    # Minimum sepLR [sigma]
  _ratesepLR     = []    # sepLR reduction rate [sigma/h]
  
  ### Integrated luminosity parameters: 
  
  _xsec          = 0.0   # Cross section (inelastic) for pile-up [mb]
  _xsecburn      = 0.0   # Cross section (total) for burn-off [mb]
  _days          = 0.0   # Dedicated time to physics [days yr-1]
  _eff           = 0.0   # Efficiency [1]
  _turnar        = 0.0   # Turn-around time [h]
  
  def __init__(self, name):
    
    print "*", "{0:36} {1:10}".format("Case requested", "-"), name
    self.loadcommon()
    
    if name != "Default":
      exec("self.load"+name+"()")
      self.checkconfig()

  def loadcommon(self):
    
    self._momeV         = 7000.0e9
    self._ppb           = 2.2e11
    self._dpp           = 0.0001074259105 # Update from Stefania, Feb 2018
    self._sigs          = 0.09
    self._epsn          = [2.5e-6, 2.5e-6]
    self._beta          = [[0.15, 0.15], [0.15, 0.15], [3.0, 3.0]] # (updated for Elias and WP2 24 Oct), old: [[0.20, 0.20], [0.20, 0.20], [3.5, 3.5]]
    self._betamindict   = {}
    self._tau_obs       = [None, None] # old: [None, 40.0] (at Elias and WP2 24 Oct, updated to None after them)
    self._tau_cc        = None
    self._kappa         = 1.0
    self._kappac        = 0.1
    self._t1            = 0.0
    self._t2            = 0.0
    
    self._incroncc      = [False, False, False]
    self._incronccrate  = 0.0
    self._niterlev      = 1
    self._IBSRF         = True
    
    self._adaptivexsing = [False, False, False]
    
    self._opticsfile    = "opt_round.madx"
    self._circ          = 26658.8832
    self._rho           = 2804
    self._alfmom        = 3.225e-4
    self._tuneb         = 64.31
    self._tunes         = 0.002342
    self._Nbunch        = 2760 # (updated for Elias and WP2 24 Oct), old: 2748
    self._ipnames       = ["IP1", "IP5", "IP8"]
    self._nip           = len(self._ipnames)
    self._nbunch        = [2748, 2748, 2572] # (updated for Elias and WP2 24 Oct), old: [2736, 2736, 2524] 
    self._xplane        = [0, 0, 0] # old: [0, 0, 1] (for Elias and WP2 24 Oct, updated to 0 after them)
    self._sepLR         = [10.5, 10.5, 21.8] # old: [12.6, 12.6, 12.0] (updated for Elias and WP2 24 Oct, updated to 12.5 after them), older: [12.5, 12.5, 12.0]  -- 21.8 gives around 230 urad = 2*(250 ext - 135 int); 72.9 gives around 770 = 2*(250 + 135)
    self._sepLRsteptime = [False for i in range(self._nip)]
    self._sepLRstep     = [0.0   for i in range(self._nip)]
    self._wcc           = 400.0e6
    self._oncc          = [380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])), 380e-6/(self._sepLR[1]*np.sqrt(self._epsn[1]/(self._momeV/cst.pmass)/self._beta[1][0])), 0.0] # = [0.7657, 0.7657, 0.0], older: = [0.7585, 0.7585, 0.0] (for Mark), oldest: [0.74, 0.74, 0.0] (Elias and WP2 24 Oct)
    self._tcc1          = [[0.0,0.0], [0.0,0.0], [0.0,0.0]]
    self._tcc2          = [[0.0,0.0], [0.0,0.0], [0.0,0.0]]
    
    self._levtech       = ["levlumi" for i in range(self._nip)]
    self._levvar        = ["Beta", "Beta", "ParSep"]
    self._levlumi       = [5.00e34, 5.00e34, 2.00e33] # [7.50e34, 7.50e34, 2.00e33] for ultimate
    self._levppus       = [1.30, 1.30, None]
    self._ck            = [False for i in range(self._nip)]
    self._constbetar    = [True  for i in range(self._nip)]
    self._sepLRconst    = [False for i in range(self._nip)]
    self._longdens      = "qGaussian"
    self._p             = 0.98
    self._step          = 600.0
    self._optimumfill   = True
    self._maxfill       = 20.0*3600.0
    self._penstp        = True
    self._timepenstp    = 0.0
    self._updatepenstp  = False
    
    self._constlong     = True
    self._ppblong       = 0.0
    self._redlong       = 0.0
    self._minlong       = self._sigs
    
    self._redsepLR      = [False for i in range(self._nip)]
    self._minsepLR      = [None for i in range(self._nip)]
    self._ratesepLR     = [None for i in range(self._nip)]
  
    self._xsec          = 81.0
    self._xsecburn      = 111.0
    self._days          = 160.0
    self._eff           = 0.5
    self._turnar        = 145./60. # old: 3.12 (updated for Elias and WP2 24 Oct), 150./60. for ultimate

  ##################################################
  ###################  THESIS  #####################
  ##################################################

  def loadApr19Thesis_LHCdesign2005(self):
    self._momeV         = 7000.0e9
    self._ppb           = 1.15e11
    self._dpp           = 1.2e-4 # assumed
    self._sigs          = 0.0755
    self._epsn          = [3.75e-6, 3.75e-6]
    self._beta          = [[0.55, 0.55], [0.55, 0.55], [1.0, 1.0]]
    self._betamindict   = {}
    self._tau_obs       = [None, None]
    self._tau_cc        = None
    self._kappa         = 1.0
    self._kappac        = 0.1
    self._t1            = 0.0
    self._t2            = 0.0
    
    self._incroncc      = [False, False, False]
    self._incronccrate  = 0.0
    self._niterlev      = 1
    self._IBSRF         = True
    
    self._adaptivexsing = [False, False, False]
    
    self._opticsfile    = "opt_round.madx" # HLLHC
    self._circ          = 26658.8832
    self._rho           = 2803.95
    self._alfmom        = 3.225e-4
    self._tuneb         = 64.31
    self._tunes         = 0.002342 # HLLHC
    self._Nbunch        = 2808
    self._ipnames       = ["IP1", "IP5", "IP8"]
    self._nip           = len(self._ipnames)
    self._nbunch        = [2808, 2808, 2572] # assumed
    self._xplane        = [0, 0, 0]
    self._sepLR         = [2*142.5e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]), 2*142.5e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[1][0]), 2*200e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])]
    self._sepLRsteptime = [False for i in range(self._nip)]
    self._sepLRstep     = [0.0   for i in range(self._nip)]
    self._wcc           = 400.0e6
    self._oncc          = [0.0, 0.0, 0.0]
    self._tcc1          = [[0.0,0.0], [0.0,0.0], [0.0,0.0]]
    self._tcc2          = [[0.0,0.0], [0.0,0.0], [0.0,0.0]]
    
    self._levtech       = ["levlumi" for i in range(self._nip)]
    self._levvar        = ["Beta", "Beta", "ParSep"]
    self._levlumi       = [1.0e40, 1.0e40, 1.0e40]
    self._levppus       = [None, None, None]
    self._ck            = [False for i in range(self._nip)]
    self._constbetar    = [True  for i in range(self._nip)]
    self._sepLRconst    = [False for i in range(self._nip)]
    self._longdens      = "Gaussian"
    self._p             = 1.00
    self._step          = 600.0
    self._optimumfill   = True
    self._maxfill       = 20.0*3600.0
    self._penstp        = True
    self._timepenstp    = 0.0
    self._updatepenstp  = False
    
    self._constlong     = True
    self._ppblong       = 0.0
    self._redlong       = 0.0
    self._minlong       = self._sigs
    
    self._redsepLR      = [False for i in range(self._nip)]
    self._minsepLR      = [None for i in range(self._nip)]
    self._ratesepLR     = [None for i in range(self._nip)]
  
    self._xsec          = 60.0
    self._xsecburn      = 100.0
    self._days          = 160.0
    self._eff           = 0.5
    self._turnar        = 145./60. # old: 3.12 (updated for Elias and WP2 24 Oct), 150./60. for ultimate

    
  def loadApr19Thesis_BaseStanNoLevel_nom(self):
    self._levlumi       = [1.0e40, 1.0e40, 1.0]
    self._p             = 1.0
    self._step          = 180.0
    self._optimumfill   = True
    self._maxfill       = 24.0*3600.0
    pass
    
  def loadEne19Thesis_RedSigs_nom(self):
    self._constlong     = True
    self._ppblong       = 2.1e11
    self._redlong       = -0.01
    self._minlong       = 0.08
    pass
      
  #### Ago18Thesis_BaseStan_nom ###
  
  def loadAgo18Thesis_BaseStan_nom(self):
    pass

  # Ago18Thesis_BaseStan_onccXpXX_nom
  
  def loadAgo18Thesis_BaseStan_oncc0p00_nom(self):
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p05_nom(self):
    self._oncc[0]       = 0.05
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p10_nom(self):
    self._oncc[0]       = 0.10
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p15_nom(self):
    self._oncc[0]       = 0.15
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p20_nom(self):
    self._oncc[0]       = 0.20
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p25_nom(self):
    self._oncc[0]       = 0.25
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p30_nom(self):
    self._oncc[0]       = 0.30
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p35_nom(self):
    self._oncc[0]       = 0.35
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p40_nom(self):
    self._oncc[0]       = 0.40
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p45_nom(self):
    self._oncc[0]       = 0.45
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p50_nom(self):
    self._oncc[0]       = 0.50
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p55_nom(self):
    self._oncc[0]       = 0.55
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p60_nom(self):
    self._oncc[0]       = 0.60
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p65_nom(self):
    self._oncc[0]       = 0.65
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p70_nom(self):
    self._oncc[0]       = 0.70
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p75_nom(self):
    self._oncc[0]       = 0.75
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p80_nom(self):
    self._oncc[0]       = 0.80
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p85_nom(self):
    self._oncc[0]       = 0.85
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p90_nom(self):
    self._oncc[0]       = 0.90
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass

  def loadAgo18Thesis_BaseStan_oncc0p95_nom(self):
    self._oncc[0]       = 0.95
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p00_nom(self):
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p05_nom(self):
    self._oncc[0]       = 1.05
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p10_nom(self):
    self._oncc[0]       = 1.10
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p15_nom(self):
    self._oncc[0]       = 1.15
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p20_nom(self):
    self._oncc[0]       = 1.20
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  # Ago18Thesis_BaseStan_ppbXpXe11_nom
  
  def loadAgo18Thesis_BaseStan_ppb1p5e11_nom(self):
    self._ppb           = 1.5e11
    pass
  
  def loadAgo18Thesis_BaseStan_ppb1p6e11_nom(self):
    self._ppb           = 1.6e11
    pass
  
  def loadAgo18Thesis_BaseStan_ppb1p7e11_nom(self):
    self._ppb           = 1.7e11
    pass
  
  def loadAgo18Thesis_BaseStan_ppb1p8e11_nom(self):
    self._ppb           = 1.8e11
    pass
  
  def loadAgo18Thesis_BaseStan_ppb1p9e11_nom(self):
    self._ppb           = 1.9e11
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p0e11_nom(self):
    self._ppb           = 2.0e11
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p1e11_nom(self):
    self._ppb           = 2.1e11
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p2e11_nom(self):
    self._ppb           = 2.2e11
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p3e11_nom(self):
    self._ppb           = 2.3e11
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p4e11_nom(self):
    self._ppb           = 2.4e11
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p5e11_nom(self):
    self._ppb           = 2.5e11
    pass
  
  # Ago18Thesis_BaseStan_epsnXpXem6_nom
  
  def loadAgo18Thesis_BaseStan_epsn1p5em6_nom(self):
    self._epsn[0]       = 1.5e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn1p6em6_nom(self):
    self._epsn[0]       = 1.6e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn1p7em6_nom(self):
    self._epsn[0]       = 1.7e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn1p8em6_nom(self):
    self._epsn[0]       = 1.8e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn1p9em6_nom(self):
    self._epsn[0]       = 1.9e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p0em6_nom(self):
    self._epsn[0]       = 2.0e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p1em6_nom(self):
    self._epsn[0]       = 2.1e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p2em6_nom(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p3em6_nom(self):
    self._epsn[0]       = 2.3e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p4em6_nom(self):
    self._epsn[0]       = 2.4e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p5em6_nom(self):
    self._epsn[0]       = 2.5e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p6em6_nom(self):
    self._epsn[0]       = 2.6e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p7em6_nom(self):
    self._epsn[0]       = 2.7e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p8em6_nom(self):
    self._epsn[0]       = 2.8e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p9em6_nom(self):
    self._epsn[0]       = 2.9e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_epsn3p0em6_nom(self):
    self._epsn[0]       = 3.0e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  # Ago18Thesis_BaseStan_beta0pXXX_nom
  
  def loadAgo18Thesis_BaseStan_beta0p050_nom(self):
    self._beta[0]       = [0.050, 0.050]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p075_nom(self):
    self._beta[0]       = [0.075, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p100_nom(self):
    self._beta[0]       = [0.100, 0.100]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass

  def loadAgo18Thesis_BaseStan_beta0p125_nom(self):
    self._beta[0]       = [0.125, 0.125]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p150_nom(self):
    self._beta[0]       = [0.150, 0.150]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p175_nom(self):
    self._beta[0]       = [0.175, 0.175]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p200_nom(self):
    self._beta[0]       = [0.200, 0.200]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p225_nom(self):
    self._beta[0]       = [0.225, 0.225]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p250_nom(self):
    self._beta[0]       = [0.250, 0.250]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p275_nom(self):
    self._beta[0]       = [0.275, 0.275]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 1.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p300_nom(self):
    self._beta[0]       = [0.300, 0.300]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 1.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  # Ago18Thesis_BaseStan_sigs0pXXX_nom
  
  def loadAgo18Thesis_BaseStan_sigs0p070_nom(self):
    self._sigs          = 0.070
    self._minlong       = self._sigs
    pass

  def loadAgo18Thesis_BaseStan_sigs0p075_nom(self):
    self._sigs          = 0.075
    self._minlong       = self._sigs
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p080_nom(self):
    self._sigs          = 0.080
    self._minlong       = self._sigs
    pass

  def loadAgo18Thesis_BaseStan_sigs0p085_nom(self):
    self._sigs          = 0.085
    self._minlong       = self._sigs
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p090_nom(self):
    self._sigs          = 0.090
    self._minlong       = self._sigs
    pass

  def loadAgo18Thesis_BaseStan_sigs0p095_nom(self):
    self._sigs          = 0.095
    self._minlong       = self._sigs
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p100_nom(self):
    self._sigs          = 0.100
    self._minlong       = self._sigs
    pass

  def loadAgo18Thesis_BaseStan_sigs0p105_nom(self):
    self._sigs          = 0.105
    self._minlong       = self._sigs
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p110_nom(self):
    self._sigs          = 0.110
    self._minlong       = self._sigs
    pass

  def loadAgo18Thesis_BaseStan_sigs0p115_nom(self):
    self._sigs          = 0.115
    self._minlong       = self._sigs
    pass
   
  def loadAgo18Thesis_BaseStan_sigs0p120_nom(self):
    self._sigs          = 0.120
    self._minlong       = self._sigs
    pass

  def loadAgo18Thesis_BaseStan_sigs0p125_nom(self):
    self._sigs          = 0.125
    self._minlong       = self._sigs
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p130_nom(self):
    self._sigs          = 0.130
    self._minlong       = self._sigs
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p135_nom(self):
    self._sigs          = 0.135
    self._minlong       = self._sigs
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p140_nom(self):
    self._sigs          = 0.140
    self._minlong       = self._sigs
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p145_nom(self):
    self._sigs          = 0.145
    self._minlong       = self._sigs
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p150_nom(self):
    self._sigs          = 0.150
    self._minlong       = self._sigs
    pass
  
  # Ago18Thesis_BaseStan_xsecburn81_nom
  
  def loadAgo18Thesis_BaseStan_xsecburn81_nom(self):
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseStan_onccXpXX_xsecburn81_nom
  
  def loadAgo18Thesis_BaseStan_oncc0p00_xsecburn81_nom(self):
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p00_xsecburn81_nom(self):
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseStan_taucc0pXXXem6_nom

  def loadAgo18Thesis_BaseStan_taucc0p000em6_nom(self):
    self._tau_cc        = 1./( (0.001e-9/self._epsn[0]) * (0.15/380e-6**2)) # 0.000e-6/2.5e-6 = 0.0%
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p025em6_nom(self):
    self._tau_cc        = 1./( (0.025e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.025e-6/2.5e-6 = 1.0%
    pass

  def loadAgo18Thesis_BaseStan_taucc0p050em6_nom(self):
    self._tau_cc        = 1./( (0.050e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.050e-6/2.5e-6 = 2.0%
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p075em6_nom(self):
    self._tau_cc        = 1./( (0.075e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.075e-6/2.5e-6 = 3.0%
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p100em6_nom(self):
    self._tau_cc        = 1./( (0.100e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.100e-6/2.5e-6 = 4.0%
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p115em6_nom(self):
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p125em6_nom(self):
    self._tau_cc        = 1./( (0.125e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.125e-6/2.5e-6 = 5.0%
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p150em6_nom(self):
    self._tau_cc        = 1./( (0.150e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.150e-6/2.5e-6 = 6.0%
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p175em6_nom(self):
    self._tau_cc        = 1./( (0.175e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.175e-6/2.5e-6 = 7.0%
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p200em6_nom(self):
    self._tau_cc        = 1./( (0.200e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.200e-6/2.5e-6 = 8.0%
    pass
  
  # Ago18Thesis_BaseStan_t1sXXXem12_t2sXXXem12_nom
  
  def loadAgo18Thesis_BaseStan_t1m100em12_t2m100em12_nom(self):
    self._t1            = -100e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m090em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m080em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m070em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m060em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m050em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m040em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m030em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m020em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m010em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p000em12_nom(self):
    self._t1            = -100e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p010em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p020em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p030em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p040em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p050em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p060em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p070em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p080em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p090em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p100em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m100em12_nom(self):
    self._t1            =  -90e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m090em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m080em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m070em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m060em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m050em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m040em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m030em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m020em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m010em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p000em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p010em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p020em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p030em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p040em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p050em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p060em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p070em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p080em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p090em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p100em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m100em12_nom(self):
    self._t1            =  -80e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m090em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m080em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m070em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m060em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m050em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m040em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m030em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m020em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m010em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p000em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p010em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p020em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p030em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p040em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p050em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p060em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p070em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p080em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p090em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p100em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m100em12_nom(self):
    self._t1            =  -70e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m090em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m080em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m070em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m060em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m050em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m040em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m030em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m020em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m010em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p000em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p010em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p020em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p030em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p040em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p050em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p060em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p070em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p080em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p090em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p100em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m100em12_nom(self):
    self._t1            =  -60e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m090em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m080em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m070em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m060em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m050em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m040em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m030em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m020em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m010em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p000em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p010em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p020em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p030em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p040em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p050em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p060em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p070em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p080em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p090em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p100em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m100em12_nom(self):
    self._t1            =  -50e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m090em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m080em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m070em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m060em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m050em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m040em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m030em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m020em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m010em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p000em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p010em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p020em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p030em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p040em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p050em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p060em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p070em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p080em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p090em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p100em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m100em12_nom(self):
    self._t1            =  -40e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m090em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m080em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m070em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m060em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m050em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m040em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m030em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m020em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m010em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p000em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p010em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p020em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p030em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p040em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p050em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p060em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p070em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p080em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p090em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p100em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m100em12_nom(self):
    self._t1            =  -30e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m090em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m080em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m070em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m060em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m050em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m040em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m030em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m020em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m010em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p000em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p010em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p020em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p030em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p040em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p050em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p060em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p070em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p080em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p090em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p100em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m100em12_nom(self):
    self._t1            =  -20e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m090em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m080em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m070em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m060em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m050em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m040em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m030em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m020em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m010em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p000em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p010em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p020em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p030em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p040em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p050em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p060em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p070em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p080em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p090em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p100em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m100em12_nom(self):
    self._t1            =  -10e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m090em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m080em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m070em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m060em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m050em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m040em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m030em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m020em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m010em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p000em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p010em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p020em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p030em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p040em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p050em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p060em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p070em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p080em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p090em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p100em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m100em12_nom(self):
    self._t1            =    0e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m090em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m080em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m070em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m060em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m050em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m040em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m030em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m020em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m010em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p000em12_nom(self):
    self._t1            =    0e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p010em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p020em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p030em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p040em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p050em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p060em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p070em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p080em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p090em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p100em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m100em12_nom(self):
    self._t1            =   10e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m090em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m080em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m070em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m060em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m050em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m040em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m030em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m020em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m010em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p000em12_nom(self):
    self._t1            =   10e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p010em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p020em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p030em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p040em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p050em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p060em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p070em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p080em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p090em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p100em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m100em12_nom(self):
    self._t1            =   20e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m090em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m080em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m070em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m060em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m050em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m040em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m030em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m020em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m010em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p000em12_nom(self):
    self._t1            =   20e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p010em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p020em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p030em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p040em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p050em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p060em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p070em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p080em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p090em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p100em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m100em12_nom(self):
    self._t1            =   30e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m090em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m080em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m070em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m060em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m050em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m040em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m030em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m020em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m010em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p000em12_nom(self):
    self._t1            =   30e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p010em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p020em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p030em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p040em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p050em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p060em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p070em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p080em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p090em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p100em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m100em12_nom(self):
    self._t1            =   40e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m090em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m080em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m070em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m060em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m050em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m040em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m030em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m020em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m010em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p000em12_nom(self):
    self._t1            =   40e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p010em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p020em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p030em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p040em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p050em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p060em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p070em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p080em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p090em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p100em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m100em12_nom(self):
    self._t1            =   50e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m090em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m080em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m070em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m060em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m050em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m040em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m030em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m020em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m010em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p000em12_nom(self):
    self._t1            =   50e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p010em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p020em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p030em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p040em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p050em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p060em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p070em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p080em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p090em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p100em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m100em12_nom(self):
    self._t1            =   60e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m090em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m080em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m070em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m060em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m050em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m040em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m030em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m020em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m010em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p000em12_nom(self):
    self._t1            =   60e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p010em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p020em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p030em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p040em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p050em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p060em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p070em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p080em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p090em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p100em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m100em12_nom(self):
    self._t1            =   70e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m090em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m080em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m070em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m060em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m050em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m040em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m030em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m020em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m010em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p000em12_nom(self):
    self._t1            =   70e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p010em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p020em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p030em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p040em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p050em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p060em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p070em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p080em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p090em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p100em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m100em12_nom(self):
    self._t1            =   80e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m090em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m080em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m070em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m060em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m050em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m040em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m030em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m020em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m010em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p000em12_nom(self):
    self._t1            =   80e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p010em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p020em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p030em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p040em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p050em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p060em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p070em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p080em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p090em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p100em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m100em12_nom(self):
    self._t1            =   90e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m090em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m080em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m070em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m060em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m050em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m040em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m030em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m020em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m010em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p000em12_nom(self):
    self._t1            =   90e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p010em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p020em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p030em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p040em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p050em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p060em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p070em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p080em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p090em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p100em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m100em12_nom(self):
    self._t1            =  100e-12
    self._t2            = -100e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m090em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m080em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m070em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m060em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m050em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m040em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m030em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m020em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m010em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p000em12_nom(self):
    self._t1            =  100e-12
    self._t2            =    0e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p010em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   10e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p020em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   20e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p030em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   30e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p040em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   40e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p050em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   50e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p060em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   60e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p070em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   70e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p080em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   80e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p090em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   90e-12
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p100em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  100e-12
    pass
   
  # Ago18Thesis_BaseStan_pXpXX_nom
  
  def loadAgo18Thesis_BaseStan_p0p90_nom(self):
    self._p             = 0.90
    pass
  
  def loadAgo18Thesis_BaseStan_p0p91_nom(self):
    self._p             = 0.91
    pass
  
  def loadAgo18Thesis_BaseStan_p0p92_nom(self):
    self._p             = 0.92
    pass
  
  def loadAgo18Thesis_BaseStan_p0p93_nom(self):
    self._p             = 0.93
    pass
  
  def loadAgo18Thesis_BaseStan_p0p94_nom(self):
    self._p             = 0.94
    pass
  
  def loadAgo18Thesis_BaseStan_p0p95_nom(self):
    self._p             = 0.95
    pass
  
  def loadAgo18Thesis_BaseStan_p0p96_nom(self):
    self._p             = 0.96
    pass
  
  def loadAgo18Thesis_BaseStan_p0p97_nom(self):
    self._p             = 0.97
    pass
  
  def loadAgo18Thesis_BaseStan_p0p98_nom(self):
    self._p             = 0.98
    pass
  
  def loadAgo18Thesis_BaseStan_p0p99_nom(self):
    self._p             = 0.99
    pass

  def loadAgo18Thesis_BaseStan_p0p9999_nom(self):
    self._p             = 0.9999
    pass
  
  def loadAgo18Thesis_BaseStan_p1p00_nom(self):
    self._p             = 1.00
    pass

  # Ago18Thesis_BaseStan_timepenstpXX_nom
  
  def loadAgo18Thesis_BaseStan_timepenstp00_nom(self):
    self._timepenstp    = 0.
    self._updatepenstp  = True
    pass

  def loadAgo18Thesis_BaseStan_timepenstp01_nom(self):
    self._timepenstp    = 1.
    self._updatepenstp  = True
    pass

  def loadAgo18Thesis_BaseStan_timepenstp02_nom(self):
    self._timepenstp    = 2.
    self._updatepenstp  = True
    pass

  def loadAgo18Thesis_BaseStan_timepenstp03_nom(self):
    self._timepenstp    = 3.
    self._updatepenstp  = True
    pass

  def loadAgo18Thesis_BaseStan_timepenstp04_nom(self):
    self._timepenstp    = 4.
    self._updatepenstp  = True
    pass

  def loadAgo18Thesis_BaseStan_timepenstp05_nom(self):
    self._timepenstp    = 5.
    self._updatepenstp  = True
    pass

  def loadAgo18Thesis_BaseStan_timepenstp06_nom(self):
    self._timepenstp    = 6.
    self._updatepenstp  = True
    pass

  def loadAgo18Thesis_BaseStan_timepenstp07_nom(self):
    self._timepenstp    = 7.
    self._updatepenstp  = True
    pass

  def loadAgo18Thesis_BaseStan_timepenstp08_nom(self):
    self._timepenstp    = 8.
    self._updatepenstp  = True
    pass

  def loadAgo18Thesis_BaseStan_timepenstp09_nom(self):
    self._timepenstp    = 9.
    self._updatepenstp  = True
    pass

  def loadAgo18Thesis_BaseStan_timepenstp10_nom(self):
    self._timepenstp    = 10.
    self._updatepenstp  = True
    pass

  # Ago18Thesis_BaseStan_turnarXXX_nom

  def loadAgo18Thesis_BaseStan_turnar120_nom(self):
    self._turnar        = 120./60.
    pass

  def loadAgo18Thesis_BaseStan_turnar125_nom(self):
    self._turnar        = 125./60.
    pass

  def loadAgo18Thesis_BaseStan_turnar130_nom(self):
    self._turnar        = 130./60.
    pass

  def loadAgo18Thesis_BaseStan_turnar135_nom(self):
    self._turnar        = 135./60.
    pass

  def loadAgo18Thesis_BaseStan_turnar140_nom(self):
    self._turnar        = 140./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar145_nom(self):
    self._turnar        = 145./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar150_nom(self):
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar155_nom(self):
    self._turnar        = 155./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar160_nom(self):
    self._turnar        = 160./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar165_nom(self):
    self._turnar        = 165./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar170_nom(self):
    self._turnar        = 170./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar175_nom(self):
    self._turnar        = 175./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar180_nom(self):
    self._turnar        = 180./60.
    pass
  
  # Ago18Thesis_BaseStan_phiXXX_nom
  
  def loadAgo18Thesis_BaseStan_phi300_nom(self):
    self._sepLR[0]      = 300e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_phi325_nom(self):
    self._sepLR[0]      = 325e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass

  def loadAgo18Thesis_BaseStan_phi350_nom(self):
    self._sepLR[0]      = 350e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_phi375_nom(self):
    self._sepLR[0]      = 375e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_phi400_nom(self):
    self._sepLR[0]      = 400e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_phi425_nom(self):
    self._sepLR[0]      = 425e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_phi450_nom(self):
    self._sepLR[0]      = 450e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_phi475_nom(self):
    self._sepLR[0]      = 475e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_phi500_nom(self):
    self._sepLR[0]      = 500e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_phi525_nom(self):
    self._sepLR[0]      = 525e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_phi550_nom(self):
    self._sepLR[0]      = 550e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseStan_phi575_nom(self):
    self._sepLR[0]      = 575e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass

  def loadAgo18Thesis_BaseStan_phi600_nom(self):
    self._sepLR[0]      = 600e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  # Ago18Thesis_BaseStan_Gaussian_nom
  
  def loadAgo18Thesis_BaseStan_Gaussian_nom(self):
    self._longdens      = "Gaussian"
    pass
  
  # Ago18Thesis_BaseStan_Gaussian_dpp_nom
  
  def loadAgo18Thesis_BaseStan_Gaussian_dpp_nom(self):
    self._dpp           = 0.0001102712121
    self._longdens      = "Gaussian"
    pass
  
  # Ago18Thesis_BaseStan_Gaussian_t1sXXXem12_t2sXXXem12_nom
  
  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m100em12_nom(self):
    self._t1            = -100e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m090em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m080em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m070em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m060em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m050em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m040em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m030em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m020em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m010em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p000em12_nom(self):
    self._t1            = -100e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p010em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p020em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p030em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p040em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p050em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p060em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p070em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p080em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p090em12_nom(self):
    self._t1            = -100e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p100em12_nom(self):
    self._t1            = -100e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m100em12_nom(self):
    self._t1            =  -90e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m090em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m080em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m070em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m060em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m050em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m040em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m030em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m020em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m010em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p000em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p010em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p020em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p030em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p040em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p050em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p060em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p070em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p080em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p090em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p100em12_nom(self):
    self._t1            =  -90e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m100em12_nom(self):
    self._t1            =  -80e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m090em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m080em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m070em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m060em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m050em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m040em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m030em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m020em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m010em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p000em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p010em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p020em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p030em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p040em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p050em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p060em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p070em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p080em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p090em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p100em12_nom(self):
    self._t1            =  -80e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m100em12_nom(self):
    self._t1            =  -70e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m090em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m080em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m070em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m060em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m050em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m040em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m030em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m020em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m010em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p000em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p010em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p020em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p030em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p040em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p050em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p060em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p070em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p080em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p090em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p100em12_nom(self):
    self._t1            =  -70e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m100em12_nom(self):
    self._t1            =  -60e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m090em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m080em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m070em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m060em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m050em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m040em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m030em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m020em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m010em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p000em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p010em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p020em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p030em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p040em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p050em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p060em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p070em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p080em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p090em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p100em12_nom(self):
    self._t1            =  -60e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m100em12_nom(self):
    self._t1            =  -50e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m090em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m080em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m070em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m060em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m050em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m040em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m030em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m020em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m010em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p000em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p010em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p020em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p030em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p040em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p050em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p060em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p070em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p080em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p090em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p100em12_nom(self):
    self._t1            =  -50e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m100em12_nom(self):
    self._t1            =  -40e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m090em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m080em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m070em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m060em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m050em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m040em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m030em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m020em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m010em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p000em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p010em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p020em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p030em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p040em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p050em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p060em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p070em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p080em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p090em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p100em12_nom(self):
    self._t1            =  -40e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m100em12_nom(self):
    self._t1            =  -30e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m090em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m080em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m070em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m060em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m050em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m040em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m030em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m020em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m010em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p000em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p010em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p020em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p030em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p040em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p050em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p060em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p070em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p080em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p090em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p100em12_nom(self):
    self._t1            =  -30e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m100em12_nom(self):
    self._t1            =  -20e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m090em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m080em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m070em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m060em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m050em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m040em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m030em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m020em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m010em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p000em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p010em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p020em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p030em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p040em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p050em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p060em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p070em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p080em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p090em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p100em12_nom(self):
    self._t1            =  -20e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m100em12_nom(self):
    self._t1            =  -10e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m090em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m080em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m070em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m060em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m050em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m040em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m030em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m020em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m010em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p000em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p010em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p020em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p030em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p040em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p050em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p060em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p070em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p080em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p090em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p100em12_nom(self):
    self._t1            =  -10e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m100em12_nom(self):
    self._t1            =    0e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m090em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m080em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m070em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m060em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m050em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m040em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m030em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m020em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m010em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p000em12_nom(self):
    self._t1            =    0e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p010em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p020em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p030em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p040em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p050em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p060em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p070em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p080em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p090em12_nom(self):
    self._t1            =    0e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p100em12_nom(self):
    self._t1            =    0e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m100em12_nom(self):
    self._t1            =   10e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m090em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m080em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m070em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m060em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m050em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m040em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m030em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m020em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m010em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p000em12_nom(self):
    self._t1            =   10e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p010em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p020em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p030em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p040em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p050em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p060em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p070em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p080em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p090em12_nom(self):
    self._t1            =   10e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p100em12_nom(self):
    self._t1            =   10e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m100em12_nom(self):
    self._t1            =   20e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m090em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m080em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m070em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m060em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m050em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m040em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m030em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m020em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m010em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p000em12_nom(self):
    self._t1            =   20e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p010em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p020em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p030em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p040em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p050em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p060em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p070em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p080em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p090em12_nom(self):
    self._t1            =   20e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p100em12_nom(self):
    self._t1            =   20e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m100em12_nom(self):
    self._t1            =   30e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m090em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m080em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m070em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m060em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m050em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m040em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m030em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m020em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m010em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p000em12_nom(self):
    self._t1            =   30e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p010em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p020em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p030em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p040em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p050em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p060em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p070em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p080em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p090em12_nom(self):
    self._t1            =   30e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p100em12_nom(self):
    self._t1            =   30e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m100em12_nom(self):
    self._t1            =   40e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m090em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m080em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m070em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m060em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m050em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m040em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m030em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m020em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m010em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p000em12_nom(self):
    self._t1            =   40e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p010em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p020em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p030em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p040em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p050em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p060em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p070em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p080em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p090em12_nom(self):
    self._t1            =   40e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p100em12_nom(self):
    self._t1            =   40e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m100em12_nom(self):
    self._t1            =   50e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m090em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m080em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m070em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m060em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m050em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m040em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m030em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m020em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m010em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p000em12_nom(self):
    self._t1            =   50e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p010em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p020em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p030em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p040em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p050em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p060em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p070em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p080em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p090em12_nom(self):
    self._t1            =   50e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p100em12_nom(self):
    self._t1            =   50e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m100em12_nom(self):
    self._t1            =   60e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m090em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m080em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m070em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m060em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m050em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m040em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m030em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m020em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m010em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p000em12_nom(self):
    self._t1            =   60e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p010em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p020em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p030em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p040em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p050em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p060em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p070em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p080em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p090em12_nom(self):
    self._t1            =   60e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p100em12_nom(self):
    self._t1            =   60e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m100em12_nom(self):
    self._t1            =   70e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m090em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m080em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m070em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m060em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m050em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m040em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m030em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m020em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m010em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p000em12_nom(self):
    self._t1            =   70e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p010em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p020em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p030em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p040em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p050em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p060em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p070em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p080em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p090em12_nom(self):
    self._t1            =   70e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p100em12_nom(self):
    self._t1            =   70e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m100em12_nom(self):
    self._t1            =   80e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m090em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m080em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m070em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m060em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m050em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m040em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m030em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m020em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m010em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p000em12_nom(self):
    self._t1            =   80e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p010em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p020em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p030em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p040em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p050em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p060em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p070em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p080em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p090em12_nom(self):
    self._t1            =   80e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p100em12_nom(self):
    self._t1            =   80e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m100em12_nom(self):
    self._t1            =   90e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m090em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m080em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m070em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m060em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m050em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m040em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m030em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m020em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m010em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p000em12_nom(self):
    self._t1            =   90e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p010em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p020em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p030em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p040em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p050em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p060em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p070em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p080em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p090em12_nom(self):
    self._t1            =   90e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p100em12_nom(self):
    self._t1            =   90e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m100em12_nom(self):
    self._t1            =  100e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m090em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m080em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m070em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m060em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m050em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m040em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m030em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m020em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m010em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p000em12_nom(self):
    self._t1            =  100e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p010em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p020em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p030em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p040em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p050em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p060em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p070em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p080em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p090em12_nom(self):
    self._t1            =  100e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p100em12_nom(self):
    self._t1            =  100e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    pass

  ### Ago18Thesis_BaseStan_ult ###
  
  def loadAgo18Thesis_BaseStan_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStan_onccXpXX_ult
  
  def loadAgo18Thesis_BaseStan_oncc0p00_ult(self):
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p05_ult(self):
    self._oncc[0]       = 0.05
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p10_ult(self):
    self._oncc[0]       = 0.10
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p15_ult(self):
    self._oncc[0]       = 0.15
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p20_ult(self):
    self._oncc[0]       = 0.20
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p25_ult(self):
    self._oncc[0]       = 0.25
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p30_ult(self):
    self._oncc[0]       = 0.30
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p35_ult(self):
    self._oncc[0]       = 0.35
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p40_ult(self):
    self._oncc[0]       = 0.40
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p45_ult(self):
    self._oncc[0]       = 0.45
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p50_ult(self):
    self._oncc[0]       = 0.50
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p55_ult(self):
    self._oncc[0]       = 0.55
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p60_ult(self):
    self._oncc[0]       = 0.60
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p65_ult(self):
    self._oncc[0]       = 0.65
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p70_ult(self):
    self._oncc[0]       = 0.70
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p75_ult(self):
    self._oncc[0]       = 0.75
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p80_ult(self):
    self._oncc[0]       = 0.80
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p85_ult(self):
    self._oncc[0]       = 0.85
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p90_ult(self):
    self._oncc[0]       = 0.90
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc0p95_ult(self):
    self._oncc[0]       = 0.95
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p00_ult(self):
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p05_ult(self):
    self._oncc[0]       = 1.05
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p10_ult(self):
    self._oncc[0]       = 1.10
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p15_ult(self):
    self._oncc[0]       = 1.15
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p20_ult(self):
    self._oncc[0]       = 1.20
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStan_ppbXeXe11_ult
  
  def loadAgo18Thesis_BaseStan_ppb1p5e11_ult(self):
    self._ppb           = 1.5e11
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_ppb1p6e11_ult(self):
    self._ppb           = 1.6e11
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_ppb1p7e11_ult(self):
    self._ppb           = 1.7e11
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_ppb1p8e11_ult(self):
    self._ppb           = 1.8e11
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_ppb1p9e11_ult(self):
    self._ppb           = 1.9e11
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p0e11_ult(self):
    self._ppb           = 2.0e11
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p1e11_ult(self):
    self._ppb           = 2.1e11
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p2e11_ult(self):
    self._ppb           = 2.2e11
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p3e11_ult(self):
    self._ppb           = 2.3e11
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p4e11_ult(self):
    self._ppb           = 2.4e11
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_ppb2p5e11_ult(self):
    self._ppb           = 2.5e11
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStan_epsnXpXem6_ult
  
  def loadAgo18Thesis_BaseStan_epsn1p5em6_ult(self):
    self._epsn[0]       = 1.5e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn1p6em6_ult(self):
    self._epsn[0]       = 1.6e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn1p7em6_ult(self):
    self._epsn[0]       = 1.7e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn1p8em6_ult(self):
    self._epsn[0]       = 1.8e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn1p9em6_ult(self):
    self._epsn[0]       = 1.9e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p0em6_ult(self):
    self._epsn[0]       = 2.0e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p1em6_ult(self):
    self._epsn[0]       = 2.1e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p2em6_ult(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p3em6_ult(self):
    self._epsn[0]       = 2.3e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p4em6_ult(self):
    self._epsn[0]       = 2.4e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p5em6_ult(self):
    self._epsn[0]       = 2.5e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p6em6_ult(self):
    self._epsn[0]       = 2.6e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p7em6_ult(self):
    self._epsn[0]       = 2.7e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p8em6_ult(self):
    self._epsn[0]       = 2.8e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn2p9em6_ult(self):
    self._epsn[0]       = 2.9e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_epsn3p0em6_ult(self):
    self._epsn[0]       = 3.0e-6
    self._epsn[1]       = self._epsn[0]
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStan_beta0pXXX_nult
  
  def loadAgo18Thesis_BaseStan_beta0p050_ult(self):
    self._beta[0]       = [0.050, 0.050]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p075_ult(self):
    self._beta[0]       = [0.075, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p100_ult(self):
    self._beta[0]       = [0.100, 0.100]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p125_ult(self):
    self._beta[0]       = [0.125, 0.125]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p150_ult(self):
    self._beta[0]       = [0.150, 0.150]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p175_ult(self):
    self._beta[0]       = [0.175, 0.175]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p200_ult(self):
    self._beta[0]       = [0.200, 0.200]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p225_ult(self):
    self._beta[0]       = [0.225, 0.225]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p250_ult(self):
    self._beta[0]       = [0.250, 0.250]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p275_ult(self):
    self._beta[0]       = [0.275, 0.275]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 1.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_beta0p300_ult(self):
    self._beta[0]       = [0.300, 0.300]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._oncc[0]       = 1.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStan_sigs0pXXX_ult
  
  def loadAgo18Thesis_BaseStan_sigs0p070_ult(self):
    self._sigs          = 0.070
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_sigs0p075_ult(self):
    self._sigs          = 0.075
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p080_ult(self):
    self._sigs          = 0.080
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_sigs0p085_ult(self):
    self._sigs          = 0.085
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p090_ult(self):
    self._sigs          = 0.090
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_sigs0p095_ult(self):
    self._sigs          = 0.095
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p100_ult(self):
    self._sigs          = 0.100
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_sigs0p105_ult(self):
    self._sigs          = 0.105
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p110_ult(self):
    self._sigs          = 0.110
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_sigs0p115_ult(self):
    self._sigs          = 0.115
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
   
  def loadAgo18Thesis_BaseStan_sigs0p120_ult(self):
    self._sigs          = 0.120
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_sigs0p125_ult(self):
    self._sigs          = 0.125
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p130_ult(self):
    self._sigs          = 0.130
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p135_ult(self):
    self._sigs          = 0.135
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p140_ult(self):
    self._sigs          = 0.140
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p145_ult(self):
    self._sigs          = 0.145
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_sigs0p150_ult(self):
    self._sigs          = 0.150
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_Ago18Thesis_BaseStan_xsecburn81_ult

  def loadAgo18Thesis_BaseStan_xsecburn81_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseStan_onccXpXX_xsecburn81_ult
  
  def loadAgo18Thesis_BaseStan_oncc0p00_xsecburn81_ult(self):
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_BaseStan_oncc1p00_xsecburn81_ult(self):
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_Ago18Thesis_BaseStan_taucc0pXXXem6_ult
  
  def loadAgo18Thesis_BaseStan_taucc0p000em6_ult(self):
    self._tau_cc        = 1./( (0.001e-9/self._epsn[0]) * (0.15/380e-6**2)) # 0.000e-6/2.5e-6 = 0.0%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p025em6_ult(self):
    self._tau_cc        = 1./( (0.025e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.025e-6/2.5e-6 = 1.0%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_taucc0p050em6_ult(self):
    self._tau_cc        = 1./( (0.050e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.050e-6/2.5e-6 = 2.0%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p075em6_ult(self):
    self._tau_cc        = 1./( (0.075e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.075e-6/2.5e-6 = 3.0%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p100em6_ult(self):
    self._tau_cc        = 1./( (0.100e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.100e-6/2.5e-6 = 4.0%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p115em6_ult(self):
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p125em6_ult(self):
    self._tau_cc        = 1./( (0.125e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.125e-6/2.5e-6 = 5.0%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p150em6_ult(self):
    self._tau_cc        = 1./( (0.150e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.150e-6/2.5e-6 = 6.0%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p175em6_ult(self):
    self._tau_cc        = 1./( (0.175e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.175e-6/2.5e-6 = 7.0%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_taucc0p200em6_ult(self):
    self._tau_cc        = 1./( (0.200e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.200e-6/2.5e-6 = 8.0%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStan_t1sXXXem12_t2sXXXem12_ult
  
  def loadAgo18Thesis_BaseStan_t1m100em12_t2m100em12_ult(self):
    self._t1            = -100e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m090em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m080em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m070em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m060em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m050em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m040em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m030em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m020em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2m010em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p000em12_ult(self):
    self._t1            = -100e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p010em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p020em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p030em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p040em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p050em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p060em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p070em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p080em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p090em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m100em12_t2p100em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m100em12_ult(self):
    self._t1            =  -90e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m090em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m080em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m070em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m060em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m050em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m040em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m030em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m020em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2m010em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p000em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p010em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p020em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p030em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p040em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p050em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p060em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p070em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p080em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p090em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m090em12_t2p100em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m100em12_ult(self):
    self._t1            =  -80e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m090em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m080em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m070em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m060em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m050em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m040em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m030em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m020em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2m010em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p000em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p010em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p020em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p030em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p040em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p050em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p060em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p070em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p080em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p090em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m080em12_t2p100em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m100em12_ult(self):
    self._t1            =  -70e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m090em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m080em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m070em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m060em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m050em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m040em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m030em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m020em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2m010em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p000em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p010em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p020em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p030em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p040em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p050em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p060em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p070em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p080em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p090em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m070em12_t2p100em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m100em12_ult(self):
    self._t1            =  -60e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m090em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m080em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m070em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m060em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m050em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m040em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m030em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m020em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2m010em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p000em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p010em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p020em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p030em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p040em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p050em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p060em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p070em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p080em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p090em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m060em12_t2p100em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m100em12_ult(self):
    self._t1            =  -50e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m090em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m080em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m070em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m060em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m050em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m040em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m030em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m020em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2m010em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p000em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p010em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p020em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p030em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p040em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p050em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p060em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p070em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p080em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p090em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m050em12_t2p100em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m100em12_ult(self):
    self._t1            =  -40e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m090em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m080em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m070em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m060em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m050em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m040em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m030em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m020em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2m010em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p000em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p010em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p020em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p030em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p040em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p050em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p060em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p070em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p080em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p090em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m040em12_t2p100em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m100em12_ult(self):
    self._t1            =  -30e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m090em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m080em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m070em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m060em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m050em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m040em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m030em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m020em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2m010em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p000em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p010em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p020em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p030em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p040em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p050em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p060em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p070em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p080em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p090em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m030em12_t2p100em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m100em12_ult(self):
    self._t1            =  -20e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m090em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m080em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m070em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m060em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m050em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m040em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m030em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m020em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2m010em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p000em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p010em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p020em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p030em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p040em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p050em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p060em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p070em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p080em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p090em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m020em12_t2p100em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m100em12_ult(self):
    self._t1            =  -10e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m090em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m080em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m070em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m060em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m050em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m040em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m030em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m020em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2m010em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p000em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p010em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p020em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p030em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p040em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p050em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p060em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p070em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p080em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p090em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1m010em12_t2p100em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m100em12_ult(self):
    self._t1            =    0e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m090em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m080em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m070em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m060em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m050em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m040em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m030em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m020em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2m010em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p000em12_ult(self):
    self._t1            =    0e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p010em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p020em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p030em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p040em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p050em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p060em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p070em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p080em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p090em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p000em12_t2p100em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m100em12_ult(self):
    self._t1            =   10e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m090em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m080em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m070em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m060em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m050em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m040em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m030em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m020em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2m010em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p000em12_ult(self):
    self._t1            =   10e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p010em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p020em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p030em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p040em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p050em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p060em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p070em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p080em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p090em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p010em12_t2p100em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m100em12_ult(self):
    self._t1            =   20e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m090em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m080em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m070em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m060em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m050em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m040em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m030em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m020em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2m010em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p000em12_ult(self):
    self._t1            =   20e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p010em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p020em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p030em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p040em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p050em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p060em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p070em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p080em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p090em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p020em12_t2p100em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m100em12_ult(self):
    self._t1            =   30e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m090em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m080em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m070em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m060em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m050em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m040em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m030em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m020em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2m010em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p000em12_ult(self):
    self._t1            =   30e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p010em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p020em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p030em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p040em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p050em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p060em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p070em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p080em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p090em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p030em12_t2p100em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m100em12_ult(self):
    self._t1            =   40e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m090em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m080em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m070em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m060em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m050em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m040em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m030em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m020em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2m010em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p000em12_ult(self):
    self._t1            =   40e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p010em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p020em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p030em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p040em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p050em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p060em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p070em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p080em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p090em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p040em12_t2p100em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m100em12_ult(self):
    self._t1            =   50e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m090em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m080em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m070em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m060em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m050em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m040em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m030em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m020em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2m010em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p000em12_ult(self):
    self._t1            =   50e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p010em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p020em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p030em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p040em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p050em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p060em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p070em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p080em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p090em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p050em12_t2p100em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m100em12_ult(self):
    self._t1            =   60e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m090em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m080em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m070em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m060em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m050em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m040em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m030em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m020em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2m010em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p000em12_ult(self):
    self._t1            =   60e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p010em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p020em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p030em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p040em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p050em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p060em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p070em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p080em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p090em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p060em12_t2p100em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m100em12_ult(self):
    self._t1            =   70e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m090em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m080em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m070em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m060em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m050em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m040em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m030em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m020em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2m010em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p000em12_ult(self):
    self._t1            =   70e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p010em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p020em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p030em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p040em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p050em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p060em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p070em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p080em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p090em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p070em12_t2p100em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m100em12_ult(self):
    self._t1            =   80e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m090em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m080em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m070em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m060em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m050em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m040em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m030em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m020em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2m010em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p000em12_ult(self):
    self._t1            =   80e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p010em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p020em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p030em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p040em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p050em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p060em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p070em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p080em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p090em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p080em12_t2p100em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m100em12_ult(self):
    self._t1            =   90e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m090em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m080em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m070em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m060em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m050em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m040em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m030em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m020em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2m010em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p000em12_ult(self):
    self._t1            =   90e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p010em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p020em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p030em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p040em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p050em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p060em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p070em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p080em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p090em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p090em12_t2p100em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m100em12_ult(self):
    self._t1            =  100e-12
    self._t2            = -100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m090em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m080em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m070em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m060em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m050em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m040em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m030em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m020em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2m010em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p000em12_ult(self):
    self._t1            =  100e-12
    self._t2            =    0e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p010em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   10e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p020em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   20e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p030em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   30e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p040em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   40e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p050em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   50e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p060em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   60e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p070em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   70e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p080em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   80e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p090em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   90e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_t1p100em12_t2p100em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  100e-12
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  # Ago18Thesis_BaseStan_pXpXX_ult

  def loadAgo18Thesis_BaseStan_p0p90_ult(self):
    self._p             = 0.90
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_p0p91_ult(self):
    self._p             = 0.91
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_p0p92_ult(self):
    self._p             = 0.92
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_p0p93_ult(self):
    self._p             = 0.93
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_p0p94_ult(self):
    self._p             = 0.94
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_p0p95_ult(self):
    self._p             = 0.95
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_p0p96_ult(self):
    self._p             = 0.96
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_p0p97_ult(self):
    self._p             = 0.97
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_p0p98_ult(self):
    self._p             = 0.98
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_p0p99_ult(self):
    self._p             = 0.99
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_p0p9999_ult(self):
    self._p             = 0.9999
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_p1p00_ult(self):
    self._p             = 1.00
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  # Ago18Thesis_BaseStan_timepenstpXX_ult
  
  def loadAgo18Thesis_BaseStan_timepenstp00_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._timepenstp    = 0.
    self._updatepenstp  = True
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_timepenstp01_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._timepenstp    = 1.
    self._updatepenstp  = True
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_timepenstp02_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._timepenstp    = 2.
    self._updatepenstp  = True
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_timepenstp03_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._timepenstp    = 3.
    self._updatepenstp  = True
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_timepenstp04_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._timepenstp    = 4.
    self._updatepenstp  = True
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_timepenstp05_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._timepenstp    = 5.
    self._updatepenstp  = True
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_timepenstp06_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._timepenstp    = 6.
    self._updatepenstp  = True
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_timepenstp07_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._timepenstp    = 7.
    self._updatepenstp  = True
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_timepenstp08_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._timepenstp    = 8.
    self._updatepenstp  = True
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_timepenstp09_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._timepenstp    = 9.
    self._updatepenstp  = True
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_timepenstp10_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._timepenstp    = 10.
    self._updatepenstp  = True
    self._turnar        = 150./60.
    pass

  # Ago18Thesis_BaseStan_turnarXXX_ult

  def loadAgo18Thesis_BaseStan_turnar120_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 120./60.
    pass

  def loadAgo18Thesis_BaseStan_turnar125_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 125./60.
    pass

  def loadAgo18Thesis_BaseStan_turnar130_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 130./60.
    pass

  def loadAgo18Thesis_BaseStan_turnar135_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 135./60.
    pass

  def loadAgo18Thesis_BaseStan_turnar140_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 140./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar145_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 145./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar150_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar155_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar160_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 160./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar165_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 165./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar170_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 170./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar175_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 175./60.
    pass
  
  def loadAgo18Thesis_BaseStan_turnar180_ult(self):
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 180./60.
    pass
  
# Ago18Thesis_BaseStan_phiXXX_ult
  
  def loadAgo18Thesis_BaseStan_phi300_ult(self):
    self._sepLR[0]      = 300e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_phi325_ult(self):
    self._sepLR[0]      = 325e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_phi350_ult(self):
    self._sepLR[0]      = 350e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_phi375_ult(self):
    self._sepLR[0]      = 375e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_phi400_ult(self):
    self._sepLR[0]      = 400e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_phi425_ult(self):
    self._sepLR[0]      = 425e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_phi450_ult(self):
    self._sepLR[0]      = 450e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_phi475_ult(self):
    self._sepLR[0]      = 475e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_phi500_ult(self):
    self._sepLR[0]      = 500e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_phi525_ult(self):
    self._sepLR[0]      = 525e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_phi550_ult(self):
    self._sepLR[0]      = 550e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStan_phi575_ult(self):
    self._sepLR[0]      = 575e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_phi600_ult(self):
    self._sepLR[0]      = 600e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStan_Gaussian_ult
  
  def loadAgo18Thesis_BaseStan_Gaussian_ult(self):
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStan_Gaussian_dpp_ult
  
  def loadAgo18Thesis_BaseStan_Gaussian_dpp_ult(self):
    self._dpp           = 0.0001102712121
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStan_Gaussian_t1sXXXem12_t2sXXXem12_ult
  
  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m100em12_ult(self):
    self._t1            = -100e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m090em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m080em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m070em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m060em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m050em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m040em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m030em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m020em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2m010em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p000em12_ult(self):
    self._t1            = -100e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p010em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p020em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p030em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p040em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p050em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p060em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p070em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p080em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p090em12_ult(self):
    self._t1            = -100e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m100em12_t2p100em12_ult(self):
    self._t1            = -100e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m100em12_ult(self):
    self._t1            =  -90e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m090em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m080em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m070em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m060em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m050em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m040em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m030em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m020em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2m010em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p000em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p010em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p020em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p030em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p040em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p050em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p060em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p070em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p080em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p090em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m090em12_t2p100em12_ult(self):
    self._t1            =  -90e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m100em12_ult(self):
    self._t1            =  -80e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m090em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m080em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m070em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m060em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m050em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m040em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m030em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m020em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2m010em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p000em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p010em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p020em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p030em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p040em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p050em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p060em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p070em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p080em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p090em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m080em12_t2p100em12_ult(self):
    self._t1            =  -80e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m100em12_ult(self):
    self._t1            =  -70e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m090em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m080em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m070em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m060em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m050em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m040em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m030em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m020em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2m010em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p000em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p010em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p020em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p030em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p040em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p050em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p060em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p070em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p080em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p090em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m070em12_t2p100em12_ult(self):
    self._t1            =  -70e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m100em12_ult(self):
    self._t1            =  -60e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m090em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m080em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m070em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m060em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m050em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m040em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m030em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m020em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2m010em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p000em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p010em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p020em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p030em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p040em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p050em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p060em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p070em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p080em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p090em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m060em12_t2p100em12_ult(self):
    self._t1            =  -60e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m100em12_ult(self):
    self._t1            =  -50e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m090em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m080em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m070em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m060em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m050em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m040em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m030em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m020em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2m010em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p000em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p010em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p020em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p030em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p040em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p050em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p060em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p070em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p080em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p090em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m050em12_t2p100em12_ult(self):
    self._t1            =  -50e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m100em12_ult(self):
    self._t1            =  -40e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m090em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m080em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m070em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m060em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m050em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m040em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m030em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m020em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2m010em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p000em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p010em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p020em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p030em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p040em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p050em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p060em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p070em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p080em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p090em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m040em12_t2p100em12_ult(self):
    self._t1            =  -40e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m100em12_ult(self):
    self._t1            =  -30e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m090em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m080em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m070em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m060em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m050em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m040em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m030em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m020em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2m010em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p000em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p010em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p020em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p030em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p040em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p050em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p060em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p070em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p080em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p090em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m030em12_t2p100em12_ult(self):
    self._t1            =  -30e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m100em12_ult(self):
    self._t1            =  -20e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m090em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m080em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m070em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m060em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m050em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m040em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m030em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m020em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2m010em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p000em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p010em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p020em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p030em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p040em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p050em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p060em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p070em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p080em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p090em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m020em12_t2p100em12_ult(self):
    self._t1            =  -20e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m100em12_ult(self):
    self._t1            =  -10e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m090em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m080em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m070em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m060em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m050em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m040em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m030em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m020em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2m010em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p000em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p010em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p020em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p030em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p040em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p050em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p060em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p070em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p080em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p090em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1m010em12_t2p100em12_ult(self):
    self._t1            =  -10e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m100em12_ult(self):
    self._t1            =    0e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m090em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m080em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m070em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m060em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m050em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m040em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m030em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m020em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2m010em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p000em12_ult(self):
    self._t1            =    0e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p010em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p020em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p030em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p040em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p050em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p060em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p070em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p080em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p090em12_ult(self):
    self._t1            =    0e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p000em12_t2p100em12_ult(self):
    self._t1            =    0e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m100em12_ult(self):
    self._t1            =   10e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m090em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m080em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m070em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m060em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m050em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m040em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m030em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m020em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2m010em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p000em12_ult(self):
    self._t1            =   10e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p010em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p020em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p030em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p040em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p050em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p060em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p070em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p080em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p090em12_ult(self):
    self._t1            =   10e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p010em12_t2p100em12_ult(self):
    self._t1            =   10e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m100em12_ult(self):
    self._t1            =   20e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m090em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m080em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m070em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m060em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m050em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m040em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m030em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m020em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2m010em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p000em12_ult(self):
    self._t1            =   20e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p010em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p020em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p030em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p040em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p050em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p060em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p070em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p080em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p090em12_ult(self):
    self._t1            =   20e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p020em12_t2p100em12_ult(self):
    self._t1            =   20e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m100em12_ult(self):
    self._t1            =   30e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m090em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m080em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m070em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m060em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m050em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m040em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m030em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m020em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2m010em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p000em12_ult(self):
    self._t1            =   30e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p010em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p020em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p030em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p040em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p050em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p060em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p070em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p080em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p090em12_ult(self):
    self._t1            =   30e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p030em12_t2p100em12_ult(self):
    self._t1            =   30e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m100em12_ult(self):
    self._t1            =   40e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m090em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m080em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m070em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m060em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m050em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m040em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m030em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m020em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2m010em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p000em12_ult(self):
    self._t1            =   40e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p010em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p020em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p030em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p040em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p050em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p060em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p070em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p080em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p090em12_ult(self):
    self._t1            =   40e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p040em12_t2p100em12_ult(self):
    self._t1            =   40e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m100em12_ult(self):
    self._t1            =   50e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m090em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m080em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m070em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m060em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m050em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m040em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m030em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m020em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2m010em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p000em12_ult(self):
    self._t1            =   50e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p010em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p020em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p030em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p040em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p050em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p060em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p070em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p080em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p090em12_ult(self):
    self._t1            =   50e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p050em12_t2p100em12_ult(self):
    self._t1            =   50e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m100em12_ult(self):
    self._t1            =   60e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m090em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m080em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m070em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m060em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m050em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m040em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m030em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m020em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2m010em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p000em12_ult(self):
    self._t1            =   60e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p010em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p020em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p030em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p040em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p050em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p060em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p070em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p080em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p090em12_ult(self):
    self._t1            =   60e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p060em12_t2p100em12_ult(self):
    self._t1            =   60e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m100em12_ult(self):
    self._t1            =   70e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m090em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m080em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m070em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m060em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m050em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m040em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m030em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m020em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2m010em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p000em12_ult(self):
    self._t1            =   70e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p010em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p020em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p030em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p040em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p050em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p060em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p070em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p080em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p090em12_ult(self):
    self._t1            =   70e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p070em12_t2p100em12_ult(self):
    self._t1            =   70e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m100em12_ult(self):
    self._t1            =   80e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m090em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m080em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m070em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m060em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m050em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m040em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m030em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m020em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2m010em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p000em12_ult(self):
    self._t1            =   80e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p010em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p020em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p030em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p040em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p050em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p060em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p070em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p080em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p090em12_ult(self):
    self._t1            =   80e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p080em12_t2p100em12_ult(self):
    self._t1            =   80e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m100em12_ult(self):
    self._t1            =   90e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m090em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m080em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m070em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m060em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m050em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m040em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m030em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m020em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2m010em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p000em12_ult(self):
    self._t1            =   90e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p010em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p020em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p030em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p040em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p050em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p060em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p070em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p080em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p090em12_ult(self):
    self._t1            =   90e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p090em12_t2p100em12_ult(self):
    self._t1            =   90e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m100em12_ult(self):
    self._t1            =  100e-12
    self._t2            = -100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m090em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m080em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m070em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m060em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m050em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m040em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m030em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m020em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2m010em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  -10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p000em12_ult(self):
    self._t1            =  100e-12
    self._t2            =    0e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p010em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   10e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p020em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   20e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p030em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   30e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p040em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   40e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p050em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   50e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p060em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   60e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p070em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   70e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p080em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   80e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p090em12_ult(self):
    self._t1            =  100e-12
    self._t2            =   90e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  def loadAgo18Thesis_BaseStan_Gaussian_t1p100em12_t2p100em12_ult(self):
    self._t1            =  100e-12
    self._t2            =  100e-12
    self._longdens      = "Gaussian"
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  ##################################################
  
  ### Ago18Thesis_BaseBCMS_nom ###
  
  def loadAgo18Thesis_BaseBCMS_nom(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    pass
  
  # Ago18Thesis_BaseBCMS_onccXpXX_nom
  
  def loadAgo18Thesis_BaseBCMS_oncc0p00_nom(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_BaseBCMS_oncc1p00_nom(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  # Ago18Thesis_BaseBCMS_xsecburn81_nom
  
  def loadAgo18Thesis_BaseBCMS_xsecburn81_nom(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseBCMS_onccXpXX_xsecburn81_nom

  def loadAgo18Thesis_BaseBCMS_oncc0p00_xsecburn81_nom(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_BaseBCMS_oncc1p00_xsecburn81_nom(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseBCMS_taucc0pXXXem6_nom
  
  def loadAgo18Thesis_BaseBCMS_taucc0p115em6_nom(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    pass
  
  # Ago18Thesis_BaseBCMS_turnarXXX_nom

  def loadAgo18Thesis_BaseBCMS_turnar130_nom(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._turnar        = (145.-15.)/60.
    pass
  
  ### Ago18Thesis_BaseBCMS_ult ###
  
  def loadAgo18Thesis_BaseBCMS_ult(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseBCMS_onccXpXX_ult
  
  def loadAgo18Thesis_BaseBCMS_oncc0p00_ult(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseBCMS_oncc1p00_ult(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseBCMS_xsecburn81_ult
  
  def loadAgo18Thesis_BaseBCMS_xsecburn81_ult(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseBCMS_onccXpXX_xsecburn81_ult
  
  def loadAgo18Thesis_BaseBCMS_oncc0p00_xsecburn81_ult(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_BaseBCMS_oncc1p00_xsecburn81_ult(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseBCMS_taucc0pXXXem6_ult
  
  def loadAgo18Thesis_BaseBCMS_taucc0p115em6_ult(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseBCMS_turnarXXX_ult

  def loadAgo18Thesis_BaseBCMS_turnar135_ult(self):
    self._Nbunch        = 2748
    self._nbunch[0]     = 2736
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2374
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = (150.-15.)/60.
    pass
  
  ##################################################
  
  ### Ago18Thesis_BaseStanWire_nom ###
  
  def loadAgo18Thesis_BaseStanWire_nom(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  # Ago18Thesis_BaseStanWire_onccXpXX_nom
  
  def loadAgo18Thesis_BaseStanWire_oncc0p00_nom(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass

  def loadAgo18Thesis_BaseStanWire_oncc1p00_nom(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    pass
  
  # Ago18Thesis_BaseStanWire_xsecburn81_nom
  
  def loadAgo18Thesis_BaseStanWire_xsecburn81_nom(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseStanWire_onccXpXX_xsecburn81_nom
  
  def loadAgo18Thesis_BaseStanWire_oncc0p00_xsecburn81_nom(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_BaseStanWire_oncc1p00_xsecburn81_nom(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseStanWire_taucc0pXXXem6_nom
  
  def loadAgo18Thesis_BaseStanWire_taucc0p115em6_nom(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    pass
  
  # Ago18Thesis_BaseStanWire_turnarXXX_nom
  
  def loadAgo18Thesis_BaseStanWire_turnar130_nom(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = (145.-15.)/60.
    pass
  
  ### Ago18Thesis_BaseStanWire_ult ###
  
  def loadAgo18Thesis_BaseStanWire_ult(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStanWire_onccXpXX_ult
  
  def loadAgo18Thesis_BaseStanWire_oncc0p00_ult(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_BaseStanWire_oncc1p00_ult(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStanWire_xsecburn81_ult
  
  def loadAgo18Thesis_BaseStanWire_xsecburn81_ult(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseStanWire_onccXpXX_xsecburn81_ult
  
  def loadAgo18Thesis_BaseStanWire_oncc0p00_xsecburn81_ult(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_BaseStanWire_oncc1p00_xsecburn81_ult(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseStanWire_taucc0pXXXem6_ult
  
  def loadAgo18Thesis_BaseStanWire_taucc0p115em6_ult(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStanWire_turnarXXX_ult
  
  def loadAgo18Thesis_BaseStanWire_turnar135_ult(self):
    self._beta[0]       = [0.13, 0.13]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 8.5
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = (150.-15.)/60.
    pass
 
  ##################################################
  
  ### Ago18Thesis_BaseStanAdaptiveXsing_nom ###
  
  def loadAgo18Thesis_BaseStanAdaptiveXsing_nom(self):
    self._sepLR[0]      = 308e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 308 urad
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 308e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))  # although the max. CC angle is up to 380 urad
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._adaptivexsing[0] = True
    self._adaptivexsing[1] = self._adaptivexsing[0]
    #self._adaptivexsing[2] remains the same
    pass

  # Ago18Thesis_BaseStanAdaptiveXsing_xsecburn81_nom
  
  def loadAgo18Thesis_BaseStanAdaptiveXsing_xsecburn81_nom(self):
    self._sepLR[0]      = 308e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 308 urad
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 308e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))  # although the max. CC angle is up to 380 urad
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._adaptivexsing[0] = True
    self._adaptivexsing[1] = self._adaptivexsing[0]
    #self._adaptivexsing[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseStanAdaptiveXsing_taucc0pXXXem6_nom
  
  def loadAgo18Thesis_BaseStanAdaptiveXsing_taucc0p115em6_nom(self):
    self._sepLR[0]      = 308e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 308 urad
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 308e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))  # although the max. CC angle is up to 380 urad
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._adaptivexsing[0] = True
    self._adaptivexsing[1] = self._adaptivexsing[0]
    #self._adaptivexsing[2] remains the same
    pass
  
  # Ago18Thesis_BaseStanAdaptiveXsing_turnarXXX_nom

  def loadAgo18Thesis_BaseStanAdaptiveXsing_turnar130_nom(self):
    self._sepLR[0]      = 308e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 308 urad
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 308e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))  # although the max. CC angle is up to 380 urad
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._adaptivexsing[0] = True
    self._adaptivexsing[1] = self._adaptivexsing[0]
    #self._adaptivexsing[2] remains the same
    self._turnar        = (145.-15.)/60.
    pass
  
  ### Ago18Thesis_BaseStanAdaptiveXsing_ult ###
  
  def loadAgo18Thesis_BaseStanAdaptiveXsing_ult(self):
    self._sepLR[0]      = 308e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 308 urad
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 308e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))  # although the max. CC angle is up to 380 urad
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._adaptivexsing[0] = True
    self._adaptivexsing[1] = self._adaptivexsing[0]
    #self._adaptivexsing[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStanAdaptiveXsing_xsecburn81_ult
  
  def loadAgo18Thesis_BaseStanAdaptiveXsing_xsecburn81_ult(self):
    self._sepLR[0]      = 308e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 308 urad
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 308e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))  # although the max. CC angle is up to 380 urad
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._adaptivexsing[0] = True
    self._adaptivexsing[1] = self._adaptivexsing[0]
    #self._adaptivexsing[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_BaseStanAdaptiveXsing_taucc0pXXXem6_ult
  
  def loadAgo18Thesis_BaseStanAdaptiveXsing_taucc0p115em6_ult(self):
    self._sepLR[0]      = 308e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 308 urad
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 308e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))  # although the max. CC angle is up to 380 urad
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._adaptivexsing[0] = True
    self._adaptivexsing[1] = self._adaptivexsing[0]
    #self._adaptivexsing[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_BaseStanAdaptiveXsing_turnarXXX_ult

  def loadAgo18Thesis_BaseStanAdaptiveXsing_turnar135_ult(self):
    self._sepLR[0]      = 308e-6/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 308 urad
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 308e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))  # although the max. CC angle is up to 380 urad
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._adaptivexsing[0] = True
    self._adaptivexsing[1] = self._adaptivexsing[0]
    #self._adaptivexsing[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = (150.-15.)/60.
    pass
  
  ##################################################
  
  ### Ago18Thesis_UltEnergyBaseStan500urad
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan500urad_onccXpXX
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_oncc0p00_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_oncc1p00_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan500urad_xsecburn81_nom
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_xsecburn81_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_UltEnergyBaseStan500urad_onccXpXX_xsecburn81_nom
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_oncc0p00_xsecburn81_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_oncc1p00_xsecburn81_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_UltEnergyBaseStan500urad_taucc0pXXXem6_nom
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_taucc0p115em6_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan500urad_turnarXXX_nom
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_turnar135_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = (150.-15.)/60.
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_turnar145_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = 145./60.
    pass
  
  ### Ago18Thesis_UltEnergyBaseStan500urad_ult ###
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan500urad_onccXpXX_ult
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_oncc0p00_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_oncc1p00_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan500urad_xsecburn81_ult
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_xsecburn81_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_UltEnergyBaseStan500urad_onccXpXX_xsecburn81_ult

  def loadAgo18Thesis_UltEnergyBaseStan500urad_oncc0p00_xsecburn81_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_oncc1p00_xsecburn81_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_UltEnergyBaseStan500urad_taucc0pXXXem6_ult
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_taucc0p115em6_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan500urad_turnarXXX_ult

  def loadAgo18Thesis_UltEnergyBaseStan500urad_turnar140_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = (155.-15.)/60.
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan500urad_turnar150_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._sepLR[0]      = 0.000496282662778985/np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]) # to give exactly 496 urad (approx 500 urad), the angle used at 7 TeV
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  ##################################################
  
  ### Ago18Thesis_UltEnergyBaseStan
  
  def loadAgo18Thesis_UltEnergyBaseStan_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan_onccXpXX
  
  def loadAgo18Thesis_UltEnergyBaseStan_oncc0p00_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan_oncc1p00_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan_xsecburn81_nom
  
  def loadAgo18Thesis_UltEnergyBaseStan_xsecburn81_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_UltEnergyBaseStan_onccXpXX_xsecburn81_nom
  
  def loadAgo18Thesis_UltEnergyBaseStan_oncc0p00_xsecburn81_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan_oncc1p00_xsecburn81_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_UltEnergyBaseStan_taucc0pXXXem6_nom
  
  def loadAgo18Thesis_UltEnergyBaseStan_taucc0p115em6_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan_turnarXXX_nom
  
  def loadAgo18Thesis_UltEnergyBaseStan_turnar135_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = (150.-15.)/60.
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan_turnar145_nom(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454 # From Stefania, Feb 2018
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = 145./60.
    pass
  
  ### Ago18Thesis_UltEnergyBaseStan_ult ###
  
  def loadAgo18Thesis_UltEnergyBaseStan_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan_onccXpXX_ult
  
  def loadAgo18Thesis_UltEnergyBaseStan_oncc0p00_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan_oncc1p00_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    # self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan_xsecburn81_ult
  
  def loadAgo18Thesis_UltEnergyBaseStan_xsecburn81_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_UltEnergyBaseStan_onccXpXX_xsecburn81_ult

  def loadAgo18Thesis_UltEnergyBaseStan_oncc0p00_xsecburn81_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan_oncc1p00_xsecburn81_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_UltEnergyBaseStan_taucc0pXXXem6_ult
  
  def loadAgo18Thesis_UltEnergyBaseStan_taucc0p115em6_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 155./60.
    pass
  
  # Ago18Thesis_UltEnergyBaseStan_turnarXXX_ult

  def loadAgo18Thesis_UltEnergyBaseStan_turnar140_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = (155.-15.)/60.
    pass
  
  def loadAgo18Thesis_UltEnergyBaseStan_turnar150_ult(self):
    self._momeV         = 7500.0e9
    self._dpp           = 0.0001037829454
    self._oncc[0]       = 354e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0])) # = 0.7133
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  ##################################################

  ### Ago18Thesis_BaseStanLHCb_betaXpX_phiXXX_levlumiXpXe34 ###
 
  def loadAgo18Thesis_BaseStanLHCb_beta1p4_phi230_levlumi2p0e33_nom(self):
    self._beta[2]       = [1.4, 1.4]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e33
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta1p4_phi230_levlumi1p0e34_nom(self):
    self._beta[2]       = [1.4, 1.4]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.0e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta1p4_phi230_levlumi1p5e34_nom(self):
    self._beta[2]       = [1.4, 1.4]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.5e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta1p4_phi230_levlumi2p0e34_nom(self):
    self._beta[2]       = [1.4, 1.4]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e34
    pass

  def loadAgo18Thesis_BaseStanLHCb_beta2p0_phi230_levlumi2p0e33_nom(self):
    self._beta[2]       = [2.0, 2.0]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e33
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta2p0_phi230_levlumi1p0e34_nom(self):
    self._beta[2]       = [2.0, 2.0]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.0e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta2p0_phi230_levlumi1p5e34_nom(self):
    self._beta[2]       = [2.0, 2.0]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.5e34
    pass
 
  def loadAgo18Thesis_BaseStanLHCb_beta2p0_phi230_levlumi2p0e34_nom(self):
    self._beta[2]       = [2.0, 2.0]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e34
    pass
    
  def loadAgo18Thesis_BaseStanLHCb_beta3p0_phi230_levlumi2p0e33_nom(self):
    self._beta[2]       = [3.0, 3.0]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e33
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta3p0_phi230_levlumi1p0e34_nom(self):
    self._beta[2]       = [3.0, 3.0]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.0e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta3p0_phi230_levlumi1p5e34_nom(self):
    self._beta[2]       = [3.0, 3.0]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.5e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta3p0_phi230_levlumi2p0e34_nom(self):
    self._beta[2]       = [3.0, 3.0]
    self._sepLR[2]      = 2*(250e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta1p4_phi770_levlumi2p0e33_nom(self):
    self._beta[2]       = [1.4, 1.4]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e33
    pass
  
  
  def loadAgo18Thesis_BaseStanLHCb_beta1p4_phi770_levlumi1p0e34_nom(self):
    self._beta[2]       = [1.4, 1.4]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.0e34
    pass
    
  def loadAgo18Thesis_BaseStanLHCb_beta1p4_phi770_levlumi1p5e34_nom(self):
    self._beta[2]       = [1.4, 1.4]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.5e34
    pass

  def loadAgo18Thesis_BaseStanLHCb_beta1p4_phi770_levlumi2p0e34_nom(self):
    self._beta[2]       = [1.4, 1.4]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e34
    pass

  def loadAgo18Thesis_BaseStanLHCb_beta2p0_phi770_levlumi2p0e33_nom(self):
    self._beta[2]       = [2.0, 2.0]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e33
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta2p0_phi770_levlumi1p0e34_nom(self):
    self._beta[2]       = [2.0, 2.0]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.0e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta2p0_phi770_levlumi1p5e34_nom(self):
    self._beta[2]       = [2.0, 2.0]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.5e34
    pass
 
  def loadAgo18Thesis_BaseStanLHCb_beta2p0_phi770_levlumi2p0e34_nom(self):
    self._beta[2]       = [2.0, 2.0]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e34
    pass

  def loadAgo18Thesis_BaseStanLHCb_beta3p0_phi770_levlumi2p0e33_nom(self):
    self._beta[2]       = [3.0, 3.0]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e33
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta3p0_phi770_levlumi1p0e34_nom(self):
    self._beta[2]       = [3.0, 3.0]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.0e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta3p0_phi770_levlumi1p5e34_nom(self):
    self._beta[2]       = [3.0, 3.0]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.5e34
    pass
 
  def loadAgo18Thesis_BaseStanLHCb_beta3p0_phi770_levlumi2p0e34_nom(self):
    self._beta[2]       = [3.0, 3.0]
    self._sepLR[2]      = 2*(250e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e34
    pass

  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi130_levlumi2p0e33_nom(self):
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*(200e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e33
    pass

  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi130_levlumi1p0e34_nom(self):
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*(200e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.0e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi130_levlumi1p5e34_nom(self):
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*(200e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.5e34
    pass

  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi130_levlumi2p0e34_nom(self):
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*(200e-6 - 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e34
    pass

  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi570_levlumi2p0e33_nom(self):
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*(150e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e33
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi570_levlumi1p0e34_nom(self):
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*(150e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.0e34
    pass

  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi570_levlumi1p5e34_nom(self):
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*(150e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 1.5e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi570_levlumi2p0e34_nom(self):
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*(150e-6 + 135e-6) /np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[2][0])
    self._levlumi[2]    = 2.0e34
    pass

  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi419V_levlumi2p0e33_nom(self):
    self._xplane[2]     = 1
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*np.sqrt((160e-6)**2 + (135e-6)**2) /np.sqrt(self._epsn[1]/(self._momeV/cst.pmass)/self._beta[2][1])
    self._levlumi[2]    = 2.0e33
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi419V_levlumi1p0e34_nom(self):
    self._xplane[2]     = 1
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*np.sqrt((160e-6)**2 + (135e-6)**2) /np.sqrt(self._epsn[1]/(self._momeV/cst.pmass)/self._beta[2][1])
    self._levlumi[2]    = 1.0e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi419V_levlumi1p5e34_nom(self):
    self._xplane[2]     = 1
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*np.sqrt((160e-6)**2 + (135e-6)**2) /np.sqrt(self._epsn[1]/(self._momeV/cst.pmass)/self._beta[2][1])
    self._levlumi[2]    = 1.5e34
    pass
  
  def loadAgo18Thesis_BaseStanLHCb_beta1p5_phi419V_levlumi2p0e34_nom(self):
    self._xplane[2]     = 1
    self._beta[2]       = [1.5, 1.5]
    self._sepLR[2]      = 2*np.sqrt((160e-6)**2 + (135e-6)**2) /np.sqrt(self._epsn[1]/(self._momeV/cst.pmass)/self._beta[2][1])
    self._levlumi[2]    = 2.0e34
    pass
  
  ##################################################

  ### Ago18Thesis_Flat_nom ###

  def loadAgo18Thesis_Flat_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass

  # Ago18Thesis_Flat_onccXpXX_nom
  
  def loadAgo18Thesis_Flat_oncc0p00_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_Flat_oncc1p00_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  # Ago18Thesis_Flat_xsecburn81_nom
  
  def loadAgo18Thesis_Flat_xsecburn81_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_Flat_onccXpXX_xsecburn81_nom

  def loadAgo18Thesis_Flat_oncc0p00_xsecburn81_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_Flat_oncc1p00_xsecburn81_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_Flat_taucc0pXXXem6_nom
  
  def loadAgo18Thesis_Flat_taucc0p115em6_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    pass
  
  # Ago18Thesis_Flat_turnarXXX_nom
  
  def loadAgo18Thesis_Flat_turnar130_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = (145.-15.)/60.
    pass
  
  ### Ago18Thesis_Flat_ult ###
  
  def loadAgo18Thesis_Flat_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  # Ago18Thesis_Flat_onccXpXX_ult
  
  def loadAgo18Thesis_Flat_oncc0p00_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
    
  def loadAgo18Thesis_Flat_oncc1p00_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  # Ago18Thesis_Flat_xsecburn81_ult
  
  def loadAgo18Thesis_Flat_xsecburn81_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass

  # Ago18Thesis_Flat_onccXpXX_xsecburn81_ult
  
  def loadAgo18Thesis_Flat_oncc0p00_xsecburn81_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass

  def loadAgo18Thesis_Flat_oncc1p00_xsecburn81_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass

  # Ago18Thesis_Flat_taucc0pXXXem6_ult
  
  def loadAgo18Thesis_Flat_taucc0p115em6_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_Flat_turnarXXX_ult
  
  def loadAgo18Thesis_Flat_turnar135_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = (150.-15.)/60.
    pass

  ##################################################
  
  ### Ago18Thesis_8b4e_nom ###
  
  def loadAgo18Thesis_8b4e_nom(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 140. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    pass
  
  # Ago18Thesis_8b4e_onccXpXX_nom
  
  def loadAgo18Thesis_8b4e_oncc0p00_nom(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 140. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    pass
  
  def loadAgo18Thesis_8b4e_oncc1p00_nom(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 140. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    pass
  
  # Ago18Thesis_8b4e_xsecburn81_nom
  
  def loadAgo18Thesis_8b4e_xsecburn81_nom(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 140. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_8b4e_onccXpXX_xsecburn81_nom
  
  def loadAgo18Thesis_8b4e_oncc0p00_xsecburn81_nom(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 140. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_8b4e_oncc1p00_xsecburn81_nom(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 140. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_8b4e_taucc0pXXXem6_nom
  
  def loadAgo18Thesis_8b4e_taucc0p115em6_nom(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.2e-6 = 5.2272727%
    self._levlumi[0]    = 140. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    pass
  
  # Ago18Thesis_8b4e_turnarXXX_nom
  
  def loadAgo18Thesis_8b4e_turnar130_nom(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 140. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._turnar        = (145.-15.)/60.
    pass
  
  ### Ago18Thesis_8b4e_ult ###
  
  def loadAgo18Thesis_8b4e_ult(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 200. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_8b4e_onccXpXX_ult
  
  def loadAgo18Thesis_8b4e_oncc0p00_ult(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 200. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_8b4e_oncc1p00_ult(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 200. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_8b4e_xsecburn81_ult
  
  def loadAgo18Thesis_8b4e_xsecburn81_ult(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 200. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass

  # Ago18Thesis_8b4e_onccXpXX_xsecburn81_ult
  
  def loadAgo18Thesis_8b4e_oncc0p00_xsecburn81_ult(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 200. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_8b4e_oncc1p00_xsecburn81_ult(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 200. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass

  # Ago18Thesis_8b4e_taucc0pXXXem6_ult
  
  def loadAgo18Thesis_8b4e_taucc0p115em6_ult(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.2e-6 = 5.2272727%
    self._levlumi[0]    = 200. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_8b4e_turnarXXX_ult
  
  def loadAgo18Thesis_8b4e_turnar135_ult(self):
    self._epsn[0]       = 2.2e-6
    self._epsn[1]       = self._epsn[0]
    self._Nbunch        = 1972
    self._nbunch[0]     = 1967
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 1886
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 200. /(self._xsec*1e-27)*self._nbunch[0]*(cst.clight/self._circ)
    self._levlumi[1]    = self._levlumi[0]
    self._levlumi[2]    = 5.601/(self._xsec*1e-27)*self._nbunch[2]*(cst.clight/self._circ)
    self._turnar        = (150.-15.)/60.
    pass
  
  ##################################################
  
  ### Ago18Thesis_200MHzFlatGauss_nom ###

  def loadAgo18Thesis_200MHzFlatGauss_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  # Ago18Thesis_200MHzFlatGauss_onccXpXX_nom
  
  def loadAgo18Thesis_200MHzFlatGauss_oncc0p00_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  def loadAgo18Thesis_200MHzFlatGauss_oncc1p00_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  # Ago18Thesis_200MHzFlatGauss_xsecburn81_nom
  
  def loadAgo18Thesis_200MHzFlatGauss_xsecburn81_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_200MHzFlatGauss_onccXpXX_xsecburn81_nom
  
  def loadAgo18Thesis_200MHzFlatGauss_oncc0p00_xsecburn81_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_200MHzFlatGauss_oncc1p00_xsecburn81_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_200MHzFlatGauss_taucc0pXXXem6_nom
  
  def loadAgo18Thesis_200MHzFlatGauss_taucc0p115em6_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    pass
  
  # Ago18Thesis_200MHzFlatGauss_turnar130_nom
  
  def loadAgo18Thesis_200MHzFlatGauss_turnar130_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = (145.-15.)/60.
    pass
  
  ### Ago18Thesis_200MHzFlatGauss_ult ###
  
  def loadAgo18Thesis_200MHzFlatGauss_ult(self):
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  # Ago18Thesis_200MHzFlatGauss_onccXpXX_ult
  
  def loadAgo18Thesis_200MHzFlatGauss_oncc0p00_ult(self):
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_200MHzFlatGauss_oncc1p00_ult(self):
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass

  # Ago18Thesis_200MHzFlatGauss_xsecburn81_ult
  
  def loadAgo18Thesis_200MHzFlatGauss_xsecburn81_ult(self):
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass

  # Ago18Thesis_200MHzFlatGauss_onccXpXX_xsecburn81_ult
  
  def loadAgo18Thesis_200MHzFlatGauss_oncc0p00_xsecburn81_ult(self):
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass
  
  def loadAgo18Thesis_200MHzFlatGauss_oncc1p00_xsecburn81_ult(self):
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 1.00
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass

  # Ago18Thesis_200MHzFlatGauss_taucc0pXXXem6_ult
  
  def loadAgo18Thesis_200MHzFlatGauss_taucc0p115em6_ult(self):
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_200MHzFlatGauss_turnar135_ult
  
  def loadAgo18Thesis_200MHzFlatGauss_turnar135_ult(self):
    self._dpp           = 1.0e-4
    self._longdens      = "Gaussian"
    self._sigs          = 0.15
    self._constlong     = False
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 380e-6/(self._sepLR[0]*np.sqrt(self._epsn[0]/(self._momeV/cst.pmass)/self._beta[0][0]))
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = (150.-15.)/60.
    pass

  ##################################################

  ### Ago18Thesis_NoCCFlat_nom ###
  
  def loadAgo18Thesis_NoCCFlat_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    pass
  
  # Ago18Thesis_NoCCFlat_xsecburn81_nom
  
  def loadAgo18Thesis_NoCCFlat_xsecburn81_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_NoCCFlat_turnarXXX_nom
  
  def loadAgo18Thesis_NoCCFlat_turnar130_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._turnar        = (145.-15.)/60.
    pass
  
  # Ago18Thesis_NoCCFlat_levppusXpXX_nom
  
  def loadAgo18Thesis_NoCCFlat_levppus1p00_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 1.00                   # NOT VALID
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus1p10_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 1.10                   # NOT VALID
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus1p20_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 1.20                   # NOT VALID
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus1p30_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 1.30
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus1p40_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 1.40
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus1p50_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 1.50
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus1p60_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 1.60
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus1p70_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 1.70
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus1p80_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 1.80
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus1p90_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 1.90                   # NOT REACHED
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    pass
  
  ### Ago18Thesis_NoCCFlat_ult ###
  
  def loadAgo18Thesis_NoCCFlat_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_NoCCFlat_xsecburn81_ult
  
  def loadAgo18Thesis_NoCCFlat_xsecburn81_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = 150./60.
    self._xsecburn      = 81.0
    pass

  # Ago18Thesis_NoCCFlat_turnarXXX_ult
  
  def loadAgo18Thesis_NoCCFlat_turnar135_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._turnar        = (150.-15.)/60.
    pass
  
  ### Ago18Thesis_NoCCFlat_levppusXpXX_ult
  
  def loadAgo18Thesis_NoCCFlat_levppus1p90_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 1.90                   # NOT VALID?
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus2p00_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 2.00                   # NOT VALID?
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus2p10_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 2.10                   # NOT VALID?
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus2p20_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 2.20                   # NOT VALID?
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus2p30_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 2.30                   # NOT VALID?
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus2p40_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 2.40                   # VALID?
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus2p50_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 2.50                   # VALID?
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus2p60_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 2.60                   # VALID?
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus2p70_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 2.70
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._turnar        = 150./60.
    pass
  
  def loadAgo18Thesis_NoCCFlat_levppus2p80_ult(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.315, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 12.6333
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.0
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levppus[0]    = 2.80                   # NOT REACHED?
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._turnar        = 150./60.
    pass
  
  ##################################################

  ### Ago18Thesis_CKaSF_nom ###
  
  def loadAgo18Thesis_CKaSF_nom(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.15, 0.15]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 12.5
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 1.0 # = 12./12.
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.075
    self._minlong       = self._sigs
    self._longdens      = "Gaussian"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 3.0
    pass

  # Ago18Thesis_CKaSF_taucc0pXXXem6_nom
  
  def loadAgo18Thesis_CKaSF_taucc0p115em6_nom(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.15, 0.15]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 12.5
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 1.0 # = 12./12.
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.075
    self._minlong       = self._sigs
    self._longdens      = "Gaussian"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 3.0
    pass

  # Ago18Thesis_CKaSF_turnarXXX_nom
  
  def loadAgo18Thesis_CKaSF_turnar145_nom(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.15, 0.15]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 12.5
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 1.0 # = 12./12.
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.075
    self._minlong       = self._sigs
    self._longdens      = "Gaussian"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 145./60.
    pass
  
  # Ago18Thesis_CKaSF_taucc0pXXXem6_turnarXXX_nom
  
  def loadAgo18Thesis_CKaSF_taucc0p115em6_turnar145_nom(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.15, 0.15]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 12.5
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 1.0 # = 12./12.
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.075
    self._minlong       = self._sigs
    self._longdens      = "Gaussian"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 145./60.
    pass
  
  ### Ago18Thesis_CKaSF_ult ###
  
  def loadAgo18Thesis_CKaSF_ult(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.15, 0.15]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 12.5
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 1.0 # = 12./12.
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.075
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._longdens      = "Gaussian"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 3.0
    pass
  
  # Ago18Thesis_CKaSF_taucc0pXXXem6_ult
  
  def loadAgo18Thesis_CKaSF_taucc0p115em6_ult(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.15, 0.15]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 12.5
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 1.0 # = 12./12.
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.075
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._longdens      = "Gaussian"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 3.0
    pass
  
  # Ago18Thesis_CKaSF_turnarXXX_ult
  
  def loadAgo18Thesis_CKaSF_turnar150_ult(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.15, 0.15]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 12.5
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 1.0 # = 12./12.
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.075
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._longdens      = "Gaussian"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_CKaSF_taucc0pXXXem6_turnarXXX_ult
  
  def loadAgo18Thesis_CKaSF_taucc0p115em6_turnar150_ult(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.15, 0.15]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 12.5
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 1.0 # = 12./12.
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.075
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._longdens      = "Gaussian"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 150./60.
    pass
  
  ##################################################
  
  ### Ago18Thesis_CKcSF_nom ###
  
  def loadAgo18Thesis_CKcSF_nom(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.30, 0.10]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 15.0
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 0.881 # approx = 9./10.2
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.61
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 3.0
    pass

  # Ago18Thesis_CKcSF_taucc0pXXXem6_nom
  
  def loadAgo18Thesis_CKcSF_taucc0p115em6_nom(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.30, 0.10]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 15.0
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 0.881 # approx = 9./10.2
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.61
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 3.0
    pass

  # Ago18Thesis_CKcSF_turnarXXX_nom
  
  def loadAgo18Thesis_CKcSF_turnar145_nom(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.30, 0.10]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 15.0
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 0.881 # approx = 9./10.2
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.61
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 145./60.
    pass
  
  # Ago18Thesis_CKcSF_taucc0pXXXem6_turnarXXX_nom
  
  def loadAgo18Thesis_CKcSF_taucc0p115em6_turnar145_nom(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.30, 0.10]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 15.0
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 0.881 # approx = 9./10.2
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.61
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 145./60.
    pass
  
  ### Ago18Thesis_CKcSF_ult ###
  
  def loadAgo18Thesis_CKcSF_ult(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.30, 0.10]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 15.0
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 0.916 # approx = 9.35/10.2
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 1.00
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 3.0
    pass
  
  # Ago18Thesis_CKcSF_taucc0pXXXem6_ult
  
  def loadAgo18Thesis_CKcSF_taucc0p115em6_ult(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.30, 0.10]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 15.0
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 0.916 # approx = 9.35/10.2
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 1.00
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 3.0
    pass
  
  # Ago18Thesis_CKcSF_turnarXXX_ult
  
  def loadAgo18Thesis_CKcSF_turnar150_ult(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.30, 0.10]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 15.0
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 0.916 # approx = 9.35/10.2
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 1.00
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 150./60.
    pass
  
  # Ago18Thesis_CKcSF_taucc0pXXXem6_turnarXXX_ult
  
  def loadAgo18Thesis_CKcSF_taucc0p115em6_turnar150_ult(self):
    self._Nbunch        = 2808
    self._nbunch[0]     = 2808
    self._nbunch[1]     = self._nbunch[0]
    self._nbunch[2]     = 2524
    self._beta[0]       = [0.30, 0.10]
    self._beta[1]       = self._beta[0]
    self._beta[2]       = [3.5, 3.5]
    self._sepLR[0]      = 15.0
    self._sepLR[1]      = self._sepLR[1]
    self._sepLR[2]      = 12.0
    self._oncc[0]       = 0.916 # approx = 9.35/10.2
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._IBSRF         = False
    self._dpp           = 1.2e-4
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levlumi[0]    = 7.5e34
    self._levlumi[1]    = self._levlumi[0]
    #self._levlumi[2] remains the same
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 1.00
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    self._xsec          = 85
    self._xsecburn      = 100
    self._turnar        = 150./60.
    pass
  
  ##################################################
  
  ### FlatCK ###
  
  def loadAgo18Thesis_FlatCK_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.9036
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.61
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    pass
  
  # Ago18Thesis_FlatCK_xsecburn81_nom
  
  def loadAgo18Thesis_FlatCK_xsecburn81_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.9036
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.61
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    self._xsecburn      = 81.0
    pass
  
  # Ago18Thesis_FlatCK_taucc0pXXXem6_nom
  
  def loadAgo18Thesis_FlatCK_taucc0p115em6_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.9036
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._tau_cc        = 1./( (0.115e-6/self._epsn[0]) * (0.15/380e-6**2)) # 0.115e-6/2.5e-6 = 4.6%
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.61
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    pass

  # Ago18Thesis_FlatCK_turnarXXX_nom
  
  def loadAgo18Thesis_FlatCK_turnar130_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.9036
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.61
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    self._turnar        = (145.-15.)/60.
    pass

  # Ago18Thesis_FlatCK_levppusXpXX_nom  # oncc optimized accordingly
  
  def loadAgo18Thesis_FlatCK_levppus0p50_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.9999876
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.50                   # NOT VALID?
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    pass
  
  def loadAgo18Thesis_FlatCK_levppus0p55_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.9805
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.55
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    pass
  
  def loadAgo18Thesis_FlatCK_levppus0p65_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.8632
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.65
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    pass
  
  def loadAgo18Thesis_FlatCK_levppus0p70_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.8193
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.70
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    pass
  
  def loadAgo18Thesis_FlatCK_levppus0p75_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.7805
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.75
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    pass
  
  def loadAgo18Thesis_FlatCK_levppus0p80_nom(self):
    self._opticsfile    = "opt_flat.madx"
    self._beta[0]       = [0.18, 0.075]
    self._beta[1]       = self._beta[0]
    #self._beta[2] remains the same
    self._sepLR[0]      = 11.43
    self._sepLR[1]      = self._sepLR[0]
    #self._sepLR[2] remains the same
    self._oncc[0]       = 0.7452
    self._oncc[1]       = self._oncc[0]
    #self._oncc[2] remains the same
    self._sigs          = 0.10
    self._minlong       = self._sigs
    self._levtech[0]    = ["levppus", "levlumi"]
    self._levtech[1]    = self._levtech[0]
    #self._levtech[2] remains the same
    self._levvar[0]     = ["PhiCK", "PhiCR"]
    self._levvar[1]     = self._levvar[0]
    #self._levvar[2] remains the same
    self._levppus[0]    = 0.80
    self._levppus[1]    = self._levppus[0]
    #self._levppus[2] remains the same
    self._ck[0]         = True
    self._ck[1]         = self._ck[0]
    #self._ck[2] remains the same
    self._niterlev      = 1
    self._longdens      = "RF800"
    pass
  





  
 
  
  #####################################################################################
  
  ### Check that list length of those parameters that are list coincides with the number of bunches. Not complete: some parameters are missing to be checked
  
  def checkconfig(self):
    
    if (len(self._nbunch) != self._nip):
      print "ERROR: wrong bunch definition:", self._nip, self._nbunch
    
    if (len(self._beta) != self._nip):
      print "ERROR: wrong beta definition:", self._nip, self._beta
    
    if (len(self._sepLR) != self._nip):
      print "ERROR: wrong sepLR definition:", self._nip, self._sepLR
    
    if (len(self._redsepLR) != self._nip):
      print "ERROR: wrong _redsepLR definition:", self._nip, self._redsepLR
      
    if (len(self._xplane) != self._nip):
      print "ERROR: wrong xplane definition:", self._nip, self._xplane
    
    if (len(self._levlumi) != self._nip):
      print "ERROR: wrong Level_Lumi definition:", self._nip, self._levlumi
    
    ### Values for this parameters are logical in the newes versions, but this keeps compatibility with older versions, where "0" or "1" was used instead:
    
    #self._longdens            = auxfuncclass.StringBinaryToLogical(self._longdens)  # Depracated since there are three options available
    self._constlong = auxfuncclass.StringBinaryToLogical(self._constlong)
    
    for i in range(len(self._constbetar)): self._constbetar[i] = auxfuncclass.StringBinaryToLogical(self._constbetar[i])

    if type(self._sepLRconst) == str:
      self._sepLRconst = [self._sepLRconst for x in range(self._nip)]
      
    if type(self._levtech) == str:
      print True
      self._levtech = [[self._levtech] for x in range(self._nip)]
    elif type(self._levtech) == list:
      for i in range(len(self._levtech)):
        if type(self._levtech[i]) == str:
          self._levtech[i] = [self._levtech[i]]
          
    if type(self._levvar) == str:
      self._levvar = [[self._levvar] for x in range(self._nip)]
    elif type(self._levvar) == list:
      for i in range(len(self._levvar)):
        if type(self._levvar[i]) == str:
          self._levvar[i] = [self._levvar[i]]

    if len(self._levtech) != len(self._levvar):
      print "ERROR: length of _levvar must match that of _levtech", self._levvar, self._levtech
