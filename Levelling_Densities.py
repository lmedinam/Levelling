import sys
import numpy as np
from metaclass import *
from Levelling_Others import Constants as cst
from Levelling_Others import AuxFunc
from Levelling_Others import qGaussianAux
from scipy.integrate import quad, dblquad

class Densities:
  
  def __init__(self):
    pass
  
  ### INTEGRANDS ###
  
  ## No CC ##
  
  def integrand_noCC_Gaussian(self, s, phi, sigs, betc, sigc, betp, t1, t2):
    """ Integrand in the longitudinal coordinate for a collission without crab cavities, for bunches with Gaussian longitudinal density. """
    sigs2  = sigs**2
    sigc2  = sigc**2
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    const  = 1./(np.sqrt(np.pi)*sigs)*np.exp(-(ct2-ct1)**2/4./sigs2)
    result = const*np.exp(
      - s*cosPHh * (s*cosPHh + (ct2-ct1)) / sigs2
      - ((s*sinPHh)**2) / sigc2/hgc2
      ) / hgc/hgp
    return result
  
  def integrand_noCC_parsep_Gaussian(self, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2):
    """ Integrand in the longitudinal coordinate for a collission without crab cavities, for bunches with Gaussian longitudinal density. Reduction factor from separation included. """
    sigs2  = sigs**2
    sigc2  = sigc**2
    sigp2  = sigp**2
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    const  = 1./(np.sqrt(np.pi)*sigs)*np.exp(-(ct2-ct1)**2/4./sigs2)
    result = const * np.exp(-parsep**2/4./sigp2/hgp2 ) * np.exp(
      - s*cosPHh * (s*cosPHh + (ct2-ct1)) / sigs2
      - ((s*sinPHh)**2) / sigc2/hgc2
      ) / hgc/hgp
    return result
  
  def integrand_noCC_RF800(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2):
    """ TIME OFFSETS DOES NOT WORK HERE. Integrand in the longitudinal and time coordinates for a collission without crab cavities, for bunches with flat (for RF800) longitudinal density. """
    sigs2  = sigs**2
    sigs4  = sigs**4
    sigc2  = sigc**2
    FactRMSGauss2 = cst.FactRMSGauss**2
    FactRMSGauss4 = cst.FactRMSGauss**4
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    const  = 2.0/(cst.NormRMSGauss**2 * FactRMSGauss2*sigs2) 
    result = const*np.exp(
        - ct**4 / FactRMSGauss4/sigs4
        - 6 * ct**2 * (s*cosPHh)**2 / FactRMSGauss4/sigs4
        - (s*cosPHh)**4 / FactRMSGauss4/sigs4
        - ((s*sinPHh)**2) / sigc2/hgc2
      ) / hgc/hgp
    return result
  
  def integrand_noCC_parsep_RF800(self, ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2):
    """ TIME OFFSETS DOES NOT WORK HERE. Integrand in the longitudinal and time coordinates for a collission without crab cavities, for bunches with flat (for RF800) longitudinal density. """
    sigs2  = sigs**2
    sigs4  = sigs**4
    sigc2  = sigc**2
    sigp2  = sigp**2
    FactRMSGauss2 = cst.FactRMSGauss**2
    FactRMSGauss4 = cst.FactRMSGauss**4
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    const  = 2.0/(cst.NormRMSGauss**2 * FactRMSGauss2*sigs2) 
    result = const * np.exp(-parsep**2/4./sigp2/hgp2 ) * np.exp(
        - ct**4 / FactRMSGauss4/sigs4
        - 6 * ct**2 * (s*cosPHh)**2 / FactRMSGauss4/sigs4
        - (s*cosPHh)**4 / FactRMSGauss4/sigs4
        - ((s*sinPHh)**2) / sigc2/hgc2
      ) / hgc/hgp
    return result

  def integrand_noCC_qGaussian(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2):
    """ Integrand in the crossing transversal coordinate, longitudinal and time coordinates for a collission without crab cavities, for bunches with flat (qGaussian) longitudinal density. """
    sigs2  = sigs**2
    sigs10 = sigs**10
    sigc2  = sigc**2
    FactFWHMGauss2  = cst.FactFWHMGauss**2
    FactFWHMGauss10 = cst.FactFWHMGauss**10
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    const  = 2./(cst.NormFWHMGauss**2 * FactFWHMGauss2*sigs2) * (2.**10/ FactFWHMGauss10/sigs10)
    result = const*(
         ((ct+ct1) - s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct1) - s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
      )**(5/2.) * np.exp(
        - ((s*sinPHh)**2) / sigc2/hgc2
      ) / hgc/hgp
    return result
  
  def integrandwithx_noCC_qGaussian(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, x):
    """ Integrand in the crossing transversal coordinate, longitudinal and time coordinates for a collission without crab cavities, for bunches with flat (qGaussian) longitudinal density. """
    sigs2  = sigs**2
    sigs10 = sigs**10
    sigc2  = sigc**2
    FactFWHMGauss2  = cst.FactFWHMGauss**2
    FactFWHMGauss10 = cst.FactFWHMGauss**10
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    const  = 2./(cst.NormFWHMGauss**2 * FactFWHMGauss2*sigs2) * (2.**10/ FactFWHMGauss10/sigs10) * (cosPHh/np.sqrt(np.pi)/sigc) # The constant changes since the integral in x is not performed
    result = const*(
         ((ct+ct1) - s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct1) - s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
      )**(5/2.) * np.exp(
        - ((s*sinPHh)**2) / sigc2/hgc2
        - ((x*cosPHh)**2) / sigc2/hgc2   # contributions from x integral
      ) / hgc2/hgp                       # mind the absence of sqrt in x, since the x integral was not done
    return result
 
  def integrandwithy_noCC_qGaussian(self, ct, s, phi, sigs, betc, sigc, betp, sigp, t1, t2, y):
    """ Integrand in the crossing transversal coordinate, longitudinal and time coordinates for a collission without crab cavities, for bunches with flat (qGaussian) longitudinal density. """
    sigs2  = sigs**2
    sigs10 = sigs**10
    sigc2  = sigc**2
    sigp2  = sigp**2
    FactFWHMGauss2  = cst.FactFWHMGauss**2
    FactFWHMGauss10 = cst.FactFWHMGauss**10
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    const  = 2./(cst.NormFWHMGauss**2 * FactFWHMGauss2*sigs2) * (2.**10/ FactFWHMGauss10/sigs10) * (1./np.sqrt(np.pi)/sigp) # The constant changes since the integral in y is not performed
    result = const*(
         ((ct+ct1) - s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct1) - s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
      )**(5/2.) * np.exp(
        - ((s*sinPHh)**2) / sigc2/hgc2
        - (y**2) / sigp2/hgp2           # contributions from y integral
      ) / hgc/hgp2                      # mind the absence of sqrt in y, since the y integral was not done
    return result
  
  def integrand_noCC_parsep_qGaussian(self, ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2):
    """ Integrand in the crossing transversal coordinate, longitudinal and time coordinates for a collission without crab cavities, for bunches with flat (qGaussian) longitudinal density. Reduction factor from separation included. """
    sigs2  = sigs**2
    sigs10 = sigs**10
    sigc2  = sigc**2
    sigp2  = sigp**2
    FactFWHMGauss2  = cst.FactFWHMGauss**2
    FactFWHMGauss10 = cst.FactFWHMGauss**10
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    const  = 2./(cst.NormFWHMGauss**2 * FactFWHMGauss2*sigs2) * (2.**10/ FactFWHMGauss10/sigs10)
    result = const * np.exp(-parsep**2/4./sigp2/hgp2 ) * (
         ((ct+ct1) - s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct1) - s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
      )**(5/2.) * np.exp(
        - ((s*sinPHh)**2) / sigc2/hgc2
      ) / hgc/hgp
    return result

  def integrandwithx_noCC_parsep_qGaussian(self, ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, x):
    """ Integrand in the crossing transversal coordinate, longitudinal and time coordinates for a collission without crab cavities, for bunches with flat (qGaussian) longitudinal density. Reduction factor from separation included. """
    sigs2  = sigs**2
    sigs10 = sigs**10
    sigc2  = sigc**2
    sigp2  = sigp**2
    FactFWHMGauss2  = cst.FactFWHMGauss**2
    FactFWHMGauss10 = cst.FactFWHMGauss**10
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    const  = 2./(cst.NormFWHMGauss**2 * FactFWHMGauss2*sigs2) * (2.**10/ FactFWHMGauss10/sigs10) * (cosPHh/np.sqrt(np.pi)/sigc) # The constant changes since the integral in x is not performed
    result = const * np.exp(-parsep**2/4./sigp2/hgp2 ) * (
         ((ct+ct1) - s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct1) - s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
      )**(5/2.) * np.exp(
        - ((s*sinPHh)**2) / sigc2/hgc2
        - ((x*cosPHh)**2) / sigc2/hgc2   # contributions from x integral
      ) / hgc2/hgp                       # mind the absence of sqrt in x, since the x integral was not done
    return result
  
  def integrandwithy_noCC_parsep_qGaussian(self, ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, y):
    """ Integrand in the crossing transversal coordinate, longitudinal and time coordinates for a collission without crab cavities, for bunches with flat (qGaussian) longitudinal density. Reduction factor from separation included. """
    sigs2  = sigs**2
    sigs10 = sigs**10
    sigc2  = sigc**2
    sigp2  = sigp**2
    FactFWHMGauss2  = cst.FactFWHMGauss**2
    FactFWHMGauss10 = cst.FactFWHMGauss**10
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    const  = 2./(cst.NormFWHMGauss**2 * FactFWHMGauss2*sigs2) * (2.**10/ FactFWHMGauss10/sigs10) * (1./np.sqrt(np.pi)/sigp) # The constant changes since the integral in y is not performed
    result = const * np.exp(-parsep**2/4./sigp2/hgp2 ) * (
         ((ct+ct1) - s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct1) - s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
      )**(5/2.) * np.exp(
        - ((s*sinPHh)**2) / sigc2/hgc2
        - (y**2) / sigp2/hgp2           # contributions from y integral
      ) / hgc/hgp2                      # mind the absence of sqrt in y, since the y integral was not done
    return result
  
  ## CC ##
  
  def integrand_CC_Gaussian(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, printcomp = None):
    """ Integrand in the longitudinal and time coordinates for a collission with crab cavities, for bunches with either Gaussian or flat longitudinal density. """
    kCR  = wcc/cst.clight*2*np.pi
    sigs2  = sigs**2
    sigc2  = sigc**2
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    sinCRh  = np.sin(phiCR/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    ct1c = cst.clight*t1c
    ct2c = cst.clight*t2c
    const  = 1./(np.pi*sigs2)
    if printcomp != None:
      print "t1 =", t1, ", t2 =", t2
      print "ct1 =", ct1, ", ct2 =", ct2
      print "ct =", ct, "s =", s
      print "ft(t)    =", - (ct**2 + ct*(ct2+ct1) + 0.5*(ct2**2+ct1**2)) / sigs2
      print "fs(s)    =", - s*cosPHh * (s*cosPHh + (ct2-ct1)) / sigs2
      print "Kxy(s,t) =", - (( (1./kCR)*sinCRh * np.sin(kCR*(s+0.5*(ct2c-ct1c))) * np.cos(kCR*(ct+0.5*(ct2c+ct1c))) - s*sinPHh )**2) / sigc2/hgc2
    result = const*np.exp(
        - (ct**2 + ct*(ct2+ct1) + 0.5*(ct2**2+ct1**2)) / sigs2
        - s*cosPHh * (s*cosPHh + (ct2-ct1)) / sigs2
        - (( (1./kCR)*sinCRh * np.sin(kCR*(s+0.5*(ct2c-ct1c))) * np.cos(kCR*(ct+0.5*(ct2c+ct1c))) - s*sinPHh )**2) / sigc2/hgc2
      ) / hgc/hgp
    return result
  
  def integrand_CC_parsep_Gaussian(self, ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, wcc, phiCR, t1c, t2c, printcomp = None):
    """ CHECK!!!!!!!! Integrand in the longitudinal and time coordinates for a collission with crab cavities, for bunches with either Gaussian or flat longitudinal density. """
    kCR  = wcc/cst.clight*2*np.pi
    sigs2  = sigs**2
    sigc2  = sigc**2
    sigp2  = sigp**2
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    sinCRh  = np.sin(phiCR/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    ct1c = cst.clight*t1c
    ct2c = cst.clight*t2c
    const  = 1./(np.pi*sigs2)
    if printcomp != None:
      print "t1 =", t1, ", t2 =", t2
      print "ct1 =", ct1, ", ct2 =", ct2
      print "ct =", ct, "s =", s
      print "ft(t)    =", - (ct**2 + ct*(ct2+ct1) + 0.5*(ct2**2+ct1**2)) / sigs2
      print "fs(s)    =", - s*cosPHh * (s*cosPHh + (ct2-ct1)) / sigs2
      print "Kxy(s,t) =", - (( (1./kCR)*sinCRh * np.sin(kCR*(s+0.5*(ct2c-ct1c))) * np.cos(kCR*(ct+0.5*(ct2c+ct1c))) - s*sinPHh )**2) / sigc2/hgc2
    result = const * np.exp(-parsep**2/4./sigp2/hgp2 ) * np.exp(
        - (ct**2 + ct*(ct2+ct1) + 0.5*(ct2**2+ct1**2)) / sigs2
        - s*cosPHh * (s*cosPHh + (ct2-ct1)) / sigs2
        - (( (1./kCR)*sinCRh * np.sin(kCR*(s+0.5*(ct2c-ct1c))) * np.cos(kCR*(ct+0.5*(ct2c+ct1c))) - s*sinPHh )**2) / sigc2/hgc2
      ) / hgc/hgp
    return result
  
  def integrandwithx_CC_Gaussian(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x):
    """ TIME OFFSETS DOES NOT WORK HERE. Integrand in the horizontal, longitudinal and time coordinates for a collission with crab cavities, for bunches with Gaussian longitudinal density. """
    # kCR  = wcc/cst.clight*2*np.pi
    # sigs2  = sigs**2
    # sigc2  = sigc**2
    # hgc2 = 1. + (s/betc)**2
    # hgp2 = 1. + (s/betp)**2
    # hgc = np.sqrt(hgc2)
    # hgp = np.sqrt(hgp2)
    # cosPHh  = np.cos(phi/2.)
    # sinPHh  = np.sin(phi/2.)
    # sinCRh  = np.sin(phiCR/2.)
    # ct1 = cst.clight*t1
    # ct2 = cst.clight*t2
    # ct1c = cst.clight*t1c
    # ct2c = cst.clight*t2c
    # const  = 1./(np.pi*sigs2) * cosPHh/(sigc*np.sqrt(np.pi)) # The constant changes since the integral in x is not performed
    # result = const*np.exp(
    #     - (ct**2 + ct*(ct2+ct1) + 0.5*(ct2**2+ct1**2)) / sigs2
    #     - s*cosPHh * (s*cosPHh + (ct2-ct1)) / sigs2
    #     - (( (1./kCR)*sinCRh * np.sin(kCR*(s+0.5*(ct2c-ct1c))) * np.cos(kCR*(ct+0.5*(ct2c+ct1c))) - s*sinPHh )**2) / sigc2/hgc2
    #     - (x*cosPHh)**2 / sigc2/hgc2 - (x*sinPHh)**2 / sigs2                        # contributions from x integral
    #     + x*ct*sinPHh/sigs2                                                         # contributions from x integral
    #     + (2*x*cosPHh) * np.cos(kCR*s) * np.sin(kCR*ct) * sinCRh / kCR / sigc2/hgc2 # contributions from x
    #   ) / hgc2/hgp
    # return result

    # OLD Formulae -- there seems to be an error in the previous one, so I am using the old one...
    sigs2  = sigs**2
    sigc2  = sigc**2
    k  = wcc/cst.clight*2*np.pi
    k2 = k**2
    cos  = np.cos(phi/2.)
    sin  = np.sin(phi/2.)
    cos2 = cos**2
    sin2 = sin**2
    sinCR   = np.sin(phiCR/2.)
    sinCR2  = sinCR**2
    sigmac2 = sigc2*(1+s**2/betc**2)
    const  = cos*np.sqrt(np.pi)/(np.pi**2 *sigc*sigs2) # The constant changes since the integral in x is not performed
    result = const*np.exp(
      -ct**2/sigs2
      -s**2*cos2/sigs2
      -1.0/(4*k2*sigmac2)*(
          # -4*np.cos(k*s)**2 * np.sin(k*ct)**2 * sinCR2  # this term would be the result from the integral in x, it is substituted with the terms in x to be integrated.
          -8*s*k*np.sin(k*s)*np.cos(k*ct) * sin*sinCR
          +2 * sinCR2
          -np.cos(2*k*(s-ct)) * sinCR2
          -np.cos(2*k*(s+ct)) * sinCR2
          +4*k2*s**2 * sin2 )
      -x**2*(cos2/sigmac2+sin2/sigs2)                   # contributions from x integral
      +x*ct*sin/sigs2                                   # contributions from x integral -- a missing factor of 2???
      +2*x*cos*np.cos(k*s)*np.sin(k*ct)*sinCR/k/sigmac2 # contributions from x
      ) / (1+s**2/betc**2) / np.sqrt(1+s**2/betp**2)  # mind the absence of sqrt in x, since the x integral was not done
    return result
  
  def integrandwithy_CC_Gaussian(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, sigp, t1p, t2p, y):
    """ TIME OFFSETS DOES NOT WORK HERE. Integrand in the vertical, longitudinal and time coordinates for a collission with crab cavities, for bunches with Gaussian longitudinal density. """
    kCR  = wcc/cst.clight*2*np.pi
    sigs2  = sigs**2
    sigc2  = sigc**2
    sigp2  = sigp**2
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    sinCRh  = np.sin(phiCR/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    ct1c = cst.clight*t1c
    ct2c = cst.clight*t2c
    const  = 1./(np.pi*sigs2) * 1./(sigp*np.sqrt(np.pi)) # The constant changes since the integral in y is not performed
    result = const*np.exp(
        - (ct**2 + ct*(ct2+ct1) + 0.5*(ct2**2+ct1**2)) / sigs2
        - s*cosPHh * (s*cosPHh + (ct2-ct1)) / sigs2
        - (( (1./kCR)*sinCRh * np.sin(kCR*(s+0.5*(ct2c-ct1c))) * np.cos(kCR*(ct+0.5*(ct2c+ct1c))) - s*sinPHh )**2) / sigc2/hgc2
        -y**2 / sigp2/hgp2
      ) / hgc/hgp2
    return result
  
    # # OLD Formulae:
    # sigs2  = sigs**2
    # sigc2  = sigc**2
    # sigp2  = sigp**2
    # k  = wcc/cst.clight*2*np.pi
    # k2 = k**2
    # cos  = np.cos(phi/2.)
    # sin  = np.sin(phi/2.)
    # cos2 = cos**2
    # sin2 = sin**2
    # sinCR   = np.sin(phiCR/2.)
    # sinCR2  = sinCR**2
    # sigmac2 = sigc2*(1+s**2/betc**2)
    # sigmap2 = sigp2*(1+s**2/betp**2)
    # const  = np.sqrt(np.pi)/(np.pi**2 *sigp*sigs2) # The constant changes since the integral in y is not performed
    # result = const*np.exp(
    #   -y**2/sigmap2
    #   -ct**2/sigs2
    #   -s**2*cos2/sigs2
    #   -1.0/(4*k2*sigmac2)*(
    #       -4*np.cos(k*s)**2 * np.sin(k*ct)**2 * sinCR2
    #       -8*s*k*np.sin(k*s)*np.cos(k*ct) * sin*sinCR
    #       +2 * sinCR2
    #       -np.cos(2*k*(s-ct)) * sinCR2
    #       -np.cos(2*k*(s+ct)) * sinCR2
    #       +4*k2*s**2 * sin2 )
    #   ) / np.sqrt(1+s**2/betc**2) / (1+s**2/betp**2)  # mind the absence of sqrt in y, since the y integral was not done
    # return result
  
  def integrand_CC_RF800(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c):
    """ TIME OFFSETS DOES NOT WORK HERE. Integrand in the longitudinal and time coordinates for a collission with crab cavities, for bunches with flat (for RF800) longitudinal density. """
    kCR  = wcc/cst.clight*2*np.pi
    sigs2  = sigs**2
    sigs4  = sigs**4
    sigc2  = sigc**2
    FactRMSGauss2 = cst.FactRMSGauss**2
    FactRMSGauss4 = cst.FactRMSGauss**4
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    sinCRh  = np.sin(phiCR/2.)
    const  = 2.0/(cst.NormRMSGauss**2 * FactRMSGauss2*sigs2) 
    result = const*np.exp(
        - ct**4 / FactRMSGauss4/sigs4
        - 6 * ct**2 * (s*cosPHh)**2 / FactRMSGauss4/sigs4
        - (s*cosPHh)**4 / FactRMSGauss4/sigs4
        - (( (1./kCR)*sinCRh * np.sin(kCR*(s)) * np.cos(kCR*(ct)) - s*sinPHh )**2) / sigc2/hgc2
      ) / hgc/hgp
    return result

  def integrand_CC_parsep_RF800(self, ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, wcc, phiCR, t1c, t2c):
    """ TIME OFFSETS DOES NOT WORK HERE. Integrand in the longitudinal and time coordinates for a collission with crab cavities, for bunches with flat (for RF800) longitudinal density. """
    kCR  = wcc/cst.clight*2*np.pi
    sigs2  = sigs**2
    sigs4  = sigs**4
    sigc2  = sigc**2
    sigp2  = sigp**2
    FactRMSGauss2 = cst.FactRMSGauss**2
    FactRMSGauss4 = cst.FactRMSGauss**4
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    sinCRh  = np.sin(phiCR/2.)
    const  = 2.0/(cst.NormRMSGauss**2 * FactRMSGauss2*sigs2) 
    result = const * np.exp(-parsep**2/4./sigp2/hgp2 ) * np.exp(
        - ct**4 / FactRMSGauss4/sigs4
        - 6 * ct**2 * (s*cosPHh)**2 / FactRMSGauss4/sigs4
        - (s*cosPHh)**4 / FactRMSGauss4/sigs4
        - (( (1./kCR)*sinCRh * np.sin(kCR*(s)) * np.cos(kCR*(ct)) - s*sinPHh )**2) / sigc2/hgc2
      ) / hgc/hgp
    return result
  
  def integrand_CC_qGaussian(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c):
    """ Integrand in the longitudinal and time coordinates for a collission with crab cavities, for bunches with flat (qGaussian) longitudinal density. """
    kCR  = wcc/cst.clight*2*np.pi
    sigs2  = sigs**2
    sigs10 = sigs**10
    sigc2  = sigc**2
    FactFWHMGauss2  = cst.FactFWHMGauss**2
    FactFWHMGauss10 = cst.FactFWHMGauss**10
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    sinCRh  = np.sin(phiCR/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    ct1c = cst.clight*t1c
    ct2c = cst.clight*t2c
    const  = 2./(cst.NormFWHMGauss**2 * FactFWHMGauss2*sigs2) * (2.**10/ FactFWHMGauss10/sigs10)
    result = const*(
         ((ct+ct1) - s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct1) - s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
      )**(5/2.) * np.exp(
        - (( (1./kCR)*sinCRh * np.sin(kCR*(s+0.5*(ct2c-ct1c))) * np.cos(kCR*(ct+0.5*(ct2c+ct1c))) - s*sinPHh )**2) / sigc2/hgc2
      ) / hgc/hgp
    return result
  
  def integrandwithx_CC_qGaussian(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x):
    """ Case added on 14/Jun/18 -- t1c and t2c don't work!!! """
    kCR  = wcc/cst.clight*2*np.pi
    sigs2  = sigs**2
    sigs10 = sigs**10
    sigc2  = sigc**2
    FactFWHMGauss2  = cst.FactFWHMGauss**2
    FactFWHMGauss10 = cst.FactFWHMGauss**10
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    sinCRh  = np.sin(phiCR/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    ct1c = cst.clight*t1c
    ct2c = cst.clight*t2c
    const  = 2./(cst.NormFWHMGauss**2 * FactFWHMGauss2*sigs2) * (2.**10/ FactFWHMGauss10/sigs10) * (cosPHh/np.sqrt(np.pi)/sigc) # The constant changes since the integral in x is not performed
    result = const*(
         ((ct+ct1) - (s*cosPHh + x*sinPHh) - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct1) - (s*cosPHh + x*sinPHh) + 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + (s*cosPHh - x*sinPHh) - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + (s*cosPHh - x*sinPHh) + 0.5*cst.FactFWHMGauss*sigs)
      )**(5/2.) * np.exp(
        - (( (1./kCR)*sinCRh * np.sin(kCR*(s+0.5*(ct2c-ct1c))) * np.cos(kCR*(ct+0.5*(ct2c+ct1c))) - s*sinPHh )**2) / sigc2/hgc2
        - ((x*cosPHh)**2) / sigc2/hgc2   # contributions from x integral
        + (2*x*cosPHh) * np.cos(kCR*s) * np.sin(kCR*ct) * sinCRh / kCR / sigc2/hgc2   # contributions from x integral   # ct1c y ct2 are missing from been added here, but it should be fine if they are zero!!!
      ) / hgc2/hgp                       # mind the absence of sqrt in x, since the x integral was not done
    return result
  
  def integrandwithy_CC_qGaussian(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, sigp, t1p, t2p, y):
    """ Case added on 14/Jun/18  -- t1c and t2c and t1p and t2p don't work!!!  """
    kCR  = wcc/cst.clight*2*np.pi
    sigs2  = sigs**2
    sigs10 = sigs**10
    sigc2  = sigc**2
    sigp2  = sigp**2
    FactFWHMGauss2  = cst.FactFWHMGauss**2
    FactFWHMGauss10 = cst.FactFWHMGauss**10
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    sinCRh  = np.sin(phiCR/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    ct1c = cst.clight*t1c
    ct2c = cst.clight*t2c
    const  = 2./(cst.NormFWHMGauss**2 * FactFWHMGauss2*sigs2) * (2.**10/ FactFWHMGauss10/sigs10) * (1./np.sqrt(np.pi)/sigp) # The constant changes since the integral in y is not performed
    result = const*(
         ((ct+ct1) - s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct1) - s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
      )**(5/2.) * np.exp(
        - (( (1./kCR)*sinCRh * np.sin(kCR*(s+0.5*(ct2c-ct1c))) * np.cos(kCR*(ct+0.5*(ct2c+ct1c))) - s*sinPHh )**2) / sigc2/hgc2
        - (y**2) / sigp2/hgp2           # contributions from y integral
      ) / hgc/hgp2                      # mind the absence of sqrt in y, since the y integral was not done
    return result
  
  #def integrand_CC_parsep_qGaussian is missing for completeness, but it is not a case that is needed
  
  ## CC and CK ##
  
  def integrand_CC_CK_Gaussian(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p):
    """ Integrand in the longitudinal and time coordinates for a collission with crab cavities and crab kissing, for bunches with Gaussian longitudinal density. """
    kCR  = wcc/cst.clight*2*np.pi
    kCK  = kCR
    sigs2  = sigs**2
    sigc2  = sigc**2
    sigp2  = sigp**2
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    sinCRh  = np.sin(phiCR/2.)
    sinCKh  = np.sin(phiCK/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    ct1c = cst.clight*t1c
    ct2c = cst.clight*t2c
    ct1p = cst.clight*t1p
    ct2p = cst.clight*t2p
    const  = 1./(np.pi*sigs2)
    result = const*np.exp(
        - (ct**2 + ct*(ct2+ct1) + 0.5*(ct2**2+ct1**2)) / sigs2
        - s*cosPHh * (s*cosPHh + (ct2-ct1)) / sigs2
        - (( (1./kCR)*sinCRh * np.sin(kCR*(s+0.5*(ct2c-ct1c))) * np.cos(kCR*(ct+0.5*(ct2c+ct1c))) - s*sinPHh )**2) / sigc2/hgc2
        - (( (1./kCK)*sinCKh * np.cos(kCK*(s+0.5*(ct2p-ct1p))) * np.sin(kCK*(ct+0.5*(ct2p+ct1p)))            )**2) / sigp2/hgp2
      ) / hgc/hgp
    return result
  
  def integrand_CC_CK_RF800(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p):
    """ TIME OFFSETS DOES NOT WORK HERE. Integrand in the longitudinal and time coordinates for a collission with crab cavities and crab kissing, for bunches with flat (for RF800) longitudinal density. """
    kCR  = wcc/cst.clight*2*np.pi
    kCK  = kCR
    sigs2  = sigs**2
    sigs4  = sigs**4
    sigc2  = sigc**2
    sigp2  = sigp**2
    FactRMSGauss2 = cst.FactRMSGauss**2
    FactRMSGauss4 = cst.FactRMSGauss**4
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    sinCRh  = np.sin(phiCR/2.)
    sinCKh  = np.sin(phiCK/2.)
    const  = 2.0/(cst.NormRMSGauss**2 * FactRMSGauss2*sigs2) 
    result = const*np.exp(
        - ct**4 / FactRMSGauss4/sigs4
        - 6 * ct**2 * (s*cosPHh)**2 / FactRMSGauss4/sigs4
        - (s*cosPHh)**4 / FactRMSGauss4/sigs4
        - (( (1./kCR)*sinCRh * np.sin(kCR*(s)) * np.cos(kCR*(ct)) - s*sinPHh )**2) / sigc2/hgc2
        - (( (1./kCK)*sinCKh * np.cos(kCK*(s)) * np.sin(kCK*(ct))            )**2) / sigp2/hgp2
      ) / hgc/hgp
    return result
  
  def integrand_CC_CK_qGaussian(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p):
    """ Integrand in the longitudinal and time coordinates for a collission without crab cavities, for bunches with flat (qGaussian) longitudinal density. """
    kCR  = wcc/cst.clight*2*np.pi
    kCK  = kCR
    sigs2  = sigs**2
    sigs10 = sigs**10
    sigc2  = sigc**2
    sigp2  = sigp**2
    FactFWHMGauss2  = cst.FactFWHMGauss**2
    FactFWHMGauss10 = cst.FactFWHMGauss**10
    hgc2 = 1. + (s/betc)**2
    hgp2 = 1. + (s/betp)**2
    hgc = np.sqrt(hgc2)
    hgp = np.sqrt(hgp2)
    cosPHh  = np.cos(phi/2.)
    sinPHh  = np.sin(phi/2.)
    sinCRh  = np.sin(phiCR/2.)
    sinCKh  = np.sin(phiCK/2.)
    ct1 = cst.clight*t1
    ct2 = cst.clight*t2
    ct1c = cst.clight*t1c
    ct2c = cst.clight*t2c
    ct1p = cst.clight*t1p
    ct2p = cst.clight*t2p
    const  = 2./(cst.NormFWHMGauss**2 * FactFWHMGauss2*sigs2) * (2.**10/ FactFWHMGauss10/sigs10)
    result = const*(
         ((ct+ct1) - s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct1) - s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh - 0.5*cst.FactFWHMGauss*sigs)
        *((ct+ct2) + s*cosPHh + 0.5*cst.FactFWHMGauss*sigs)
      )**(5/2.) * np.exp(
        - (( (1./kCR)*sinCRh * np.sin(kCR*(s+0.5*(ct2c-ct1c))) * np.cos(kCR*(ct+0.5*(ct2c+ct1c))) - s*sinPHh )**2) / sigc2/hgc2
        - (( (1./kCK)*sinCKh * np.cos(kCK*(s+0.5*(ct2p-ct1p))) * np.sin(kCK*(ct+0.5*(ct2p+ct1p)))            )**2) / sigp2/hgp2
      ) / hgc/hgp
    return result
  
  ### EVALUATION OF INTEGRAND AND INTEGRALS ###
  
  def Evaluate_integrand(self, ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep, printcomp = None):
    if ip == 2:
      if oncc == 0.0:
        if longdens == "Gaussian":
          result = self.integrand_CC_parsep_Gaussian(ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, wcc,   0.0, 0.0, 0.0, printcomp)
        elif longdens == "RF800":
          result = self.integrand_noCC_parsep_RF800( ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2)
        elif longdens == "qGaussian":
          result = self.integrand_noCC_parsep_qGaussian(ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2)
        else:
          sys.exit("\n[!] ERROR.\n")
      else:
        sys.exit("\n[!] SETTINGS HAVE TO BE ONCC = 0 FOR IP8.\n")
    else:
      if ck == False:
        if oncc == 0.0:
          if longdens == "Gaussian":
            result = self.integrand_CC_Gaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc,   0.0, 0.0, 0.0, printcomp)
          elif longdens == "RF800":
            result = self.integrand_noCC_RF800( ct, s, phi, sigs, betc, sigc, betp, t1, t2)
          elif longdens == "qGaussian":
            result = self.integrand_noCC_qGaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
        else:
          if longdens == "Gaussian":
            result = self.integrand_CC_Gaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, printcomp)
          elif longdens == "RF800":
            result = self.integrand_CC_RF800(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c)
          elif longdens == "qGaussian":
            result = self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          result = self.integrand_CC_CK_Gaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p)
        elif longdens == "RF800":
          result = self.integrand_CC_CK_RF800(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p)
        elif longdens == "qGaussian":
          result = self.integrand_CC_CK_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p)
        else:
          sys.exit("\n[!] Longitudinal density unknown.\n")
    return result # [1/m^2]
  
  def IntCT(self, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep):
    if ip == 2:
      if oncc == 0.0:
        if longdens == "Gaussian":
          result = quad(lambda ct:    self.integrand_CC_parsep_Gaussian( ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, wcc,   0.0, 0.0, 0.0), -np.inf, np.inf)
        elif longdens == "RF800":
            result = quad(lambda ct:    self.integrand_noCC_parsep_RF800(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       -np.inf, np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          if  s >= qGAux.szero:
            result1 = quad(lambda ct:    self.integrand_noCC_parsep_qGaussian( ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.Zct(s), qGAux.Cct(s))
            result4 = quad(lambda ct:    self.integrand_noCC_parsep_qGaussian( ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.Bct(s), qGAux.Zct(s))
            result = (result1[0]+result4[0],)
          else:
            result2 = quad(lambda ct:    self.integrand_noCC_parsep_qGaussian( ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.Zct(s), qGAux.Act(s))
            result3 = quad(lambda ct:    self.integrand_noCC_parsep_qGaussian( ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.Dct(s), qGAux.Zct(s))
            result = (result2[0]+result3[0],)
        else:
          sys.exit("\n[!] ERROR.\n")
      else:
        sys.exit("\n[!] SETTINGS HAVE TO BE ONCC = 0 FOR IP8.\n")
    else:
      if ck == False:
        if oncc == 0.0:
          if longdens == "Gaussian":
            result = quad(lambda ct:    self.integrand_CC_Gaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc,   0.0, 0.0, 0.0), -np.inf, np.inf)
          elif longdens == "RF800":
            result = quad(lambda ct:    self.integrand_noCC_RF800(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       -np.inf, np.inf)
          elif longdens == "qGaussian":
           # print 'you are here in intCT'
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            if  s >= qGAux.szero:
              result1 = quad(lambda ct:    self.integrand_noCC_qGaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.Zct(s), qGAux.Cct(s))
              result4 = quad(lambda ct:    self.integrand_noCC_qGaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.Bct(s), qGAux.Zct(s))
              result = (result1[0]+result4[0],)
            else:
              result2 = quad(lambda ct:    self.integrand_noCC_qGaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.Zct(s), qGAux.Act(s))
              result3 = quad(lambda ct:    self.integrand_noCC_qGaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.Dct(s), qGAux.Zct(s))
              result = (result2[0]+result3[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
        else:
          if longdens == "Gaussian":
            result = quad(lambda ct:    self.integrand_CC_Gaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf)
          elif longdens == "RF800":
            result = quad(lambda ct:    self.integrand_CC_RF800(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            if  s >= qGAux.szero:
              result1 = quad(lambda ct:    self.integrand_CC_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.Zct(s), qGAux.Cct(s))
              result4 = quad(lambda ct:    self.integrand_CC_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.Bct(s), qGAux.Zct(s))
              result = (result1[0]+result4[0],)
            else:
              result2 = quad(lambda ct:    self.integrand_CC_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.Zct(s), qGAux.Act(s))
              result3 = quad(lambda ct:    self.integrand_CC_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.Dct(s), qGAux.Zct(s))
              result = (result2[0]+result3[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          result = quad(lambda ct:    self.integrand_CC_CK_Gaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf)
        elif longdens == "RF800":
          result = quad(lambda ct:    self.integrand_CC_CK_RF800(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          if  s >= qGAux.szero:
            result1 = quad(lambda ct:    self.integrand_CC_CK_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.Zct(s), qGAux.Cct(s))
            result4 = quad(lambda ct:    self.integrand_CC_CK_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.Bct(s), qGAux.Zct(s))
            result = (result1[0]+result4[0],)
          else:
            result2 = quad(lambda ct:    self.integrand_CC_CK_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.Zct(s), qGAux.Act(s))
            result3 = quad(lambda ct:    self.integrand_CC_CK_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.Dct(s), qGAux.Zct(s))
            result = (result2[0]+result3[0],)
        else:
          sys.exit("\n[!] Longitudinal density unknown.\n")
    return result[0] # [1/m]
  
  def IntS(self, ct, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep):
    if ip == 2:
      if oncc == 0.0:
        if longdens == "Gaussian":
          result = quad(lambda s:    self.integrand_CC_parsep_Gaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, wcc,   0.0, 0.0, 0.0), -np.inf, np.inf)
        elif longdens == "RF800":
            result = quad(lambda s:    self.integrand_noCC_parsep_RF800(   ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       -np.inf, np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          if ct >= qGAux.ctzero:
            result1 = quad(lambda s:    self.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.Zs(ct), qGAux.Cs(ct))
            result4 = quad(lambda s:    self.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.As(ct), qGAux.Zs(ct))
            result = (result1[0]+result4[0],)
          else:
            result2 = quad(lambda s:    self.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.Zs(ct), qGAux.Bs(ct))
            result3 = quad(lambda s:    self.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                       qGAux.Ds(ct), qGAux.Zs(ct))
            result = (result2[0]+result3[0],)
        else:
          sys.exit("\n[!] ERROR.\n")
      else:
        sys.exit("\n[!] SETTINGS HAVE TO BE ONCC = 0 FOR IP8.\n")
    else:
      if ck == False:
        if oncc == 0.0:
          if longdens == "Gaussian":
            result = quad(lambda s:    self.integrand_CC_Gaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc,   0.0, 0.0, 0.0), -np.inf, np.inf)
          elif longdens == "RF800":
            result = quad(lambda s:    self.integrand_noCC_RF800(   ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       -np.inf, np.inf)
          elif longdens == "qGaussian":
            #print 'you are here in intS'
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            if ct >= qGAux.ctzero:
              result1 = quad(lambda s:    self.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.Zs(ct), qGAux.Cs(ct))
              result4 = quad(lambda s:    self.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.As(ct), qGAux.Zs(ct))
              result = (result1[0]+result4[0],)
            else:
              result2 = quad(lambda s:    self.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.Zs(ct), qGAux.Bs(ct))
              result3 = quad(lambda s:    self.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                       qGAux.Ds(ct), qGAux.Zs(ct))
              result = (result2[0]+result3[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
        else:
          if longdens == "Gaussian":
            result = quad(lambda s:    self.integrand_CC_Gaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf)
          elif longdens == "RF800":
            result = quad(lambda s:    self.integrand_CC_RF800(     ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            if ct >= qGAux.ctzero:
              result1 = quad(lambda s:    self.integrand_CC_qGaussian(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.Zs(ct), qGAux.Cs(ct))
              result4 = quad(lambda s:    self.integrand_CC_qGaussian(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.As(ct), qGAux.Zs(ct))
              result = (result1[0]+result4[0],)
            else:
              result2 = quad(lambda s:    self.integrand_CC_qGaussian(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.Zs(ct), qGAux.Bs(ct))
              result3 = quad(lambda s:    self.integrand_CC_qGaussian(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.Ds(ct), qGAux.Zs(ct))
              result = (result2[0]+result3[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          result = quad(lambda s:    self.integrand_CC_CK_Gaussian( ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf)
        elif longdens == "RF800":
          result = quad(lambda s:    self.integrand_CC_CK_RF800(    ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf)
        elif longdens == "qGaussian":
          if ct >= qGAux.ctzero:
            result1 = quad(lambda s:    self.integrand_CC_CK_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.Zs(ct), qGAux.Cs(ct))
            result4 = quad(lambda s:    self.integrand_CC_CK_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.As(ct), qGAux.Zs(ct))
            result = (result1[0]+result4[0],)
          else:
            result2 = quad(lambda s:    self.integrand_CC_CK_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.Zs(ct), qGAux.Bs(ct))
            result3 = quad(lambda s:    self.integrand_CC_CK_qGaussian(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.Ds(ct), qGAux.Zs(ct))
            result = (result2[0]+result3[0],)
        else:
          sys.exit("\n[!] Longitudinal density unknown.\n")
    return result[0] # [1/m], to be converted into 1/s afterwards
  
  def IntAll(self, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep):
    if ip == 2:
      if oncc == 0.0:
        if longdens == "Gaussian":
          result =    quad(lambda s    :    self.integrand_noCC_parsep_Gaussian(    s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                     -np.inf, np.inf)
        elif longdens == "RF800":
            result = dblquad(lambda s, ct:    self.integrand_noCC_parsep_RF800(   ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                     -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          result1 = dblquad(lambda s, ct:    self.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                     qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs )
          result2 = dblquad(lambda s, ct:    self.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                     qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs )
          result3 = dblquad(lambda s, ct:    self.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                     qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs )
          result4 = dblquad(lambda s, ct:    self.integrand_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2),                     qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs )
          result = (result1[0]+result2[0]+result3[0]+result4[0],)
        else:
          sys.exit("\n[!] ERROR.\n")
      else:
        sys.exit("\n[!] SETTINGS HAVE TO BE ONCC = 0 FOR IP8.\n")
    else:
      if ck == False:
        if oncc == 0.0:
          if longdens == "Gaussian":
            result =    quad(lambda s    :    self.integrand_noCC_Gaussian(    s, phi, sigs, betc, sigc, betp, t1, t2),                     -np.inf, np.inf)
          elif longdens == "RF800":
            result = dblquad(lambda s, ct:    self.integrand_noCC_RF800(   ct, s, phi, sigs, betc, sigc, betp, t1, t2),                     -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            result1 = dblquad(lambda s, ct:    self.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                     qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs )
            result2 = dblquad(lambda s, ct:    self.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                     qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs )
            result3 = dblquad(lambda s, ct:    self.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                     qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs )
            result4 = dblquad(lambda s, ct:    self.integrand_noCC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2),                     qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs )
            result = (result1[0]+result2[0]+result3[0]+result4[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
        else:
          if longdens == "Gaussian":
            result = dblquad(lambda s, ct:    self.integrand_CC_Gaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
          elif longdens == "RF800":
            result = dblquad(lambda s, ct:    self.integrand_CC_RF800(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            result1 = dblquad(lambda s, ct:    self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs )
            result2 = dblquad(lambda s, ct:    self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs )
            result3 = dblquad(lambda s, ct:    self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs )
            result4 = dblquad(lambda s, ct:    self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs )
            result = (result1[0]+result2[0]+result3[0]+result4[0],)
            # # Old: (has been reduced by a factor of 2)
            #result = dblquad(lambda ct, s: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), 0.0, 0.5*cst.FactFWHMGauss*sigs/np.cos(phi/2.), lambda u: -0.5*cst.FactFWHMGauss*sigs+u*np.cos(phi/2.), lambda u: 0.5*cst.FactFWHMGauss*sigs-u*np.cos(phi/2.)) 
            # 
            # # New: (extensive):
            # 
            # # print "ctmax =", qGAux.ctmax
            # # print "ctmin =", qGAux.ctmin
            # # print "smax  =", qGAux.smax
            # # print "smin  =", qGAux.smin
            # 
            # print "\nCT S\n"
            # resultcts1 = dblquad(lambda ct, s: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.szero, qGAux.smax,  qGAux.Zct, qGAux.Cct)
            # resultcts2 = dblquad(lambda ct, s: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.smin,  qGAux.szero, qGAux.Zct, qGAux.Act)
            # resultcts3 = dblquad(lambda ct, s: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.smin,  qGAux.szero, qGAux.Dct, qGAux.Zct)
            # resultcts4 = dblquad(lambda ct, s: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.szero, qGAux.smax,  qGAux.Bct, qGAux.Zct)
            # print "resultcts1 =", resultcts1
            # print "resultcts2 =", resultcts2
            # print "resultcts3 =", resultcts3
            # print "resultcts4 =", resultcts4
            # print ""
            # 
            # resultcts41 = dblquad(lambda ct, s: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.szero, qGAux.smax,  qGAux.Bct, qGAux.Cct)
            # resultcts32 = dblquad(lambda ct, s: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.smin,  qGAux.szero, qGAux.Dct, qGAux.Act)
            # print "resultcts41 =", resultcts41
            # print "resultcts32 =", resultcts32
            # print ""
            # 
            # print "\nS CT\n"
            # resultsct1 = dblquad(lambda s, ct: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs)
            # resultsct2 = dblquad(lambda s, ct: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs)
            # resultsct3 = dblquad(lambda s, ct: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs)
            # resultsct4 = dblquad(lambda s, ct: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs)
            # print "resultsct1 =", resultsct1
            # print "resultsct2 =", resultsct2
            # print "resultsct3 =", resultsct3
            # print "resultsct4 =", resultsct4
            # print ""
            # 
            # resultsct41 = dblquad(lambda s, ct: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.szero,  qGAux.smax,   qGAux.As,  qGAux.Cs )
            # resultsct32 = dblquad(lambda s, ct: self.integrand_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c), qGAux.smin,   qGAux.szero,  qGAux.Ds,  qGAux.Bs )
            # print "resultsct41 =", resultsct41
            # print "resultsct32 =", resultsct32
            # print ""
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          result = dblquad(lambda s, ct:    self.integrand_CC_CK_Gaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
        elif longdens == "RF800":
          result = dblquad(lambda s, ct:    self.integrand_CC_CK_RF800(   ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), -np.inf, np.inf,     lambda u: -np.inf, lambda u: np.inf)
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          result1 = dblquad(lambda s, ct:    self.integrand_CC_CK_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs )
          result2 = dblquad(lambda s, ct:    self.integrand_CC_CK_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs )
          result3 = dblquad(lambda s, ct:    self.integrand_CC_CK_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs )
          result4 = dblquad(lambda s, ct:    self.integrand_CC_CK_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p), qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs )
          result = (result1[0]+result2[0]+result3[0]+result4[0],)
        else:
          sys.exit("\n[!] Longitudinal density unknown.\n")
    return result[0] # [1]
  
  def IntAll_withx(self, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep, x):
    if ip == 2:
      if oncc == 0.0:
        if longdens == "Gaussian":
          sys.exit("Not ready")
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          result1 = dblquad(lambda s, ct:    self.integrandwithx_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, x),                     qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs )
          result2 = dblquad(lambda s, ct:    self.integrandwithx_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, x),                     qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs )
          result3 = dblquad(lambda s, ct:    self.integrandwithx_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, x),                     qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs )
          result4 = dblquad(lambda s, ct:    self.integrandwithx_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, x),                     qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs )
          result = (result1[0]+result2[0]+result3[0]+result4[0],)
        else:
          sys.exit("\n[!] SETTINGS HAVE TO BE GAUSSIAN/QGAUSSIAN FOR IP8.\n")
      else:
        sys.exit("\n[!] SETTINGS HAVE TO BE ONCC = 0 FOR IP8.\n")
    else:
      if ck == False:
        if oncc == 0.0:
          if longdens == "Gaussian":
            result = dblquad(self.integrandwithx_CC_Gaussian, -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, 0.0, t1c, t2c, x)) # READY
          elif longdens == "RF800":
            # result = dblquad(self.integrandwithx_noCC_RF800,    -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, t1c, t2c, t1p, t2p, x))
            sys.exit("Not ready")
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2, x)
            result1 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, 0.0, t1c, t2c, x),  qGAux.ctzero, qGAux.ctmax,  qGAux.Zswithx,  qGAux.Cswithx ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result2 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, 0.0, t1c, t2c, x),  qGAux.ctmin,  qGAux.ctzero, qGAux.Zswithx,  qGAux.Bswithx ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result3 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, 0.0, t1c, t2c, x),  qGAux.ctmin,  qGAux.ctzero, qGAux.Dswithx,  qGAux.Zswithx ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result4 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, 0.0, t1c, t2c, x),  qGAux.ctzero, qGAux.ctmax,  qGAux.Aswithx,  qGAux.Zswithx ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result = (result1[0]+result2[0]+result3[0]+result4[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
        else:
          if longdens == "Gaussian":
            result = dblquad(self.integrandwithx_CC_Gaussian, -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x)) # READY
          elif longdens == "RF800":
            # result = dblquad(self.integrandwithx_CC_RF800,    -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x))
            sys.exit("Not ready")
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2, x)
            #print x, phi, ' ', 
            #print '\n', x, qGAux.As(1e-9), qGAux.Aswithx(1e-9), -cst.FactFWHMGauss*sigs/2., - x*np.sin(phi/2.), 1e-9, t1, t2
            result1 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x),  qGAux.ctzero, qGAux.ctmax,  qGAux.Zswithx,  qGAux.Cswithx ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result2 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x),  qGAux.ctmin,  qGAux.ctzero, qGAux.Zswithx,  qGAux.Bswithx ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result3 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x),  qGAux.ctmin,  qGAux.ctzero, qGAux.Dswithx,  qGAux.Zswithx ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result4 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x),  qGAux.ctzero, qGAux.ctmax,  qGAux.Aswithx,  qGAux.Zswithx ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            #result1 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x),  qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            #result2 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x),  qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            #result3 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x),  qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            #result4 = dblquad(lambda s, ct:    self.integrandwithx_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x),  qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs ) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            #print result1[0], result2[0], result3[0], result4[0], '->',
            result = (result1[0]+result2[0]+result3[0]+result4[0],)
            #print result
            #sys.exit("Not ready")
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          # result = dblquad(self.integrandwithx_CC_CK_Gaussian, -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, x))
          #sys.exit("Not ready")
          result = (np.NaN,)
        elif longdens == "RF800":
          # result = dblquad(self.integrandwithx_CC_CK_RF800,    -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, x))
          #sys.exit("Not ready")
          result = (np.NaN,)
        elif longdens == "qGaussian":
          # result = dblquad(self.integrand_CC_CK_qGaussian,   -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p)) # wrong
          #sys.exit("Not ready")
          result = (np.NaN,)
        else:
          sys.exit("\n[!] Longitudinal density unknown.\n")
    return result[0] # [1]
            
  def IntAll_withy(self, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep, y):
    if ip == 2:
      if oncc == 0.0:
        if longdens == "Gaussian":
          sys.exit("Not ready")
        elif longdens == "qGaussian":
          qGAux = qGaussianAux(sigs, phi, t1, t2)
          result1 = dblquad(lambda s, ct:    self.integrandwithy_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, y),                     qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs )
          result2 = dblquad(lambda s, ct:    self.integrandwithy_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, y),                     qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs )
          result3 = dblquad(lambda s, ct:    self.integrandwithy_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, y),                     qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs )
          result4 = dblquad(lambda s, ct:    self.integrandwithy_noCC_parsep_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, sigp, parsep, t1, t2, y),                     qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs )
          result = (result1[0]+result2[0]+result3[0]+result4[0],)
        else:
          sys.exit("\n[!] SETTINGS HAVE TO BE GAUSSIAN/QGAUSSIAN FOR IP8.\n")
      else:
        sys.exit("\n[!] SETTINGS HAVE TO BE ONCC = 0 FOR IP8.\n")
    else:
      if ck == False:
        if oncc == 0.0:
          if longdens == "Gaussian":
            result = dblquad(self.integrandwithy_CC_Gaussian, -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, 0.0, t1c, t2c, sigp, t1p, t2p, y)) # READY
          elif longdens == "RF800":
            # result = dblquad(self.integrandwithy_noCC_RF800,    -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, t1c, t2c, t1p, t2p, sigp, y))
            sys.exit("Not ready")
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            result1 = dblquad(lambda s, ct:    self.integrandwithy_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, 0.0, t1c, t2c, sigp, t1p, t2p, y),  qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result2 = dblquad(lambda s, ct:    self.integrandwithy_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, 0.0, t1c, t2c, sigp, t1p, t2p, y),  qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result3 = dblquad(lambda s, ct:    self.integrandwithy_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, 0.0, t1c, t2c, sigp, t1p, t2p, y),  qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result4 = dblquad(lambda s, ct:    self.integrandwithy_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, 0.0, t1c, t2c, sigp, t1p, t2p, y),  qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result = (result1[0]+result2[0]+result3[0]+result4[0],)
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
        else:
          if longdens == "Gaussian":
            result = dblquad(self.integrandwithy_CC_Gaussian, -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, sigp, t1p, t2p, y)) # READY
          elif longdens == "RF800":
            # result = dblquad(self.integrandwithy_CC_RF800,    -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, sigp, t1p, t2p, y))
            sys.exit("Not ready")
          elif longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            result1 = dblquad(lambda s, ct:    self.integrandwithy_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, sigp, t1p, t2p, y),  qGAux.ctzero, qGAux.ctmax,  qGAux.Zs,  qGAux.Cs) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result2 = dblquad(lambda s, ct:    self.integrandwithy_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, sigp, t1p, t2p, y),  qGAux.ctmin,  qGAux.ctzero, qGAux.Zs,  qGAux.Bs) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result3 = dblquad(lambda s, ct:    self.integrandwithy_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, sigp, t1p, t2p, y),  qGAux.ctmin,  qGAux.ctzero, qGAux.Ds,  qGAux.Zs) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result4 = dblquad(lambda s, ct:    self.integrandwithy_CC_qGaussian(  ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, sigp, t1p, t2p, y),  qGAux.ctzero, qGAux.ctmax,  qGAux.As,  qGAux.Zs) # EDIT 14/Jun/18 - t1c and t2c have to be zero!!!
            result = (result1[0]+result2[0]+result3[0]+result4[0],)
            #sys.exit("Not ready")
          else:
            sys.exit("\n[!] Longitudinal density unknown.\n")
      else:
        if longdens == "Gaussian":
          # result = dblquad(self.integrandwithy_CC_CK_Gaussian, -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, y))
          #sys.exit("Not ready")
          result = (np.NaN,)
        elif longdens == "RF800":
          # result = dblquad(self.integrandwithy_CC_CK_RF800,    -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, y))
          #sys.exit("Not ready")
          result = (np.NaN,)
        elif longdens == "qGaussian":
          # result = dblquad(self.integrandwithy_CC_CK_qGaussian,   -np.inf, np.inf, lambda u: -np.inf, lambda u: np.inf, args = (phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, y))
          #sys.exit("Not ready")
          result = (np.NaN,)
        else:
          sys.exit("\n[!] Longitudinal density unknown.\n")
    return result[0] # [1]
  
  ### LINE AND SURFACE DISTRIBUTIONS ###
  
  def GetDistribution(self, resultsname, ip, step):
    """ Computes and writes to external files the PU s-, t-, and 2D (s,t)-distributions in given ranges. """
    #densities = Densities()
    auxfunc = AuxFunc()
    
    ip = int(ip)
    
    print "\n"
    # Read desired row (steo) from the input table
    table = twiss(resultsname)
    if step == "FIRST" or step == "LAST":
      idx  = table.TYPESTEP.index(step)
    elif step == "END":
      try:    idx = table.TYPESTEP.index(step)
      except: idx = table.nts
    else:
      idx = int(step)
    print idx
    
    t1 = table.t1
    t2 = table.t2
    wcc      = table.wcc
    gamma    = table.gamma
    betarel  = table.betarel
    longdens = table.longdens # auxfunc.StringBinaryToLogical(table.longdens) # depreacted since there are three options of longitudinal distributions
    sigs     = table.SIGS[idx]
    if ip == 0:
      t1c = table.tcc1x0
      t2c = table.tcc2x0
      t1p = table.tcc1y0
      t2p = table.tcc2y0
      xplane   = table.xplane0
      ck       = auxfunc.StringBinaryToLogical(table.ck0)
      if xplane == 0:
        betc   = table.BETX0[idx]
        betp   = table.BETY0[idx]
        epsnc  = table.EPSXN[idx]
        epsnp  = table.EPSYN[idx]
      elif xplane == 1:
        betc   = table.BETY0[idx]
        betp   = table.BETX0[idx]
        epsnc  = table.EPSYN[idx]
        epsnp  = table.EPSXN[idx]
      oncc     = table.ONCC0[idx]
      phi      = table.PHI0[idx]
      phiCR    = table.PHICR0[idx]
      phiCK    = table.PHICK0[idx]
      pu       = table.PU0[idx]
      parsep   = None
    elif ip == 2:
      t1c = table.tcc1x2
      t2c = table.tcc2x2
      t1p = table.tcc1y2
      t2p = table.tcc2y2
      xplane   = table.xplane2
      ck       = auxfunc.StringBinaryToLogical(table.ck2)
      if xplane == 0:
        betc   = table.BETX2[idx]
        betp   = table.BETY2[idx]
        epsnc  = table.EPSXN[idx]
        epsnp  = table.EPSYN[idx]
      elif xplane == 1:
        betc   = table.BETY2[idx]
        betp   = table.BETX2[idx]
        epsnc  = table.EPSYN[idx]
        epsnp  = table.EPSXN[idx]
      oncc     = table.ONCC2[idx]
      phi      = table.PHI2[idx]
      phiCR    = table.PHICR2[idx]
      phiCK    = table.PHICK2[idx]
      pu       = table.PU2[idx]
      parsep   = table.PARSEP2[idx]
    else:
      sys.exit('\nINVALID IP\n')
    sigc     = np.sqrt( epsnc*betc/(gamma*betarel) )
    sigp     = np.sqrt( epsnp*betp/(gamma*betarel) )
    
    nstep0 = 200
    nstep1 = 200
    
    print 'ip =', ip
    #print 'longdens =', longdens
    if ip == 2:
      xrang0 = 100e-6
      yrang0 = 100e-6
      srang0 = 250e-3 # pu: 0.25 mm-1
      trang0 = 1.0e-9 # pu: 0.75 ns-1   
    else:
      xrang0 = 100e-6
      yrang0 = 100e-6
      srang0 = 250e-3
      trang0 = 1.0e-9
    #print '\n'
    #print 't1 =', t1
    #print 't2 =', t2
    #print 't1c =', t1c
    #print 't2c =', t2c
    #print 't1p =', t1p
    #print 't2p =', t2p
    #print '\n'
    
    # Norm
    norm = self.IntAll(phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) / pu  #  [1/events]
    #print norm
    
    #if ip == 2: proj = ['x'] #, 'y'] #'y', 's', 't'] #, 't', 'st'] # ['x', 'y', 's', 't', 'st', 'xt', 'yt']:
    #else:       proj = ['s', 't', 'st']  # 'x', 'y', 
    
    proj = ['x', 'y', 's', 't', 'st']
    
    if 'x' in proj:
      # Pile-up x-density
      print "\nPile-up x-density"
      FILEnamex = resultsname + '.' + str(ip) + "." + 'x'
      FILEnamex = FILEnamex + "." + str(step)
      FILEdensx = open(FILEnamex, "w")
      print >> FILEdensx, "* {0:12} {1:14}".format("X",   "PUX") # x in m, pux in events/m (usually rewritten in cm vs. events/mm)
      print >> FILEdensx, "$ {0:12} {1:14}".format("%le", "%le")
      xrang  = xrang0
      xnstep = nstep0
      for x in np.arange(-xrang, xrang*(1.+2./xnstep), 2*xrang/xnstep):
        pux = self.IntAll_withx(phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep, x) / norm
        print >> FILEdensx, "{0:14} {1:14}".format("%1.6e" %x, "%1.6e" %pux)
        #print x, pux
      print "\n> Writting '" + FILEnamex + "'...", "\n"
      FILEdensx.close()
      
    if 'y' in proj:
      # Pile-up y-density
      print "\nPile-up y-density"
      FILEnamey = resultsname + '.' + str(ip) + "." + 'y'
      FILEnamey = FILEnamey + "." + str(step)
      FILEdensy = open(FILEnamey, "w")
      print >> FILEdensy, "* {0:12} {1:14}".format("Y",   "PUY") # y in m, puy in events/m (usually rewritten in cm vs. events/mm)
      print >> FILEdensy, "$ {0:12} {1:14}".format("%le", "%le")
      yrang  = yrang0
      ynstep = nstep0
      for y in np.arange(-yrang, yrang*(1.+2./ynstep), 2*yrang/ynstep):
        puy = self.IntAll_withy(phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep, y) / norm
        print >> FILEdensy, "{0:14} {1:14}".format("%1.6e" %y, "%1.6e" %puy)
        #print y, puy
      print "\n> Writting '" + FILEnamey + "'...", "\n"
      FILEdensy.close()
    
    if 's' in proj:
      # Pile-up s-density
      print "\nPile-up s-density"
      FILEnames = resultsname + '.' + str(ip) + "." + 's'
      FILEnames = FILEnames + "." + str(step)
      FILEdenss = open(FILEnames, "w")
      print >> FILEdenss, "* {0:12} {1:14}".format("S",   "PUS") # s in m, pus in events/m (usually rewritten in m vs. events/mm)
      print >> FILEdenss, "$ {0:12} {1:14}".format("%le", "%le")
      srang  = srang0
      snstep = nstep0
      for s in np.arange(-srang, srang*(1.+2./snstep), 2*srang/snstep):
        pus = self.IntCT(s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) / norm
        print >> FILEdenss, "{0:14} {1:14}".format("%1.8f" %s, "%1.6e" %pus)
      print "\n> Writting '" + FILEnames + "'...", "\n"
      FILEdenss.close()
      
    if 't' in proj:
      # Pile-up t-density
      print "\nPile-up t-density"
      FILEnamet = resultsname + '.' + str(ip) + "." + 't'
      FILEnamet = FILEnamet + "." + str(step)
      FILEdenst = open(FILEnamet, "w")
      print >> FILEdenst, "* {0:12} {1:14}".format("T",   "PUT") # t in s, put in events/s (usually rewritten in ns vs. events/ps)
      print >> FILEdenst, "$ {0:12} {1:14}".format("%le", "%le")
      trang  = trang0
      tnstep = nstep0
      for t in np.arange(-trang, trang*(1.+2./tnstep), 2*trang/tnstep):
        ct = cst.clight*t
        put = self.IntS(ct, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) / norm * cst.clight
        print >> FILEdenst, "{0:14} {1:14}".format("%1.6e" %t, "%1.6e" %put)
      print "\n> Writting '" + FILEnamet + "'...", "\n"
      FILEdenst.close()
      
    if 'st' in proj:
      # Pile-up 2D-density (s, t)
      print "\nPile-up st-density"
      FILEname = resultsname + '.' + str(ip) + "." + 'st'
      FILEname = FILEname + "." + str(step)
      FILEdens = open(FILEname, "w")
      print >> FILEdens, "* {0:12} {1:14} {2:14}".format("S",   "T",   "PU") # s in m, t in s, pu in events/m/s
      print >> FILEdens, "$ {0:12} {1:14} {2:14}".format("%le", "%le", "%le")
      srang  = srang0
      snstep = nstep1
      trang  = trang0
      tnstep = nstep1
      for s in np.arange(-srang, srang*(1.+2./snstep), 2*srang/snstep):
        for t in np.arange(-trang, trang*(1.+2./tnstep), 2*trang/tnstep):
          ct = cst.clight*t
          if longdens == "qGaussian":
            qGAux = qGaussianAux(sigs, phi, t1, t2)
            #print "(s >= qGAux.szero and (ct <= qGAux.Cct(s) and ct >= qGAux.Bct(s)) =", s, ">=", qGAux.szero, "and", "(", ct, "<=", qGAux.Cct(s), "and", ct, ">=", qGAux.Bct(s), ") =", s >= qGAux.szero, "and", "(", ct <= qGAux.Cct(s), "and", ct >= qGAux.Bct(s), ")", "->", s >= qGAux.szero and (ct <= qGAux.Cct(s) and ct >= qGAux.Bct(s))
            #print "(s <  qGAux.szero and (ct <= qGAux.Act(s) and ct >= qGAux.Dct(s)) =", s, "< ", qGAux.szero, "and", "(", ct, "<=", qGAux.Act(s), "and", ct, ">=", qGAux.Dct(s), ") =", s <  qGAux.szero, "and", "(", ct <= qGAux.Act(s), "and", ct >= qGAux.Dct(s), ")", "->", s <  qGAux.szero and (ct <= qGAux.Act(s) and ct >= qGAux.Dct(s))
            #if longdens != "qGaussian" or (s >= qGAux.szero and (ct <= qGAux.Cct(s) and ct >= qGAux.Bct(s))) or (s <  qGAux.szero and (ct <= qGAux.Act(s) and ct >= qGAux.Dct(s))):
            if (s >= qGAux.szero and (ct <= qGAux.Cct(s) and ct >= qGAux.Bct(s))) or (s <  qGAux.szero and (ct <= qGAux.Act(s) and ct >= qGAux.Dct(s))):
              pu = self.Evaluate_integrand(ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) / norm * cst.clight
              #print s, ct, pu
            else:
              pu = np.NaN
              #print s, ct, "X"
          else:
            pu = self.Evaluate_integrand(ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, phiCK, sigp, t1p, t2p, ck, oncc, longdens, ip, parsep) / norm * cst.clight
          print >> FILEdens, "{0:14} {1:14} {2:14}".format("%1.8f" %s, "%1.6e" %t, "%1.6e" %pu)
        print >> FILEdens, ""
      print "\n> Writting '" + FILEname + "'...", "\n"
      FILEdens.close()
      
    if 'xt' in proj:
      # Pile-up 2D-density (x, t)
      print "\nPile-up xt-density"
      FILEname = resultsname + '.' + str(ip) + "." + 'xt'
      FILEname = FILEname + "." + str(step)
      FILEdens = open(FILEname, "w")
      print >> FILEdens, "* {0:12} {1:14} {2:14}".format("X",   "T",   "PU") # x in m, t in s, pu in events/m/s
      print >> FILEdens, "$ {0:12} {1:14} {2:14}".format("%le", "%le", "%le")
      xrang  = xrang0
      xnstep = nstep1
      trang  = trang0
      tnstep = nstep1
      #for x in np.arange(-xrang, xrang*(1.+2./xnstep), 2*xrang/xnstep):
      #  for t in np.arange(-trang, trang*(1.+2./tnstep), 2*trang/tnstep):
      #    #print x, t
      #    ct = cst.clight*t
      #    pu = quad(lambda s: self.integrandwithx_CC_Gaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, x),  -np.inf, np.inf)[0] / norm * cst.clight, 
      #    print >> FILEdens, "{0:14} {1:14} {2:14}".format("%1.8f" %x, "%1.6e" %t, "%1.6e" %pu)
      #  print >> FILEdens, ""
      #print "\n> Writting '" + FILEname + "'...", "\n"
      FILEdens.close()
      
    if 'yt' in proj:
      # Pile-up 2D-density (y, t)
      print "\nPile-up yt-density"
      FILEname = resultsname + '.' + str(ip) + "." + 'yt'
      FILEname = FILEname + "." + str(step)
      FILEdens = open(FILEname, "w")
      print >> FILEdens, "* {0:12} {1:14} {2:14}".format("Y",   "T",   "PU") # y in m, t in s, pu in events/m/s
      print >> FILEdens, "$ {0:12} {1:14} {2:14}".format("%le", "%le", "%le")
      yrang  = yrang0
      ynstep = nstep1
      trang  = trang0
      tnstep = nstep1
      #for y in np.arange(-yrang, yrang*(1.+2./ynstep), 2*yrang/ynstep):
      #  for t in np.arange(-trang, trang*(1.+2./tnstep), 2*trang/tnstep):
      #    #print y, t
      #    ct = cst.clight*t
      #    pu = quad(lambda s: self.integrandwithy_CC_Gaussian(ct, s, phi, sigs, betc, sigc, betp, t1, t2, wcc, phiCR, t1c, t2c, sigp, t1p, t2p, y),  -np.inf, np.inf)[0] / norm * cst.clight, 
      #    print >> FILEdens, "{0:14} {1:14} {2:14}".format("%1.8f" %y, "%1.6e" %t, "%1.6e" %pu)
      #  print >> FILEdens, ""
      #print "\n> Writting '" + FILEname + "'...", "\n"
      FILEdens.close()

