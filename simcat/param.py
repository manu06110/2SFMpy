import numpy as np

# Maximum mass
logMmax = 13.
dlogM = 0.01

# stellar mass for SF gal (Davidzon et al. 2017)
z_SF = np.array([0.35,0.65,0.95,1.3,1.75,2.25,2.75,3.25,3.75,5.0])

logMknee_SF = np.array([10.26,10.40,10.35,10.42,10.40,10.45,10.39,10.83,10.77,11.30])+0.24
Phiknee1_SF = np.array([2.410,1.661,1.739,1.542,1.156,0.441,0.441,0.086,0.052,0.003])*1.e-3
Phiknee2_SF = np.array([1.30,0.86,0.95,0.49,0.46,0.38,0.13,0.0,0.0,0.0])*1.e-3
alpha1_SF = np.array([-1.29,-1.32,-1.29,-1.21,-1.24,-1.50,-1.52,-1.78,-1.84,-2.12])
alpha2_SF = np.array([1.01,0.84,0.81,1.11,0.90,0.59,1.05,0.0,0.0,0.0])

# stellar mass for QUI gal (Davidzon et al. 2017)
z_QUI = [0.35,0.65,0.95,1.3,1.75,2.25,2.75,3.25,3.75]
  
logMknee_QUI = np.array([10.83,10.83,10.75,10.56,10.54,10.69,10.24,10.10,10.10])+0.24
Phiknee1_QUI = np.array([0.098,0.012,1.724,0.757,0.251,0.068,0.028,0.010,0.004])*1.e-3
Phiknee2_QUI = np.array([1.58,1.44,0.,0.,0.,0.,0.,0.,0.])*1.e-3
alpha1_QUI = np.array([-1.3,-1.46,-0.07,0.53,0.93,0.17,1.15,1.15,1.15])
alpha2_QUI = np.array([-0.39,-0.21,0.,0.,0.,0.,0.,0.,0.])



# ;;SFR Schreiber
# m0 = 0.5
# eM0 = 0.07
# a0 = 1.5
# ea0 = 0.15
# a1 = 0.3
# ea1 = 0.08
# m1 = 0.36
# em1 = 0.3
# a2 = 2.5
# ea2 = 0.6

# sigma_MS = 0.31 ;;dex
# esigma_MS = 0.02
# sigma_SB = sigma_MS
# fSB = 0.033
# efSB = 0.015
# B_SB = 5.3
# eB_SB = 0.4
# x0 = 0.87
# ex0 = 0.04

# ;; Kennicutt
# Kuv = 1.65*1.7e-10     
# ;;1.70*1d-10 
# Kir = 1.65*1.09e-10     
# ;;1.09*1d-10

# ;; Heinis2013
# alpha_H13 = 0.72
# ealpha_H13 = 0.08
# IRX0 = 1.32
# eIRX0 = 0.04

# ;;Pannella2015
# alpha_P15 = 1.6
# Auv0_P15 = -13.5
# scat_P15 = 0.5