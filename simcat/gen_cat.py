import numpy as np
import init
import pandas as pd

from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from gen_func import *

# DVLP modules
import pdb
import matplotlib.pyplot as plt

def compute_dVdz(z, lumdist):
	# from Hogg 1999
	# dV [Mpc/sr]
	omegaM = cosmo.Om0
	omegaL = 1. - omegaM
	H0 = cosmo.H0

	E = np.sqrt(omegaM * (1.+z)**3. + omegaL)
	DH = 3000./(H0/100.)
	Da = lumdist/(1.+z)**2.
	dVdz = DH * (1.+z)**2. * Da**2. / E

	dVdz_dim = np.array(dVdz/u.s*u.km/u.Mpc**3.)

	return dVdz_dim

# Store variables from the init file
catName = init.cat_name
fieldSize = init.fieldSize
logMcut = init.logMcut
zMax = init.zmax

# Field size en sr
Omega = (2.*np.pi/360.)**2*fieldSize

# Create the redshift grid
dlog1plusz = 0.01
log1pluszMax = np.log10(zMax + 1)
Nlog1plusz = log1pluszMax/dlog1plusz
log1plusz = np.linspace(0.01,log1pluszMax+0.01,Nlog1plusz)
z = 10**log1plusz - 1.
dz = np.diff(z)
z = z[0:-1]

# Create the luminosity distance grid
lumdist = cosmo.luminosity_distance(z)

# Compute the dV
dVdz = compute_dVdz(z, lumdist)
dV = dVdz * dz

logM = []
Ngal = 0
zgal = []
SF = []
for z_i, dz_i, dV_i in zip(z, dz, dV):
	
	# Generate the stellar mass SF
	logM_i, N_i = gen_massSF(logMcut, z_i, dV_i, Omega)
	zgal_i = np.random.uniform(z_i, z_i+dz_i, N_i)
	SF_i = np.ones(N_i)

	logM.append(logM_i)
	zgal.append(zgal_i)
	SF.append(SF_i)
	Ngal += N_i


	# Generate the stellar mass QUI
	logM_i, N_i = gen_massQUI(logMcut, z_i, dV_i, Omega)
	zgal_i = np.random.uniform(z_i, z_i+dz_i, N_i)
	SF_i = np.zeros(N_i)

	logM.append(logM_i)
	zgal.append(zgal_i)
	SF.append(SF_i)
	Ngal += N_i

# Generate SFR 


logM = np.hstack(logM)
zgal = np.hstack(zgal)
SF = np.hstack(SF)
ID = np.linspace(1,Ngal,Ngal, dtype = int)

result_dict = {'ID' : ID, 'z' : zgal, 'logM' : logM, 'type' : SF}
df = pd.DataFrame(result_dict, columns = ['ID', 'z', 'logM', 'type'])
df.to_csv('../catalogues/'+catName+'.csv', index = False)