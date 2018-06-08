import numpy as np
import param

from scipy import interpolate

# DVLP modules
import pdb
import matplotlib.pyplot as plt

#-------------------------------------------------------
# Generate the stellar masses for SF gal (Davidzon+2017)
#-------------------------------------------------------
def gen_massSF(logMmin, z, dV, Omega):

	logMmax = param.logMmax
	dlogM = param.dlogM
	NlogM = int((logMmax - logMmin) / dlogM + 1)

	logM_vec = np.linspace(logMmin, logMmax, NlogM)
	Mvec = 10**logM_vec

	if z < 6:
		Mknee = np.interp(z, param.z_SF, 10.**param.logMknee_SF)
		Phiknee1 = np.interp(z, param.z_SF, param.Phiknee1_SF)
		Phiknee2 = np.interp(z, param.z_SF, param.Phiknee2_SF)
		a1 = np.interp(z, param.z_SF, param.alpha1_SF)
		a2 = np.interp(z, param.z_SF, param.alpha2_SF)
		
		Mrat = Mvec/Mknee 
		Phi = np.exp(-Mrat) * (Phiknee1*(Mrat)**a1 + Phiknee2*(Mrat)**a2)/Mknee*Mvec*np.log(10)

	else:
		# For z>4, we use a different redshift evolution of the parameters as suggested in
		# Sargent+2012
		alpha = 1.3*(1+(z-2.))**0.19               
		Mknee = 10.**(np.log10(10.**11.20)-0.3*(z-3.5))
		phi_knee = 0.00470*(1.+z)**(-2.4)

		Phi = phi_knee*(Mvec/Mknee)**(-alpha)*np.exp(-Mvec/Mknee)*(Mvec/Mknee)*np.log(10)
	
	PDF = Phi/np.sum(Phi)/dlogM
	cdf = np.cumsum(PDF[::-1]* dlogM)[::-1] 

	Ngal = int(cdf[0]*np.sum(Phi)*dlogM * dV * Omega)

	X = np.random.uniform(0.,1.,Ngal)

	logM = interpolate.griddata(cdf,logM_vec,X)

	return logM, Ngal

#-------------------------------------------------------
# Generate the stellar masses for QUI gal (Davidzon+2017)
#-------------------------------------------------------
def gen_massQUI(logMmin, z, dV, Omega):

	logMmax = param.logMmax
	dlogM = param.dlogM
	NlogM = int((logMmax - logMmin) / dlogM + 1)

	logM_vec = np.linspace(logMmin, logMmax, NlogM)
	Mvec = 10**logM_vec

	Mknee = np.interp(z, param.z_QUI, 10.**param.logMknee_QUI)
	Phiknee1 = np.interp(z, param.z_QUI, param.Phiknee1_QUI)
	Phiknee2 = np.interp(z, param.z_QUI, param.Phiknee2_QUI)
	a1 = np.interp(z, param.z_QUI, param.alpha1_QUI)
	a2 = np.interp(z, param.z_QUI, param.alpha2_QUI)
	
	Mrat = Mvec/Mknee 
	Phi = np.exp(-Mrat) * (Phiknee1*(Mrat)**a1 + Phiknee2*(Mrat)**a2)/Mknee*Mvec*np.log(10)
	
	PDF = Phi/np.sum(Phi)/dlogM
	cdf = np.cumsum(PDF[::-1]* dlogM)[::-1]

	Ngal = int(cdf[0]*np.sum(Phi)*dlogM * dV * Omega)

	X = np.random.uniform(0.,1.,Ngal)

	logM = interpolate.griddata(cdf,logM_vec,X)

	return logM, Ngal





