import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pdb

from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u

field_size = 10.0 #"

df = pd.read_csv('./My_First_Catalogue.csv')







pdb.set_trace()

o = np.where(df.type == 1.0)[0]
dfSF = df.iloc[o]

Rsb = dfSF.SFR/dfSF.SFRms

bins = np.linspace(-5.,5.,500)
plt.hist(np.log10(Rsb), bins = bins, )

o = np.where(df.type == 2.0)[0]
dfSB = df.iloc[o]
Rsb = dfSB.SFR/dfSB.SFRms

plt.hist(np.log10(Rsb), bins = bins)

o = np.where(df.type == 0.0)[0]
dfqui = df.iloc[o]
Rsb = dfqui.SFR/dfqui.SFRms

plt.hist(np.log10(Rsb), bins = bins)
plt.yscale('log')
plt.axis([-1.5,2.0,1e1,1e4])
plt.show()



# o = np.where((df.z > 0.6) & (df.z<0.7))[0]
# dfz = df.iloc[o]

# bins = np.linspace(8.,12.5,100)

# o = np.where(dfz.type == 1)[0]
# dfzSF = dfz.iloc[o]

# Dc_min = cosmo.luminosity_distance(min(dfzSF.z))/(1.+min(dfzSF.z))
# Dc_max = cosmo.luminosity_distance(max(dfzSF.z))/(1.+max(dfzSF.z))
# Vslice = field_size*(np.pi/180.)**2 / 3.*(Dc_max**3-Dc_min**3)
  
# pdfSF, bins = np.histogram(dfzSF.logM, bins = bins)
# dM = np.diff(bins)
# mfSF = pdfSF/dM/Vslice 

# o = np.where(dfz.type == 0)[0]
# dfzQUI = dfz.iloc[o]

# Dc_min = cosmo.luminosity_distance(min(dfzQUI.z))/(1.+min(dfzQUI.z))
# Dc_max = cosmo.luminosity_distance(max(dfzQUI.z))/(1.+max(dfzQUI.z))
# Vslice = field_size*(np.pi/180.)**2 / 3.*(Dc_max**3-Dc_min**3)

# pdfQUI, bins = np.histogram(dfzQUI.logM, bins = bins)
# dM = np.diff(bins)
# mfQUI = pdfQUI/dM/Vslice 

# logMvec = np.linspace(8.0,13.0,100)
# Mvec = 10**logMvec

# p = np.array([10.40+0.24, -1.32, 1.661*1e-3, 0.84, 0.86*1e-3])
# Mrat = Mvec/10**p[0]
# PhiSF = np.exp(-Mrat) * (p[2]*(Mrat)**p[1] + p[4]*(Mrat)**p[3])/10**p[0]*Mvec*np.log(10)

# p = np.array([10.83+0.24, -1.46, 0.012*1e-3, -0.21, 1.44*1e-3])
# Mrat = Mvec/10**p[0]
# PhiQUI = np.exp(-Mrat) * (p[2]*(Mrat)**p[1] + p[4]*(Mrat)**p[3])/10**p[0]*Mvec*np.log(10)

# p = np.array([10.83+0.24, -1.46, 0.015*1e-3, -0.21, 1.44*1e-3])
# Mrat = Mvec/10**p[0]
# PhiQUI = np.exp(-Mrat) * (p[2]*(Mrat)**p[1] + p[4]*(Mrat)**p[3])/10**p[0]*Mvec*np.log(10)


# plt.plot(logMvec, PhiSF, 'r-')
# plt.plot(logMvec, PhiQUI, 'r--')
# plt.plot(bins[0:-1]+dM/2., mfSF, 'b.')
# plt.plot(bins[0:-1]+dM/2., mfQUI, 'bx')
# plt.yscale('log')
# plt.axis([8.0,12.5,1e-6,1e-1])
# plt.show()

pdb.set_trace()
