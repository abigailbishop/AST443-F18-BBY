"""

BBY
AST 443
Lab0 4.5 
4.5 Spectral Lineup Plot

"""

import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal

# load data
data = np.loadtxt('spectral_lines_lineup.txt',skiprows=1)
pixel = data[:,0]
wavelength = data[:,1]

# perform linear regression
linearfit = np.polynomial.polynomial.polyfit(pixel,wavelength, deg=1)
x=np.linspace(0,550,10**3)
def line(x,b,m):
    return b+m*x
lin_reg = line(x,*linearfit)
m = round(Decimal(linearfit[1]),2)
b = round(Decimal(linearfit[0]),2)
fit_label = 'linear fit, y = '+ str(m)+ 'x + '+ str(b)
# m = -0.98 A/pix

# plot data and fit
plt.figure(1)
plt.title(r'Wavelength Calibration')
plt.errorbar(pixel,wavelength,fmt='x',label='data')
plt.plot(x,lin_reg,label=fit_label)
plt.ylabel(r'Wavelength, angstroms')
plt.xlabel(r'Pixels')
plt.legend(loc='best')
plt.savefig('4.3.6_wavelengthCalibration.pdf',format='pdf',dpi=1000,bbox_inches='tight')
