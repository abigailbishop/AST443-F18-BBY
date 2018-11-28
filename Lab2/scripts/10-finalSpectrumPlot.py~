"""

BBY
AST 443
Lab0 4.5 
4.5 Calibrated Arclamp Spectrum Plot

"""

import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal

# load data
data = np.loadtxt('calibrated_arcspec.txt')
x = np.flip(data[:,0],axis=0)
y = data[:,1]

ymax_idx = y.argmax(axis=0)
shift = x[ymax_idx]

x = x - shift + 5852.49

# plot data and fit
plt.figure(1)
plt.title(r'Arclamp Spectrum')
plt.errorbar(x,y,label='data')
plt.xlabel(r'Wavelength, angstroms')
plt.ylabel(r'Relative Brightness')
plt.annotate(r'$5852.49{\AA}$', xy = (5852.49,0.12),xytext=(-60,-20),
             textcoords='offset points', 
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
plt.annotate(r'$6143.06{\AA}$', xy = (6143.06,0.055),xytext=(-60,20),
             textcoords='offset points', 
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
plt.annotate(r'$6266.49{\AA}$', xy = (6266.49,0.043),xytext=(0,20),
             textcoords='offset points', 
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
plt.savefig('4.5.8_calibratedSpectrum.pdf',format='pdf',dpi=1000,bbox_inches='tight')
