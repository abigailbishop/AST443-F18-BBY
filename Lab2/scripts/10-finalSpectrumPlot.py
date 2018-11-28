"""

BBY
AST 443
Lab2
Final Spectrum Plot

"""

import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal

# load data
data = np.loadtxt('../specCalFinal/nebula_avg.cal.txt')
#x = np.flip(data[:,0],axis=0)
x = data[:,0]
y = data[:,1]

# plot data
plt.figure(1)
plt.title(r'NGC 7662 Spectrum')
plt.errorbar(x,y,label='data')
plt.xlabel(r'Wavelength, angstroms')
plt.ylabel(r'Relative Brightness')

plt.annotate(r'$H_{\gamma}$', xy = (4340,4e-11),xytext=(4340,1.6e-10),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.annotate(r'[O III]', xy = (4360,1.2e-11),xytext=(4385.5,6e-11),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.annotate(r'He II', xy = (4680,6e-11),xytext=(4642,1.3e-10),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.annotate(r'[Ar IV]', xy = (4710,5e-12),xytext=(4726,1.3e-10),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.annotate(r'[Ar IV]', xy = (4736.53,8.2e-12),xytext=(4763,6e-11),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.annotate(r'$H_{\beta}$', xy = (4860,1e-10),xytext=(4885.42,1.8e-10),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.annotate(r'[O III]', xy = (4957.51,4e-10),xytext=(4913,4.7e-10),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.annotate(r'[O III]', xy = (5007,1.15e-9),xytext=(4935.57,1e-9),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))


plt.savefig('../finalSpectrumPlot.pdf',format='pdf',dpi=1000,bbox_inches='tight')

