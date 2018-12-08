import numpy as np
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math
from scipy.stats import norm
from scipy import asarray as ar,exp
from scipy.optimize import curve_fit
from astropy.stats import gaussian_sigma_to_fwhm
from operator import itemgetter

# Load Constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

AU2km = 1.496e+8 
alpha_litVal = 1919.0*(1/3600.0)*(np.pi/180) # NASA, 1919 seconds of arc
#alpha_litVal = 0.0087
d_litVal = alpha_litVal * AU2km
d_prelim = 1.540135e+06
alpha_prelim = 1.029502e-02
pcov1 = 1.38319182e-10

# Extract Data
fname1 = '../52c-Sun_prelimPlotdata.txt'
sunData = np.loadtxt(fname1)
B_sun = sunData[:,0]
sig_Bsun = sunData[:,1]
V_sun = sunData[:,2]
sig_Vsun = sunData[:,3]

fname2 = '../52c-Satellite_prelimPlotdata.txt'
satData = np.loadtxt(fname2)
V_sat = satData[:,2]
sig_Vsat = satData[:,3]

# Normalize
Vsun_norm = []
sig_Vsun_norm = []
for i in range(len(V_sun)):
    Vsun_norm.append(V_sun[i] / V_sat[i])
    sig_Vsun_norm.append(Vsun_norm[i]*np.sqrt((sig_Vsun[i]/V_sun[i])**2
                       + (sig_Vsat[i]/V_sat[i])**2))

# Define a sinc function
def sincFunc(BLambda, alpha):
    return abs( np.sinc(BLambda * alpha) )

# Fit Normalized Data and Calculate Diameter
popt, pcov = curve_fit(sincFunc, B_sun, Vsun_norm, p0=[1.0e+6/AU2km],
                       sigma=sig_Vsun_norm, absolute_sigma=True)

alpha = popt[0]
#sigAlpha = pcov[0]
sigAlpha = abs(alpha - alpha_prelim)
d = AU2km * alpha # small angle approx to find diameter
#sigd = AU2km * sigAlpha # prop. uncertainty
sigd = abs(d - d_prelim) # take diff between norm and prenorm for uncer.

# Literature Agreement
d_litAgree1 = abs(d - d_litVal)/sigd
d_litAgree2 = abs(d_prelim - d_litVal)/sigd
alpha_litAgree1 = abs(alpha - alpha_litVal)/sigAlpha
alpha_litAgree2 = abs(alpha_prelim - alpha_litVal)/sigAlpha

print('alpha_lit = 1919 arcsec, NASA')
print('diameter lit = {:e} km, calc from alpha'.format(d_litVal))
print('DIAMETER FROM SAT NORM')
print('Xi^2 = ', pcov[0])
print('alpha = {:e} pm {:e} radians'.format(alpha, sigAlpha))
print('Agreement = {:e} sigma'.format(alpha_litAgree1))
print('d =  {:e} pm {:e} km'.format(d, sigd))
print('Agreement = {:e} sigma'.format(d_litAgree1))
print(' ')
print('DIAMETER SANS SAT NORM')
print('Xi^2 = ', pcov1)
print('alpha = {:e} pm {:e} radians'.format(alpha_prelim, sigAlpha))
print('Agreement = {:e} sigma'.format(alpha_litAgree2))
print('d =  {:e} pm {:e} km'.format(d_prelim, sigd))
print('Agreement = {:e} sigma'.format(d_litAgree2))

# Plot
fittedB = np.arange(min(B_sun), max(B_sun),
                    step=( (max(B_sun)-min(B_sun))/100. ) )
fittedV = []
trueV = []
for b in fittedB:
    fittedV.append(sincFunc(b, *popt))
    trueV.append(sincFunc(b, alpha_litVal))

plt.errorbar(B_sun, Vsun_norm, xerr=sig_Bsun, yerr=sig_Vsun_norm,
             fmt='.', label='Sat. Normalized Data')
plt.plot(fittedB, fittedV, label='Fitted Sinc Function')
plt.plot(fittedB, trueV,'--', label='Expected Sinc Function')
plt.xlabel(r'$B_{\lambda}$')
plt.ylabel(r'Visibility, $V_0(B_{\lambda})$')
plt.minorticks_on()
plt.legend()
plt.title('Normalized Sun Interferometer Visibility')
plt.savefig(info['images'] + 'normVis-Sun.pdf', ppi=300)
