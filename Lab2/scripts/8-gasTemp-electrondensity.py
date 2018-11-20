import numpy as np

outputs = open('../finalValues.txt', 'w')
outputs.write("Calculations from our data of the Snowball Nebula's electron\n")
outputs.write("  density and gas temperature.\n")


# Load Lab constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Load data
fname = info['specCal'] + 'fluxMeas_knownLines.txt'
data = np.loadtxt(fname,delimiter=',',skiprows=1)
centers = data[:,0]
centers = centers.tolist()
flux = data[:,1]
flux = flux.tolist()

# Extract wavelength and flux measurements
f5007 = [centers[0:3],flux[0:3]]
f4959 = [centers[3:6],flux[3:6]]
f4861 = [centers[6:9],flux[6:9]]
f4740 = [centers[9:12],flux[9:12]]
f4711 = [centers[12:15],flux[12:15]]
f4685 = [centers[15:18],flux[15:18]]
f4363 = [centers[18:21],flux[18:21]]
f4341 = [centers[21:],flux[21:]]

# Initialize values of temps we're looking at
gastemps = []
temps = np.linspace(1.e4, 1.5e4, num=50000)
rhss = []
lhss = []
n = 1.1e3     # electron density
for temp in temps:
    rhss.append(
        (7.90 * np.exp( (3.29e4) / temp ) / 
        (1 + ( (4.5e-4) * n * temp**(-0.5) ) ) )
    )

# Function to determine the index and value closest to a given value
def nearest(array, value):
    array = np.asarray(array)
    i = (np.abs(array - value)).argmin()
    return [array[i], i]

# Determine gas temp for each permutation of our measured values
for low in f4363[1]:
    for mid in f4959[1]:
        for high in f5007[1]:
            lhs = (mid + high) / low
            lhss.append(lhs)
            temp = nearest(rhss, lhs)
            gastemps.append(temps[temp[1]])
#print(gastemps)

gastemp = np.mean(gastemps)
gastemp_uncert = np.std(gastemps) / np.sqrt(len(gastemps))
outputs.write('OIII line ratios: %.6f +- %.6f\n' % 
               (np.mean(lhss), (np.std(lhss) / np.sqrt(len(lhss))) ) )
outputs.write('Gas Temp: %.6f +- %.6f\n' % (gastemp, gastemp_uncert))


# Get the Argon line fraction
Ar_fracs = []
for small in f4711[1]:
    for big in f4740[1]:
        Ar_fracs.append(small / big)
outputs.write('Argon lines ratio: %.6f +- %.6f\n' % 
        (np.mean(Ar_fracs) , (np.std(Ar_fracs) / np.sqrt(len(Ar_fracs))) ) )

# From that info we used image ../ArIV_doublet_density to get electron density
outputs.write('Electron density via Argon: 1100 +- 100\n')
