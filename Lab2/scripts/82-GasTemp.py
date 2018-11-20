import numpy as np

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
f5007 = [[centers[0:3]],[flux[0:3]]]
f4959 = [[centers[3:6]],[flux[3:6]]]
f4861 = [[centers[6:9]],[flux[6:9]]]
f4740 = [[centers[9:12]],[flux[9:12]]]
f4711 = [[centers[12:15]],[flux[12:15]]]
f4685 = [[centers[15:18]],[flux[15:18]]]
f4363 = [[centers[18:21]],[flux[18:21]]]
f4341 = [[centers[21:]],[flux[21:]]]
'''
f5007 = [[5006.67, 5006.67, 5006.67], 
         [5.592e-9, 5.612e-9, 5.63e-9]] 
f4959 = [[4957.89, 4957.89, 4957.89], 
         [1.346e-9, 1.131e-9, 1.364e-9]] 
f4861 = [[4858.82, 4858.84, 4858.81], 
         [3.56e-10, 3.66e-10, 3.56e-10]] 
f4740 = [[4736.64, 4736.59, 4736.64], 
         [1.83e-11, 1.63e-11, 1.82e-11]] 
f4711 = [[4707.92, 4707.87, 4707.89], 
         [2.13e-11, 2.00e-11, 2.26e-11]] 
f4685 = [[4681.91, 4681.91, 4681.91], 
         [1.70e-10, 1.75e-10, 1.74e-10]] 
f4363 = [[4362.34, 4362.27, 4362.21], 
         [4.81e-11, 4.06e-11, 4.30e-11]] 
f4341 = [[4339.97, 4339.92, 4339.94], 
         [9.39e-11, 9.50e-11, 9.88e-11]] 
'''
# Initialize values of temps we're looking at
gastemps = []
temps = np.linspace(1.e4, 1.5e4, num=50000)
rhss = []
for temp in temps:
    rhss.append(
        (7.90 * np.exp( (3.29e4) / temp ) / 
        (1 + ( (4.5e-4) * (1.e3) * temp**(-0.5) ) ) )
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
            temp = nearest(rhss, lhs)
            gastemps.append(temps[temp[1]])

gastemp = np.mean(gastemps)
gastemp_uncert = np.std(gastemps) / np.sqrt(len(gastemps))
print('Gas Temp: %.6f +- %.6f' % (gastemp, gastemp_uncert))
