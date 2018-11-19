import numpy as np

outputs = open('../finalValues.txt', 'w')
outputs.write("Calculations from our data of the Snowball Nebula's electron\n")
outputs.write("  density and gas temperature.\n")

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
