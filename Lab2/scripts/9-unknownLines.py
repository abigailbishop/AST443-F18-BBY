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
fname = '../flux_measurements_unknownLines.txt'
data = np.loadtxt(fname,delimiter=',',skiprows=1)
centers = data[:,0]
centers = centers.tolist()
flux = data[:,1]
flux = flux.tolist()

N = 22
centersMean = np.zeros(N)
sig_centers = np.zeros(N)
fluxMean = np.zeros(N)
sig_flux = np.zeros(N)

# Extract wavelength and flux measurements, uncertainty on mean (std/sqrt(3))
i=0
while i < N:
    j=i*3
    centersMean[i] = np.mean(centers[j:j+3])
    sig_centers[i] = np.std(centers[j:j+3])/np.sqrt(3)
    fluxMean[i] = np.mean(flux[j:j+3])
    sig_flux[i] = np.std(flux[j:j+3])/np.sqrt(3)
    i=i+1

# save to txt file
fname_out = '../unknownLinesAvg.txt'
np.savetxt(fname_out,np.c_[centersMean,sig_centers,fluxMean,sig_flux],
           delimiter=',',header='wavelength,sig_wl,flux,sig_flux')
