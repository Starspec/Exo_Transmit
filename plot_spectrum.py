import numpy as np
from scipy.signal import savgol_filter as smoothing
import matplotlib.pyplot as plt
import pdb

# Constants
fname = 'default.dat'

# Read spectrum
spec = np.genfromtxt(f'Spectra/{fname}', skip_header=2)

# Convert to nm
spec[:,0] *= 1e9

# Smooth the spectrum for clarity
smooth_fluxes = smoothing(spec[:,1], 25, 3)
smooth_spec = np.zeros_like(spec)
smooth_spec[:,0] = spec[:,0]
smooth_spec[:,1] = smooth_fluxes.flatten()[:]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(spec[:,0], spec[:,1], 'r-', alpha=0.2)
ax.plot(smooth_spec[:,0], smooth_spec[:,1],'k-')
ax.set_ylabel('Pressure (bars)')
ax.set_yscale('log')
ax.set_xlabel('Wavelength (nm)')

plt.show()
