import sncosmo
import matplotlib.pyplot as plt
import numpy as np

'''
Building a very simple model. In practise, this could be used to build a model of a supernova.
In this case, it has been used to create a linear model.
phase - a time array
disp - Range of wavelengths
flux - The elements within the flux array represent how the spectrum is 'shaped' as time progresses.
     - flux = [t1, t2] where t1 = [flux at wavelength 1,flux at wavelength 2, flux at wavelength 3] etc

Note: len(flux) = len(phase) and len(flux[0])=len(disp)
'''

time_points = 7
wavelength_points = 4

phase = np.linspace(2455454.88645 - 50., 2455454.88645 + 50., time_points)
disp = np.linspace(1000., 9000., wavelength_points)

# Below, the flux array is made and manipulated across different times and wavelengths. These can be changed, or
# inputted manually, changing the final plots.
x = np.arange(wavelength_points)
flux = []
for i in range(1, time_points + 1):
    flux.append(x * i)
flux = np.array(flux)

for i in range(wavelength_points):
    for z in flux:
        z[i] += i

# Set the created arrays as model
model = sncosmo.TimeSeriesSource(phase, disp, flux)

# Plot the graph of observed luminosity through different filters over time
timearray = np.arange(2455454.88645 - 120., 2455454.88645 + 120, 0.01)

# Registering ptf48r filter
import os
import sys
sys.path.append(os.path.abspath('..'))
import Bandpass

# Plotting relative flux against time
filters = ['ptf48r', 'bessellr', 'bessellv', 'bessellb' ]
colors = ['black','red','green','blue']
for i in range(len(filters)):
    plt.scatter(timearray, model.bandflux(filters[i], timearray), color = colors[i%4], label= str(filters[i]))
plt.xlabel('Time/ days')
plt.ylabel('Relative Flux')
plt.legend()
plt.show()

# A 3D plot with axes representing time, wavelength and observed luminosity
coords = []
for i in range(len(flux)):
    for j in range(len(flux[0])):
        coords.append([phase[i], disp[j], flux[i][j]])
x = []
y = []
z = []
for i in coords:
    x.append(i[0])
    y.append(i[1])
    z.append(i[2])

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z)
ax.set_xlabel('Time /days')
ax.set_ylabel('Wavelength /A')
ax.set_zlabel('Relative Intensity')
plt.show()
