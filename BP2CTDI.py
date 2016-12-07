"""
Program to calculate the CTDI in air as a function of radial distance from isocenter using a known beam profile.
"""

import numpy as np
import math
from numpy import sin, cos, tan, pi, arctan2, sqrt, arccos
import scipy
from scipy import interpolate
import matplotlib.pyplot as plt

BP = np.array([1.0, 0.84211, 0.55526, 0.24868, 0.15, 0.10921, 0.08289 ]) # Measured beam profile values (normalized)
d = np.array([0., 5., 10., 15., 20., 25., 30.])  # Measurement positions

theta = arctan2(d, 57.0) # Angular position of the beam profile measured values relative to central ray

r = np.sqrt(np.square(d)+57**2.0)  # radial distance from focal spot to measurement point

P = BP*np.square(r)/57.0**2.0 # Beam Profile corrected to constant 57 cm distance from FS

profile = interpolate.interp1d(theta, P,  kind = "cubic") # cubic spline fit of BP vs theta

ctdi = np.zeros(d.size)

for x in  np.arange(0,d.size):
    for a in np.arange(0,360): # 360 angular tube positions to calculate the CTDI
        Tx = -57 * sin(a * np.pi / 180)  # Focal Spot x coordinates (0 degrees is overhead so use sin instead of cos)
        Ty = 57 * cos(a * np.pi / 180)  # Focal Spot y coordinates
        dp = np.dot([-Tx,-Ty],[d[x]-Tx,-Ty])
        Phi = arccos(dp/(np.sqrt((d[x]-Tx)**2+Ty**2)*np.sqrt(Tx**2+Ty**2)))
        ctdi[x] = ctdi[x] + profile(Phi)*(57**2)/np.sqrt((d[x]-Tx)**2+Ty**2)**2

plt.plot(d,ctdi, 'o', d, BP, '-')
plt.show()

print (d)
print (BP)
print (ctdi)

