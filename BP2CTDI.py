"""
Program to calculate the CTDI in air as a function of radial distance from isocenter using a known beam profile.
"""

import numpy as np
from numpy import sin, cos, tan, pi, arctan2, sqrt, arccos
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

fs2iso = 57.0 # Focal spot to isocenter distance in cm
BP = np.array([1.0, 0.84211, 0.55526, 0.24868, 0.15, 0.10921, 0.08289 ]) # Measured beam profile values on x axis (normalized)
d = np.array([0., 5., 10., 15., 20., 25., 30.])  # Measurement positions

theta = arctan2(d, fs2iso) # Angular position of the beam profile values relative to central ray

r = np.sqrt(np.square(d)+fs2iso**2.0)  # radial distance from focal spot to measurement point

P = BP*np.square(r)/fs2iso**2.0 # Beam Profile corrected to constant fs2iso cm distance from FS

profile = interpolate.PchipInterpolator(theta, P, extrapolate = True) # cubic spline fit of BP vs theta

ctdi = np.zeros(d.size)

for x in  np.arange(0,d.size): # for each measured BP value, calculate ctdi at that same location
    for a in np.arange(0,360): # 360 angular tube positions to calculate the CTDI
        Tx = -fs2iso * sin(a * np.pi / 180)  # Focal Spot x coordinates (0 degrees is overhead so use sin instead of cos)
        Ty = fs2iso * cos(a * np.pi / 180)  # Focal Spot y coordinates
        dp = np.dot([-Tx,-Ty],[d[x]-Tx,-Ty])
        r= np.sqrt((d[x]-Tx)**2+Ty**2)
        cosTheta=dp/(fs2iso*r)
        if (cosTheta >1):
            cosTheta = 1
        Phi = arccos(cosTheta)
        ctdi[x] += profile(Phi)*(fs2iso**2)/np.sqrt((d[x]-Tx)**2+Ty**2)**2 # accumulate inverse square corrected contribution

maxval = ctdi.max()
ctdi = ctdi/maxval

dSmooth = np.linspace(0,d.max(), 100)
BPSmooth = interpolate.spline(d, BP, dSmooth,kind = 'smoothest')
ctdiSmooth = interpolate.spline(d, ctdi, dSmooth,kind = 'smoothest')


red_patch = mpatches.Patch(color='red', label='Beam Profile')
green_patch = mpatches.Patch(color='green', label='CTDI')
plt.legend(handles=[green_patch, red_patch])


plt.plot(d,ctdi, 'o', ms=10)
plt.plot(d, BP, '+', ms=10)
plt.plot(dSmooth, BPSmooth, 'r', label='Beam Profile')
plt.plot(dSmooth, ctdiSmooth, 'g', label='CTDI')

plt.show()

print (d)
print (BP)
print (ctdi)

