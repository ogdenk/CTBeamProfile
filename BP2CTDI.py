"""
Program to calculate the CTDI in air as a function of radial distance from isocenter using a known beam profile.
"""

import numpy
import math
from numpy import sin, cos, tan, pi, arctan2, sqrt, arccos

BP = [1.0, 0.84211, 0.55526, 0.24868, 0.15, 0.10921, 0.08289 ] # Measured beam profile values (normalized)
d = [0., 5., 10., 15., 20., 25., 30.]  # Measurement positions

theta = arctan2(d, 57.0) # Angular position of the beam profile measured values relative to central ray

r = numpy.sqrt(numpy.square(d)+57**2.0)  # radial distance from focal spot to measurement point

print (r)

