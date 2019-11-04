"""Basic mathematical functions not included in numpy or math"""

import numpy as np

def csc(x):
    """Return the cosecant of x"""
    return 1/np.sin(x)

def sec(x):
    """Return the secant of x"""
    return 1/np.cos(x)

def cot(x):
    """Return the cotangent of x"""
    return 1/np.tan(x)

def arccsc(x):
    """Return the inverse cosecant of x"""
    return np.arcsin(1/x)

def arcsec(x):
    """Return the inverse secant of x"""
    return np.arccos(1/x)

def arccot(x):
    """Return the inverse cotangent of x"""
    return np.arctan(1/x)
