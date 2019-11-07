#import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

n = len(np.arange(5, 11))
phi = len(np.arange(0.3, 0.9, 0.1))
psi = len(np.arange(0.8, 2, 0.1))
Lambda = len(np.arange(0.45,0.55,0.05))
AR = len(np.arange(1.0, 2.0, 0.2))
dho = len(np.arange(1, 1.2, 0.02))
a1 = len(np.arange(0, 20, 2))

num = AR*phi**2*psi**2*dho*a1*n

print(num)
