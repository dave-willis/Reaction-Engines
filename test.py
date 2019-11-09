#import itertools
import numpy as np

n = len(np.arange(5, 11))
phi = len(np.arange(0.2, 0.9, 0.3))
psi = len(np.arange(0.6, 2, 0.3))
Lambda = len(np.arange(0.2,0.9,0.3))
AR = len(np.arange(0.9, 2.0, 0.5))
dho = len(np.arange(1, 1.2, 0.1))
a1 = len(np.arange(0, 20, 2))

num = AR*phi**2*psi**2*dho**2*Lambda**2

print(num)
