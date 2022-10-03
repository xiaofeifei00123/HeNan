# %%
<<<<<<< HEAD
import numpy as np
from numba import jit
import time

# %%
# %time
# @np.vectorize
@jit(nopython=True)
def myloop(x,y):
    c = 0
    # c = x*y
    for i in x:
        for j in y:
        # print(i)
            c = i*2+j+c
    return c
    # for i in x:
        # for j in y:
            # c = c+i*j 
    # return x*y

x = np.array(np.arange(1, 10000, 0.1))
y = np.array(np.arange(1000))
start = time.process_time()
cc = myloop(x,y)
# print(cc)
end = time.process_time()
print('different is %6.3f' % (end - start))
# cc
=======
import xarray as xr
import numpy as np
from waveletAnalysis import get_data
import matplotlib.pyplot as plt

# %%
da = get_data()
# da.distance
xn = da.values
xn
# %%
xk = np.fft.fft(xn)
N = len(xk)
# N = 5
# N = da.distance.values
# n = np.linspace(0, N-1, N)
n = da.distance.values

f = n/N

sx = xk**2/N

plt.plot(f, np.abs(sx))




>>>>>>> gravity
