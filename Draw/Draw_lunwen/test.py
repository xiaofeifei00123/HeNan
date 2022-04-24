# %%
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




