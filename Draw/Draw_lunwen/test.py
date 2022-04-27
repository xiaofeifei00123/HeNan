# %%
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

x = np.array(np.arange(10000))
y = np.array(np.arange(2000))
start = time.process_time()
cc = myloop(x,y)
# print(cc)
end = time.process_time()
print('different is %6.3f' % (end - start))
# cc