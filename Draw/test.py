# %%
import numpy as np
import xarray as xr
import pandas as pd
# %%
aa = np.array([[1,np.nan], 
               [2,2], 
               [3, 3],
               [np.nan, 4]])
aa
# cc = pd.DataFrame(aa)
# cc
dd = xr.DataArray(aa,
                  coords={
                      'type':['lat','lon'],
                      'id':[0,1, 2, 3 ],
                  },
                  dims=['id', 'type'])
dd
# %%
dd.mean(dim='type')
# d
