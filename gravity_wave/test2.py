# %%
from geopy.distance import distance  # 根据经纬度计算两点距离
from gravity_wave.draw_rain_distribution import Draw, Rain
import matplotlib.pyplot as plt

# %%

# cm = 1/2.54
# fig = plt.figure(figsize=(8*cm, 8*cm), dpi=300)
# ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# dr = Draw(fig, ax)
# # %%
# loc1 = (dr.station['C']['lat'], dr.station['C']['lon'])
# loc2 = (dr.station['B']['lat'], dr.station['B']['lon'])
# distance(loc1, loc2)
dr = Rain()
for i in dr.station.values():
    print(i)