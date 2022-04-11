# %%
import matplotlib

# %%
print("hello")
# %%
from matplotlib import font_manager

for font in font_manager.fontManager.ttflist:
     print(font.name,) # 查看字体名


