# %%
import cartopy.crs as ccrs
import re

# %%
def get_information(flnm):
    """根据namelist.wps文件，获取地图的基本信息

    Args:
        flnm ([type]): [description]

    Returns:
        [type]: [description]
    """

    ## 设置正则表达式信息
    pattern = {}
    pattern['dx'] = 'dx\s*=\s*\d*,'
    pattern['dy'] = 'dy\s*=\s*\d*,'
    pattern['max_dom'] = 'max_dom\s*=\s*\d\s*,'
    pattern[
        'parent_grid_ratio'] = 'parent_grid_ratio\s*=\s*\d,\s*\d,\s*\d,\s*\d,'
    pattern['j_parent_start'] = 'j_parent_start\s*=\s*\d,\s*\d*,\s*\d*,\s*\d*,'
    pattern['i_parent_start'] = 'i_parent_start\s*=\s*\d,\s*\d*,\s*\d*,\s*\d*,'
    pattern['e_sn'] = 'e_sn\s*=\s*\d*,\s*\d*,\s*\d*,\s*\d*'
    pattern['e_we'] = 'e_we\s*=\s*\d*,\s*\d*,\s*\d*,\s*\d*'
    pattern['ref_lat'] = 'ref_lat\s*=\s*\d*.?\d*,'
    pattern['ref_lon'] = 'ref_lon\s*=\s*\d*.?\d*,'
    pattern['true_lat1'] = 'truelat1\s*=\s*\d*.?\d*,'
    pattern['true_lat2'] = 'truelat2\s*=\s*\d*.?\d*,'

    f = open(flnm)
    fr = f.read()

    def get_var(var, pattern=pattern, fr=fr):
        """处理正则表达式得到的数据"""
        ff1 = re.search(pattern[var], fr, flags=0)
        str_f1 = ff1.group(0)

        str1 = str_f1.replace('=', ',')
        aa = str1.split(',')
        bb = []
        for i in aa[1:]:
            if i != '':
                bb.append(i.strip())
        return bb

    dic_return = {}
    aa = get_var('parent_grid_ratio')

    var_list = [
        'dx',
        'dy',
        'max_dom',
        'parent_grid_ratio',
        'j_parent_start',
        'i_parent_start',
        'e_sn',
        'e_we',
        'ref_lat',
        'ref_lon',
        'true_lat1',
        'true_lat2',
    ]

    for i in var_list:
        aa = get_var(i)
        if i in [
                'parent_grid_ratio',
                'j_parent_start',
                'i_parent_start',
                'e_we',
                'e_sn',
        ]:
            bb = aa
            bb = [float(i) for i in bb]
        else:
            bb = float(aa[0])
        dic_return[i] = bb

    return dic_return


def get_proj():
    
    file_folder = "./"
    file_name = "namelist.wps"
    flnm = file_folder + file_name

    info = get_information(flnm)  # 获取namelist.wps文件信息
    # print(info['ref_lat'])
    # ax = create_map(info)  # 在domain1区域内，添加地理信息，创建底图

    

    """创建一个包含青藏高原区域的Lambert投影的底图

    Returns:
        ax: 坐标图对象
    """
    ## --创建画图空间

    ref_lat = info['ref_lat']
    ref_lon = info['ref_lon']
    true_lat1 = info['true_lat1']
    true_lat2 = info['true_lat2']
    false_easting = (info['e_we'][0] - 1) / 2 * info['dx']
    false_northing = (info['e_sn'][0] - 1) / 2 * info['dy']

    # false_easting = (info['e_we'][2] - 1) / 2 * info['dx']/9
    # false_northing = (info['e_sn'][2] - 1) / 2 * info['dy']/9
    # print(info['e_we'])
    
    extent_area = [0, false_easting * 2, 0, false_northing * 2]
                  
    proj_lambert = ccrs.LambertConformal(
        central_longitude=ref_lon,
        central_latitude=ref_lat,
        standard_parallels=(true_lat1, true_lat2),
        cutoff=-30,
        false_easting=false_easting,
        false_northing=false_northing,
    )
    return proj_lambert, extent_area


if __name__ == '__main__':
    aa = get_proj()
    print(aa)
    
# %%
