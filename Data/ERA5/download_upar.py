#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:

-----------------------------------------
Author           :Forxd
Version          :1.0
Time：2021/09/09/ 10:30
'''
import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'pressure_level': [
            '1', '2', '3',
            '5', '7', '10',
            '20', '30', '50',
            '70', '100', '125',
            '150', '175', '200',
            '225', '250', '300',
            '350', '400', '450',
            '500', '550', '600',
            '650', '700', '750',
            '775', '800', '825',
            '850', '875', '900',
            '925', '950', '975',
            '1000',
        ],
        'variable': [
            'divergence', 'geopotential', 'ozone_mass_mixing_ratio',
            'potential_vorticity', 'relative_humidity', 'specific_humidity',
            'specific_rain_water_content', 'specific_snow_water_content', 'temperature',
            'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity',
            'vorticity',
        ],
        'year': '2021',
        'month': '07',
        'day': [
            '16', '17', '18',
            '19', '20', '21',
            '22', '23',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'area': [
            60, 90, 10,
            150,
        ],
    },
    'upar.nc')
