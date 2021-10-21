import cdsapi, os
from pathlib import Path


c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
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
        'year': [
            '2005',
        ],
        'month': [
            '08','09', '10',
        ],
        'area': [
            35, -130, 10,
            10,
        ],
        'variable': [
            'specific_humidity', 'temperature',
        ],
        'time': '00:00',
    },
    'era5_pl.nc')

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable':['mean_sea_level_pressure','land_sea_mask','sea_surface_temperature'],
        'year': [
            '2005',
        ],
        'month': [
            '08','09', '10',
        ],
        'time': '00:00',
        'area': [
            35, -130, 10,
            10,
        ],
    },
    'era5_mslp.nc')


# Merge with CDO and remove
abspath=str(Path().absolute())
savepath="/".join(abspath.split("/")[:-1])+"/Data/"
print(savepath)
os.system("cdo -O merge era5_mslp.nc era5_pl.nc %sera5_tc_lev.nc"%savepath)
os.system("cdo -O invertlev %sera5_tc_lev.nc %sera5_tc_ilev.nc"%(savepath,savepath))
rms=["era5_mslp.nc", "era5_pl.nc"]
for i in rms: os.remove(i)


