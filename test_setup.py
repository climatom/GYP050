{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "85709738",
   "metadata": {},
   "outputs": [],
   "source": [
    "import os, cdsapi\n",
    "import zipfile\n",
    "#import xarray as xa\n",
    "import numpy as np, pandas as pd\n",
    "import matplotlib.pyplot as plt"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "fd18b76e",
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "2021-10-18 11:20:52,128 INFO Welcome to the CDS\n",
      "2021-10-18 11:20:52,130 INFO Sending request to https://cds.climate.copernicus.eu/api/v2/resources/reanalysis-era5-pressure-levels\n",
      "2021-10-18 11:20:52,151 INFO Request is queued\n",
      "2021-10-18 11:39:11,068 INFO Request is completed\n",
      "2021-10-18 11:39:11,069 INFO Downloading https://download-0011.copernicus-climate.eu/cache-compute-0011/cache/data1/adaptor.mars.internal-1634553460.889749-23572-13-ae2003ce-f3cc-4e86-ac4c-b52ea686680c.nc to test.nc (2M)\n",
      "2021-10-18 11:39:11,487 INFO Download rate 4.8M/s   \n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "Result(content_length=2086232,content_type=application/x-netcdf,location=https://download-0011.copernicus-climate.eu/cache-compute-0011/cache/data1/adaptor.mars.internal-1634553460.889749-23572-13-ae2003ce-f3cc-4e86-ac4c-b52ea686680c.nc)"
      ]
     },
     "execution_count": 2,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "#!/usr/bin/env python\n",
    "c = cdsapi.Client()\n",
    "c.retrieve(\"reanalysis-era5-pressure-levels\",\n",
    "{\n",
    "\"variable\": \"temperature\",\n",
    "\"pressure_level\": \"1000\",\n",
    "\"product_type\": \"reanalysis\",\n",
    "\"year\": \"2020\",\n",
    "\"month\": \"05\",\n",
    "\"day\": \"01\",\n",
    "\"time\": \"12:00\",\n",
    "\"format\": \"netcdf\"\n",
    "}, \"test.nc\")"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.8"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
