#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter (Modified by Benjamin Heutte, December 2020)

Function creating the unbeach velocity for the NEMO data (C-grid)
"""


import xarray as xr
import numpy as np

data_dir = 'D:/Barents_Run_Data/ORCA0083-N006/' % res

datasetM = xr.open_dataset(data_dir + 'domain/mesh_hgr.nc', decode_cf=False)
datasetU = xr.open_dataset(data_dir + 'ORCA0083-N06_20090105d05U.nc' % (res, day), decode_cf=False)
datasetV = xr.open_dataset(data_dir + 'ORCA0083-N06_20090105d05V.nc' % (res, day), decode_cf=False)

dataArrayLonF = datasetM.glamf
dataArrayLatF = datasetM.gphif
dataArrayTime = datasetM.time
dataArrayTime.attrs['time_origin'] = '1950-JAN-01 00:00:00'
dataArrayTime.attrs['units'] = 'seconds since 1950-01-01 00:00:00'

U = np.array(datasetU.uo)
V = np.array(datasetV.vo)
U[np.isnan(U)] = 0
V[np.isnan(V)] = 0
dataArrayLonU = datasetU.nav_lon
dataArrayLatU = datasetU.nav_lat
dataArrayLonV = datasetV.nav_lon
dataArrayLatV = datasetV.nav_lat

unBeachU = np.zeros(U.shape[2:])
unBeachV = np.zeros(V.shape[2:])

for j in range(1, U.shape[2]-2):
    for i in range(1, U.shape[3]-2):
        if U[0, 0, j+1, i] == 0 and U[0, 0, j+1, i+1] == 0 and V[0, 0, j, i+1] == 0 and V[0, 0, j+1, i+1] == 0:
            if abs(U[0, 0, j+1, i-1]) > 1e-10:
                unBeachU[j+1, i] = -1
            if abs(U[0, 0, j+1, i+2]) > 1e-10:
                unBeachU[j+1, i+1] = 1
            if abs(V[0, 0, j-1, i+1]) > 1e-10:
                unBeachV[j, i+1] = -1
            if abs(V[0, 0, j+2, i+1]) > 1e-10:
                unBeachV[j+1, i+1] = 1

coordsU = {'glamu': dataArrayLonU,
           'gphiu': dataArrayLatU}
dataArrayUnBeachU = xr.DataArray(unBeachU, name='unBeachU', coords=coordsU, dims=('y', 'x'))
coordsV = {'glamv': dataArrayLonV,
           'gphiv': dataArrayLatV}
dataArrayUnBeachV = xr.DataArray(unBeachV, name='unBeachV', coords=coordsV, dims=('y', 'x'))

dataset = xr.Dataset()
dataset[dataArrayLonF.name] = dataArrayLonF
dataset[dataArrayLatF.name] = dataArrayLatF
dataset[dataArrayTime.name] = dataArrayTime
dataset[dataArrayUnBeachU.name] = dataArrayUnBeachU
dataset[dataArrayUnBeachV.name] = dataArrayUnBeachV
dataset.to_netcdf(path='D:/Barents_Run_Data/ORCA0083-N006/ORCA0083-N06_unbeaching_vel.nc' % res, engine='scipy')
