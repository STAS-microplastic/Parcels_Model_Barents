#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter

Function creating the unbeach velocity for the CMEMS data (A-grid)
"""


import xarray as xr
import numpy as np


data_dir = '/projects/0/topios/hydrodynamic_data/CMEMS/NORTHWESTSHELF_REANALYSIS_PHYS_004_009/MetO-NWS-REAN-PHYS-daily-CUR/2000/01/'
datasetM = xr.open_dataset(data_dir + 'metoffice_foam1_amm7_NWS_RFVL_dm20000101.nc')

dataArrayLonF = datasetM.lon
dataArrayLatF = datasetM.lat

U = np.array(datasetM.vozocrtx)
V = np.array(datasetM.vomecrty)
U[np.isnan(U)] = 0
V[np.isnan(V)] = 0

unBeachU = np.zeros(U.shape[2:])
unBeachV = np.zeros(V.shape[2:])


def island(U, V, j, i):
    if U[0, 0, j, i] == 0 and U[0, 0, j, i+1] == 0 and U[0, 0, j+1, i+1] == 0 and U[0, 0, j+1, i+1] == 0 and\
       V[0, 0, j, i] == 0 and V[0, 0, j, i+1] == 0 and V[0, 0, j+1, i+1] == 0 and V[0, 0, j+1, i+1] == 0:
        return True
    else:
        return False


for j in range(1, U.shape[2]-2):
    for i in range(1, U.shape[3]-2):
        if island(U, V, j, i):
            if not island(U, V, j, i-1):
                unBeachU[j, i] = -1
                unBeachU[j+1, i] = -1
            if not island(U, V, j, i+1):
                unBeachU[j, i+1] = 1
                unBeachU[j+1, i+1] = 1
            if not island(U, V, j-1, i):
                unBeachV[j, i] = -1
                unBeachV[j, i+1] = -1
            if not island(U, V, j+1, i):
                unBeachV[j+1, i] = 1
                unBeachV[j+1, i+1] = 1
            if not island(U, V, j, i-1) and not island(U, V, j+1, i) and island(U, V, j+1, i-1):
                print('Watch out: one cell width land [%d %d]: %g %g' %
                      (j, i, dataArrayLonF[i], dataArrayLatF[j]))
            if not island(U, V, j, i+1) and not island(U, V, j+1, i) and island(U, V, j+1, i+1):
                print('Watch out: one cell width land [%d %d]: %g %g' %
                      (j, i, dataArrayLonF[i], dataArrayLatF[j]))
            if not island(U, V, j, i-1) and not island(U, V, j-1, i) and island(U, V, j-1, i-1):
                print('Watch out: one cell width land [%d %d]: %g %g' %
                      (j, i, dataArrayLonF[i], dataArrayLatF[j]))
            if not island(U, V, j, i+1) and not island(U, V, j-1, i) and island(U, V, j-1, i+1):
                print('Watch out: one cell width land [%d %d]: %g %g' %
                      (j, i, dataArrayLonF[i], dataArrayLatF[j]))

coords = {'lon': dataArrayLonF,
          'lat': dataArrayLatF}
dataArrayUnBeachU = xr.DataArray(unBeachU, name='unBeachU', dims=('lat', 'lon'))
dataArrayUnBeachV = xr.DataArray(unBeachV, name='unBeachV', dims=('lat', 'lon'))

dataset = xr.Dataset()
dataset[dataArrayLonF.name] = dataArrayLonF
dataset[dataArrayLatF.name] = dataArrayLatF
dataset[dataArrayUnBeachU.name] = dataArrayUnBeachU
dataset[dataArrayUnBeachV.name] = dataArrayUnBeachV
dataset.to_netcdf(path='cmems_NWS_rean_004_009_unbeaching_vel.nc', engine='scipy')
