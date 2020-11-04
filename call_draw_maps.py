#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter

Master postprocessing file
"""


from draw_maps import draw


path = '../ncfiles/'

draw(path+'nemo0083.nc')
draw(path+'nemo025.nc')
draw(path+'nemo0083_cmems.nc')
draw(path+'nemo0083_stokes.nc')
draw(path+'nemo0083_diff.nc')
draw(path+'nemo0083_3D.nc', part_3d=True)
