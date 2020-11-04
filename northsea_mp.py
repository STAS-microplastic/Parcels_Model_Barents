#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter

Simulation master file
"""

from run_northsea_mp import run_northsea_mp

path = './ncfiles/'

run_northsea_mp(path+'nemo0083.nc',
                nemo_res='0083', cmems=False, stokes=False, diffusion=0)
run_northsea_mp(path+'nemo025.nc',
                nemo_res='025', cmems=False, stokes=False, diffusion=0)
run_northsea_mp(path+'nemo0083_cmems.nc',
                nemo_res='0083', cmems=True, stokes=False, diffusion=0)
run_northsea_mp(path+'nemo0083_stokes.nc',
                nemo_res='0083', cmems=False, stokes=True, diffusion=0)
run_northsea_mp(path+'nemo0083_diff.nc',
                nemo_res='0083', cmems=False, stokes=False, diffusion=1)
