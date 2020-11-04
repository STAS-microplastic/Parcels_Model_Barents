#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter

Function running floating MP particles released in Thames and Rhiver estuaries
"""


from parcels import FieldSet, Field, NestedField, VectorField, ParticleFile, ParticleSet, JITParticle, Variable
from parcels import ErrorCode
import numpy as np
from glob import glob
import time as timelib
from datetime import timedelta as delta
from northsea_mp_kernels import *


def get_nemo_fieldset(res='0083', run3D=False):
    data_dir = '/projects/0/topios/hydrodynamic_data/NEMO-MEDUSA/ORCA%s-N006/' % res
    ufiles = sorted(glob(data_dir+'means/ORCA%s-N06_200?????d05U.nc' % res))
    vfiles = sorted(glob(data_dir+'means/ORCA%s-N06_200?????d05V.nc' % res))
    mesh_mask = data_dir + 'domain/coordinates.nc'

    if run3D:
        wfiles = sorted(glob(data_dir+'means/ORCA%s-N06_200?????d05W.nc' % res))
        filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
                     'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles},
                     'W': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles}}
        variables = {'U': 'uo',
                     'V': 'vo',
                     'W': 'wo'}
        dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                      'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                      'W': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}}
    else:
        filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
                     'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles}}

        variables = {'U': 'uo',
                     'V': 'vo'}
        dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'},
                      'V': {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}}

    fieldset = FieldSet.from_nemo(filenames, variables, dimensions)
    fieldset.nemo_res = res
    if run3D:
        def compute(fieldset):
            fieldset.W.data[:, 0, :, :] = 0.

        fieldset.compute_on_defer = compute

    fieldset.cmems = False
    return fieldset


def set_cmems(fieldset):
    data_dir = '/projects/0/topios/hydrodynamic_data/CMEMS/NORTHWESTSHELF_REANALYSIS_PHYS_004_009/MetO-NWS-REAN-PHYS-daily-CUR/'
    fnames = []
    years = range(2000, 2005)
    for y in years:
        basepath = data_dir + str(y) + '/*/' + 'metoffice_foam1_amm7_NWS_RFVL_dm*.nc'
        fnames += sorted(glob(str(basepath)))
    dimensionsU = {'lon': 'lon', 'lat': 'lat', 'time': 'time'}
    dimensionsV = {'lon': 'lon', 'lat': 'lat', 'time': 'time'}
    indices = {'lon': range(1, 296), 'lat': range(1, 374)}  # cmems puts nan values at its borders
    Ucmems = Field.from_netcdf(fnames, ('Ucmems', 'vozocrtx'), dimensionsU, fieldtype='U', indices=indices, allow_time_extrapolation=False)
    Vcmems = Field.from_netcdf(fnames, ('Vcmems', 'vomecrty'), dimensionsV, fieldtype='V', indices=indices, allow_time_extrapolation=False, grid=Ucmems.grid, dataFiles=Ucmems.dataFiles)
    fieldset.add_field(Ucmems)
    fieldset.add_field(Vcmems)
    fieldset.Ucmems.vmax = 5
    fieldset.Vcmems.vmax = 5

    fieldset.Unemo = fieldset.U
    fieldset.Unemo.name = 'Unemo'
    fieldset.Vnemo = fieldset.V
    fieldset.Vnemo.name = 'Vnemo'

    U = NestedField('U', [fieldset.Ucmems, fieldset.Unemo])
    V = NestedField('V', [fieldset.Vcmems, fieldset.Vnemo])
    fieldset.U = U
    fieldset.V = V

    fieldset.cmems = True


def set_unbeaching(fieldset):
    files = '/home/philippe/data/ORCA%s-N06_unbeaching_vel.nc' % fieldset.nemo_res
    filenames = files
    variables = {'Unemo_unbeach': 'unBeachU',
                 'Vnemo_unbeach': 'unBeachV'}
    dimensions = {'lon': 'glamf', 'lat': 'gphif'}
    fieldsetUnBeach = FieldSet.from_nemo(filenames, variables, dimensions, tracer_interp_method='cgrid_velocity')
    fieldset.add_field(fieldsetUnBeach.Unemo_unbeach)
    fieldset.add_field(fieldsetUnBeach.Vnemo_unbeach)

    if fieldset.cmems:
        fname = '/home/philippe/data/cmems_NWS_rean_004_009_unbeaching_vel.nc'
        dimensionsU = {'lon': 'lon', 'lat': 'lat'}
        Ucmems_unbeach = Field.from_netcdf(fname, ('Ucmems_unbeach', 'unBeachU'), dimensionsU, fieldtype='U')
        dimensionsV = {'lon': 'lon', 'lat': 'lat'}
        Vcmems_unbeach = Field.from_netcdf(fname, ('Vcmems_unbeach', 'unBeachV'), dimensionsV, fieldtype='V')
        fieldset.add_field(Ucmems_unbeach)
        fieldset.add_field(Vcmems_unbeach)

        UVnemo_unbeach = VectorField('UVnemo_unbeach', fieldset.Unemo_unbeach, fieldset.Vnemo_unbeach)
        UVcmems_unbeach = VectorField('UVcmems_unbeach', fieldset.Ucmems_unbeach, fieldset.Vcmems_unbeach)
        UVunbeach = NestedField('UVunbeach', [UVcmems_unbeach, UVnemo_unbeach])
        fieldset.add_vector_field(UVunbeach)
    else:
        UVunbeach = VectorField('UVunbeach', fieldset.Unemo_unbeach, fieldset.Vnemo_unbeach)
        fieldset.add_vector_field(UVunbeach)


def set_stokes(fieldset):
    data_dir = '/projects/0/topios/hydrodynamic_data/WWIII/'
    fnames = []
    years = range(2000, 2005)
    for y in years:
        basepath = data_dir + str(y) + '/ww3.*_uss.nc'
        fnames += sorted(glob(str(basepath)))
    dimensionsU = {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}
    dimensionsV = {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}
    Uuss = Field.from_netcdf(fnames, ('Uuss', 'uuss'), dimensionsU, fieldtype='U', allow_time_extrapolation=False)
    Vuss = Field.from_netcdf(fnames, ('Vuss', 'vuss'), dimensionsV, fieldtype='V', allow_time_extrapolation=False, grid=Uuss.grid, dataFiles=Uuss.dataFiles)
    fieldset.add_field(Uuss)
    fieldset.add_field(Vuss)
    fieldset.Uuss.vmax = 5
    fieldset.Vuss.vmax = 5
    uv_uss = VectorField('UVuss', fieldset.Uuss, fieldset.Vuss)
    fieldset.add_vector_field(uv_uss)


def set_diffusion(fieldset, diffusivity):
    fname = '/home/philippe/data/ORCA0083-N06_meshSize.nc'
    dimensions = {'lon': 'glamt', 'lat': 'gphit'}
    meshSize = Field.from_netcdf(fname, 'meshSize', dimensions, interp_method='nearest')
    fieldset.add_field(meshSize)
    fieldset.add_field(Field('Kh_zonal', data=diffusivity*np.ones(meshSize.data.shape),
                             grid=meshSize.grid, mesh='spherical'))
    fieldset.add_field(Field('Kh_meridional', data=diffusivity*np.ones(meshSize.data.shape),
                             grid=meshSize.grid, mesh='spherical'))


def get_particle_set(fieldset, run3D=False):

    class PlasticParticle(JITParticle):
        age = Variable('age', dtype=np.float32, initial=0.)
        # beached : 0 sea, 1 beached, 2 after non-beach dyn, 3 after beach dyn, 4 please unbeach
        beached = Variable('beached', dtype=np.int32, initial=0.)
        unbeachCount = Variable('unbeachCount', dtype=np.int32, initial=0.)

    # meshgrid containing 11x11 points uniformly distributed in a [0,1]x[0,1] quad
    vec = np.linspace(0, 1, 11)
    xsi, eta = np.meshgrid(vec, vec)

    # get particles in Rhine estuary
    lonCorners = [2.96824026, 3.22713804, 3.26175451, 3.002671]
    latCorners = [51.60693741, 51.58454132, 51.73711395, 51.759758]
    lon_r = (1-xsi)*(1-eta) * lonCorners[0] + xsi*(1-eta) * lonCorners[1] + \
        xsi*eta * lonCorners[2] + (1-xsi)*eta * lonCorners[3]
    lat_r = (1-xsi)*(1-eta) * latCorners[0] + xsi*(1-eta) * latCorners[1] + \
        xsi*eta * latCorners[2] + (1-xsi)*eta * latCorners[3]

    # get particles in Thames estuary
    lonCorners = [1.37941658, 1.63887346, 1.67183721, 1.41217935]
    latCorners = [51.58309555, 51.56196213, 51.71636581, 51.73773575]
    lon_t = (1-xsi)*(1-eta) * lonCorners[0] + xsi*(1-eta) * lonCorners[1] + \
        xsi*eta * lonCorners[2] + (1-xsi)*eta * lonCorners[3]
    lat_t = (1-xsi)*(1-eta) * latCorners[0] + xsi*(1-eta) * latCorners[1] + \
        xsi*eta * latCorners[2] + (1-xsi)*eta * latCorners[3]

    lons = np.concatenate((lon_r.flatten(), lon_t.flatten()))
    lats = np.concatenate((lat_r.flatten(), lat_t.flatten()))

    if run3D:
        seeding_depths = np.array([.1, .5, 1.])

        depths = np.repeat(seeding_depths, len(lons))
        lons = np.tile(lons, len(seeding_depths))
        lats = np.tile(lats, len(seeding_depths))

    # gather all particles, released every day for one year
    times = np.arange(np.datetime64('2000-01-05'), np.datetime64('2001-01-05'))

    if run3D:
        return ParticleSet.from_list(fieldset, PlasticParticle,
                                     lon=np.tile(lons, [len(times)]),
                                     lat=np.tile(lats, [len(times)]),
                                     depth=np.tile(depths, [len(times)]),
                                     time=np.repeat(times, len(lons)))
    else:
        return ParticleSet.from_list(fieldset, PlasticParticle,
                                     lon=np.tile(lons, [len(times)]),
                                     lat=np.tile(lats, [len(times)]),
                                     time=np.repeat(times, len(lons)))


def run_northsea_mp(outfile, nemo_res='0083', cmems=False, stokes=False, diffusion=0, run3D=False):
    fieldset = get_nemo_fieldset(nemo_res, run3D)
    if cmems:
        set_cmems(fieldset)
    if stokes:
        set_stokes(fieldset)
    if diffusion > 0:
        set_diffusion(fieldset, diffusion)

    set_unbeaching(fieldset)
    pset = get_particle_set(fieldset, run3D)

    kernel = pset.Kernel(AdvectionRK4_3D) if run3D else pset.Kernel(AdvectionRK4)
    BeachTesting = BeachTesting_3D if run3D else BeachTesting_2D
    kernel += pset.Kernel(BeachTesting) + pset.Kernel(UnBeaching)
    if stokes:
        kernel += pset.Kernel(StokesDrag) + pset.Kernel(BeachTesting)
    if diffusion > 0:
        kernel += pset.Kernel(BrownianMotion2D) + pset.Kernel(BeachTesting)
    kernel += pset.Kernel(Ageing)

    pfile = ParticleFile(outfile, pset)
    pfile.write(pset, pset[0].time)

    tic = timelib.time()
    ndays = 365*4+100
    for d in range(ndays/2):
        day = 2 * d
        print('running %d / %d [time %g s]: %d particles ' % (day, ndays, timelib.time()-tic, len(pset)))
        pset.execute(kernel, runtime=delta(days=2), dt=900, verbose_progress=False, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
        pfile.write(pset, pset[0].time)
