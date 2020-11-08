from parcels import FieldSet, Field, NestedField, VectorField, ParticleFile, ParticleSet, JITParticle, Variable
from parcels import ErrorCode
import numpy as np
from glob import glob
import time as timelib
from datetime import timedelta as delta
from barents_mp_kernels import *
from datetime import timedelta, datetime

def get_nemo_fieldset(run3D=False):
    data_dir = '/Volumes/Barents/OceanParcels/Data/ORCA0083/'
    ufiles = sorted(glob(data_dir+'ORCA0083-N06_20??????d05U.nc'))
    vfiles = sorted(glob(data_dir+'ORCA0083-N06_20??????d05V.nc'))
    mesh_mask = data_dir + 'mesh_hgr.nc'

    if run3D:
        wfiles = sorted(glob(data_dir+'ORCA0083-N06_200?????d05W.nc'))
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

    fieldset = FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation = True)
    if run3D:
        def compute(fieldset):
            fieldset.W.data[:, 0, :, :] = 0.

        fieldset.compute_on_defer = compute

    return fieldset

def set_unbeaching(fieldset):
    files = '/Users/fionarussell/STAS_PROJECT/Barents/ORCA0083-N06_unbeaching_vel.nc'
    filenames = files
    variables = {'Unemo_unbeach': 'unBeachU',
                 'Vnemo_unbeach': 'unBeachV'}
    dimensions = {'lon': 'glamf', 'lat': 'gphif'}
    fieldsetUnBeach = FieldSet.from_nemo(filenames, variables, dimensions, tracer_interp_method='cgrid_velocity')
    fieldset.add_field(fieldsetUnBeach.Unemo_unbeach)
    fieldset.add_field(fieldsetUnBeach.Vnemo_unbeach)

    UVunbeach = VectorField('UVunbeach', fieldset.Unemo_unbeach, fieldset.Vnemo_unbeach)
    fieldset.add_vector_field(UVunbeach)
    
def get_particle_set(fieldset, run3D=False):

    class PlasticParticle(JITParticle):
        age = Variable('age', dtype=np.float32, initial=0.)
        # beached : 0 sea, 1 beached, 2 after non-beach dyn, 3 after beach dyn, 4 please unbeach
        beached = Variable('beached', dtype=np.int32, initial=0.)
        unbeachCount = Variable('unbeachCount', dtype=np.int32, initial=0.)

    # meshgrid containing 11x11 points uniformly distributed in a [0,1]x[0,1] quad
    vec = np.linspace(0, 1, 11)
    xsi, eta = np.meshgrid(vec, vec)

    # get particles in Norwegian sea
    lonCorners = [23.96824026, 24.22713804, 24.26175451, 24.002671]
    latCorners = [71.60693741, 71.58454132, 71.73711395, 71.759758]
    lon_r = (1-xsi)*(1-eta) * lonCorners[0] + xsi*(1-eta) * lonCorners[1] + \
        xsi*eta * lonCorners[2] + (1-xsi)*eta * lonCorners[3]
    lat_r = (1-xsi)*(1-eta) * latCorners[0] + xsi*(1-eta) * latCorners[1] + \
        xsi*eta * latCorners[2] + (1-xsi)*eta * latCorners[3]

    # get particles in Barrents sea
    lonCorners = [38.37941658, 38.63887346, 38.67183721, 38.41217935]
    latCorners = [71.58309555, 71.56196213, 71.71636581, 71.73773575]
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
    times = np.arange(np.datetime64('2011-01-05'), np.datetime64('2011-01-06'))
    # repeat release at the same location for every 80 days
    repeatdt=delta(days=80)

    if run3D:
        return ParticleSet.from_list(fieldset, PlasticParticle,
                                     lat=np.tile(lats, [len(times)]),
                                     lon=np.tile(lons, [len(times)]),
                                     depth=np.tile(depths, [len(times)]),
                                     time=np.repeat(times, len(lons)))
    else:
        return ParticleSet.from_list(fieldset, PlasticParticle,
                                     lon=np.tile(lons, [len(times)]),
                                     lat=np.tile(lats, [len(times)]),
                                     repeatdt=repeatdt) # addition of the repeat argument
                                     

def run_barents_mp(outfile, run3D=False):
    fieldset = get_nemo_fieldset(run3D)
    
    set_unbeaching(fieldset)
    pset = get_particle_set(fieldset, run3D)

    kernel = pset.Kernel(AdvectionRK4_3D) if run3D else pset.Kernel(AdvectionRK4)
    BeachTesting = BeachTesting_3D if run3D else BeachTesting_2D
    kernel += pset.Kernel(BeachTesting) + pset.Kernel(UnBeaching) 
    kernel += pset.Kernel(Ageing)

    pfile = ParticleFile(outfile, pset)
    pfile.write(pset, pset[0].time)

    tic = timelib.time()
    ndays = 84
    
    pset.show(savefile='before')
    
    for d in range(int(ndays/2)):
        day = 2 * d
        print('running %d / %d [time %g s]: %d particles ' % (day, ndays, timelib.time()-tic, len(pset)))
        pset.execute(kernel, runtime=delta(days=2), dt=900, verbose_progress=False, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
        pfile.write(pset, pset[0].time)    
    
    pset.show(savefile='after_84j_repeat')