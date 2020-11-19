from parcels import FieldSet, Field, NestedField, VectorField, ParticleFile, ParticleSet, JITParticle, Variable
from parcels import ErrorCode
import numpy as np
from glob import glob
import time as timelib
from datetime import timedelta as delta
from barents_mp_kernels import *
from datetime import timedelta, datetime

def get_nemo_fieldset(run3D=False):
    data_dir = 'D:/Barents_Run_Data/ORCA0083-N006/'
    ufiles = sorted(glob(data_dir+'ORCA0083-N06_2011????d05U.nc'))
    vfiles = sorted(glob(data_dir+'ORCA0083-N06_2011????d05V.nc'))
    mesh_mask = data_dir + 'domain/mesh_hgr.nc'

    filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'data': ufiles},
                 'V': {'lon': mesh_mask, 'lat': mesh_mask, 'data': vfiles}}

    variables = {'U': 'uo',
                 'V': 'vo'}
    dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'},
                  'V': {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}}

    fieldset = FieldSet.from_nemo(filenames, variables, dimensions)

    return fieldset

def set_unbeaching(fieldset):
    files = 'D:/Barents_Run_Data/ORCA0083-N006/ORCA0083-N06_unbeaching_vel.nc'
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

    # get particles in Norwegian sea
    npart = 100  # number of particles to be released
    lon = np.linspace(21.4 , 19.6, npart, dtype=np.float32)
    lat = np.linspace(71.3 , 72.3, npart, dtype=np.float32)
    repeatdt = delta(days=60)  # release from the same set of locations every 3 hours

   
    return ParticleSet.from_list(fieldset, PlasticParticle, lon=lon, lat=lat, repeatdt=repeatdt)

                                     
def run_barents_mp(outfile, run3D=False):
    fieldset = get_nemo_fieldset(run3D)
    
    set_unbeaching(fieldset)
    pset = get_particle_set(fieldset, run3D)

    kernel = pset.Kernel(AdvectionRK4)
    BeachTesting = BeachTesting_2D
    kernel += pset.Kernel(BeachTesting) + pset.Kernel(UnBeaching) 
    kernel += pset.Kernel(Ageing)

    pfile = ParticleFile(outfile, pset)
    pfile.write(pset, pset[0].time)

    tic = timelib.time()
    ndays = 365*4
    
    pset.show(savefile='before')
    
    for d in range(int(ndays/2)):
        day = 2 * d
        print('running %d / %d [time %g s]: %d particles ' % (day, ndays, timelib.time()-tic, len(pset)))
        pset.execute(kernel, runtime=delta(days=2), dt=900, verbose_progress=False, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
        pfile.write(pset, pset[0].time)   
        
    pset.show(savefile='after_4y')
    