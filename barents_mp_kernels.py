#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter (Modified by Benjamin Heutte, December 2020)

Kernels defining floating MP particle dynamics
"""

from parcels import rng as random
import math


def AdvectionRK4(particle, fieldset, time):
    if particle.beached == 0:
        (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)

        (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
        lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)

        (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)

        (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle.beached = 2


def AdvectionRK4_3D(particle, fieldset, time):
    if particle.beached == 0:
        (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        lon1 = particle.lon + u1*.5*particle.dt
        lat1 = particle.lat + v1*.5*particle.dt
        dep1 = particle.depth + w1*.5*particle.dt
        (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1]
        lon2 = particle.lon + u2*.5*particle.dt
        lat2 = particle.lat + v2*.5*particle.dt
        dep2 = particle.depth + w2*.5*particle.dt
        (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2]
        lon3 = particle.lon + u3*particle.dt
        lat3 = particle.lat + v3*particle.dt
        dep3 = particle.depth + w3*particle.dt
        (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle.depth += (w1 + 2*w2 + 2*w3 + w4) / 6. * particle.dt
        particle.beached = 2


def StokesDrag(particle, fieldset, time):
    if particle.lat < 80 and particle.beached == 0:
        (u_uss, v_uss) = fieldset.UVuss[time, particle.depth, particle.lat, particle.lon]
        particle.lon += u_uss * particle.dt
        particle.lat += v_uss * particle.dt
        particle.beached = 3


def BrownianMotion2D(particle, fieldset, time):
    if particle.beached == 0:
        kh_meridional = fieldset.Kh_meridional[time, particle.depth, particle.lat, particle.lon]
        kh_zonal = fieldset.Kh_zonal[time, particle.depth, particle.lat, particle.lon]
        dx = fieldset.meshSize[time, particle.depth, particle.lat, particle.lon]
        dx0 = 1000

        particle.lat += random.uniform(-1., 1.) * math.sqrt(2*math.fabs(particle.dt) * kh_meridional * math.pow(dx/dx0, 1.33))
        particle.lon += random.uniform(-1., 1.) * math.sqrt(2*math.fabs(particle.dt) * kh_zonal      * math.pow(dx/dx0, 1.33))
        particle.beached = 3


def BeachTesting_2D(particle, fieldset, time):
    if particle.beached == 2 or particle.beached == 3:
        (u, v) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        if fabs(u) < 1e-14 and fabs(v) < 1e-14:
            if particle.beached == 2:
                particle.beached = 4
            else:
                particle.beached = 1
        else:
            particle.beached = 0


def BeachTesting_3D(particle, fieldset, time):
    if particle.beached == 2 or particle.beached == 3:
        (u, v, w) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        if fabs(u) < 1e-14 and fabs(v) < 1e-14:
            if particle.beached == 2:
                particle.beached = 4
            else:
                particle.beached = 1
        else:
            particle.beached = 0


def UnBeaching(particle, fieldset, time):
    if particle.beached == 4:
        (ub, vb) = fieldset.UVunbeach[time, particle.depth, particle.lat, particle.lon]
        particle.lon += ub * particle.dt
        particle.lat += vb * particle.dt
        particle.beached = 0
        particle.unbeachCount += 1


def Ageing(particle, fieldset, time):
    particle.age += particle.dt


def DeleteParticle(particle, fieldset, time):
    print("Particle [%d] lost !! (%g %g %g %g)" % (particle.id, particle.lon, particle.lat, particle.depth, particle.time))
    particle.delete()
