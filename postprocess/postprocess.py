#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter

Postprocessing class and useful functions
"""


import numpy as np
from netCDF4 import Dataset
import numpy.ctypeslib as npct
from ctypes import c_int, c_float, POINTER
import subprocess
from os import path, getenv
from struct import calcsize


class ParticleData(object):
    def __init__(self, pfile=None, part_3d=False):
        if pfile:
            dataset = Dataset(pfile)
            self.id = dataset.variables['trajectory'][:]
            self.time = dataset.variables['time'][:]
            self.age = dataset.variables['age'][:] / 86400.
            self.lon = dataset.variables['lon'][:]
            self.lat = dataset.variables['lat'][:]
            self.depth = dataset.variables['z'][:] if part_3d is True else None
            self.beached = dataset.variables['beached'][:]
            self.age.fill_value = 9999
            self.age = self.age.filled()
            self.npart = self.lon.shape[0]
            dataset.close()
        else:
            self.id = None
            self.time = None
            self.age = None
            self.lon = None
            self.lat = None
            self.depth = None
            self.beached = None
            self.npart = -1

    def keepParticles(self, pind, keep='part'):
        if keep == 'part':
            self.lon = self.lon[pind[0], :]
            self.lat = self.lat[pind[0], :]
            if self.depth is not None:
                self.depth = self.depth[pind[0], :]
            self.id = self.id[pind[0], :]
            self.time = self.time[pind[0], :]
            self.age = self.age[pind[0], :]
            self.beached = self.beached[pind[0], :]
            self.npart = self.lon.shape[0]
        else:
            self.lon = self.lon[:, pind[1]]
            self.lat = self.lat[:, pind[1]]
            if self.depth is not None:
                self.depth = self.depth[:, pind[1]]
            self.id = self.id[:, pind[1]]
            self.time = self.time[:, pind[1]]
            self.age = self.age[:, pind[1]]
            self.beached = self.beached[:, pind[1]]

    def copy(self):
        p = ParticleData()
        p.id = self.id.copy()
        p.time = self.time.copy()
        p.age = self.age.copy()
        p.lon = self.lon.copy()
        p.lat = self.lat.copy()
        p.depth = self.depth.copy() if self.depth is not None else None
        p.beached = self.beached.copy()
        p.npart = self.npart
        return p


class Utils(object):
    def __init__(self):
        src_file = 'utils.c'
        lib_file = 'utils.so'
        log_file = 'utils.log'
        c_compiler = GNUCompiler()
        c_compiler.compile(src_file, lib_file, log_file)
        self.lib = npct.load_library(lib_file, '.')

    def locate(self, plon, plat, glon, glat):
        pxi = np.floor((plon - glon[0]) / (glon[1]-glon[0]))
        pxi = pxi.astype(c_int)
        pyi = np.floor((plat - glat[0]) / (glat[1]-glat[0]))
        pyi = pyi.astype(c_int)
        if not pxi.flags.c_contiguous:
            pxi = pxi.copy()
            pyi = pyi.copy()
        return pxi, pyi

    def locate_depth(self, pdepth, gdepth):
        pzi = np.floor((pdepth - gdepth[0]) / (gdepth[1]-gdepth[0]))
        pzi = np.where(pzi > len(gdepth)-1, len(gdepth)-1, pzi)
        pzi = pzi.astype(c_int)
        if not pzi.flags.c_contiguous:
            pzi = pzi.copy()
        return pzi

    def histogram(self, pxi, pyi, shape):
        hist = np.zeros(shape, dtype=c_int)
        cfunc = self.lib.histogram
        for tab in [pxi, pyi, hist]:
            if not tab.flags.c_contiguous:
                print 'error'
                exit(-1)

        cfunc(pxi.ctypes.data_as(POINTER(POINTER(c_int))),
              pyi.ctypes.data_as(POINTER(POINTER(c_int))),
              hist.ctypes.data_as(POINTER(POINTER(c_int))),
              pxi.shape[0], pxi.shape[1],
              shape[0], shape[1])
        return hist

    def mean_age(self, pxi, pyi, page, shape, map_count):
        age = np.zeros(shape, dtype=np.float32)
        if not page.dtype == np.float32:
            page = page.astype(np.float32)
        page = page.copy()
        cfunc = self.lib.cumulative_age
        for tab in [pxi, pyi, page, age]:
            if not tab.flags.c_contiguous:
                print 'error'
                exit(-1)

        cfunc(pxi.ctypes.data_as(POINTER(POINTER(c_int))),
              pyi.ctypes.data_as(POINTER(POINTER(c_int))),
              page.ctypes.data_as(POINTER(POINTER(c_float))),
              age.ctypes.data_as(POINTER(POINTER(c_float))),
              pxi.shape[0], pxi.shape[1],
              shape[0], shape[1])

        map_count0 = map_count.copy()
        map_count0 = np.where(map_count0 == 0, -1, map_count0)
        mean_age = age / map_count0
        mean_age = np.where(map_count0 < 0, -1, mean_age)
        return mean_age

    def depth_fraction(self, pxi, pyi, pdepth, depth_lim, shape, map_count):
        depth_count = np.zeros(shape, dtype=np.float32)
        if not pdepth.dtype == np.float32:
            pdepth = pdepth.astype(np.float32)
        pdepth = pdepth.copy()
        cfunc = self.lib.cumulative_depth_lim
        for tab in [pxi, pyi, pdepth, depth_count]:
            if not tab.flags.c_contiguous:
                print 'error'
                exit(-1)

        cfunc.argtypes = [POINTER(POINTER(c_int)),
                          POINTER(POINTER(c_int)),
                          POINTER(POINTER(c_float)),
                          POINTER(POINTER(c_float)),
                          c_int, c_int, c_int, c_int, c_float]

        cfunc(pxi.ctypes.data_as(POINTER(POINTER(c_int))),
              pyi.ctypes.data_as(POINTER(POINTER(c_int))),
              pdepth.ctypes.data_as(POINTER(POINTER(c_float))),
              depth_count.ctypes.data_as(POINTER(POINTER(c_float))),
              pxi.shape[0], pxi.shape[1],
              shape[0], shape[1], depth_lim)

        map_count0 = map_count.copy()
        map_count0 = np.where(map_count0 == 0, -1, map_count0)
        depth_frac = depth_count / map_count0
        depth_frac = np.where(map_count0 < 0, -1, depth_frac)
        return 100*depth_frac

    def touched(self, pxi, pyi, shape):
        touch = np.zeros(shape, dtype=np.float32)
        cfunc = self.lib.touched

        cfunc(pxi.ctypes.data_as(POINTER(POINTER(c_int))),
              pyi.ctypes.data_as(POINTER(POINTER(c_int))),
              touch.ctypes.data_as(POINTER(POINTER(c_float))),
              pxi.shape[0], pxi.shape[1],
              shape[0], shape[1])

        touch /= pxi.shape[0]
        return 100*touch

    def zone_concentration(self, zones, page, exportdt, pxi, pyi):
        page_int = page.copy().astype(c_int) / exportdt
        nstep = int(np.max(page_int))+1
        nzone = int(np.max(zones))
        zone_concentration = np.zeros((nzone, nstep), dtype=np.float32)
        zones = zones.astype(c_int)
        cfunc = self.lib.zone_concentration
        for tab in [pxi, pyi, page_int, zones]:
            if not tab.flags.c_contiguous:
                print 'error'
                exit(-1)

        cfunc(pxi.ctypes.data_as(POINTER(POINTER(c_int))),
              pyi.ctypes.data_as(POINTER(POINTER(c_int))),
              page_int.ctypes.data_as(POINTER(POINTER(c_int))),
              zones.ctypes.data_as(POINTER(POINTER(c_int))),
              zone_concentration.ctypes.data_as(POINTER(POINTER(c_float))),
              pxi.shape[0], pxi.shape[1],
              zones.shape[0], zones.shape[1],
              nzone, nstep)

        zone_concentration = zone_concentration / np.sum(zone_concentration, axis=0)
        return zone_concentration


def right_of(lon, lat, lon1, lat1, lon2, lat2):
    a = (lat2-lat1+0.)/(lon2-lon1)
    b = lat1-a*lon1
    cond = lat-a*lon-b < 0
    cond.shape
    return cond


def NWcontinental_shelf_zones(glon, glat):
    zones = np.zeros(glon.shape)

    # Danish Strait
    c1 = zones == 0
    c1 *= np.logical_not(right_of(glon, glat, 7.8, 58, 8.6, 56))
    c1 *= glon < 14
    c1 *= glat > 54
    c1 *= glat < 60
    zones[c1] = 1

    # North Sea
    c1 = zones == 0
    c1 *= glat < 58.5
    c1 *= np.logical_or(np.logical_not(right_of(glon, glat, -5, 58.5, 0, 50)),
                        np.logical_and(np.logical_not(right_of(glon, glat, -5.5, 50, -4, 48.2)),
                                       np.logical_and(right_of(glon, glat, -5.5, 50, 0, 52),
                                                      glat > 48)))
    c1 *= right_of(glon, glat, 8, 58, 13, 50)
    zones[c1] = 2

    # Norwegian Sea
    c1 = zones == 0
    c1 *= glat >= 58.5
    c1 *= glat < 73
    c1 *= right_of(glon, glat, -5, 58, 10, 74)
    c1 *= glon <= 27.8
    c1 *= np.logical_not(right_of(glon, glat, 10, 58, 28, 70))
    zones[c1] = 3

    # Arctic
    c1 = zones == 0
    c1 *= np.logical_or(glat > 70, np.logical_and(glon > 27.8, glat > 62))
    c1 *= np.logical_or(glat > 80, glon > -30)
    zones[c1] = 4

    # Atlantic
    c1 = zones == 0
    c1 *= glon < 7
    zones[c1] = 5

    # Baltic Sea
    c1 = zones == 0
    c1 *= glon > 10
    c1 *= glon < 30
    c1 *= glat > 53
    c1 *= glat < 80
    zones[c1] = 6

    return zones


def get_package_dir():
    return path.abspath(path.dirname(__file__))


class Compiler(object):
    """A compiler object for creating and loading shared libraries.

    :arg cc: C compiler executable (uses environment variable ``CC`` if not provided).
    :arg cppargs: A list of arguments to the C compiler (optional).
    :arg ldargs: A list of arguments to the linker (optional)."""

    def __init__(self, cc=None, cppargs=None, ldargs=None):
        if cppargs is None:
            cppargs = []
        if ldargs is None:
            ldargs = []

        self._cc = getenv('CC') if cc is None else cc
        self._cppargs = cppargs
        self._ldargs = ldargs

    def compile(self, src, obj, log):
        if path.isfile(obj) and path.getmtime(obj) > path.getmtime(src):
            print('%s does already exist and is newer than %s. Ignoring compilation' % (obj, src))
            return
        cc = [self._cc] + self._cppargs + ['-o', obj, src] + self._ldargs
        with open(log, 'w') as logfile:
            logfile.write("Compiling: %s\n" % " ".join(cc))
            try:
                subprocess.check_call(cc, stdout=logfile, stderr=logfile)
            except OSError:
                err = """OSError during compilation
Please check if compiler exists: %s""" % self._cc
                raise RuntimeError(err)
            except subprocess.CalledProcessError:
                with open(log, 'r') as logfile2:
                    err = """Error during compilation:
Compilation command: %s
Source file: %s
Log file: %s

Log output: %s""" % (" ".join(cc), src, logfile.name, logfile2.read())
                raise RuntimeError(err)
        print("Compiled %s ==> %s\n" % (src, obj))


class GNUCompiler(Compiler):
    """A compiler object for the GNU Linux toolchain.

    :arg cppargs: A list of arguments to pass to the C compiler
         (optional).
    :arg ldargs: A list of arguments to pass to the linker (optional)."""
    def __init__(self, cppargs=None, ldargs=None):
        if cppargs is None:
            cppargs = []
        if ldargs is None:
            ldargs = []

        opt_flags = ['-g', '-O3']
        arch_flag = ['-m64' if calcsize("P") is 8 else '-m32']
        cppargs = ['-Wall', '-fPIC', '-I%s' % path.join(get_package_dir(), 'include')] + opt_flags + cppargs
        cppargs += arch_flag
        ldargs = ['-shared'] + ldargs + arch_flag
        super(GNUCompiler, self).__init__("gcc", cppargs=cppargs, ldargs=ldargs)
