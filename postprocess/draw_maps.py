#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter

Function postprocessing floating MP particles released
in Thames and Rhiver estuaries
"""


import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import AutoMinorLocator, LinearLocator, FixedFormatter
import numpy as np
import postprocess as postpro
from postprocess import NWcontinental_shelf_zones
from parcels import RectilinearZGrid, Field
from os import path


def draw(filename, part_3d=False):
    p = postpro.ParticleData(filename, part_3d=part_3d)
    basename = path.splitext(path.basename(filename))[0]
    print('Total number of particles: %d' % p.npart)

    pland = np.logical_and(abs(p.lon[:, 0]-p.lon[:, -1]) < 1e-4,
                           abs(p.lat[:, 0]-p.lat[:, -1]) < 1e-4)
    pind = np.where(np.logical_not(pland))
    p.keepParticles(pind)
    wetParticles = p.npart
    print('Number of initially wet particles: %d' % wetParticles)

    p.origin = np.where(p.lon[:, 0] > 2, 0, 1)

    pind = np.logical_and(p.age > 0*365, p.age <= 4*365)
    p.keepParticles(pind, 'time')

    glon = np.arange(-120, 120, 1)
    glat = np.arange(35, 90, .5)
    glon_m, glat_m = np.meshgrid(glon, glat)

    utils = postpro.Utils()

    pxi, pyi = utils.locate(p.lon, p.lat, glon, glat)
    # map_count = utils.histogram(pxi, pyi, (len(glat), len(glon)))
    # map_mean_age = utils.mean_age(pxi, pyi, p.age, (len(glat), len(glon)), map_count)
    map_touched = utils.touched(pxi, pyi, (len(glat), len(glon)))

    grid = RectilinearZGrid(glon, glat, mesh='spherical')
    cellSizeField = Field('cellsize', np.zeros((len(glat), len(glon))), grid=grid)
    cellSizeField.data[0, :] = cellSizeField.cell_areas()
    pind = np.logical_and(p.age > 0*365, p.age <= 2*365)
    p.keepParticles(pind, 'time')
    pxi, pyi = utils.locate(p.lon, p.lat, glon, glat)
    map_density_4th_year = utils.histogram(pxi, pyi, (len(glat), len(glon)))
    map_density_4th_year = map_density_4th_year.astype(np.float32) / cellSizeField.data[0, :]

    def plot(data, title, cmap_name='plasma_r', log=False, vmin=None, vmax=None, under=None, over=None, show=True, fname=''):

        if fname:
            fig = plt.figure(figsize=(14, 8.5), dpi=150, facecolor='w', edgecolor='k')
        else:
            fig = plt.figure(figsize=(14, 8.5), dpi=60, facecolor='w', edgecolor='k')
        ax = fig.add_axes([.05, .05, .9, .9])

        m = Basemap(width=3.4e6, height=2.5e6,
                    resolution='l', projection='stere',
                    lat_ts=70, lat_0=75, lon_0=40.)

        m.drawcoastlines(zorder=3)
        m.drawparallels(np.arange(-90, 91, 2), labels=[True, False, False, False], fontsize=15, zorder=3)
        m.drawmeridians(np.arange(-180, 181, 10), labels=[False, False, False, True], fontsize=15, zorder=3)
        m.fillcontinents(color='blanchedalmond')

        legend_text = title
        xshift = 0 if len(legend_text) > 6 else .095
        ax.text(0.92+xshift, 1.05, legend_text,
                transform=ax.transAxes,
                verticalalignment='top',
                fontsize=24)

        if not under and not over:
            extend = 'neither'
            under = 'white'
            over = 'white'
        elif not under:
            extend = 'max'
            under = 'white'
        elif not over:
            extend = 'min'
            over = 'white'
        else:
            extend = 'both'
        cmap = plt.get_cmap(cmap_name)
        cmap.set_under(under)
        cmap.set_over(over)

        xs, ys = m(glon_m, glat_m)

        if vmin is None:
            vmin = np.min(data)
        if vmax is None:
            vmax = np.max(data)
        if under == 'white':
            data = np.ma.masked_where(data < vmin+1e-15, data)
        print('min/max: %f %f' %( np.min(data), np.max(data)))

        if log:
            m.pcolormesh(xs, ys, data, cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax), zorder=2)
        else:
            m.pcolormesh(xs, ys, data, cmap=cmap, vmin=vmin, vmax=vmax, zorder=2)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        cbar = plt.colorbar(extend=extend, cax=cax)
        cbar.ax.tick_params(labelsize=15)

        if show:
            plt.show()
        if fname:
            plt.savefig(fname)
        plt.close(fig)

    def plot_vertProfile(dens, dens_lat, dens_depth, show, fname):
        if fname:
            fig = plt.figure(figsize=(14, 8.5), dpi=150, facecolor='w', edgecolor='k')
        else:
            fig = plt.figure(figsize=(14, 8.5), dpi=60, facecolor='w', edgecolor='k')

        ax1 = fig.add_axes([.07, .06, .83, .68])
        cmap = plt.get_cmap('plasma')
        cmap.set_over('green')
        im1 = ax1.pcolormesh(dens_lat, -dens_depth, dens, norm=colors.LogNorm(vmin=1, vmax=8e5), cmap=cmap)
        ax1.yaxis.set_major_locator(LinearLocator(numticks=5))
        ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax1.yaxis.set_major_formatter(FixedFormatter(list(np.arange(2000, 10, -500)) + ["0 m"]))
        ax1.tick_params(which='major', length=7)
        ax1.tick_params(which='minor', length=5)
        for label in ax1.get_ymajorticklabels():
            label.set_fontsize(20)
        for label in ax1.get_xmajorticklabels():
            label.set_fontsize(20)
        ax1.set_xticks([40+5*i for i in range(10)])
        ax1.set_xticklabels([u"%g\u00B0N" % (40+5*i) for i in range(10)])

        ax2 = fig.add_axes([.07, .78, .83, .19])
        dens2 = np.where(dens == 0, np.nan, dens)
        cmap = plt.get_cmap('jet')
        cmap.set_over('green')
        im2 = ax2.pcolormesh(dens_lat[:11, :], -dens_depth[:11, :], dens2[:11, :], vmin=1, vmax=8e5, cmap=cmap)
        ax2.yaxis.set_major_locator(LinearLocator(numticks=3))
        ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax2.yaxis.set_major_formatter(FixedFormatter(["50", "25", "0 m"]))
        ax2.tick_params(which='major', length=7)
        ax2.tick_params(which='minor', length=5)
        ax2.set_xticks(())
        for label in ax2.get_ymajorticklabels():
            label.set_fontsize(20)
        ax2.set_xticks([40+5*i for i in range(10)])
        ax2.set_xticklabels(['' for i in range(10)])

        ax3 = fig.add_axes([.91, .06, .017, .68])
        cbar = plt.colorbar(im1, cax=ax3, extend='max')
        cbar.ax.tick_params(labelsize=20)

        ax4 = fig.add_axes([.91, .78, .017, .19])
        cbar = plt.colorbar(im2, cax=ax4, extend='max', ticks=([1] + [1e5*i for i in range(1, 9)]))
        cbar.ax.set_yticklabels([1, '', '$2 \cdot 10^5$', '', '$4 \cdot 10^5$', '', '$6 \cdot 10^5$', '', '$8 \cdot 10^5$'])
        cbar.ax.tick_params(labelsize=20)

        if show:
            plt.show()
        if fname:
            plt.savefig(fname)

    # title = 'Age (days)'
    # plot(map_mean_age, title, cmap_name='plasma_r', log=False, vmin=0, vmax=3*365, show=False, fname=basename+'_mean_age.png')
    # title = 'Age (days)'
    # plot(map_mean_age, title, cmap_name='plasma_r', log=False, vmin=0, vmax=365, show=False, fname=basename+'_mean_age_short.png')

    title = '# Part / km$^2$'
    plot(map_density_4th_year, title, cmap_name='hot_r', log=True, vmin=1e-10, vmax=1e-5, under='white', over='blue', show=False, fname=basename+'_4th_yr_density_log.png')
    plot(map_density_4th_year, title, cmap_name='hot_r', log=False, vmin=1e-7, vmax=1e-5, under='white', over='blue', show=False, fname=basename+'_4th_yr_density_lin.png')

    title = '%'
    plot(map_touched, title, cmap_name='Spectral_r', log=True, vmin=.1, vmax=100, under='white', show=False, fname=basename+'_touched_log.png')
    plot(map_touched, title, cmap_name='Spectral_r', log=False, vmin=.1, vmax=100, under='white', show=False, fname=basename+'_touched_lin.png')

    if part_3d:
        plot_vertProfile(dens, dens_lat, dens_depth, show=False, fname=basename+'_vertProfile.png')
