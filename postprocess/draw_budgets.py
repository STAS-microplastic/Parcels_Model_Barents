#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on December 1 2018

@author: Philippe Delandmeter

Function plotting the different MP budgets per zone
"""


import numpy as np
import matplotlib.pyplot as plt

zones = {}
zones['0083'] = np.load('nemo0083_zone_concentration.npy')
zones['0083_cmems'] = np.load('nemo0083_cmems_zone_concentration.npy')
zones['025'] = np.load('nemo025_zone_concentration.npy')
zones['0083_stokes'] = np.load('nemo0083_stokes_zone_concentration.npy')
zones['0083_diff'] = np.load('nemo0083_diff_zone_concentration.npy')
zones['0083_3d'] = np.load('nemo0083_3d_zone_concentration.npy')


def plot_run(data, name, show_legend=False):
    fig = plt.figure(figsize=(14, 8.5), dpi=300, facecolor='w', edgecolor='k')
    ax = fig.add_axes([.1, .1, .8, .8])
    ax.axis([0, 600, 0, 100.5])

    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    j = [1, 0, 2, 3, 4, 5]
    color = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
    # linestyles https://matplotlib.org/gallery/lines_bars_and_markers/linestyles.html
    linestyles = ['-', '--', ':', '-.', (0, (3, 1, 1, 1)), (0, (3, 5, 1, 5, 1, 5))]
    for i in range(6):
        ax.plot(100*data[j[i]], linestyle=linestyles[i], linewidth=4, color=color[i])

    plt.xticks([365/4. * i for i in range(7)],
               ("0", ".5", "1", "1.5", "2", "2.5", "3"))
    plt.yticks([20 * i for i in range(6)],
               ("0", "20", "40", "60", "80", "100%"))
    for label in ax.get_xmajorticklabels():
        label.set_fontsize(26)
    for label in ax.get_ymajorticklabels():
        label.set_fontsize(26)

    ax.set_xlabel('years', fontsize=26)
    ax.xaxis.set_label_coords(.98, -0.01)

    zone_names = ['North Sea', 'Danish Strait',
                  'Norwegian Sea', 'Arctic', 'Atlantic', 'Baltic Sea']

    if show_legend:
        legend = ax.legend(zone_names, fontsize=30, ncol=2)
        legend.get_frame().set_linewidth(0.0)
    fig.savefig('budget_%s.png' % name)


for name, data in zones.items():
    plot_run(data, name)
