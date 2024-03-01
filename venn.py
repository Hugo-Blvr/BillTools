
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@original_author: tctianchi
@lien github : https://github.com/tctianchi/pyvenn

@edit by: Hugo
"""

from itertools import chain
try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors
import math

default_colors = [
    # r, g, b, a
    [92, 192, 98, 0.5],
    [90, 155, 212, 0.5],
    [246, 236, 86, 0.6],
    [241, 90, 96, 0.4],
    [255, 117, 0, 0.3],
    [82, 82, 190, 0.2],
]
default_colors = [
    [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
    for i in default_colors
]


def draw_ellipse(ax, x, y, w, h, a, fillcolor):
    e = patches.Ellipse(
        xy=(x, y),
        width=w,
        height=h,
        angle=a,
        color=fillcolor)
    ax.add_patch(e)


def draw_text(ax, x, y, text, fontsize=14, ha="center", va="center"):
    ax.text(
        x, y, text,
        horizontalalignment=ha,
        verticalalignment=va,
        fontsize=fontsize,
        color="black")


def venn2(labels, names, **options):
    colors = options.get('colors', [default_colors[i] for i in range(2)])
    figsize = options.get('figsize', (9, 7))
    dpi = options.get('dpi', 96)
    fontsize = options.get('fontsize', 14)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=0.7)
    ax.set_xlim(left=0.0, right=1.0)

    # body
    draw_ellipse(ax, 0.375, 0.3, 0.5, 0.5, 0.0, colors[0])
    draw_ellipse(ax, 0.625, 0.3, 0.5, 0.5, 0.0, colors[1])
    draw_text(ax, 0.74, 0.30, labels.get('10', ''), fontsize=fontsize)
    draw_text(ax, 0.26, 0.30, labels.get('01', ''), fontsize=fontsize)
    draw_text(ax, 0.50, 0.30, labels.get('11', ''), fontsize=fontsize)

    # legend
    draw_text(ax, 0.20, 0.56, names[0], fontsize=fontsize, ha="right", va="bottom")
    draw_text(ax, 0.80, 0.56, names[1], fontsize=fontsize, ha="left", va="bottom")
    leg = ax.legend(names, loc='center left', bbox_to_anchor=(0.98, 0.5), fancybox=True)
    leg.get_frame().set_alpha(0.5)

    return fig, ax


def venn3(labels, names, **options):
    colors = options.get('colors', [default_colors[i] for i in range(3)])
    figsize = options.get('figsize', (9, 9))
    dpi = options.get('dpi', 96)
    fontsize = options.get('fontsize', 14)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)

    # body
    draw_ellipse(ax, 0.333, 0.633, 0.5, 0.5, 0.0, colors[0])
    draw_ellipse(ax, 0.666, 0.633, 0.5, 0.5, 0.0, colors[1])
    draw_ellipse(ax, 0.500, 0.310, 0.5, 0.5, 0.0, colors[2])
    draw_text(ax, 0.50, 0.27, labels.get('100', ''), fontsize=fontsize)
    draw_text(ax, 0.73, 0.65, labels.get('010', ''), fontsize=fontsize)
    draw_text(ax, 0.61, 0.46, labels.get('110', ''), fontsize=fontsize)
    draw_text(ax, 0.27, 0.65, labels.get('001', ''), fontsize=fontsize)
    draw_text(ax, 0.39, 0.46, labels.get('101', ''), fontsize=fontsize)
    draw_text(ax, 0.50, 0.65, labels.get('011', ''), fontsize=fontsize)
    draw_text(ax, 0.50, 0.51, labels.get('111', ''), fontsize=fontsize)

    # legend
    draw_text(ax, 0.15, 0.87, names[0], fontsize=fontsize, ha="right", va="bottom")
    draw_text(ax, 0.85, 0.87, names[1], fontsize=fontsize, ha="left", va="bottom")
    draw_text(ax, 0.50, 0.02, names[2], fontsize=fontsize, va="top")
    leg = ax.legend(names, loc='center left', bbox_to_anchor=(0.98, 0.5), fancybox=True)
    leg.get_frame().set_alpha(0.5)

    return fig, ax


def venn4(labels, names, **options):
    colors = options.get('colors', [default_colors[i] for i in range(4)])
    figsize = options.get('figsize', (12, 12))
    dpi = options.get('dpi', 96)
    fontsize = options.get('fontsize', 14)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)

    # body
    draw_ellipse(ax, 0.350, 0.400, 0.72, 0.45, 140.0, colors[0])
    draw_ellipse(ax, 0.450, 0.500, 0.72, 0.45, 140.0, colors[1])
    draw_ellipse(ax, 0.544, 0.500, 0.72, 0.45, 40.0, colors[2])
    draw_ellipse(ax, 0.644, 0.400, 0.72, 0.45, 40.0, colors[3])
    draw_text(ax, 0.85, 0.42, labels.get('0100', ''), fontsize=fontsize)
    draw_text(ax, 0.68, 0.72, labels.get('0010', ''), fontsize=fontsize)
    draw_text(ax, 0.77, 0.59, labels.get('0110', ''), fontsize=fontsize)
    draw_text(ax, 0.32, 0.72, labels.get('0001', ''), fontsize=fontsize)
    draw_text(ax, 0.71, 0.30, labels.get('0101', ''), fontsize=fontsize)
    draw_text(ax, 0.50, 0.66, labels.get('0011', ''), fontsize=fontsize)
    draw_text(ax, 0.65, 0.50, labels.get('0111', ''), fontsize=fontsize)
    draw_text(ax, 0.14, 0.42, labels.get('1000', ''), fontsize=fontsize)
    draw_text(ax, 0.50, 0.17, labels.get('1100', ''), fontsize=fontsize)
    draw_text(ax, 0.29, 0.30, labels.get('1010', ''), fontsize=fontsize)
    draw_text(ax, 0.39, 0.24, labels.get('1110', ''), fontsize=fontsize)
    draw_text(ax, 0.23, 0.59, labels.get('1001', ''), fontsize=fontsize)
    draw_text(ax, 0.61, 0.24, labels.get('1101', ''), fontsize=fontsize)
    draw_text(ax, 0.35, 0.50, labels.get('1011', ''), fontsize=fontsize)
    draw_text(ax, 0.50, 0.38, labels.get('1111', ''), fontsize=fontsize)

    # legend
    i = names[0]
    j = names[1]
    k = names[2]
    names[0] = names[3]
    names[1] = i
    names[2] = j
    names[3] = k

    draw_text(ax, 0.13, 0.18, names[0], fontsize=fontsize, ha="right")
    draw_text(ax, 0.18, 0.83, names[1], fontsize=fontsize, ha="right", va="bottom")
    draw_text(ax, 0.82, 0.83, names[2], fontsize=fontsize, ha="left", va="bottom")
    draw_text(ax, 0.87, 0.18, names[3], fontsize=fontsize, ha="left", va="top")
    leg = ax.legend(names, loc='center left', bbox_to_anchor=(0.98, 0.5), fancybox=True)
    leg.get_frame().set_alpha(0.5)

    return fig, ax


def venn5(labels, names, **options):
    colors = options.get('colors', [default_colors[i] for i in range(5)])
    figsize = options.get('figsize', (13, 13))
    dpi = options.get('dpi', 96)
    fontsize = options.get('fontsize', 14)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)

    # body
    draw_ellipse(ax, 0.428, 0.449, 0.87, 0.50, 155.0, colors[0])
    draw_ellipse(ax, 0.469, 0.543, 0.87, 0.50, 82.0, colors[1])
    draw_ellipse(ax, 0.558, 0.523, 0.87, 0.50, 10.0, colors[2])
    draw_ellipse(ax, 0.578, 0.432, 0.87, 0.50, 118.0, colors[3])
    draw_ellipse(ax, 0.489, 0.383, 0.87, 0.50, 46.0, colors[4])
    draw_text(ax, 0.10, 0.61, labels.get('00010', ''), fontsize=fontsize)
    draw_text(ax, 0.18, 0.50, labels.get('00110', ''), fontsize=fontsize)
    draw_text(ax, 0.20, 0.31, labels.get('10010', ''), fontsize=fontsize)
    draw_text(ax, 0.21, 0.37, labels.get('10110', ''), fontsize=fontsize)
    draw_text(ax, 0.25, 0.58, labels.get('00111', ''), fontsize=fontsize)
    draw_text(ax, 0.27, 0.11, labels.get('10000', ''), fontsize=fontsize)
    draw_text(ax, 0.27, 0.70, labels.get('00011', ''), fontsize=fontsize)
    draw_text(ax, 0.28, 0.39, labels.get('10111', ''), fontsize=fontsize)
    draw_text(ax, 0.33, 0.72, labels.get('01011', ''), fontsize=fontsize)
    draw_text(ax, 0.34, 0.25, labels.get('10011', ''), fontsize=fontsize)
    draw_text(ax, 0.36, 0.66, labels.get('01111', ''), fontsize=fontsize)
    draw_text(ax, 0.39, 0.15, labels.get('10001', ''), fontsize=fontsize)
    draw_text(ax, 0.42, 0.78, labels.get('01001', ''), fontsize=fontsize)
    draw_text(ax, 0.50, 0.15, labels.get('11001', ''), fontsize=fontsize)
    draw_text(ax, 0.51, 0.22, labels.get('11011', ''), fontsize=fontsize)
    draw_text(ax, 0.51, 0.47, labels.get('11111', ''), fontsize=fontsize)
    draw_text(ax, 0.51, 0.74, labels.get('01101', ''), fontsize=fontsize)
    draw_text(ax, 0.51, 0.90, labels.get('00001', ''), fontsize=fontsize)
    draw_text(ax, 0.55, 0.13, labels.get('11000', ''), fontsize=fontsize)
    draw_text(ax, 0.64, 0.67, labels.get('11101', ''), fontsize=fontsize)
    draw_text(ax, 0.65, 0.23, labels.get('11010', ''), fontsize=fontsize)
    draw_text(ax, 0.67, 0.76, labels.get('00101', ''), fontsize=fontsize)
    draw_text(ax, 0.70, 0.71, labels.get('10101', ''), fontsize=fontsize)
    draw_text(ax, 0.72, 0.11, labels.get('01000', ''), fontsize=fontsize)
    draw_text(ax, 0.74, 0.40, labels.get('11110', ''), fontsize=fontsize)
    draw_text(ax, 0.76, 0.25, labels.get('01010', ''), fontsize=fontsize)
    draw_text(ax, 0.76, 0.55, labels.get('11100', ''), fontsize=fontsize)
    draw_text(ax, 0.78, 0.64, labels.get('10100', ''), fontsize=fontsize)
    draw_text(ax, 0.81, 0.37, labels.get('01110', ''), fontsize=fontsize)
    draw_text(ax, 0.84, 0.41, labels.get('01100', ''), fontsize=fontsize)
    draw_text(ax, 0.91, 0.58, labels.get('00100', ''), fontsize=fontsize)

    # legend
    i = names[0]
    names[0] = names[1]
    names[1] = i
    draw_text(ax, 0.02, 0.72, names[0], fontsize=fontsize, ha="right")
    draw_text(ax, 0.72, 0.94, names[1], fontsize=fontsize, va="bottom")
    draw_text(ax, 0.97, 0.74, names[2], fontsize=fontsize, ha="left")
    draw_text(ax, 0.88, 0.05, names[3], fontsize=fontsize, ha="left")
    draw_text(ax, 0.12, 0.05, names[4], fontsize=fontsize, ha="right")

    leg = ax.legend(names, loc='center left', bbox_to_anchor=(0.98, 0.5), fancybox=True)
    leg.get_frame().set_alpha(0.5)

    return fig, ax
