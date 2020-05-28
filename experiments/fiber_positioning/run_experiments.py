#!/usr/bin/env python

#
# Script that runs the experiment and outputs a plot. If the results should
# be recomputed (takes several minutes), you should remove the positioning.log
# file before running this script.
#

import matplotlib
import numpy as np
import os
import subprocess

matplotlib.use("PDF")

from matplotlib import pyplot as plt

# The list of data points
positions = list(0.001 + 0.002 * i for i in range(500))

if not os.path.isfile("positioning.log"): 
    for pos in positions:
        cmd = ['../../apps/universal/universalapp', 'hansbo_beam.ini']
        cmd.append('-grid.fibres.first.start')
        cmd.append('-0.5 {}'.format(pos))
        cmd.append('-grid.fibres.first.end')
        cmd.append('4.5 {}'.format(pos))
        
        subprocess.call(cmd)

displacements = []
with open("positioning.log", 'r') as f:
    for line in f:
        _, ydis = line.split()
        ydis = float(ydis)
        displacements.append(ydis)

# Do some truncating
min_val = -2.0
max_val =  2.0
displacements = [max(min(d, max_val), min_val) for d in displacements]

# Plot
plt.plot(positions, displacements)

# Add vertical lines at vertex positionings
ticks = np.arange(0.0, 1.1, 0.1)
ax = plt.axes()
ax.set_xticks(ticks)
ax.xaxis.grid()

plt.savefig("positioning.pdf")
