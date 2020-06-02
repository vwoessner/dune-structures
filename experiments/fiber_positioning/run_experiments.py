#!/usr/bin/env python

#
# Script that runs a batch of experiments for the positioning
# of a fibre in a cantilever.
#

import os
import subprocess
import sys
import matplotlib
import numpy as np

matplotlib.use("PDF")

from matplotlib import pyplot as plt

# Extract data from command line arguments
ini = sys.argv[1]
cells = int(sys.argv[2])
left = float(sys.argv[3])
right = float(sys.argv[4])
samples = int(sys.argv[5])
stab = sys.argv[6]

# The list of data points
h = (right - left) / samples
positions = list(left + 0.5 * h + h * i for i in range(samples))

filename = "positioning_cells{}_l{}_r{}_n{}_beta{}.log".format(cells, left, right, samples, stab)
if not os.path.isfile(filename):
    for pos in positions:
        cmd = ['../../apps/universal/universalapp', ini]
        cmd.append('-grid.fibres.first.start')
        cmd.append('-0.5 {}'.format(pos))
        cmd.append('-grid.fibres.first.end')
        cmd.append('4.5 {}'.format(pos))
        cmd.append('-fibre_operator.stabilization_parameter')
        cmd.append(stab)
        cmd.append('-filelogger.filename')
        cmd.append(filename)
        cmd.append('-grid.N')
        cmd.append('{} {}'.format(4*cells, cells))

        subprocess.call(cmd)

displacements = []
with open(filename, 'r') as f:
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
#ticks = np.arange(0.0, 1.1, 0.1)
ax = plt.axes()
#ax.set_xticks(ticks)
ax.xaxis.grid()

ax.set_xlabel("Beam positioning in y-direction")
ax.set_ylabel("y-Displacement")

plt.savefig("{}.pdf".format(filename))
