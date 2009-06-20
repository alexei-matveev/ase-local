# -*- coding: utf-8 -*-
import os
# import numpy as np
# from numpy.linalg import eigh, solve

from ase.optimize import Optimizer
# from ase.data import atomic_numbers
from gxfile import gxread, gxwrite


class GxOptimizer(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=None):
        """An interface to geometry optimizer of ParaGauss.
        """
        
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

#       if maxstep is not None:
#           if maxstep > 1.0:
#               raise ValueError('You are using a much too large value for ' +
#                                'the maximum step size: %.1f Å' % maxstep)
#           self.maxstep = maxstep

    def initialize(self):
        #raise NotImplementedError;
        self.loop = 0

    def read(self):
        raise NotImplementedError;

#   def step(self, f):
#       atoms = self.atoms
#       r = atoms.get_positions()
#       f = f.reshape(-1)
#       self.update(r.flat, f)
#       omega, V = eigh(self.H)
#       dr = np.dot(V, np.dot(f, V) / np.fabs(omega)).reshape((-1, 3))
#       #dr = solve(self.H, f).reshape((-1, 3))
#       steplengths = (dr**2).sum(1)**0.5
#       dr /= np.maximum(steplengths / self.maxstep, 1.0).reshape(-1, 1)
#       print 'QN: step=',dr
#       atoms.set_positions(r + dr)
#       self.r0 = r.flat.copy()
#       self.f0 = f.copy()
#       self.dump((self.H, self.r0, self.f0, self.maxstep))

    def step(self, forces):
        """Called in a loop by .run() method of Optimizer class
	with forces computed for current geometry.
	"""

	# read the metadata from gxfile in columns:
	atnums, positions, isyms, inums, iconns, ivars, grads, energy = gxread('gxfile')

	# use current positions as returned by the framework,
	# in case the on-disk version is outdated (units?):
        positions = self.atoms.get_positions()

	# the energy corresponding to the current geometry (units?):
	energy = self.atoms.get_potential_energy()

	# for use in gxfile:
	self.loop += 1

	print "GxOptimizer: loop=\n", self.loop
	print "GxOptimizer: energy=\n", energy
	print "GxOptimizer: positions=\n", positions
	print "GxOptimizer: forces=\n", forces

        # write gxfile to disk, note that energy gradients == - forces (units?):
        gxwrite(atnums, positions, isyms, inums, iconns, ivars, -forces, energy, file='gxfile', loop=self.loop)

        # run external executable to update geometry:
        exitcode = os.system('optimizer.exe')

        # read the updated geometry from the gxfile:
	atnums, positions1, isyms, inums, iconns, ivars, grads, energy = gxread('gxfile')

	print "GxOptimizer: new positions=\n", positions1
	print "GxOptimizer: step=\n", positions1 - positions

        self.atoms.set_positions(positions1)

    def update(self, r, f):
        raise NotImplementedError;

    def replay_trajectory(self, traj):
        raise NotImplementedError;

# You need to add "set modeline" and eventually "set modelines=5"
# to your ~/.vimrc for this to take effect.
# Dont (accidentally) delete these lines! Unless you do it intentionally ...
# Default options for vim:sw=4:expandtab:smarttab:autoindent:syntax
