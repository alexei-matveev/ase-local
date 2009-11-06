# -*- coding: utf-8 -*-
import os
# import numpy as np
# from numpy.linalg import eigh, solve

from ase.optimize import Optimizer
# from ase.data import atomic_numbers
from ase.gxfile import gxread, gxwrite, is_dummy
# conversion factors to Anstrom and eV
# for usage in the write and read routines for the gxfile
EINEV=27.2113845
LINA=0.529177

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
        self._loop = 0
        self._converged = False

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
	atnums, positions, isyms, inums, iconns, ivars, grads, energy = gxread('gxfile', LINA, EINEV )

	# use current positions as returned by the framework,
	# in case the on-disk version is outdated (units?):
        positions = self.atoms.get_positions()

	# the energy corresponding to the current geometry (units?):
	energy = self.atoms.get_potential_energy()

	# for use in gxfile:
	self._loop += 1

	print "GxOptimizer: loop=", self._loop
	print "GxOptimizer: energy=", energy
	print "GxOptimizer: positions=\n", positions
	print "GxOptimizer: forces=\n", forces

        # write gxfile to disk, note that energy gradients == - forces (units?):
        gxwrite(atnums, positions, isyms, inums, iconns, ivars, -forces, energy, file='gxfile', loop=self._loop, lunit=LINA, eunit=EINEV)

        # run external executable to update geometry:
	# exitcode = os.system('optimizer.exe')
        tty = os.popen("optimizer.exe","r")
	for line in tty:
	    print "PGOptimizer: ", line.rstrip("\n")

	# the last line of the output should tell if optimizer thinks it is converged:
	line = line.rstrip("\n")

	if line.startswith(" optimizer_main: converged="):
	  # in regular runs the last line should tell if convergence was reached:
	  self._converged = line.endswith("T")
	else:
	  # on errors the last line could be anything, abort by faking convergence:
	  print "GxOptimizer: ERROR! UNEXPECTED OUTPUT FROM EXTRENAL OPTIMIZER!"
	  self._converged = True

	print "GxOptimizer: converged=", self._converged

        # read the updated geometry from the gxfile:
	atnums, positions1, isyms, inums, iconns, ivars, grads, energy = gxread('gxfile', LINA, EINEV )

	print "GxOptimizer: new positions=\n", positions1
	print "GxOptimizer: step=\n", positions1 - positions

        self.atoms.set_positions(positions1)

#   # this is the default class method of Optimizer (see __init__.py):
#   def converged(self, forces=None):
#       """Did the optimization converge?"""
#       if forces is None:
#           forces = self.atoms.get_forces()
#       return (forces**2).sum(axis=1).max() < self.fmax**2

    # we will use the return status of PG-optimizer:
    def converged(self, forces=None):
        """Did the optimization converge?"""
        return self._converged

    def update(self, r, f):
        raise NotImplementedError;

    def replay_trajectory(self, traj):
        raise NotImplementedError;

# You need to add "set modeline" and eventually "set modelines=5"
# to your ~/.vimrc for this to take effect.
# Dont (accidentally) delete these lines! Unless you do it intentionally ...
# Default options for vim:sw=4:expandtab:smarttab:autoindent:syntax
