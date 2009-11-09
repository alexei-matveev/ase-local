"""This module defines an ASE interface to ParaGauss.

"""
import os

import numpy as np2
from ase.gxfile import gxread, gxwrite


class ParaGauss:
    """Class for doing ParaGauss calculations.

    ParaGauss needs at least one input file, which gives to
    it all the needed parameters

    """
    def __init__(self,
          input = "input",
          writeowninput = False,
          runpg = "runpg",
          executable = "mainscf_V3.1.4-32",
          silence = True
          ):


      """
        Parameters
        ==========
      input: name of the input file wich contains all the informations
             ParaGauss needs
      writeowninput: if true the Programm does not build its own input
                     file, but reads in the one named as input parameter
                     needs this file to exist
                     so far the option false is not implemented
      runpg:
      executable: ParaGauss version, needs to be in working directory
                  or with Path in front
      """
      self.input = input
      self.writeowninput = writeowninput
      self.runpg = runpg
      self.executable = executable
      self.silence = silence

      self.converged = False

    def update(self, atoms):
        """
        decides whether and how to calculate energy and forces
        if the stored positions have not changed nothing is done
        if start or change in the settings, an initialization will
        take place
        """
        if (not self.converged or
            len(self.atnums) != len(atoms) or
            (self.atnums != atoms.get_atomic_numbers()).any()):
           self.initialize(atoms)
           self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any()):
           self.calculate(atoms)

    def initialize(self, atoms):

        self.converged = False

    def get_potential_energy(self, atoms, force_consistent=False):
        """
        makes sure energy (and forces) are up to date and afterwards
        gives energy back (energy is energy calculated with ParaGauss
        in atomic units (Hartree))
        """
        self.update(atoms)
        if self.__energy == None:
            print "WARNING: no energy available when energy wanted"
            print "There seems to be gone something wrong with the energy calculation"
            print "So I better stop here"
            sys.exit()
        return self.__energy

    def get_forces(self, atoms):
        """
        same as get_potential_energy but for forces
        units are Hartree/Bohrs
        """
        self.update(atoms)

        if self.__grads == None :
            print "WARNING: no forces available when they are wanted"
            print "There seems to be gone something wrong with the force calculation"
            print "So I better stop here"
            sys.exit()
        # note that the forces are negative of the energy gradients:
        return -self.__grads

    def get_stress(self, atoms):
            raise NotImplementedError


    def calculate(self, atoms):
        """
         calculates the energy and forces with
         ParaGauss
         uses gxfile to comunicate with the rest of the system
        """
        # read in actual positions and atomic numbers
        self.positions = atoms.get_positions().copy()
        self.atnums = atoms.get_atomic_numbers().copy()
        n = len(self.atnums)
        # there may be a gxfile from gxoptimizer
        # we must not disturb its internal coordinates
        if os.path.exists('gxfile'):
            atnums_d, xyz_d, self.isyms, inums, iconns, ivars, grads_dummy, energy_dummy = gxread('gxfile')
            if (atnums_d != self.atnums).any() :
                print "WARNING: gxfile does not fit!"
                print "Please delete or change it before restart"
                sys.exit()
        else:
            self.isyms = np2.ones(n)
            inums = np2.zeros(n)
            iconns = np2.zeros((n,3))
            ivars = np2.zeros((n,3))

        # create gxfile with actual geometry for calculation
        gxwrite(self.atnums, self.positions, self.isyms, inums, iconns, ivars, None, None, loop=1, file='gxfile' )
        # the actual calcualtion
        exetfile = self.runpg + ' ' + self.executable + ' ' + self.input
        if self.silence:
            exetfile+ =  ' >> ParaGauss.out'
        tty = os.system(exetfile)
        # reads in new energy and forces
        self.read()

        self.converged = True


    def read(self):
       # the interisting part to read in are the grads and energy, rest will be ignored afterwards
       atnums_d, xyz_d, self.isyms, inums, iconns, ivars, self.__grads, self.__energy = gxread('gxfile')

