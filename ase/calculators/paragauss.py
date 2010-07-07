"""This module defines an ASE interface to ParaGauss.

"""
import os, sys
from os.path import basename
import numpy as np2
from ase.gxfile import gxread, gxwrite
from ase.units import Bohr, Hartree


class ParaGauss:
    """Class for doing ParaGauss calculations.

    ParaGauss needs at least one input file, which gives to
    it all the needed parameters

    """
    def __init__(self,
          input = "input",
          cmdline = "runpg /users/alexei/exe/openmpi/mainscf_V3.1.4b7-64",
          silence = True,
          ):


      """
      Parameters
      ==========
      |input|         name of the input file wich contains all the informations
                      ParaGauss needs

      |cmdline|       Shell command to start ParaGauss, it will be executed in working directory.
                      A typical command line reads:

                      runpg /users/alexei/exe/openmpi/mainscf_V3.1.4

      |writeowninput| if true the Programm does not build its own input
                      file, but reads in the one named as input parameter
                      needs this file to exist
                      so far the option false is not implemented
      """
      self.input = input
      self.cmdline = cmdline
      self.silence = silence

      self.converged = False

      file = open(self.input, "r")

      self.inputstring = file.read()
      file.close()
      #print self.inputstring

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
        in atomic units (energy transformed to ASE units))
        """
        self.update(atoms)
        if self.__energy == None:
            print "ParaGauss WARNING: no energy available when energy wanted"
            print "There seems to be gone something wrong with the energy calculation"
            print "So I better stop here"
            sys.exit()
        return self.__energy * Hartree

    def get_forces(self, atoms):
        """
        same as get_potential_energy but for forces
        units are transformed
        """
        self.update(atoms)

        if self.__grads == None :
            print "ParaGauss WARNING: no forces available when they are wanted"
            print "There seems to be gone something wrong with the force calculation"
            print "So I better stop here"
            sys.exit()
        # note that the forces are negative of the energy gradients:
        return -self.__grads * Hartree / Bohr

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
        loop = 1
        # there may be a gxfile from gxoptimizer
        # we must not disturb its internal coordinates
        if os.path.exists('gxfile'):
            atnums_d, xyz_d, self.isyms, inums, iconns, ivars, grads_dummy, energy_dummy, loop = gxread('gxfile')
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
        # units of positions should be Bohrs in here, so they are changed
        gxwrite(self.atnums, self.positions/Bohr, self.isyms, inums, iconns, ivars, None, None, loop, file='gxfile' )
        input = basename(self.input)
        inputfile = open(input, "w")
        inputfile.write(self.inputstring)
        inputfile.close()
        # the actual calcualtion
        cmd = self.cmdline + ' ' + input
        if self.silence:
            cmd +=  ' >> ParaGauss.out'
        tty = os.system(cmd)
        # reads in new energy and forces
        self.read()

        self.converged = True


    def read(self):
        # the interisting part to read in are the grads and energy, rest will be ignored afterwards
        # FIXME: Somehow we need (sometimes) to pass some time here, before we can find the
        # gxfile as output. This is especially valid when running several calculations in parallel.
        #It's done this way and not with sleep as the problems seems to be related
        # to a file that should but isn't there. So it makes sense to read all the files that are there,
        # even if we don't need the output.
        os.system("ls > /dev/null")
        if os.path.exists('gxfile'):
            atnums_d, xyz_d, self.isyms, inums, iconns, ivars, self.__grads, self.__energy, loopi_d = gxread('gxfile')
            if self.__energy is not None:
                return
        else:
            print "ParaGauss ERROR: Found no gxfile to read energy or forces from"
            print "There should be at least the one I created"
            print "Therefore something very strange happened"
            print "ERROR: I quit!!"
            sys.exit(1)

        self.__energy = self.parse_output('o.' + self.input + '/output')

    def parse_output(self, output): # not actually using |self|
        """
        Currently only returns the SCF energy.

        Energy lines in SCF section look like this:

        ^  e_sum  =            -1.521590696368  [      0.000000000000]
        """

        import re

        pattern = re.compile(r'\s*e_sum\s*=\s*(\S+)')

        # in case we dont find anything:
        e_sum = None

        lines = open(output,'r')
        for line in lines:
            match = pattern.search(line)
            if match is not None:
                # print line
                e_sum = float(match.group(1))
        # print 'e_sum=',e_sum
        # FIXME: should we somehow close() the file? It happens that close() is not in scope.
        return e_sum

