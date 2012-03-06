from __future__ import with_statement
"""This module defines an ASE interface to ParaGauss.

"""
import os, sys
from os.path import basename
from os.path import isfile
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
                 optimizer = None,
                 copy_input = "always"
                 ):


        """
        Parameters
        ==========

        |input|

            name of the input  file wich contains all the informations
            ParaGauss needs

        |cmdline|

            Shell command  to start ParaGauss, it will  be executed in
            working directory.  A typical command line reads:

                runpg /users/alexei/exe/openmpi/mainscf_V3.1.4

        |silence|

            if  True (is  as default)  ParaGauss stdout  will go  to a
            separate file if False it would go to the normal stdout

        |optimizer|

            If optimizer input is  needed for a ParaGauss single point
            calculation the programm  takes the content from optimizer
            and provides  it as  optimizer.input in the  directory the
            calculation runs

        |copy_input|

            Allows three different modes:

            always

                (is  the  default) will  create  new  input file  from
                storage  each  time  a quantum  chemistry  calculation
                starts

            never

                will never create an input file

            inexistent

                will create  a new input file for  a quantum chemistry
                calculation if it finds that the file does not exist

                Both always  and inexistent will fetch  the input file
                they  will  create  lateron  in  the  current  working
                directory during initalization
        """

        self.input = input
        self.cmdline = cmdline
        self.silence = silence
        assert (copy_input in ["always", "never", "inexistent"])
        self.copy_input = copy_input

        self.converged = False

        # store  metadata here, it  might be  needed (even  in another
        # directory)
        self.data = {}

        if not self.copy_input == "never":
            with open(self.input, "r") as file:
                self.inputstring = file.read()

        if optimizer == None:
            self.optimizer = None
        else:
            with open(optimizer, "r") as file:
                self.optimizer = file.read()
        # print self.inputstring

        self.atnums = None

        # there may be  a gxfile from gxoptimizer we  must not disturb
        # its internal coordinates
        if os.path.exists('gxfile'):
            self.atnums, __, self.data["isyms"], self.data["inums"], self.data["iconns"], self.data["ivars"], \
                self.data["additional"], __, __, loop = gxread('gxfile')


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
            print >> sys.stderr, "ERROR: (ParaGauss) no energy available"
            print >> sys.stderr, "Aborting."
            raise Exception("ParaGauss: no energy available")

        return self.__energy * Hartree

    def get_forces(self, atoms):
        """
        same as get_potential_energy but for forces
        units are transformed
        """

        self.update(atoms)

        if self.__grads == None :
            print >> sys.stderr, "ERROR: (ParaGauss) no forces available!"
            print >> sys.stderr, "Try enabling geometry optimization, by setting OPERATIONS_GEO_OPT = true"
            print >> sys.stderr, "and setting MAX_GEO_ITERATION = 0"
            print >> sys.stderr, "Aborting."
            raise Exception("ParaGauss: no forces available")

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
        atnums = atoms.get_atomic_numbers().copy()
        if (self.atnums == None):
            self.atnums = atnums
        elif (atnums != self.atnums).any() :
            print >> sys.stderr, "ERROR: (ParaGauss) gxfile does not fit!"
            print >> sys.stderr, "ERROR: (ParaGauss) gxfile contains wrong atoms!"
            print >> sys.stderr, "Please delete or change it before restart"
            raise Exception("gxfile does not fit, delete or adjust!")

        n = len(self.atnums)
        loop = 1
        # there may be a gxfile from another source
        # make sure it contains the same meta data than our source:
        t_gx = {}
        if os.path.exists('gxfile'):
            atnums, __, t_gx["isyms"], t_gx["inums"], t_gx["iconns"], t_gx["ivars"], t_gx["additional"], __, __, loop = gxread('gxfile')
            for dat in self.data.keys():
                if (np2.asarray(self.data[dat]) != np2.array(t_gx[dat])).any():
                    print >> sys.stderr, "ERROR: (ParaGauss) gxfile does not fit!"
                    print >> sys.stderr, "ERROR: (ParaGauss) gxfile contains wrong " + dat +" !"
                    print >> sys.stderr, "Please delete or change it before restart"
                    raise Exception("gxfile does not fit, delete or adjust!")

            if (np2.array(atnums) != self.atnums).any():
                print >> sys.stderr, "ERROR: (ParaGauss) gxfile does not fit!"
                print >> sys.stderr, "ERROR: (ParaGauss) gxfile contains wrong atoms!"
                print >> sys.stderr, "Please delete or change it before restart"
                raise Exception("gxfile does not fit, delete or adjust!")

        # Needs not to know size of system at init, but soon they will be needed
        if "isyms" not in self.data:
            if "isyms" in t_gx:
                self.data.update(t_gx)
            else:
                def dummy_or_not(at):
                    if at == 0:
                        return 0
                    else:
                        return 1

                self.data["isyms"] = np2.array([dummy_or_not(at) for at in atnums])
                self.data["inums"] = np2.array(range(1,n+1))
                self.data["iconns"] = np2.zeros((n,3))
                self.data["ivars"] = np2.zeros((n,3))
                self.data["additional"] = None
        # create gxfile with actual geometry for calculation
        # units of positions should be Bohrs in here, so they are changed
        gxwrite(self.atnums, self.positions/Bohr, self.data["isyms"], self.data["inums"], self.data["iconns"],\
                     self.data["ivars"], self.data["additional"], None, None, loop, file='gxfile' )
        input = basename(self.input)

        copy_inp = (self.copy_input == "always")

        if self.copy_input == "inexistent":
            copy_inp = (not isfile(input))

        if copy_inp:
            inputfile = open(input, "w")
            inputfile.write(self.inputstring)
            inputfile.close()

        if not self.optimizer == None:
            optifile = open("optimizer.input", "w")
            optifile.write(self.optimizer)
            optifile.close()

        # the actual calcualtion
        cmd = self.cmdline + ' ' + input
        if self.silence:
            cmd +=  ' > ParaGauss.out'
        tty = os.system(cmd)
        # reads in new energy and forces
        self.read()

        self.report(atoms, "ParaGauss.xyz")

        self.converged = True


    def report(self, atoms, file):
        #
        # Report the energy (and the geometry currently calculated on)
        # after a finshed calculation in ASE units
        #
        symbols = atoms.get_chemical_symbols()
        natoms = len(symbols)
        f = open(file, "w")
        f.write('%d\nE = %22.15f eV\n' % (natoms, self.__energy * Hartree))
        for s, (x, y, z) in zip(symbols, atoms.get_positions()):
            f.write('%-2s %22.15f %22.15f %22.15f\n' % (s, x, y, z))
        f.close()


    def read(self):
        # the interisting part to read in are the grads and energy, rest will be ignored afterwards
        # FIXME: Somehow we need (sometimes) to pass some time here, before we can find the
        # gxfile as output. This is especially valid when running several calculations in parallel.
        #It's done this way and not with sleep as the problems seems to be related
        # to a file that should but isn't there. So it makes sense to read all the files that are there,
        # even if we don't need the output.
        os.system("ls > /dev/null")

        if isfile('o.' + basename(self.input) + '/trace_output'):
            # If there are some trace files keep content of them, as if several calcualtions
            # have been performed after another, only for the last iteration the trace would
            # be available
            f = open('o.' + basename(self.input) + '/trace_output', "r")
            keep = f.read()
            f.close()
            f = open("keep_traces", "a")
            f.write(keep)
            f.close()

        if os.path.exists('gxfile'):
            __, __, __, __, __, __,__, self.__grads, self.__energy, loopi_d = gxread('gxfile')
            if self.__energy is not None:
                return
        else:
            print "ParaGauss ERROR: Found no gxfile to read energy or forces from"
            print "There should be at least the one I created"
            print "Therefore something very strange happened"
            print "ERROR: I quit!!"
            sys.exit(1)

        self.__energy = self.parse_output('o.' + basename(self.input) + '/output')

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

