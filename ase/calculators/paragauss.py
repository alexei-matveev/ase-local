from __future__ import with_statement
"""This module defines an ASE interface to ParaGauss.

"""
import os, sys
from os.path import basename
from os.path import isfile
import numpy as np2
from numpy import array, any
from ase.gxfile import gxread, gxwrite
from ase.units import Bohr, Hartree

from general import Calculator

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

        # I am  getting tired of  this voodoo, FIXME: factor  this out
        # into a function:
        if input.startswith("i."):
            # e.g. i.h2o:
            input_base_name = input[2:]
        elif input.endswith(".scm") or input.endswith(".nml"):
            input_base_name = input[:-4]
        else:
            # e.g. input, or anything else:
            input_base_name = input

        if input_base_name == "input":
            self.output = "output"
        else:
            self.output = "o." + input_base_name + "/output"

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


    def calculate (self, atoms):
        """
         Calculates the energy and forces with ParaGauss.  Uses gxfile
         to comunicate with the rest of the system.
        """
        # read in actual positions and atomic numbers
        self.positions = atoms.get_positions().copy()
        atnums = atoms.get_atomic_numbers().copy()
        if (self.atnums == None):
            self.atnums = atnums

        if len (atnums) != len (self.atnums) or any (array (atnums) != array (self.atnums)):
            print >> sys.stderr, """
            ERROR: (ParaGauss) gxfile does not fit!  Gxfile contains
            wrong atoms!  Please delete or change it before restart.
            """
            raise Exception ("gxfile does not fit, delete or adjust!")

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

        self.__energy = self.parse_output(self.output)

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


class PG_nml():
  def __init__( self, nml_name, nml_keys={}, nml_data=[] ):
    # A PG namelist consists of
    # 1) a title (string)
    # 2) some key-value pairs (dictionary)
    # 3) an optional data-appendix (list, arbitrary type)
    self.name = nml_name
    self.keys = nml_keys
    self.data = nml_data
  def write( self ):
    # Leave whitespace for better readability
    inputtext =  ['']
    # Namelist header
    inputtext += [ '  &' + self.name ]
    # Namelist entries
    for key, val in self.keys.iteritems():
      # Should allow every type
      inputtext += ['    ' + key + ' = '+''.join(str([val])).replace('[','').replace(',','').replace(']','').replace('\'','') ]
    # Conclude Namelist
    inputtext += [ '  /' + self.name ]
    # Data entries
    for line in self.data:
      inputtext += ['    '+''.join(str([line])).replace('[','').replace(',','').replace(']','') ]
    return inputtext
  #
class PG_annotation():
  def __init__( self, nml_text='#' ):
    # A PG annotiation consists of
    # 1) a simple text (string)
    self.text = nml_text
  def write( self ):
    return [self.text]
  #
class PG(Calculator):
  #
  # Alternative calculator, in a more ASE-like style to set up calculations quickly
  #
  def __init__( self
              , exe = "/home/soini/PROGRAMMING/paragauss_mac_working/bin/runpg /home/soini/PROGRAMMING/paragauss_mac_working/mainscf_4test_suite"
              ,**kwargs ):
    #
    #
    self.__cmd__ = exe
    self.flag_keys   = {'uks'        : False
		       ,'jexact'     : False
		       ,'saveread_ks': True
		       ,'timers'     : False
		       }
    self.int_keys    = {'max_scf'    : 20
		       ,'mix_beg'    : 5
		       ,'ndiis'      : 5
		       ,'nrad'       : 150
		       ,'nang'       : 291
		       }
    self.real_keys   = {'e_conv'     : 1.0e-8
		       ,'d_conv'     : 1.0e-6
		       ,'scale_crit' : 1.0
		       ,'mix_fix'    : 0.25
		       ,'smear_val'  : 0.0
		       }
    self.str_keys    = {'task'       : '"Gradients"'
                       ,'sym'        : '"C1"'
                       ,'rel'        : '"FALSE"'
                       ,'xc'         : '"PBE"'
                       ,'mix_scf'    : '"diis"'
                       ,'smear'      : '"FALSE"'
		       }
    self.list_keys   = {'ea'         : [1]
		       ,'basis'      : {}
		       }
    #
    for key, val in kwargs.iteritems():
      if self.flag_keys.has_key( key ):
        self.flag_keys[key] = val
      elif self.int_keys.has_key( key ):
        self.int_keys[key] = val
      elif self.real_keys.has_key( key ):
        self.real_keys[key] = val
      elif self.str_keys.has_key( key ):
        self.str_keys[key] = '"' + val + '"'
      elif self.list_keys.has_key( key ):
        self.list_keys[key] = val
      else:
        self.__pg_error__( 'Illegal keyword '+key )
    #
    self.folder = 'o.input/'
    #
    self.atoms = None
    self.mixing = self.__check__( self.str_keys['mix_scf'].lower(), ['"diis"', '"chargefit"'], 'mix_scf' )
    self.smear = self.__check__( self.str_keys['smear'].lower(), ['"false"', '"gauss"', '"fermi"', '"sinus"'], 'smear' ).replace('"','')
    self.scale_crit = 1.0
    #
    self.__e_tot = None
    self.__forces = None
    #
    os.system('rm -rf trace.log')
    #
  def set_atoms(self, atoms):
    if (atoms != self.atoms):
      self.__got_output = False
    self.atoms = atoms.copy()
    #
  def get_potential_energy(self, atoms):
    from ase.units import Hartree
    self.update(atoms)
    if self.__got_output:
      energies = []
      for line in self.read( 'e_sum' ):
        if '[' in line:
          energies += [float(line[2])]
      try:
        self.__e_tot = energies[-1]*Hartree
      except:
        self.__pg_error__( 'No energies found in output.\n            Check Paragauss.out!                ' )
    else:
      self.__pg_error__( 'Failed while trying to retrieve energy from output' )
    return self.__e_tot
    #
  def get_forces(self, atoms):
    from ase.units import Hartree, Bohr
    from numpy import zeros
    self.update(atoms)
    if self.__got_output:
      grads = zeros( (atoms.get_number_of_atoms(), 3) )
      i_run = 0
      print >> sys.stdout, " "
      print >> sys.stdout, " Retrieving Gradients from output: "
      for line in self.read( 'Equal Center:' ):
        grads[i_run,0] = line[2]
        grads[i_run,1] = line[3]
        grads[i_run,2] = line[4]
        print >> sys.stdout, grads[i_run,0:2]
        i_run += 1
      self.__forces = - grads * Hartree / Bohr
      print >> sys.stdout, " "
    else:
      self.__pg_error__( 'Failed while trying to retrieve forces from output' )
    return self.__forces
    #
  def get_stress(self, atoms):
    raise NotImplementedError
    #
  def update( self, atoms ):
    from os import path
    self.__scale_convergency_criteria__()
    if self.atoms != atoms or not self.__got_output:
      self.set_atoms( atoms )
      self.__write_input__( self.atoms )
      self.calculate()
      # We should have an output at this point
      if path.exists(self.folder+'output'):
        self.__got_output = True
      else:
        self.__pg_error__( 'No output in file output or o.input/output. Something went terribly wrong' )
    #
  def calculate( self ):
    from os import path
    os.system('rm -rf')
    cmd = self.__cmd__ + ' ' + 'input'
    cmd +=  ' > ParaGauss.out'
    tty = os.system(cmd)
    if path.exists( 'o.input' ):
      self.folder = 'o.input/'
    else:
      self.folder = ''
    if self.real_keys['scale_crit'] > 1.0:
      os.system('echo "      Automatically setting convergency criteria to:  e_conv = '+str(self.real_keys['e_conv']*self.scale_crit)+'    d_conv = '+str(self.real_keys['d_conv']*self.scale_crit)+' " >> trace.log')
    os.system('cat '+self.folder+'trace_output >> trace.log')
    #
  def read( self, arg ):
    resultlines = []
    for line in open(self.folder+'output'):
      if arg in line:
	resultlines += [line.split()]
    return resultlines
    #
  def __write_input__( self, atoms ):
    #
    self.p_ua, self.n_ua = self.check_sym( atoms )
    self.n_el, self.b_ua = self.check_bas( self.p_ua, atoms, self.list_keys['basis'] )
    #
    nmls = self.define_nmls( atoms )
    #
    input = open('input','w')
    for nml in nmls:
      for line in nml.write():
        input.write("%s\n" % line )
    input.close()
    #
  def __write_namelist__( self, namelist, entries ):
    inputtext = [ '', ' &'+namelist ]
    for key, val in entries.iteritems():
      inputtext.append( '   '+key+' = '+str(val) )
    inputtext.append( ' /'+namelist )
    return inputtext
    #
    #
  def __check__( self, entry, allowed=[True,False], entryname='abcdefg' ):
    if entry not in allowed:
      print 'ERROR: only'
      for x in allowed:
	print str(x)
      print '       allowed for '+str(entryname)
      sys.exit()
    return entry
    #
  def __scale_convergency_criteria__( self ):
    from ase.units import Hartree, Bohr
    if self.real_keys['scale_crit'] > 1.0:
      if self.__forces != None:
	# Scale according to forces. Upper bound is 0.1, lower bound is e_conv and d_conv
        self.scale_crit = min([ 0.1/max([self.real_keys['e_conv'],self.real_keys['d_conv']])
			      , max( [ 1.0, self.real_keys['scale_crit'] * self.__max_force__( self.__forces * Bohr / Hartree )] )] )
      else:
	# Scale assuming forces of 1.0
	self.scale_crit = self.real_keys['scale_crit']
	#
  def __max_force__( self, grads ):
    max_force = 0.0
    for grad in grads:
      max_force = max( [max_force, (grad[0]*grad[0]+grad[1]*grad[1]+grad[2]*grad[2])**0.5 ] )
    return max_force
    #
  def __pg_error__( self, msg ):
    import sys
    print
    print ' #### ERROR: '+msg+' ####'
    print
    sys.exit()
    #
  def __point_groups__( self ):
    return ['"C1"','"C2"','"C3"','"C4"','"C5"','"C6"','"C7"','"C8"','"C9"','"C10"','"Ci"'
	   ,'"CS"','"S4"','"S6"','"S8"','"S10"','"S12"','"S14"','"S16"','"S18"','"S20"'
           ,'"D2"','"D3"','"D4"','"D5"','"D6"','"D7"','"D8"','"D9"','"D10"'
	   ,'"D2H"','"D3H"','"D4H"','"D5H"','"D6H"','"D7H"','"D8H"','"D9H"','"D10H"','"Dinh"'
           ,'"D2D"','"D3D"','"D4D"','"D5D"','"D6D"','"D7D"','"D8D"','"D9D"','"D10D"'
           ,'"C2V"','"C3V"','"C4V"','"C5V"','"C6V"','"C7V"','"C8V"','"C9V"','"C10V"','"Cinv"'
           ,'"C2H"','"C3H"','"C4H"','"C5H"','"C6H"','"C7H"','"C8H"','"C9H"','"C10H"'
	   ,'"O"','"T"','"OH"','"TH"','"TD"','"I"','"IH"']
    #
  def check_sym( self, atoms ):
    #
    self.str_keys['sym'] = self.__check__( self.str_keys['sym'].upper()
                                         , allowed=self.__point_groups__()
                                         , entryname='sym' )
    #
    if self.str_keys['sym'] == '"C1"' or atoms.get_number_of_atoms() == 1 and sum(abs(atoms.get_positions())) == 0.0:
      self.list_keys['ea'] = [1]*atoms.get_number_of_atoms()
    #
    p_ua     = [None]*len( self.list_keys['ea'] )
    i_run         = 0
    for i_ua in range(len(p_ua)):
      p_ua[i_ua] = i_run
      i_run     += self.list_keys['ea'][i_ua]
    self.__check__( i_run, [atoms.get_number_of_atoms()], entryname='p_ua' )
    n_ua     = len( self.list_keys['ea'] )
    return p_ua, n_ua
    #
  def check_bas( self, p_ua, atoms, bas_raw ):
    from os import path
    # check given basis for consistency
    el = {}
    bl = []
    for i_run in p_ua:
      symbol = atoms.get_chemical_symbols()[i_run].lower()
      # count elements
      try:
        el[symbol] += 1
      except:
        el[symbol] = 1
      #
      if bas_raw.has_key( symbol ):
        basisfilename = bas_raw[symbol]
        if path.exists( basisfilename ):
          bl += [basisfilename]
        else:
          self.__pg_error__( 'Basis set file '+basisfilename+' missing' )
      else:
        self.__pg_error__( 'Need basis set for element '+symbol )
    return len( el ), bl
    #
  def define_nmls( self, atoms ):
    from time import clock
    from os import path
    # Define input as a list of PG_nml types.
    # Where applicable, pass input through __check__ filter
    #
    # Add header for section
    head1 = PG_annotation( '\n#'+('# define calculation #').center(80,'~')+'#' )
    #
    # namelist TASKS
    tasks = PG_nml( 'tasks',
      { 'task': self.__check__( self.str_keys['task'].lower(), ['"singlepoint"','"gradients"'], 'task' ) } )
    #
    # namelist MAIN_OPTIONS
    maino = PG_nml( 'main_options'
                  , { 'spin_restricted'     : not self.flag_keys['uks']
                    , 'relativistic'        : self.__check__( self.str_keys['rel'].lower(), ['"true"','"false"','"adkh"'], 'rel' )
                    , 'perturbation_theory' : self.__check__( self.str_keys['mix_scf'].lower(), ['"diis"', '"chargefit"'], 'mix_scf' ) == '"chargefit"'
                    # override internal ParaGauss defaults !!
                    , 'integrals_on_file'   : 'False # predefined by PG-calculator' } )
    #
    # NAMELIST OUTPUT TIMING
    timers = PG_nml( 'output_timing'
                   , { 'output_timing_summary'           : self.flag_keys['timers']
                     , 'output_timing_detailedsummary'   : self.flag_keys['timers']
                     , 'output_timing_integrals'         : self.flag_keys['timers']
                     , 'output_timing_detailedintegrals' : self.flag_keys['timers']
                     , 'output_timing_scfloops'          : self.flag_keys['timers']
                     , 'output_timing_scf'               : self.flag_keys['timers']
                     , 'output_timing_detailedscf'       : self.flag_keys['timers']
                     , 'output_timing_post_scf'          : self.flag_keys['timers']
                     , 'output_timing_detailedpostscf'   : self.flag_keys['timers']
                     , 'output_timing_slaves'            : self.flag_keys['timers']
                     , 'output_timing_interrupts'        : self.flag_keys['timers'] } )
    #
    # NAMELIST RECOVER_OPTIONS
    recoo = PG_nml( 'recover_options'
                  , { 'save_ksmatrix' : self.flag_keys['saveread_ks']
                    , 'read_ksmatrix' : self.flag_keys['saveread_ks'] and path.exists('saved_ksmatrix.dat') } )
    #
    if self.mixing == '"diis"':
      # NAMELIST MIXING
      mixin = PG_nml( 'mixing'
                    , { 'chmix' : 1.0
                      , 'spmix' : 1.0
                      , 'xcmix' : 1.0
                      , 'start_after_cycle' : 1000000 } )
      #
      # NAMELIST DIIS
      diis  = PG_nml( 'diis'
                    , { 'diis_on'    : self.int_keys['ndiis'] > 0
                      , 'mmax'       : self.int_keys['ndiis']
                      , 'loop_start' : self.int_keys['mix_beg']
                      , 'threshold'  : 0.15
                      , 'cfix'       : self.real_keys['mix_fix'] } )
    else:
      mixin = PG_nml( 'mixing'
                    , { 'chmix' : self.real_keys['mix_fix']
                      , 'spmix' : 1.0
                      , 'xcmix' : 1.0
                      , 'start_after_cycle' : self.int_keys['mix_beg'] } )
      diis  = PG_nml( 'diis'   , { } )
    #
    # NAMELIST CONVERGENCE_LIST
    convl = PG_nml( 'convergence_list'
                  , { 'max_iteration': self.int_keys['max_scf']
                    , 'energy_criterion': self.real_keys['e_conv']*self.scale_crit
                    , 'density_criterion': self.real_keys['d_conv']*self.scale_crit
                    , 'energy_dev_checked': 3 } )
    #
    # NAMELIST XC_CONTROL
    xccnt = PG_nml( 'xc_control'
                  , { 'xc': self.str_keys['xc'] } )
    #
    # NAMELIST FERMI

    if self.real_keys['smear_val'] > 0.0 or self.smear != 'FALSE':
      smear = PG_nml( 'fermi'
                    , { 'fermi_'+self.smear: 'true'
                      , 'fermi_sigma': self.real_keys['smear_val'] } )
    #
    # NAMELIST ERI4C
    eri4c = PG_nml( 'eri4c'
                  , { 'j_exact': self.flag_keys['jexact'] } )
    #
    head2 = PG_annotation( '\n#'+('# define system #').center(80,'~')+'#' )
    #
    # NAMELIST SYMMETRY_GROUP
    symgr = PG_nml( 'symmetry_group'
                  , { 'point_group': self.str_keys['sym'] } )
    #
    # NAMELIST UNIQUE_ATOM_NUMBER
    uanum = PG_nml( 'unique_atom_number'
                  , { 'n_unique_atoms': self.n_ua } )
    #
    uanml = []
    for i_ua in range(self.n_ua):
      i_run = self.p_ua[i_ua]
      #
      # NAMELIST N_UNIQUE_ATOMS
      uanml += [ PG_nml( 'unique_atom # '+str(i_ua+1)
                       , { 'name'         : '"'+atoms.get_chemical_symbols()[i_run]+'"'
                         , 'z'            : str(atoms.get_atomic_numbers()[i_run])+'.0'
                         , 'n_equal_atoms': self.list_keys['ea'][i_ua] }
                       , [list(atoms.get_positions()[i_run]/Bohr)] ) ]
    #
    head3 = PG_annotation( '\n#'+('# define grid and basis #').center(80,'~')+'#' )
    #
    # NAMELIST GRID
    grid  = PG_nml( 'grid', { 'sym_reduce': True } )
    #
    ganml = []
    for i_ua in range(self.n_ua):
      # NAMELIST GRIDATOM
      ganml += [ PG_nml( 'gridatom # '+str(i_ua+1)
                       , { 'nrad': self.int_keys['nrad']
                         , 'nang': self.int_keys['nang'] } ) ]
    #
    blist = []
    for i_ua in range(self.n_ua):
      # GIVE BASIS SET FILE # no explicit namelists for now
      blist +=[ PG_annotation( '\n~'+self.b_ua[i_ua] ) ]
    #
    # FINAL LINE
    final = PG_annotation( '\n#'+('# compiled at '+str(clock())+' #').center(80,'~')+'#\n' )
    #
    return [ head1, tasks, maino, timers, recoo, mixin, diis, convl, smear, xccnt, eri4c
           , head2, symgr, uanum ]+uanml+[
             head3, grid]+ganml+blist+[ final ]

