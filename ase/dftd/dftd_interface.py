#!/usr/bin/python
import numpy as np
from math import sqrt, exp, ceil
from decimal import Decimal
from ase.dftd.dftd_native import check_interaction_group_input
from ase.dftd.dftd_native import dft_d_pbc
from ase.dftd.dftd_native import maxdist
from ase.dftd.dftd_native import minbox

# General parameters
AU_TO_ANG    = 0.52917726
K1_PARAMETER = 16.0
DF_CUTOFF_RADIUS = 30.0
#
#
def d2_pbc(atoms, functional):
    """Main function making the DFT-D2 correction available for
    isolated systems.

    """
    # Check if compiled f2py interfaces are present
    try:
      from ase.dftd.dftd_module import d2_gradients as dftd2_gradients
    except ImportError:
      dftd_autocompile()
      # Check if compilation was successful
      try:
        from ase.dftd.dftd_module import d2_gradients as dftd2_gradients
      except ImportError:
        print ' Compilation of DFT-D module failed!'
        print
    #
    #
    dispersion_correction = 0.0
    gradient_contribution = 0.0
    #
    # Obtain data about system
    atom_numbers                       = atoms.get_atomic_numbers()
    positions                          = atoms.get_positions()
    periodic_directions                = atoms.get_pbc()
    elem_cell                          = atoms.get_cell()
    interactionlist, interactionmatrix = get_interaction_controls(atoms)
    #
    # Number of atoms within a single copy
    N_atoms = len(positions)
    #
    # Check input i.e. if interactionlist and interactionmatrix are set properly
    interactionlist, interactionmatrix = check_interaction_group_input(N_atoms, interactionlist, interactionmatrix)
    #
    #
    # Start with Calculation
    def func_get_dftd2(t_vec):
	#
        dispersion_correction, gradient_contribution = dftd2_gradients(atoms.get_atomic_numbers(), atoms.get_positions(), t_vec, interactionlist, interactionmatrix, functional)
	#
        return dispersion_correction, gradient_contribution
    #
    # Initialize (defines shapes)
    dispersion_correction0 = 0.0
    gradient_contribution0 = np.zeros((N_atoms,3))
    #
    dispersion_correction, gradient_contribution = lattice_sum(func_get_dftd2, dispersion_correction0, gradient_contribution0, positions, elem_cell, periodic_directions, DF_CUTOFF_RADIUS)
    #
    # return results in a.u.
    return dispersion_correction, gradient_contribution
    #
# End of function d2_pbc
#
#
def d3_pbc(atoms, functional):
    """Main function making the DFT-D3 correction available for
    isolated systems.

    """
    # Check if compiled f2py interfaces are present
    try:
      from ase.dftd.dftd_module import d3_gradients as dftd3_gradients
    except ImportError:
      dftd_autocompile()
      # Check if compilation was successful
      try:
        from ase.dftd.dftd_module import d3_gradients as dftd3_gradients
      except ImportError:
        print ' Compilation of DFT-D module failed!'
        print
    #
    #
    dispersion_correction = 0.0
    gradient_contribution = 0.0
    #
    # Obtain data about system
    atom_numbers                       = atoms.get_atomic_numbers()
    positions                          = atoms.get_positions()
    periodic_directions                = atoms.get_pbc()
    elem_cell                          = atoms.get_cell()
    interactionlist, interactionmatrix = get_interaction_controls(atoms)
    #
    # Number of atoms within a single copy
    N_atoms = len(positions)
    #
    # Check input i.e. if interactionlist and interactionmatrix are set properly
    interactionlist, interactionmatrix = check_interaction_group_input(N_atoms, interactionlist, interactionmatrix)
    #
    # Define function for lattice summation of coordination numbers
    def func_get_cn(t_vec):
        # tmp_func contains cn, dcn2, and dcn3
        cn, dcn = get_coordination_numbers(atoms, t_vec)
        return cn, dcn
    #
    # Initialize (defines shapes)
    cn0  = np.zeros(N_atoms)
    dcn0 = np.zeros((N_atoms,N_atoms,3))

    # Calculate the coordination numbers for atoms in the unit cell
    cn, dcn = lattice_sum(func_get_cn, cn0, dcn0, positions, elem_cell, periodic_directions, 8.)
    #
    # Transform coordination number results
    #
    # Start with Calculation
    def func_get_dftd3(t_vec):
	#
        dispersion_correction, gradient_contribution = dftd3_gradients(atoms.get_atomic_numbers(), atoms.get_positions(), t_vec, interactionlist, interactionmatrix, functional, cn, dcn)
	#
        return dispersion_correction, gradient_contribution
    #
    # Initialize (defines shapes)
    dispersion_correction0 = 0.0
    gradient_contribution0 = np.zeros((N_atoms,3))
    #
    dispersion_correction, gradient_contribution = lattice_sum(func_get_dftd3, dispersion_correction0, gradient_contribution0, positions, elem_cell, periodic_directions, DF_CUTOFF_RADIUS)
    #
    # return results in a.u.
    return dispersion_correction, gradient_contribution
    #
# End of function d3_pbc
#
#
#
# Outsourced from dftd_functions.f and gdisp.f
def get_coordination_numbers(atoms, t_vec = np.array([0., 0., 0.]), cutoff_radius = 8.0):
    #
    # Initialize
    xyz  = atoms.get_positions() / AU_TO_ANG
    tvec = t_vec / AU_TO_ANG
    cn   = np.zeros(len(atoms.get_atomic_numbers()))
    dcn = np.zeros((len(atoms.get_atomic_numbers()),len(atoms.get_atomic_numbers()),3))
    #
    #
    for ind_i in range(0, len(atoms.get_atomic_numbers())):
        xn = 0.0
        for ind_iat in range(0, len(atoms.get_atomic_numbers())):
             if ind_i == ind_iat:
	        if np.sum(t_vec) == 0:
	            continue
	     #
	     dx = xyz[ind_iat,0] - xyz[ind_i,0] + tvec[0]
	     dy = xyz[ind_iat,1] - xyz[ind_i,1] + tvec[1]
	     dz = xyz[ind_iat,2] - xyz[ind_i,2] + tvec[2]
	     #
	     # Absolute distance
	     r = np.sqrt(dx * dx + dy * dy + dz * dz)
	     if r > cutoff_radius:
	         continue
	     #
	     # Covalent distance in Bohr
	     rco = get_rcov_values(atoms.get_atomic_numbers()[ind_i]) + get_rcov_values(atoms.get_atomic_numbers()[ind_iat])
	     #
	     # Ratio covalent to absolute distance
	     rr = rco / r
	     #
	     # Derivative w.r.t. dx, dy, or dz
	     rrr = 1.0 / (r * r * r)
	     #
	     # Exponential term
	     tmp1 = np.exp(- K1_PARAMETER * (rr - 1.0))
	     #
	     # Damping function
	     tmp2 = 1.0 / (tmp1 + 1.0)
	     #
	     # Intermediate for derivative
	     tmp3 = tmp1 * tmp2 * tmp2 * K1_PARAMETER * rco * rrr

	     # add
	     xn                   += tmp2
	     #
	     dcn[ind_i,ind_i,0]  -= tmp3 * dx
	     dcn[ind_i,ind_i,1]  -= tmp3 * dy
	     dcn[ind_i,ind_i,2]  -= tmp3 * dz
	     #
	     dcn[ind_i,ind_iat,0] = tmp3 * dx
	     dcn[ind_i,ind_iat,1] = tmp3 * dy
	     dcn[ind_i,ind_iat,2] = tmp3 * dz
        #
	# Transfer
        cn[ind_i] = xn
    #
    return cn, dcn
    #
#end get_coordination_numbers
#
#
def lattice_sum(func, f0, f0_prime, positions, elem_cell=np.eye(3), periodic_directions=(False, False, False), cutoff_radius=DF_CUTOFF_RADIUS):
    """Computes the lattice sum of interactions within the (0,0,0) copy
    as well as the interactions to all other copies that lie whithin
    a sphere defined by the "cutoff_radius". Thereby the interaction
    is defined by a distance-vector-dependent function "func"

        >>> from math import pi
        >>> v, g = lattice_sum(lambda x: (1, 1), [[0., 0., 0.]], periodic_directions=3*(True,))
        >>> v / (4. * pi / 3. * DF_CUTOFF_RADIUS**3)
        1.0056889511012566
    """

    # Number of atoms within a single copy
    N_atoms = len(positions)

    # Maximum distance of two atoms within a copy
    # Serves as measure for extension of a single copy
    Rij_max = maxdist(positions)

    # use this below as cutoff radius:
    cutoff = cutoff_radius / AU_TO_ANG + Rij_max

    #
    # So far we have to provide the meaningfull (non-linearly dependent) cell vectors
    # also for directions that are not periodic in order to be able to invert the 3x3
    # matrix. FIXME: Is there a better way? Also ASE sets the system in that way?
    #

    # compute the size of the box enclosing a sphere of radius |cutoff|
    # in "fractional" coordinates:
    box = minbox(elem_cell, cutoff)

    # round them in very conservative fashion, all box[i] >= 1,
    box = [ int(ceil(k)) for k in box ]

    #
    # This rounding approach appears to never restrict the summation
    # later to a single cell, rather to at least three at
    #
    #   -1, 0 , +1
    #
    # of the cell vector in each direction.
    #

    # reset the values for non-periodic directions to zero:
    for i in range(len(box)):
        if not periodic_directions[i]:
            # in this direction treat only the unit cell itself:
            box[i] = 0

    #
    f       = f0
    f_prime = f0_prime
    # Now we are ready to sum over all cells in the box
    #
    # scan over indices u, v, w
    for u in xrange(-box[0], box[0] + 1):
        for v in xrange(-box[1], box[1] + 1):
            for w in xrange(-box[2], box[2] + 1):

                # Calculate translation vector to actual copy
                t_vec = u * elem_cell[0] + v * elem_cell[1] + w * elem_cell[2]

                # Calculate the squared distance to copy, cycle if that is too long:
                t2 = np.dot(t_vec, t_vec)
                if t2 > cutoff**2: continue

                # Calculate the dispersion interaction between
                # the considered pair of copies
                f1, fprime1 = func(t_vec)
                #
                # Actualize dispersion correction and contributions to forces
                f       += f1
                f_prime += fprime1
            #
        #
    #
    return f, f_prime
# End of function lattice_sum
#
#
def get_rcov_values(z_number):
    #
    rcov_raw_data = np.array([ None ,  # 0
                0.32 ,  # H
                0.46 ,  # He
                1.20 ,  # Li
                0.94 ,  # Be
                0.77 ,  # B
                0.75 ,  # C
                0.71 ,  # N
                0.63 ,  # O
                0.64 ,  # F
                0.67 ,  # Ne
                1.40 ,  # Na
                1.25 ,  # Mg
                1.13 ,  # Al
                1.04 ,  # Si
                1.10 ,  # P
                1.02 ,  # S
                0.99 ,  # Cl
                0.96 ,  # Ar
                1.76 ,  # K
                1.54 ,  # Ca
                1.33 ,  # Sc
                1.22 ,  # V
                1.21 ,  # Ti
                1.10 ,  # Cr
                1.07 ,  # Mn
                1.04 ,  # Fe
                1.00 ,  # Co
                0.99 ,  # Ni
                1.01 ,  # Cu
                1.09 ,  # Zn
                1.12 ,  # Ga
                1.09 ,  # Ge
                1.15 ,  # As
                1.10 ,  # Se
                1.14 ,  # Br
                1.17 ,  # Kr
                1.89 ,  # Rb
                1.67 ,  # Sr
                1.47 ,  # Y
                1.39 ,  # Zr
                1.32 ,  # Nb
                1.24 ,  # Mo
                1.15 ,  # Tc
                1.13 ,  # Ru
                1.13 ,  # Rh
                1.08 ,  # Pd
                1.15 ,  # Ag
                1.23 ,  # Cd
                1.28 ,  # In
                1.26 ,  # Sn
                1.26 ,  # Sb
                1.23 ,  # Te
                1.32 ,  # I
                1.31 ,  # Xe
                2.09 ,  # Cs
                1.76 ,  # Ba
                1.62 ,  # La
                1.47 ,  # Ce
                1.58 ,  # Pr
                1.57 ,  # Nd
                1.56 ,  # Pm
                1.55 ,  # Sm
                1.51 ,  # Eu
                1.52 ,  # Gd
                1.51 ,  # Tb
                1.50 ,  # Dy
                1.49 ,  # Ho
                1.49 ,  # Er
                1.48 ,  # Tm
                1.53 ,  # Yb
                1.46 ,  # Lu
                1.37 ,  # Hf
                1.31 ,  # Ta
                1.23 ,  # W
                1.18 ,  # Re
                1.16 ,  # Os
                1.11 ,  # Ir
                1.12 ,  # Pt
                1.13 ,  # Au
                1.32 ,  # Hg
                1.30 ,  # Tl
                1.30 ,  # Pb
                1.36 ,  # Bi
                1.31 ,  # Po
                1.38 ,  # At
                1.42 ,  # Rn
                2.01 ,  # Fr
                1.81 ,  # Ra
                1.67 ,  # Ac
                1.58 ,  # Th
                1.52 ,  # Pa
                1.53 ,  # U
                1.54 ,  # Np
                1.55 ]  # Pu
                , dtype=np.float64)
    #
    rcov = rcov_raw_data[z_number] * 4.0 / (3.0 * AU_TO_ANG)
    return rcov
# end get_rcov_values
#
#
def lattice_sumX(func, positions, elem_cell=np.eye(3), periodic_directions=(False, False, False), cutoff_radius=DF_CUTOFF_RADIUS):
    """Computes the lattice sum of interactions within the (0,0,0) copy
    as well as the interactions to all other copies that lie whithin
    a sphere defined by the "cutoff_radius". Thereby the interaction
    is defined by a distance-vector-dependent function "func"
    """

    # Number of atoms within a single copy
    N_atoms = len(positions)

    # Maximum distance of two atoms within a copy
    # Serves as measure for extension of a single copy
    Rij_max = maxdist(positions)

    # use this below as cutoff radius:
    cutoff = cutoff_radius + Rij_max

    #
    # So far we have to provide the meaningfull (non-linearly dependent) cell vectors
    # also for directions that are not periodic in order to be able to invert the 3x3
    # matrix. FIXME: Is there a better way? Also ASE sets the system in that way?
    #

    # compute the size of the box enclosing a sphere of radius |cutoff|
    # in "fractional" coordinates:
    box = minbox(elem_cell, cutoff)

    # round them in very conservative fashion, all box[i] >= 1,
    box = [ int(ceil(k)) for k in box ]

    #
    # This rounding approach appears to never restrict the summation
    # later to a single cell, rather to at least three at
    #
    #   -1, 0 , +1
    #
    # of the cell vector in each direction.
    #

    # reset the values for non-periodic directions to zero:
    for i in range(len(box)):
        if not periodic_directions[i]:
            # in this direction treat only the unit cell itself:
            box[i] = 0

    #
    # Now we are ready to sum over all cells in the box
    #

    # Initialization of output data by call for the 000 cell
    f = func(np.array([0.,0.,0.]))
    #
    # Running over all residual copies in the box
    u = 0
    v = 0
    # u, v = 0 case. Summation over w, with w != 0
    for w in xrange(-box[2], 0):
        t_vec = u * elem_cell[0] + v * elem_cell[1] + w * elem_cell[2]
        t2 = np.dot(t_vec, t_vec)
        if t2 > cutoff**2: continue
        f  += func(t_vec)
    for w in xrange(1, box[2] + 1):
        t_vec = u * elem_cell[0] + v * elem_cell[1] + w * elem_cell[2]
        t2 = np.dot(t_vec, t_vec)
        if t2 > cutoff**2: continue
        f  += func(t_vec)
    #
    # u = 0 case. Summation over v and w, with v != 0
    for v in xrange(-box[1], 0):
        for w in xrange(-box[2], box[2] + 1):
            t_vec = u * elem_cell[0] + v * elem_cell[1] + w * elem_cell[2]
            t2 = np.dot(t_vec, t_vec)
            if t2 > cutoff**2: continue
            f  += func(t_vec)
        #
        for w in xrange(-box[2], box[2] + 1):
            t_vec = u * elem_cell[0] - v * elem_cell[1] + w * elem_cell[2]
            t2 = np.dot(t_vec, t_vec)
            if t2 > cutoff**2: continue
            f  += func(t_vec)
    #
    # General case. Summation over u, v, and w, with u != 0
    for u in xrange(-box[0], 0):
        for v in xrange(-box[1], box[1] + 1):
            for w in xrange(-box[2], box[2] + 1):
                t_vec = u * elem_cell[0] + v * elem_cell[1] + w * elem_cell[2]
                t2 = np.dot(t_vec, t_vec)
                if t2 > cutoff**2: continue
                f  += func(t_vec)
        #
        for v in xrange(-box[1], box[1] + 1):
            for w in xrange(-box[2], box[2] + 1):
                t_vec = -u * elem_cell[0] + v * elem_cell[1] + w * elem_cell[2]
                t2 = np.dot(t_vec, t_vec)
                if t2 > cutoff**2: continue
                f  += func(t_vec)
    #
    return f
    #
# End of function lattice_sumX
#
#
# check for interaction group input
def get_interaction_controls(atoms):
    # Check if interaction list is present
    try:
	interactionlist   = np.array(atoms.interactionlist)
    except AttributeError:
	interactionlist   = None
    #
    # Check if interaction matrix is present
    if interactionlist != None:
	try:
	    interactionmatrix   = np.array(atoms.interactionmatrix, dtype=bool)
	except AttributeError:
	    interactionmatrix   = None
    else:
        interactionmatrix   = None
    #
    return interactionlist, interactionmatrix
# End of function get_interaction_controls
#
#
# Convert xc_functional designators
def xc_name(raw_xc_name):
  conv_xc_name = '000'
  if (str.lower(raw_xc_name) == 'pbe'):
    conv_xc_name = 'pbe'
  if (str.lower(raw_xc_name) == 'ps'):
    conv_xc_name = 'pbesol'
  if (str.lower(raw_xc_name) == 'rp'):
    conv_xc_name = 'revpbe'
  if (str.lower(raw_xc_name) == 'bp' or str.lower(raw_xc_name) == 'bp86'):
    conv_xc_name = 'b-p'
  if (str.lower(raw_xc_name) == 'tpss'):
    conv_xc_name = 'tpss'
  return conv_xc_name
#
# Autocompilation subroutine
def dftd_autocompile():
  import os
  import sys
  print
  print ' DFT-D module seems not to be available! Attempt of compilation ...'
  #
  # Get list of PYTHONPATH-entries
  paths = os.environ['PYTHONPATH'].split(os.pathsep)
  #
  # Search within PYTHONPATH for a subdirectory ase/dftd and compile
  # with makefile DFTD_Makefile
  for path in paths:
    print
    print '    Looking for ase in ', path
    print
    os.system('cd '+path+'; cd ase'+os.sep+'dftd ; make -f DFTD_Makefile')
    print

