#!/usr/bin/python

from ase import Atoms
import numpy as np
from math import sqrt, exp, ceil

# General Parameters
# steepness of damping function
ALPHA = 20.0
# Default value of cutoff radius in Angstrom
DF_CUTOFF_RADIUS = 9.0

def test():
    h3 = Atoms('H3')

    x1 = ( 1, 0, 0)
    x2 = ( 0, 0, 1)
    x3 = (-1, 0, 0)

    h3.set_positions([x1,x2,x3])
    h3.set_atomic_numbers([1,1,1])
    h3.set_cell([(4.0, 0.0, 0.0), (1.0, 3.0, 0.0), (1.0, 3.0, 3.0)])
    h3.set_pbc([1,1,1])
    liste = [0, 2, 6]
    dc, dcg = dft_d_pbc(h3,interactionlist=liste)
    print dc

# Function dft_d_pbc: Main function making the D-G06 DFT-D correction available.
#                     Checks consistency of input (interactionlist and interactionmatrix).
#                     Applies the lattice summation defined in "lattice_sum" to
#                     the intra- and inter-cell evaluation of the D-G06 correction
#                     carried out in "d_g06_cell".
def dft_d_pbc(atoms, scaling_factor=0.75,interactionlist=[None,],interactionmatrix=[None,],cutoff_radius=DF_CUTOFF_RADIUS):
    #
    # Obtain data about system
    atom_numbers        = atoms.get_atomic_numbers()
    positions           = atoms.get_positions()
    periodic_directions = atoms.get_pbc()
    elem_cell           = atoms.get_cell()
    #
    # Number of atoms within a single copy
    N_atoms = len(positions)
    #
    # Check input i.e. if interactionlist and interactionmatrix are set properly
    if interactionlist == [None]:
        # Case 1: Nothing given, all atoms in the same group and interacting with each other
	interactionlist = [0,]*len(atom_numbers)
	interactionmatrix = np.ones((1,1))
    elif interactionlist[1] != [None] and interactionmatrix==[None]:
        # Case 2: List given, but interaction matrix has to be set to its default
        N_groups = max(interactionlist) + 1
        if len(interactionlist) != len(atom_numbers): raise SyntaxError, 'DFT-D: interaction list not of expected size'
        interactionmatrix = np.ones((N_groups,N_groups)) - np.eye((N_groups))
    elif interactionlist[1] != [None] and interactionmatrix != [None]:
        # Case 3: List and interaction matrix given
        N_groups = max(interactionlist) + 1
        if len(interactionlist) != len(atom_numbers): raise SyntaxError, 'DFT-D: interaction list not of expected size'
        if not (len(interactionmatrix[:,1]) == N_groups and len(interactionmatrix[1,:]) == N_groups): raise SyntaxError, 'DFT-D: interaction matrix not of expected size'
    #
    # Start with Calculation
    # Get DFT-D parameters
    params = [d_g06_parameters(ind_1) for ind_1 in atom_numbers ]
    #
    # Define function "func" as "d_g06" for the lattice summation done in "lattice_sum"
    def func(t):
        return d_g06_cell(N_atoms, params, positions, t, interactionlist, interactionmatrix, cutoff_radius)
    #
    # Call lattice summation with function "d_g06_cell"
    dispersion_correction, gradient_contribution = lattice_sum(func, positions, elem_cell, periodic_directions, cutoff_radius)
    #
    # Scale dispersion correction and forces according to used XC-functional
    dispersion_correction = dispersion_correction * scaling_factor
    gradient_contribution = gradient_contribution * scaling_factor
    #
    return dispersion_correction, gradient_contribution
# End of function dft_d_pbc
#
# Function lattice_sum: Computes the lattice sum of interactions within the (0,0,0) copy
#                       as well as the interactions to all other copies that lie whithin
#                       a sphere defined by the "cutoff_radius". Thereby the interaction
#                       is defined by a distance-vector-dependent function "func"
def lattice_sum(func, positions, elem_cell=np.eye(3), periodic_directions=3*(True,), cutoff_radius=DF_CUTOFF_RADIUS):
    """
        >>> from math import pi
	>>> v, g = lattice_sum(lambda x: (1, 1), [[0., 0., 0.]])
	>>> v / (4. * pi / 3. * DF_CUTOFF_RADIUS**3)
	1.0056889511012566
    """
    #
    # Number of atoms within a single copy
    N_atoms = len(positions)
    #
    # Initialization of output values
    f       = 0.0
    f_prime = np.zeros((3,N_atoms))
    #
    # Maximum distance of two atoms within a copy
    # Serves as measure for extension of a single copy
    Rij_max               = 0.0
    for ind_i in range(0, N_atoms - 1):
        for ind_j in range(ind_i + 1, N_atoms ):
	    # Determine pairwise distance
            dirij = positions[ind_j,:] - positions[ind_i,:]
	    Rij   = sqrt(sum(dirij**2))
	    # Check if larger than actual maximum
            if Rij_max < Rij: Rij_max = Rij
    #
    # Determine dual vectors to the unit-cell vectors
    D_vecs = np.matrix(elem_cell).I
    #
    # Search for maximum indices: max_ind = |D_vec|*|R_cut+Rij_max|
    max_ind = [0,]*3
    for ind_i in range(0,3):
        if periodic_directions[ind_i]:
            D_vec_norm = sqrt(sum(np.multiply(D_vecs[:,ind_i],D_vecs[:,ind_i])))
            max_ind[ind_i] = int(ceil(D_vec_norm * (cutoff_radius + Rij_max)))
    #
    # scan over indices u, v, w
    for ind_u in range(-max_ind[0],max_ind[0]+1):
        for ind_v in range(-max_ind[1],max_ind[1]+1):
            for ind_w in range(-max_ind[2],max_ind[2]+1):
	        # Calculate translation vector to actual copy
	        t_vec = ind_u * elem_cell[0] + ind_v * elem_cell[1] + ind_w * elem_cell[2]
		# Calculate distance to copy
		t_len = sqrt(sum(t_vec**2))
 		if t_len > cutoff_radius + Rij_max: continue
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
# Function d_g06_cell: Evaluation of the dispersion correction between two cells
#                      or within the (0,0,0) cell
def d_g06_cell(N_atoms, parameters, positions, t_vec, interactionlist, interactionmatrix, cutoff_radius):
    #
    # Initialization of output values
    dispersion_correction_cell = 0.0
    gradient_contribution_cell = np.zeros((3,N_atoms))
    #
    for ind_i in range(0, N_atoms):
        for ind_j in range(0, N_atoms):
            # Check if allowed by interaction groups
            if interactionmatrix[interactionlist[ind_i],interactionlist[ind_j]]:
                # Determine direction Atom<ind_1> -> Atom<ind_2>
                dirij = positions[ind_j,:] + t_vec - positions[ind_i,:]
                # Calculate Interatomic Distance Rij from dirij
                Rij = sqrt(sum(dirij**2))
                # Skip current loop if distance exceed cutoff radius
                if Rij > cutoff_radius: continue
                # Determine Pairwise Parameters R0ij
                R0ij = parameters[ind_i][1] + parameters[ind_j][1]
		# Skip current loop if distance too small
		if Rij < R0ij * 0.3: continue
                # Pairwise Parameter C6ij
                C6ij = sqrt(parameters[ind_i][0] * parameters[ind_j][0])
		# Pairwise interaction in terms of dimensionless variable dirij / R0ij
                Dij, dDij_ddirij = d_g06(dirij / R0ij)# * C6ij / (R0ij**6)
		Dij         = Dij * C6ij / (R0ij**6)
		dDij_ddirij = dDij_ddirij * C6ij / (R0ij**6)
		# Add contributions to correction term and its gradient
		# note the negative sign of the gradient: direction defined as i -> j
                dispersion_correction_cell    += Dij / 2.
                gradient_contribution_cell[ind_j] -= dDij_ddirij / 2.
	    #
        #	
    #
    return dispersion_correction_cell, gradient_contribution_cell
# End of function d_g06_cell
#
# Function d_g06: Pairwise evaluation of the D-G06 dispersion correction in terms of
#                 the dimensionless vector rij_vec / Rij0
def d_g06(rij_vec_Rij0):
    """
        >>> d_g06(np.array([1., 0., 0.]))
	(-0.5, array([-2., -0., -0.]))
        >>> d_g06(np.array([1., 1., 1.]))
	(-0.037037020814322689, array([ 0.07407385,  0.07407385,  0.07407385]))
        >>> d_g06(np.array([10., 0., 0.]))
        (-9.9999999999999995e-07, array([  6.00000000e-07,   0.00000000e+00,   0.00000000e+00]))
    """
    # Dimensionless distance between atoms from dimensionless vector
    Rij_Rij0 = sqrt(sum(rij_vec_Rij0**2))
    # Direction vector with length 1
    dirij = rij_vec_Rij0 / Rij_Rij0
    #
    # Exponential term in damping function 
    eij = exp(- ALPHA * (Rij_Rij0 - 1.))
    # Contribution of pair to dispersion correction
    Dij = - 1. / ((1. + eij) * Rij_Rij0**6)
    #
    # Building up derivatives
    deij_dRij_Rij0 = ALPHA * eij
    dDij_dRij_Rij0 = Dij * (deij_dRij_Rij0 / (1. + eij) - 6. / Rij_Rij0)
    #
    # Calculate gradients
    dDij_drij_vec_Rij0 = dirij * dDij_dRij_Rij0
    #
    return Dij, dDij_drij_vec_Rij0
# End of function d_g06
#
# Function d_g06_parameters: Delivers the atomic parameters C6 and R0 of the D-G06
#                            DFT-D correction
def d_g06_parameters(pse_number):
    """
        >>> d_g06_parameters(1)
        [1.4509976368338922, 1.0009999999999999]
        >>> d_g06_parameters(6)
        [18.13747046042365, 1.452]
        >>> d_g06_parameters(54)
        [310.82442234748873, 1.881]
    """
    # Element Specific Parameters
    PARS = [ None            , #
	    [ 0.14 , 1.001 ] , # H
	    [ 0.08 , 1.012 ] , # He
	    [ 1.61 , 0.825 ] , # Li
	    [ 1.61 , 1.408 ] , # Be
	    [ 3.13 , 1.485 ] , # B
	    [ 1.75 , 1.452 ] , # C
	    [ 1.23 , 1.397 ] , # N
	    [ 0.70 , 1.342 ] , # O
	    [ 0.75 , 1.287 ] , # F
	    [ 0.63 , 1.243 ] , # Ne
	    [ 5.71 , 1.144 ] , # Na
	    [ 5.71 , 1.364 ] , # Mg
	    [10.79 , 1.639 ] , # Al
	    [ 9.23 , 1.716 ] , # Si
	    [ 7.84 , 1.705 ] , # P
	    [ 5.57 , 1.683 ] , # S
	    [ 5.07 , 1.639 ] , # Cl
	    [ 4.61 , 1.595 ] , # Ar
	    [10.80 , 1.485 ] , # K
	    [10.80 , 1.474 ] , # Ca
	    [10.80 , 1.562 ] , # Sc
	    [10.80 , 1.562 ] , # V
	    [10.80 , 1.562 ] , # Ti
	    [10.80 , 1.562 ] , # Cr
	    [10.80 , 1.562 ] , # Mn
	    [10.80 , 1.562 ] , # Fe
	    [10.80 , 1.562 ] , # Co
	    [10.80 , 1.562 ] , # Ni
	    [10.80 , 1.562 ] , # Cu
	    [10.80 , 1.562 ] , # Zn
	    [16.99 , 1.650 ] , # Ga
	    [17.10 , 1.727 ] , # Ge
	    [16.37 , 1.760 ] , # As
	    [12.64 , 1.771 ] , # Se
	    [12.47 , 1.749 ] , # Br
	    [12.01 , 1.727 ] , # Kr
	    [24.67 , 1.628 ] , # Rb
	    [24.67 , 1.606 ] , # Sr
	    [24.67 , 1.639 ] , # Y
	    [24.67 , 1.639 ] , # Zr
	    [24.67 , 1.639 ] , # Nb
	    [24.67 , 1.639 ] , # Mo
	    [24.67 , 1.639 ] , # Tc
	    [24.67 , 1.639 ] , # Ru
	    [24.67 , 1.639 ] , # Rh
	    [24.67 , 1.639 ] , # Pd
	    [24.67 , 1.639 ] , # Ag
	    [24.67 , 1.639 ] , # Cd
	    [37.32 , 1.672 ] , # In
	    [38.71 , 1.804 ] , # Sn
	    [38.44 , 1.881 ] , # Sb
	    [31.74 , 1.892 ] , # Te
	    [31.50 , 1.892 ] , # I
	    [29.99 , 1.881 ] ] # Xe
    # Conversion of C6-parameter from  J * nm^6 / mol  ->  eV * Angstrom^6
    pars = PARS[pse_number]
    C_6_par = pars[0] * 10.3642688345278
    R_0_par = pars[1]
    #
    return [C_6_par, R_0_par]
# End of function d_g06_parameters
#
# python mueller_brown.py [-v]:
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    test()


# Default options for vim:sw=4:expandtab:smarttab:autoindent:syntax

