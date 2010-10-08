#!/usr/bin/python
import numpy as np
from ase.dftd.dft_d2_native import check_interaction_group_input
from ase.dftd.dftd_module import d2_gradients as dftd2_gradients
from ase.dftd.dftd_module import d3_gradients as dftd3_gradients

# General parameters
AU_TO_ANG    = 0.52917726
K1_PARAMETER = 16.0
#
#
def d3_pbc(atoms, functional):
    """Main function making the DFT-D2 correction available for
    isolated systems.

    """
    #
    # Obtain data about system
    atom_numbers                       = atoms.get_atomic_numbers()
    positions                          = atoms.get_positions()
    interactionlist, interactionmatrix = get_interaction_controls(atoms)
    #
    # Number of atoms within a single copy
    N_atoms = len(positions)
    #
    # Check input i.e. if interactionlist and interactionmatrix are set properly
    interactionlist, interactionmatrix = check_interaction_group_input(N_atoms, interactionlist, interactionmatrix)
    #
    # Start with Calculatio
    dispersion_correction, gradient_contribution = dftd3_gradients(atoms.get_atomic_numbers(), atoms.get_positions(), interactionlist, interactionmatrix, functional)
    #
    print dispersion_correction
    print gradient_contribution
    # return results in a.u.
    return dispersion_correction, gradient_contribution
    #
# End of function d3_pbc
#
#
def d2_pbc(atoms, functional):
    """Main function making the DFT-D2 correction available for
    isolated systems.

    """
    #
    # Obtain data about system
    atom_numbers                       = atoms.get_atomic_numbers()
    positions                          = atoms.get_positions()
    interactionlist, interactionmatrix = get_interaction_controls(atoms)
    #
    # Number of atoms within a single copy
    N_atoms = len(positions)
    #
    # Check input i.e. if interactionlist and interactionmatrix are set properly
    interactionlist, interactionmatrix = check_interaction_group_input(N_atoms, interactionlist, interactionmatrix)
    #
    # Start with Calculatio
    dispersion_correction, gradient_contribution = dftd2_gradients(atoms.get_atomic_numbers(), atoms.get_positions(), interactionlist, interactionmatrix, functional)
    #
    print dispersion_correction
    print gradient_contribution
    # return results in a.u.
    return dispersion_correction, gradient_contribution
    #
# End of function d2_pbc
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
                1.55 ]) # Pu
    #
    rcov = rcov_raw_data[z_number] * 4.0 / (3.0 * AU_TO_ANG)
    return rcov
# end get_rcov_values
