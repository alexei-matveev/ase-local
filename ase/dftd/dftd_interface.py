#!/usr/bin/python
import numpy as np
from ase.dftd.dft_d2_native import check_interaction_group_input
from ase.dftd.dftd_module import d2_gradients as dftd2_gradients
from ase.dftd.dftd_module import d3_gradients as dftd3_gradients

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

