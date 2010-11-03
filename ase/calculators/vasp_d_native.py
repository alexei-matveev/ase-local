"""This module defines an alternative ASE interface to VASP.
In contrast to the original VASP-interface an additional
correction can be added to the potential energy and the forces.
"""
import numpy as np
import vasp
import ase.dftd.dft_d2_native as dftd

class Vasp_d(vasp.Vasp):
    #
    def get_potential_energy(self, atoms, force_consistent=False, dft_d_cutoff_radius=25.0):
        """
	Altered version of the original get_potential_energy function
        in the class Vasp. The function obtains the DFT energy by using
	the original call for the VASP package and adds the DFT-D contribution
	to the converged SCF-energy.
	"""
#        self.update(atoms)
        #
	# Get interaction list and matrix
	interactionlist, interactionmatrix = get_interaction_controls(atoms)
	#
        # Calling original VASP-calculator for energy
        self.energy_free_or_zero = vasp.Vasp.get_potential_energy(self, atoms, force_consistent)
        #
        # Call DFT-D module: Energy and gradients
	self.dispersion_correction, self.dft_d_gradient_contribution = dftd.dft_d_pbc(atoms,
            interactionlist=interactionlist,
	    interactionmatrix=interactionmatrix,
	    cutoff_radius=dft_d_cutoff_radius)
        #
        # Print out components (Useful?)
        print
        print 'DFT total energy : ', self.energy_free_or_zero
        print 'DFT-D correction : ', self.dispersion_correction
        print
        print 'DFT-D final corrected energy: ', self.energy_free_or_zero + self.dispersion_correction
        print
        #
        # Adding correction contribution to energy
	return self.energy_free_or_zero + self.dispersion_correction

    def get_forces(self, atoms, dft_d_cutoff_radius=25.0):
        """
	Altered version of the original get_forces function in the Vasp-class.
	The function obtains the DFT forces by using the original call for the
	VASP package and adds the DFT-D contribution to the calculated forces.
	"""
#        self.update(atoms)
        #
        # Calling original VASP-calculator for forces
        self.energy_free_or_zero = vasp.Vasp.get_forces(self, atoms)
        #
	# Get interaction list and matrix
	interactionlist, interactionmatrix = get_interaction_controls(atoms)
        #
        # Call DFT-D module: Energy and gradients
	self.dispersion_correction, self.dft_d_gradient_contribution = dftd.dft_d_pbc(atoms,
            interactionlist=interactionlist,
	    interactionmatrix=interactionmatrix,
	    cutoff_radius=dft_d_cutoff_radius)
        #
        # Print out components (Useful?)
       #print
       #print 'DFT total forces   : ', self.forces[0]
       #for ind_i in range(1,len(self.forces)):
       #    print '                     ', self.forces[ind_i]
       #print 'DFT-D contribution : ', -self.dft_d_gradient_contribution[0]
       #for ind_i in range(1,len(self.forces)):
       #    print '                     ', -self.dft_d_gradient_contribution[ind_i]
       #print
        print
	print 'DFT-D total gradients:', -self.dft_d_gradient_contribution[0] + self.forces[0]
        for ind_i in range(1,len(self.forces)):
            print '                      ', -self.dft_d_gradient_contribution[ind_i] + self.forces[ind_i]
        print
        #
        # Adding correction contributions to forces
        # Note the (-) sign: DFT-D module delivers gradients, not forces
        return self.forces - self.dft_d_gradient_contribution

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

