"""This module defines an alternative ASE interface to VASP.
In contrast to the original VASP-interface an additional
correction can be added to the potential energy and the forces.
"""
import numpy as np
import vasp
import sys
from ase.dftd.dftd_interface import d2_pbc as d2_pbc
from ase.dftd.dftd_interface import d3_pbc as d3_pbc
from ase.dftd.dftd_interface import xc_name as xc_name


class Vasp_d2(vasp.Vasp):
    #
    def get_potential_energy(self, atoms, force_consistent=False, dft_d_cutoff_radius=30.0):
	"""
	Altered version of the original get_potential_energy function
        in the class Vasp. The function obtains the DFT energy by using
	the original call for the VASP package and adds the DFT-D2 contribution
	to the converged SCF-energy.
	"""
        self.update(atoms)
	# Conversion factors a.u. -> eV
	Eh__2__eV          = 27.211396132
	#
	# Get functional name as string
	functional = xc_name(str.lower(vasp.Vasp.get_xc_functional(self)))
	#
	# Calling original VASP-calculator for energy
	self.energy_free_or_zero = vasp.Vasp.get_potential_energy(self, atoms, force_consistent)
	#
	# Call DFT-D module: Energy and gradients
	self.dispersion_correction, self.dft_d_gradient_contribution = d2_pbc(atoms, functional)
	#
	# Convert to proper units
	self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
	#
	# Print out components (Useful?)
	print >> sys.stdout,  ' '
	print >> sys.stdout,  'DFT total energy  : ', self.energy_free_or_zero
	print >> sys.stdout,  'DFT-D2 correction : ', self.dispersion_correction
	print >> sys.stdout,  ' '
	print >> sys.stdout,  'DFT-D2 final corrected energy: ', self.energy_free_or_zero + self.dispersion_correction
	print >> sys.stdout,  ' '
	#
	# Adding correction contribution to energy
	return self.energy_free_or_zero + self.dispersion_correction

    def get_forces(self, atoms, dft_d_cutoff_radius=30.0):
	"""
	Altered version of the original get_forces function in the Vasp-class.
	The function obtains the DFT forces by using the original call for the
	VASP package and adds the DFT-D2 contribution to the calculated forces.
	"""
        self.update(atoms)
	# Conversion factors a.u. -> eV and a.u. -> eV/Angst
	Eh__2__eV          = 27.211396132
	Eh_rb__2__eV_Angst = 51.422086162
	#
	# Get functional name as string
	functional = xc_name(str.lower(vasp.Vasp.get_xc_functional(self)))
	#
	# Calling original VASP-calculator for forces
	self.dft_forces = vasp.Vasp.get_forces(self, atoms)
	#
	# Call DFT-D module: Energy and gradients
	self.dispersion_correction, self.dft_d_gradient_contribution = d2_pbc(atoms, functional)
	#
	# Convert to proper units
	self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
	self.dft_d_gradient_contribution = self.dft_d_gradient_contribution * Eh_rb__2__eV_Angst
	#
	print
	print 'DFT-D total gradients:', -self.dft_d_gradient_contribution[0] + self.forces[0]
	for ind_i in range(1,len(self.forces)):
	    print '                      ', -self.dft_d_gradient_contribution[ind_i] + self.forces[ind_i]
	print
	#
	# Adding correction contributions to forces
	# Note the (-) sign: DFT-D module delivers gradients, not forces
	return self.dft_forces - self.dft_d_gradient_contribution
        #
# End of class Vasp_d2
#
#
class Vasp_d3(vasp.Vasp):
    #
    def get_potential_energy(self, atoms, force_consistent=False, dft_d_cutoff_radius=30.0):
	"""
	Altered version of the original get_potential_energy function
	in the class Vasp. The function obtains the DFT energy by using
	the original call for the VASP package and adds the DFT-D3 contribution
	to the converged SCF-energy.
	"""
#        self.update(atoms)
	# Conversion factors a.u. -> eV
	Eh__2__eV          = 27.211396132
	#
	# Get functional name as string
	functional = xc_name(str.lower(vasp.Vasp.get_xc_functional(self)))
	#
	# Calling original VASP-calculator for energy
	self.energy_free_or_zero = vasp.Vasp.get_potential_energy(self, atoms, force_consistent)
	#
	# Call DFT-D module: Energy and gradients
	self.dispersion_correction, self.dft_d_gradient_contribution = d3_pbc(atoms, functional)
	#
	# Convert to proper units
	self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
	#
	# Print out components (Useful?)
        print >> sys.stdout, ' '
        print >> sys.stdout, 'DFT total energy  : ', self.energy_free_or_zero
        print >> sys.stdout, 'DFT-D3 correction : ', self.dispersion_correction
        print >> sys.stdout, ' '
	print >> sys.stdout, 'DFT-D3 final corrected energy: ', self.energy_free_or_zero + self.dispersion_correction
        print >> sys.stdout, ' '
	#
	# Adding correction contribution to energy
	return self.energy_free_or_zero + self.dispersion_correction

    def get_forces(self, atoms, dft_d_cutoff_radius=30.0):
        """
	Altered version of the original get_forces function in the Vasp-class.
	The function obtains the DFT forces by using the original call for the
	VASP package and adds the DFT-D3 contribution to the calculated forces.
	"""
#        self.update(atoms)
	# Conversion factors a.u. -> eV and a.u. -> eV/Angst
	Eh__2__eV          = 27.211396132
	Eh_rb__2__eV_Angst = 51.422086162
	#
	# Get functional name as string
	functional = xc_name(str.lower(vasp.Vasp.get_xc_functional(self)))
	#
	# Calling original VASP-calculator for forces
	self.dft_forces = vasp.Vasp.get_forces(self, atoms)
	#
	# Call DFT-D module: Energy and gradients
	self.dispersion_correction, self.dft_d_gradient_contribution = d3_pbc(atoms, functional)
	#
	# Convert to proper units
	self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
	self.dft_d_gradient_contribution = self.dft_d_gradient_contribution * Eh_rb__2__eV_Angst
	#
       #print
       #print 'DFT-D3 total gradients:', -self.dft_d_gradient_contribution[0] + self.forces[0]
       #for ind_i in range(1,len(self.forces)):
       #    print '                      ', -self.dft_d_gradient_contribution[ind_i] + self.forces[ind_i]
       #print
	#
        # Adding correction contributions to forces
        # Note the (-) sign: DFT-D module delivers gradients, not forces
	return self.dft_forces - self.dft_d_gradient_contribution
        #
# End of class Vasp_d3
