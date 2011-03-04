"""This module defines a classical ASE interface to GPAW.
"""
import numpy as np
from gpaw import GPAW
from ase.dftd.dftd_interface import d2_pbc as d2_pbc
from ase.dftd.dftd_interface import d3_pbc as d3_pbc
from ase.dftd.dftd_interface import xc_name as xc_name

class GPAW_d2(GPAW):
    #
    def get_potential_energy(self, atoms, force_consistent=False):
        """
        Altered version of the original get_potential_energy function
        in the class GPAW. The function obtains the DFT energy by using
        the original call for the GPAW package and adds the DFT-D2 contribution
        to the converged SCF-energy.
        """
        #
        # Conversion factors a.u. -> eV
        Eh__2__eV          = 27.211396132
        #
        # Calling original GPAW-calculator for energy
        self.ks_scf_energy = GPAW.get_potential_energy(self, atoms, force_consistent)
        #
        # Get functional name as string
	self.functional = xc_name(str.lower(GPAW.get_xc_functional(self)))
        #
        # Call DFT-D module: Energy and gradients
        self.dispersion_correction, self.dft_d_gradient_contribution = d2_pbc(atoms, self.functional)
        #
        # Convert to proper units
        self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
        #
#       # Print out components (Useful?)
#       print
#       print 'DFT total energy  : ', self.ks_scf_energy
#       print 'DFT-D2 correction : ', self.dispersion_correction
#       print
#       print 'DFT-D2 final corrected energy: ', self.ks_scf_energy + self.dispersion_correction
#       print
        #
        # Adding correction contribution to energy
        return self.ks_scf_energy + self.dispersion_correction
        #
    def get_forces(self, atoms):
        """
        Altered version of the original get_forces function in the GPAW-class.
        The function obtains the DFT forces by using the original call for the
        GPAW package and adds the DFT-D contribution to the calculated forces.
        """
        #
        # Conversion factors a.u. -> eV and a.u. -> eV/Angst
        Eh__2__eV          = 27.211396132
        Eh_rb__2__eV_Angst = 51.422086162
        #
        # Calling original VASP-calculator for forces
        self.dft_forces = GPAW.get_forces(self, atoms)
        #
        # Get functional name as string
	self.functional = xc_name(str.lower(GPAW.get_xc_functional(self)))
        #
        # Call DFT-D module: Energy and gradients
        self.dispersion_correction, self.dft_d_gradient_contribution = d2_pbc(atoms,self.functional)
        #
        # Convert to proper units
        self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
        self.dft_d_gradient_contribution = self.dft_d_gradient_contribution * Eh_rb__2__eV_Angst
        #
        # Adding correction contributions to forces
        # Note the (-) sign: DFT-D module delivers gradients, not forces
        return self.dft_forces - self.dft_d_gradient_contribution
        #
# End of class GPAW_d2
#
class GPAW_d3(GPAW):
    #
    def get_potential_energy(self, atoms, force_consistent=False):
        """
        Altered version of the original get_potential_energy function
        in the class GPAW. The function obtains the DFT energy by using
        the original call for the GPAW package and adds the DFT-D3 contribution
        to the converged SCF-energy.
        """
        #
        # Conversion factors a.u. -> eV
        Eh__2__eV          = 27.211396132
        #
        # Calling original GPAW-calculator for energy
        self.ks_scf_energy = GPAW.get_potential_energy(self, atoms, force_consistent)
        #
        # Get functional name as string
        self.functional = str.lower(GPAW.get_xc_functional(self))
        #
        # Call DFT-D module: Energy and gradients
        self.dispersion_correction, self.dft_d_gradient_contribution = d3_pbc(atoms, self.functional)
        #
        # Convert to proper units
        self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
        #
#       # Print out components (Useful?)
#       print
#       print 'DFT total energy  : ', self.ks_scf_energy
#       print 'DFT-D3 correction : ', self.dispersion_correction
#       print
#       print 'DFT-D3 final corrected energy: ', self.ks_scf_energy + self.dispersion_correction
#       print
        #
        # Adding correction contribution to energy
        return self.ks_scf_energy + self.dispersion_correction
        #
    def get_forces(self, atoms):
        """
        Altered version of the original get_forces function in the GPAW-class.
        The function obtains the DFT forces by using the original call for the
        GPAW package and adds the DFT-D contribution to the calculated forces.
        """
        #
        # Conversion factors a.u. -> eV and a.u. -> eV/Angst
        Eh__2__eV          = 27.211396132
        Eh_rb__2__eV_Angst = 51.422086162
        #
        # Calling original VASP-calculator for forces
        self.dft_forces = GPAW.get_forces(self, atoms)
        #
        # Get functional name as string
	self.functional = xc_name(str.lower(GPAW.get_xc_functional(self)))
        #
        # Call DFT-D module: Energy and gradients
        self.dispersion_correction, self.dft_d_gradient_contribution = d3_pbc(atoms, self.functional)
        #
        # Convert to proper units
        self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
        self.dft_d_gradient_contribution = self.dft_d_gradient_contribution * Eh_rb__2__eV_Angst
        #
        # Adding correction contributions to forces
        # Note the (-) sign: DFT-D module delivers gradients, not forces
        return self.dft_forces - self.dft_d_gradient_contribution
        #
# End of class GPAW_d3
#
#
class d3(GPAW):
    #
    def get_potential_energy(self, atoms, force_consistent=False):
        """
        Altered version of the original get_potential_energy function
        in the class GPAW. The function obtains the DFT energy by using
        the original call for the GPAW package and adds the DFT-D2 contribution
        to the converged SCF-energy.
        """
        #
        # Conversion factors a.u. -> eV
        Eh__2__eV          = 27.211396132
        #
        # Calling original GPAW-calculator for energy
        self.ks_scf_energy = 0.0 #GPAW.get_potential_energy(self, atoms, force_consistent)
        #
        # Get functional name as string
	self.functional = xc_name(str.lower(GPAW.get_xc_functional(self)))
        #
        # Call DFT-D module: Energy and gradients
        self.dispersion_correction, self.dft_d_gradient_contribution = d3_pbc(atoms, self.functional)
        #
        # Convert to proper units
        self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
        #
        # Print out components (Useful?)
        print
        print 'DFT total energy  : ', self.ks_scf_energy
        print 'DFT-D2 correction : ', self.dispersion_correction
        print
        print 'DFT-D2 final corrected energy: ', self.ks_scf_energy + self.dispersion_correction
        print
        #
        # Adding correction contribution to energy
        return self.ks_scf_energy + self.dispersion_correction
        #
    def get_forces(self, atoms):
        """
        Altered version of the original get_forces function in the GPAW-class.
        The function obtains the DFT forces by using the original call for the
        GPAW package and adds the DFT-D contribution to the calculated forces.
        """
        #
        # Conversion factors a.u. -> eV and a.u. -> eV/Angst
        Eh__2__eV          = 27.211396132
        Eh_rb__2__eV_Angst = 51.422086162
        #
        # Calling original VASP-calculator for forces
        #self.dft_forces = GPAW.get_forces(self, atoms)
        #
        # Get functional name as string
	self.functional = xc_name(str.lower(GPAW.get_xc_functional(self)))
        #
        # Call DFT-D module: Energy and gradients
        self.dispersion_correction, self.dft_d_gradient_contribution = d3_pbc(atoms,self.functional)
        #
        # Convert to proper units
        self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
        self.dft_d_gradient_contribution = self.dft_d_gradient_contribution * Eh_rb__2__eV_Angst
        #
        # Adding correction contributions to forces
        # Note the (-) sign: DFT-D module delivers gradients, not forces
        return - self.dft_d_gradient_contribution
        #
# End of class d3
#
class d2(GPAW):
    #
    def get_potential_energy(self, atoms, force_consistent=False):
        """
        Altered version of the original get_potential_energy function
        in the class GPAW. The function obtains the DFT energy by using
        the original call for the GPAW package and adds the DFT-D2 contribution
        to the converged SCF-energy.
        """
        #
        # Conversion factors a.u. -> eV
        Eh__2__eV          = 27.211396132
        #
        # Calling original GPAW-calculator for energy
        self.ks_scf_energy = 0.0 #GPAW.get_potential_energy(self, atoms, force_consistent)
        #
        # Get functional name as string
        self.functional = 'pbe' #str.lower(GPAW.get_xc_functional(self))
        #
        # Call DFT-D module: Energy and gradients
        self.dispersion_correction, self.dft_d_gradient_contribution = d2_pbc(atoms, self.functional)
        #
        # Convert to proper units
        self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
        #
        # Print out components (Useful?)
        print
        print 'DFT total energy  : ', self.ks_scf_energy
        print 'DFT-D2 correction : ', self.dispersion_correction
        print
        print 'DFT-D2 final corrected energy: ', self.ks_scf_energy + self.dispersion_correction
        print
        #
        # Adding correction contribution to energy
        return self.ks_scf_energy + self.dispersion_correction
        #
    def get_forces(self, atoms):
        """
        Altered version of the original get_forces function in the GPAW-class.
        The function obtains the DFT forces by using the original call for the
        GPAW package and adds the DFT-D contribution to the calculated forces.
        """
        #
        # Conversion factors a.u. -> eV and a.u. -> eV/Angst
        Eh__2__eV          = 27.211396132
        Eh_rb__2__eV_Angst = 51.422086162
        #
        # Calling original VASP-calculator for forces
        #self.dft_forces = GPAW.get_forces(self, atoms)
        #
        # Get functional name as string
        self.functional = 'pbe' #str.lower(GPAW.get_xc_functional(self))
        #
        # Call DFT-D module: Energy and gradients
        self.dispersion_correction, self.dft_d_gradient_contribution = d2_pbc(atoms,self.functional)
        #
        # Convert to proper units
        self.dispersion_correction       = self.dispersion_correction       * Eh__2__eV
        self.dft_d_gradient_contribution = self.dft_d_gradient_contribution * Eh_rb__2__eV_Angst
        #
        # Adding correction contributions to forces
        # Note the (-) sign: DFT-D module delivers gradients, not forces
        return - self.dft_d_gradient_contribution
        #
# End of class d2
#
