from ase.atom import Atom
from ase.atoms import Atoms
from ase.gxfile import fromgx, gxwrite
from ase.units import Bohr, Hartree
from numpy import ones, zeros

def read_gx(file = "gxfile"):
    """
    Reads in like the other input files (result is Atom object)
    """
    ats = [ Atom(*a) for a in fromgx(file)]
    atoms = Atoms(ats)
    return atoms


def write_gx(filename, atoms, with_energy = False, with_gradients = False, loop = 1, isyms = None, inums = None, iconns = None, ivars = None):
    energy = None
    if with_energy:
        energy = atoms.get_potential_energy() / Hartree

    grads = None
    if with_gradients:
        grads = -atoms.get_forces() / Hartree * Bohr

    positions = atoms.get_positions().copy()
    atnums = atoms.get_atomic_numbers().copy()
    n = len(atnums)
    if isyms == None:
        isyms = ones(n)
    if inums == None:
        inums = zeros(n)
    if iconns == None:
        iconns = zeros((n,3))
    if ivars == None:
        ivars = zeros((n,3))

    gxwrite(atnums, positions/Bohr, isyms, inums, iconns, ivars, grads, energy, loop, file=filename )

