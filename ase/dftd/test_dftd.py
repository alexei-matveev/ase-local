import os
import ase
from dftd_interface import d3_pbc as dftd3
from dftd_interface import d2_pbc as dftd2
import numpy

print 'single cell'
system = ase.read('./test/POSCAR5', format='vasp')
c, d = dftd3(system ,'pbe')
print repr(c)
print d

print 'mupltiple cell'
system = ase.read('./test/POSCAR6', format='vasp')
c, d = dftd3(system ,'pbe')
print repr(c)
print d

print 'reference'
os.system('./test/dftd3 ./test/POSCAR5.xyz -func pbe -noprint -grad')
