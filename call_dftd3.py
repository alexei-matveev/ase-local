import ase
import dftd_module
import numpy


system = ase.read('example_h2o.xyz')

#print system.get_positions()

#dftd_module.multi(5,4)

#dftd_module.dftd.division(36,6)

coords = [[1.,2.1,-0.43 ],[0.235, -4.0, 0.0001 ]]
elements = [1, 2]
#c = dftd_module.d3_energy(system.get_positions())
#c = dftd_module.d3_energy(a)
#c = dftd_module.d3_energy(numpy.transpose(numpy.array(a)))
print numpy.asarray(coords)
c = dftd_module.d3_energy(numpy.asarray(elements),numpy.asarray(coords))


print c

