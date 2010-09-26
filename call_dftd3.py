import ase
import dftd_module
import numpy


system = ase.read('example_h2o.xyz')

#print system.get_positions()

#dftd_module.multi(5,4)

#dftd_module.dftd.division(36,6)

coords=[[     .000000 ,   .000000 ,   .114079 ],
        [     .000000 ,   .780362 ,  -.456316 ],
        [     .000000 ,  -.780362 ,  -.456316 ],
        [    1.200000 ,   .000000 ,   .114079 ],
        [    1.200000 ,   .780362 ,  -.456316 ],
        [    1.200000 ,  -.780362 ,  -.456316 ] ]

elements = [8, 1, 1, 8, 1, 1]
#c = dftd_module.d3_energy(system.get_positions())
#c = dftd_module.d3_energy(a)
#c = dftd_module.d3_energy(numpy.transpose(numpy.array(a)))
print numpy.asarray(coords)
#c = dftd_module.d3_energy(numpy.asarray(elements),numpy.asarray(coords),'pbe')
c = dftd_module.d2_energy(numpy.asarray(elements),numpy.asarray(coords),'pbe')


print c

