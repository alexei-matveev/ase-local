"""Interfaces to different ASE compatible force-calculators."""

import numpy as np

from ase import _deprecate_things_from_ase_module
from ase.calculators.lj import LennardJones
from ase.calculators.emt import EMT
from ase.calculators.siesta import Siesta
from ase.calculators.dacapo import Dacapo
from ase.calculators.vasp import Vasp
from ase.calculators.singlepoint import SinglePointCalculator
from ase.calculators.test import numeric_force, numeric_forces, TestPotential
from ase.calculators.vasp_d_native import Vasp_d
from ase.calculators.vasp_d import Vasp_d3
from ase.calculators.vasp_d import Vasp_d2
from ase.calculators.aims import Aims, AimsCube
from ase.calculators.turbomole import Turbomole
from ase.calculators.paragauss import ParaGauss
from ase.calculators.exciting import Exciting
from ase.calculators.dftb import Dftb

if _deprecate_things_from_ase_module:
    from ase.utils.deprecate import Deprecate
    _locals = locals()
    for name in ['LennardJones', 'EMT', 'Siesta', 'Dacapo', 'Vasp',
                 'Aims', 'AimsCube', 'Turbomole', 'Exciting', 'Dftb',
                 'SinglePointCalculator', 'numeric_force', 'numeric_forces',
                 'TestPotential']:
        obj = _locals[name]
        _locals[name] = Deprecate(obj, name, obj.__module__, 'ase.calculators')
