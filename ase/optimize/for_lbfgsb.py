from ase.optimize.sciopt import Converged, SciPyOptimizer

from ase.optimize import Optimizer


class SciPyFminLBFGSB(SciPyOptimizer):

    """Quasi-Newton method (Broydon-Fletcher-Goldfarb-Shanno)"""
    def call_fmin(self, fmax, steps):
        output = opt.fmin_bfgs(self.f,
                               self.x0(),
                               fprime=self.fprime,
                               #args=(), 
                               gtol=fmax * 0.1, #Should never be reached
                               norm=np.inf,
                               #epsilon=1.4901161193847656e-08, 
                               maxiter=steps,
                               #full_output=1, 
                               disp=0,
                               #retall=0, 
                               callback=self.callback
                              )


