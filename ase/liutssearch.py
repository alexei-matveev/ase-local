import numpy as np
from math import sqrt
import time, operator
from ase.optimize import Dynamics
from ase import BFGS
from ase.optimize.gxoptimizer import GxOptimizer

class LiuTSsearch(Dynamics):
    def __init__(self, atoms, restart=None, logfile='-' , trajectory=None, soften = 0, factsoft = 0.7,
                 outer_optimizer=BFGS, finish_optimizer = GxOptimizer, opt_args=None, relax_max=2, liuconstr = None,
                 switchfinish = False, treat1 = True, treat2 = True):
        """Structure optimizer object.
       source of the code is the following paper:
        H.-F. Wang and Z.-P. Liu: Comprehensive Mechanism and Structure-Sensitivity of Ethanol
        Oxidation on Platinum: New Transition-State Searching Method for Resolving the Complex
        Reaction Network, JACS 130 (2008), 10996

        atoms: Atoms object
            The Atoms object to relax.
        restart: str
            Filename for restart file.  Default value is *None*.
        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.
        trajectory: Trajectory object or str
            Attach trajectory object.  If *trajectory* is a string a
            PickleTrajectory will be constructed.  Use *None* for no
            trajectory.
        outer_optimizer: optimizer for the "Broyden steps"

        relax_max: after this amount of steps (integer) an update step
                  is performed
        liuconstr: constraint object
        """

        Dynamics.__init__(self, atoms, logfile, trajectory)
        # set optimizer
        if opt_args == None:
            self.dyn = outer_optimizer( atoms )
        else:
        # set optimizer with user given variables 
            self.dyn = outer_optimizer( atoms, **opt_args)

        self.restart = restart
        self.relax_max = relax_max
        self.stepsmax = None
        self.dyn.initialize()
        self.liuconstr = liuconstr
        self.treat1 = treat1
        self.treat2 = treat2
        self.soften = soften
        self.factsoft = factsoft
        self.switchfinish = switchfinish
        self.finish_optimizer = finish_optimizer
        #self.logout = open('liu.output','w')


    def setliuconstr(self, liuconstr):
        self.liuconstr = liuconstr


    def run(self, fmax=0.05, smax = 0.05, steps=100000000):
        """Run structure transition state search algorithm from Liu.

        This method will return when the forces on all individual
        atoms are less than *fmax* and stepsize are less than *smax*
        or when the number of steps exceeds
        *steps*."""

        self.fmax = fmax
        self.smax = smax
        step = 0
        switch_finish = 0
        # next one for counting the relaxation steps
        step_relax = 1
        self.writegeometry( self.atoms.get_positions() , second = 0)
        while step < steps:
            f = self.atoms.get_forces()
            oldval = self.atoms.get_positions()
            # perform actual Broydens step proposal
            # and set steps to it
            self.dyn.step(f)
            newval = self.atoms.get_positions()
            stepdiff = newval - oldval
            # end iterating if converged, compare forces and stepsize for it
            # sets also stepsmax to maximal value of steps
            if self.converged(f, stepdiff):
                print "System is converged, calcualtion stopped!"
                return
            # write some output
            self.log(f, self.stepsmax)
            # second part of liu algorithm, update every relax_max steps
            # or if the steps don't change anything
            # else only relax (set constraint back to old value)
            # what is actual constraint is decided by the Liuc...
            if (step_relax == self.relax_max or (self.stepsmax < self.smax)):
                self.dyn.initialize()
                if (self.treat2):
                    self.__treat2_init(newval, self.soften)
                    # print "before treatment2", newval
                    #self.writegeometry( newval)
                change =  self.liuconstr.update(oldval, newval, f)
                if (self.treat2):
                    #self.writegeometry(newval)
                    self.treatment2(newval, self.soften, self.factsoft, change)
                    #self.writegeometry( newval)
                    #print "after treatment2", newval
                step_relax = 0
            else:
                self.liuconstr.relaxation( oldval, newval, f, self.treat1)
            # the positions now have to given back to the atoms object
            self.atoms.set_positions(newval)
            # some more output
            self.writegeometry( self.atoms.get_positions())
            self.call_observers()
            self.nsteps += 1
            step += 1
            step_relax += 1
            if (self.switchfinish):
                 if (self.stepsmax < self.smax):
                      switch_finish += 1
                 else:
                     switch_finish = 0
                 if (switch_finish > self.relax_max * 3 + 4 ):
                     self.switchforlast(steps - step)
        print "maximum number of steps exceeded, calculation stopped"

    def converged(self, forces=None, stepdiff = None):
        """Did the optimization converge?"""
        # this if for convergence test without steps,then
        # stepsmax should always be below smax
        self.stepsmax = self.smax / 100
        if forces is None:
            forces = self.atoms.get_forces()
        if stepdiff is not None:
           self.stepsmax = (stepdiff ** 2).sum(axis=1).max()
        return ((forces ** 2).sum(axis=1).max() < self.fmax ** 2 and self.stepsmax < self.smax )

    def log(self, forces, stmax):
        fmax2 = sqrt((forces**2).sum(axis=1).max())
        e = self.atoms.get_potential_energy()
        T = time.localtime()
        if self.logfile is not None:
            name = self.__class__.__name__
            self.logfile.write('%s: %3d  %02d:%02d:%02d %15.6f %12.4f %12.4f\n' %
                               (name, self.nsteps, T[3], T[4], T[5], e, fmax2, stmax))
            self.logfile.flush()

    def switchforlast(steps):
        self.dyn = self.finish_optimizer(atoms)
        self.dyn.initialize()
        self.dyn.run(self.fmax, steps)

    def dump(self, data):
        if rank == 0 and self.restart is not None:
            pickle.dump(data, open(self.restart, 'wb'), protocol=2)

    def load(self):
        return pickle.load(open(self.restart))

    def writegeometry(self, new, fileout = 'geodata.xyz' , second = 1):
       # writes geometry in xyz file
       print "output to geodata.xyz"
       syms = self.atoms.get_chemical_symbols()
       if fileout == None:
           writegeo = stdout.write
       else:
           if fileout.endswith('.xyz'):
               pass
           else :
               fileout = fileout + '.xyz'
           if second==0:
             writegeo = open(fileout,"w").write
           else:
             writegeo = open(fileout,"a").write
       # firsts lines give number of atoms and comment
       writegeo("%i \n" % (len(syms) ) )
       writegeo("cartesian geometry in Angstrom\n")
       for num, pos in enumerate(new):
           writegeo("%2s " % syms[num])
           writegeo("%22.12f %22.12f %22.12f \n" % (pos[0], pos[1], pos[2]) )


    def __treat2_init( self, new, soften):
        # builds up the networks for the centers given by the LiuC class
        centers = self.liuconstr.centerofnetwork()
        self.networks =  [[] for inum in range(len(centers))]
        self.memberinshell =  [[] for inum in range(len(centers))]
        for inum, center in enumerate(centers):
            self.networks[inum], self.memberinshell[inum] = self.atomicnetwork(center, new)
            # print self.memberinshell[inum][0]
            # eliminate other centers out of list of other shells
            for cen in centers:
                self.networks[inum][cen].shell = 0
        if (soften < 0):
             self.allnetwork = self.mergenetworks(self.networks, centers)


    def treatment2( self, new, soften, p, ch,  cutoff = None):
        # the actual treatment2, changes *new* positions with the help of some factors
        # to adjust the positions of the other atoms
        # print "changes", ch
        if (soften < 0):
             maxshell = 0  
             for memberinshell in self.memberinshell:
                 maxshell = max(maxshell, len(memberinshell))
             maxshell -= 1
             if (cutoff != None):
                 maxshell = min(cutoff, maxshell)
             for s in range(maxshell):
                 for i in range(len(self.allnetwork)):
                      if (self.allnetwork[i].shell == s + 1):
                          # for this treatment2 the update has to be make shell wise
                          for k, smo  in enumerate(self.allnetwork[i].nextN):
                               distold = self.allnetwork[i].dist_to_next[k]
                               print "change of", self.allnetwork[i].number, "of shell", self.allnetwork[i].shell
                               print "old value", new[self.allnetwork[i].number-1]
                               dt = self.vect_of_change(new[self.allnetwork[i].number-1], new[smo -1], distold)
                               print dt
                               new[self.allnetwork[i].number-1] += self._adjustas(dt / self.allnetwork[i].multiplicity, soften, p, s + 1)
                               print "new value", new[self.allnetwork[i].number-1]
        else:
            for inum, network in enumerate(self.networks):
                maxshell = len(self.memberinshell[inum])
                if (cutoff != None ):
                    maxshell = min(cutoff, maxshell)
                for i in range(len(network)):
                    if (network[i].shell > 0 and network[i].shell < maxshell):
                        # print "change of", network[i].number, "of shell", network[i].shell
                        # print "old value", new[network[i].number-1]
                        new[network[i].number-1] += self._adjustas(ch[inum] , soften, p, network[i].shell)
                        # print "new value", new[network[i].number-1]

    def vect_of_change(self, pos1, pos2, dold):
         # vector from atom to be changed to atom to be fix
         # lengt of vector is difference between old and new bond length
         vec = pos2 - pos1
         dnew = self.distance(pos1, pos2)
         vec *= (dnew - dold)/dnew
         return vec

    def _adjustas(self, dt, soften, p, shell):
         # decides how much the positons of the actual atom are changed
         if (soften < 1):
              return ( p ** shell ) * dt
         else:
              if ( p * shell > 1.0 ):
                  return (1.0 - p * shell) * dt
              else:
                  return 0.0

    def atomicnetwork(self, center,  pos, cutoff = None):
        '''builds an atomic network around atom *center* (position in atomics list)
           pos are the positions considerd
           for a given cutoff, the network will stop with shell = cutoff -1
        '''
        # atomic network objects are stored in netinfo class in variable atnet
        # (will be given back)
        atnet = []
        # memberinshell counts the members for a given shell, the value
        # memberinshell[0] stores the amount of all atoms (except center) negativ for better handling
        # (will be given back)
        memberinshell = [ - (len(pos) - 1) ]
        # for the radii one needs to know what atomic numbers the atoms have
        atnums  = self.atoms.get_atomic_numbers()
        # There will now be stored all the atoms with initial value for number, there distance to
        # the center atom of the network and their radius, needed further on
        for i in range(len(pos)):
            d = self.distance(pos[i], pos[center])
            atnet.append(netinfo(i+1, dist_to_center = d, radius = self.__giveradius( atnums[i])))
        # start calculating the first shell
        # memberinshell needs new value for members of shell 0
        memberinshell.append(0)
        # for i in range(len(atnet)):
        #     print atnet[i]
        ra = atnet[center].radius
        # find out if atom i is in first shell ( d(i center) < radius(i) + radius(center))
        # if true, write shell 1 in atnet of atom i, multiplicity is always 0, because there is only one
        # atom in shell "0"
        # atom center should be in shell 1 too after this calculation, but that will be changed lateron
        for i in range(len(atnet)):
            if( atnet[i].dist_to_center <= (ra + atnet[i].radius)):
                atnet[i].shell = 1
                atnet[i].multiplicity = 1
                atnet[i].nextN = [center +1]
                atnet[i].dist_to_next = [atnet[i].dist_to_center] 
                memberinshell[1] += 1
        # if every atom left is in its own shell (one big line), the maximum number of shells is reached
        # there is no use in calculating more, if cutoff exists it replaces maxshell
        maxshell = -memberinshell[0] - memberinshell[1] + 1
        if (cutoff != None ):
             maxshell = min(cutoff - 2, maxshell)
             if (maxshell < 0):
                 maxshell = 0
             print "maximum shell considerd reduced to", maxshell
        # calculate all the others shells one ofter another, consider j starts with 0, but the first shell considered
        # is shell number 2 (j+2)
        for j in range( maxshell):
             # to store the number of atoms in this shell membershell needs j+2'th number
             memberinshell.append(0)
             for i in range(len(atnet)):
                 # as all atoms start in shell 0, the atoms who have not yet any shell assigned are in shell 0
                 if(atnet[i].shell < 1):
                    # make sure the other starting values are correct for unassigned atom
                    atnet[i].multiplicity = 0
                    atnet[i].nextN = []
                    atnet[i].dist_to_next = []
                    # need to consider all atoms from actual shell -1 as potential direct neigbors for this atom
                    for k in range(len(atnet)):
                       if(atnet[k].shell == j+1 ):
                           d = self.distance(pos[atnet[i].number-1], pos[atnet[k].number-1])
                           # there exist a bond to the considerd atom if following condition is fullfilled:
                           if ( d <= (atnet[i].radius + atnet[k].radius)):
                               # if at least one bond exists, the shell of atom i has been found (j + 2 because start 
                               # with shell 2 but j = 0), as there may be several atoms in shell j +1 which are conected
                               # to atom i, the multipllicity will be increased
                               atnet[i].shell = j + 2
                               atnet[i].multiplicity += 1
                               # print i, k, atnet[i].number, atnet[k].number
                               # store for further usage the number (identifier) and distance (to) of the atom which gave
                               # the last bond
                               atnet[i].nextN.append(atnet[k].number)
                               atnet[i].dist_to_next.append(d)
                    # if there exist at least one bond, the shell has a new member
                    if (atnet[i].multiplicity > 0):
                       memberinshell[j + 2] += 1
             # there is no possibility that after an empty shell in the next shell there will be found an atom
             # also it is no use to go on if all atoms have found their shell
             # (remember memberinshell[0] = -(sum of all atoms -1))
             if (memberinshell[j+2] == 0 or (sum(memberinshell) >= 1)  ):
                  break
        # the actual code puts the center atom itsself in the first shell, this sets it back
        atnet[center].shell = 0
        memberinshell[1] -= 1
        #print 'Data from network of atom', center+1
        #for i in range(len(atnet)):
        #    print atnet[i].number, atnet[i].dist_to_center, atnet[i].shell, atnet[i].multiplicity, atnet[i].nextN, atnet[i].dist_to_next
        #print memberinshell
        # return network and shelldistribution
        return atnet, memberinshell

    def mergenetworks(self, networks, centers):
        # merge networks for single atoms to one for several atoms all together
        for k, network in enumerate(networks):
             if (k == 0):
                 # just start with the first network
                 mergednetwork = network
             else:
                 for i in range(len(mergednetwork)):
                      # if the shell for the atom is smaller in the compared network it is taken (only for shells > 0), because
                      # 0 resembles no connection or center atom), if atom has a shell in new network but not in old one, new one
                      # is taken
                      if ((network[i].shell > 0 and network[i].shell < mergednetwork[i].shell) or (mergednetwork[i].shell == 0)):
                           mergednetwork[i] = network[i]
                      if ((network[i].shell > 0 and network[i].shell == mergednetwork[i].shell)):
                          mergednetwork[i].multiplicity += network[i].multiplicity
                          mergednetwork[i].nextN = mergednetwork[i].nextN + network[i].nextN
                          mergednetwork[i].dist_to_next = mergednetwork[i].dist_to_next + network[i].dist_to_next
        for center in centers:
             mergednetwork[center].shell = 0
        print 'Data from merged network of atom'
        for i in range(len(mergednetwork)):
            print mergednetwork[i].number, mergednetwork[i].dist_to_center, mergednetwork[i].shell, mergednetwork[i].multiplicity, mergednetwork[i].nextN, mergednetwork[i].dist_to_next
        return mergednetwork

    def distance(self, atoma, atomb):
        avec = atoma - atomb
        return sqrt(np.dot(avec, avec))

    def __giveradius(self, z):
        # gives the radius back
        gradius = [ None,
             0.3200, 0.9300, 1.2300, 0.9000, 0.8200,
             0.7700, 0.7500, 0.7300, 0.7200, 0.7100,
             1.5400, 1.3600, 1.1800, 1.1100, 1.0600,
             1.0200, 0.9900, 0.9800, 2.0300, 1.7400,
             1.4400, 1.3200, 1.2200, 1.1800, 1.1700,
             1.1700, 1.1600, 1.1500, 1.1700, 1.2500,
             1.2600, 1.2200, 1.2000, 1.1600, 1.1400,
             1.1200, 2.1600, 1.9100, 1.6200, 1.4500,
             1.3400, 1.3000, 1.2700, 1.2500, 1.2500,
             1.2800, 1.3400, 1.4800, 1.4400, 1.4100,
             1.4000, 1.3600, 1.3300, 1.3100, 2.3500,
             1.9800, 1.6900, 1.6500, 1.6500, 1.6400,
             1.6300, 1.6200, 1.8500, 1.6100, 1.5900,
             1.5900, 1.5800, 1.5700, 1.5600, 1.5600,
             1.5600, 1.4400, 1.3400, 1.3000, 1.2800,
             1.2600, 1.2700, 1.3000, 1.3400, 1.4900,
             1.4800, 1.4700, 1.4600, 1.4600, 1.4500,
             1.0000, 1.0000, 1.0000, 1.0000, 1.6500,
             1.0000, 1.4200, 1.0000, 1.0000, 1.0000,
             1.0000, 1.0000, 1.0000, 1.0000, 0.8000,
             1.0000, 1.0000, 1.0000 ]
        return gradius[z]/ 0.529177

class LiuCBond:
    """Reaction coordinate of Liu algorithm is a
       bond length
      Treatment 1 and 2 are also considered

    atomA, atomB are the numbers for the atoms whose bond length is the reaction coordinate
    treat1, treat2 are treatment 1 and 2
    damp_factor dams update, corresponds to beta in the paper
    """
    def __init__(self, atomA, atomB, damp_factor = 0.8):
        # lists in python start with o element, so the n'th atom is at position n-1 !!
        self.atoma = atomA-1
        self.atomb = atomB-1
        self.beta   = damp_factor

    def relaxation(self, old, new, forces, t1):
        # positions of atoms a and b are reset to old values
        pa = old[self.atoma]
        pb = old[self.atomb]
        d = pa - pb
        p = sqrt(np.dot(d, d))
        print "relaxation to value", p
        self.adjust_positions( new, forces, p, t1)

    def update(self, old, new, forces):
        # positions of atoms a and b are set to a new value
        # but not the one the minimizator calculated
        # change in coordinates is given back for possible
        # treatment 2
        pa = old[self.atoma]
        pb = old[self.atomb]
        d = pa - pb
        p = sqrt(np.dot(d, d))
        qa = new[self.atoma].copy()
        qb = new[self.atomb].copy()
        d = qa -qb
        q = sqrt(np.dot(d, d))
        p -= self.beta * (q - p)
        print "update to value", p
        self.adjust_positions( new, forces, p, False)
        return [ new[self.atoma] - qa, new[self.atomb] -qb]

    def adjust_positions(self, new, forces, dt, t1):
        # actually move the atoms a and b in such a way,
        # that there distance is dt, set new (all atoms) value to it
        # treatment 1 may be considerd
        qa = new[self.atoma]
        qb = new[self.atomb]
        d = qa -qb
        q = sqrt(np.dot(d, d))
        qmiddle = 0.5 * (qa + qb)
        d *= 0.5
        qa = qmiddle + dt / q * d
        qb = qmiddle - dt / q * d
        if t1:
            fa = sqrt(np.dot(forces[self.atoma], forces[self.atoma]))
            fb = sqrt(np.dot(forces[self.atomb], forces[self.atomb]))
            lam = (fa - fb) / (fa + fb)
            qa += lam * (1 - dt / q) * d
            qb += lam * (1 - dt / q) * d
        new[self.atoma] = qa
        new[self.atomb] = qb

    def centerofnetwork(self):
        # the liualgorithm has to know around which atoms to build a network
        return [self.atoma , self.atomb]


class netinfo:
     def __init__(self, number, dist_to_center = None, shell = 0,
                  multiplicity = 0, nextN = [], dist_to_next = [], radius = None):
         '''Gives the info of a atomicnetwork for one (current) atom, related to a center
         the variables stored there are:
            number: the number of the atom, consider that the lists start with 0, but
                    the first atom is 1, but number can be used to get place in list back
            dist_to_center: gives the distance in units of atom posititions from the current
                            atom to the center of the network
            shell: integer number, gives the number of atoms over which one has to go (at least) to
                   reach the center, starting from the current atom
            nextN: direct neighbors of current atom (of shell -1), current atom has bond to them and they are
                   all nearer to center atom than current atom
            dist_to_next: distance to the direct neigbors, in the same order given
            multiplicity: dimension of nextN, how many direct neigbohrs are there
            radius: radius of the current atom, just to store it and have easier access
         '''
         self.number = number
         self.dist_to_center = dist_to_center
         self.shell = shell
         self.multiplicity = multiplicity
         self.nextN = nextN
         self.dist_to_next = dist_to_next
         self.radius = radius
