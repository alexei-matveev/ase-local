#!/usr/bin/python

__all__ = ["gxread", "gxwrite", "is_dummy"]

from sys import stdout
from warnings import warn

# to wrap the 3D vectors as V([x,y,z]):
from vector import Vector as V
# from numpy import array as V
# from derivatives import DerivVector as V

# Dummy atoms have to be distinguished from real  ones.
# In PG dummy atoms historically have atomic number 99.
# Within the python sources in this repo, let the dummy
# atoms to be identified by this:
DUMMY = 0

class GxAtom:
    """The record type to hold information from the gxfile

    Attributes:
    charge :: Atomic number
    pos    :: V([float,float,float]) -- postition of the atom
    id     :: int -- zero-based id used to define connectivities

    bond   :: ((int,int),int) -- ids of atoms defining the bond and bond id
    angle  :: ((int,int,int),int) -- ids of atoms defining the angle and angle id
    plane  :: ((int,int,int,int),int) -- ids of atoms defining the plane and plane id
              The first int in the last three fields is the
              id of the atom itself.

              bond[1]  -- zero based index of variable bond
              angle[1] -- zero based index of variable bond angle
              plane[1] -- zero based index of variable dihedral angle

              These fields are -1 for if the corresponding coordinate
              is fixed, the indices for symmetry-equivalent variables
              are identical

    grad   :: V([float,float,float]) -- cartesian energy gradient wrt
              atomic position. Note that grad == (-force) !

    num    :: int -- zero-based index of real atoms excluding dummies
              Dont use it for indexing into array of atoms, the array
              does include dummies. Use |id| instead.
    """
    def __str__(self):
        # dir(obj) returns a list of all object members,
        # filter out those starting with '_':
        records = filter( lambda x: x[0] != '_', dir(self) )
        out = ""
        for rec in records:
            if out: out += ", " # not in the first iteration
            out += rec + "=" + str(getattr(self,rec))
        return out
    #end def

    def _tostr_geom(self):
        """Return a line of gxfile in following format:
        ^ 1.00        -5.184294677816         2.993153927796         0.427911418660   5  18     6   1   3     4   7   0
        """

        # extract the fields from records:

        if is_dummy(self):
            # gxfile convention for dummy atoms:
            charge = 99.0
        else:
            charge = self.charge

        # conversion to floats is for the case when vector components are DerivVars:
        pos    = tuple( map(float, self.pos) )
        unique = self.unique
        id     = self.id

        # connectivities (three fields):
        conn   = self.plane[0][1:]

        # variable indices (three fields):
        vars   = (self.bond[1], self.angle[1], self.plane[1])

        #
        # WARNING: all the numerical ids are rebased to 1,
        #          so that those (invalid) -1s become 0.
        #

        unique += 1
        id     += 1
        conn    = tuple([ i + 1 for i in conn])
        vars    = tuple([ i + 1 for i in vars])

        # concatenate tuples, the first one is a singleton:
        fields = (charge,) + pos + (unique, id) + conn + vars

        return ( "%5.2f %22.12f %22.12f %22.12f %3i %3i   %3i %3i %3i   %3i %3i %3i\n" \
                     % fields )
    #end def

    def _tostr_grad(self):
        """Return a line of gxfile containg cartesian gradients:
        ^   20       -0.000000596917  -0.000000344630   0.000002806462
        """
        #
        # WARNING: all the numerical ids are rebased to 1,
        #          so that those (invalid) -1s become 0.
        #

        num  = self.num + 1

        # conversion to floats is for the case when vector components are DerivVars:
        grad = tuple( map(float, self.grad) )

        # concatenate tuples, the first one is a singleton:
        fields = (num,) + grad
        return ( "%5i      %16.12f %16.12f %16.12f\n" % fields )
    #end def

    def _parse_geom(self,line):
        """Parse a line of gxfile in following format:
        ^ 1.00        -5.184294677816         2.993153927796         0.427911418660   5  18     6   1   3     4   7   0
        """
        fields = line.split()
        # print(fields)

        # nuclear charge, may even be fractional, for dummy atoms Z=99:
        charge = int( float( fields[0] ) )

        # gxfile convention for dummy atoms:
        if charge == 99.0: charge = DUMMY

        self.charge = charge

        # exit on negative charge:
        if charge < 0:
            # this is the only usefull field in this line,
            # but there is no use for it here ...
            return

        # 3D vector of atomic position:
        self.pos    = [ float(i) for i in fields[1:4] ] # yes, three fields 1 <= xyz < 4

        # Vrap 3D vectors in vector type, see "import ... as V" at the top:
        self.pos = V(self.pos)

        #
        # WARNING: all the numerical ids are rebased to zero,
        #          so that those (invalid) zeros become -1.
        #          This is to use them for indexing into arrays/lists.
        #

        # uniqie atom id, indexing groups of symmetry equivalent atoms:
        self.unique = int( fields[4] ) - 1

        # numeric ID used to define connectiviites:
        self.id     = int( fields[5] ) - 1

        # connectivities:
        conn        = [ int(i) - 1 for i in fields[6:9] ] # three fields

        # integer lables for "variable" internal coordinates,
        # zero for "constrained" internal coordinates:
        vars        = [ int(i) - 1 for i in fields[9:12] ] # three fields

        bond        = ( self.id, conn[0] )
        angle       = ( self.id, conn[0], conn[1] )
        plane       = ( self.id, conn[0], conn[1], conn[2] )

        self.bond   = ( bond , vars[0] )
        self.angle  = ( angle, vars[1] )
        self.plane  = ( plane, vars[2] )
    #enddef

    def _parse_grad(self,line):
        """Parse a line of gxfile containg cartesian gradients:
        ^   20       -0.000000596917  -0.000000344630   0.000002806462
        """
        fields = line.split()
        self.grad = map( float, fields[1:4] )

        # Vrap 3D vectors in vector type, see "import ... as V" at the top:
        self.grad = V(self.grad)
    #enddef
#end class

# tiny helper function:
def is_dummy(atom):
    return atom.charge == DUMMY

# for private use only:
class EOF(Exception): pass

def gxread(file):
    """Read gxfile in following format,
    (first and last lines are not parts of the file, rather for field width extimation):
    #12345 1234567890123456789012 1234567890123456789012 1234567890123456789012 123 123   123 123 123   123 123 123
    ^92.00         0.000000000000         0.000000000000         0.000000000000   1   1     0   0   0     0   0   0
    ^99.00         1.000000000000         0.000000000000         0.000000000000   0   2     1   0   0     0   0   0
    ^ 8.00         0.000000000000         0.000000000000         3.354241540616   2   3     1   2   0     1   0   0
    ^ 8.00         0.000000000000         0.000000000000        -3.354241540616   2   4     1   2   3     1   0   0
    ^ ...
    ^ 1.00        -5.184294677816         2.993153927796         0.427911418660   5  18     6   1   3     4   7   0
    ^ 1.00         0.000000000000        -5.986307855591         0.427911418660   5  19     7   1   3     4   7   0
    ^ 1.00         5.184294677816        -2.993153927796        -0.427911418660   5  20     8   1   4     4   7   0
    ^ 1.00        -5.184294677816        -2.993153927796        -0.427911418660   5  21     9   1   4     4   7   0
    ^ 1.00         0.000000000000         5.986307855591        -0.427911418660   5  22    10   1   4     4   7   0
    ^  -1.0         0.000000000000         0.000000000000         0.000000000000   0   0     0   0   0     0   0   0    0
    ^     -28546.827681881085     -28546.827681881085
    ^    1        0.000000000000   0.000000000000   0.000000000000
    ^    2        0.000000000000   0.000000000000  -0.000493202624
    ^   ...
    ^   20       -0.000000596917  -0.000000344630   0.000002806462
    ^   21        0.000000000000   0.000000689261   0.000002806462
    #12345      1234567890123456 1234567890123456 1234567890123456
     """
    atoms = []

    # iterator over lines in a file:
    lines = open(file)

    #
    # First,  read in geometry:
    #
    numall = 0 # running atom counter
    num = 0 # not counting dummies
    for line in lines:
        # new empty record:
        atom = GxAtom()

        # parse a line, saving fields as records of GxAtom:
        atom._parse_geom(line)

        # exit on negative charge:
        if atom.charge < 0:
            # this is the only usefull field in this line:
            loop = - atom.charge
            break
        #end if

        # otherwise this seems to be a valid atom:
        numall += 1

        # but some of them are dummy atoms:
        if is_dummy(atom):
            atom.num = -1
        else:
            num += 1
            # sequence index, zero-based:
            atom.num = num - 1
        #end if

        # not usual, therefore print warning:
        if atom.id != numall - 1:
            warn("gxread: atom indicies are not sequential")

        # append new atom to the list:
        atoms.append(atom)
    #end for

    #
    # Second, read in energy:
    #
    try: # if the energy/forces section is present ...
        try: # if the energy line is present ...
            fields = lines.next().split()
            energy = float( fields[0] )
        except StopIteration: # is raised by .next() on EOF
            warn("gxread: gxfile does not contain energy")
            energy = None
            raise EOF, "gxread: no energies"
            # or maybe "return atoms, loop, None" right away?
        #end try

        #
        # Third, read in forces:
        #
        for atom in atoms:
            if is_dummy(atom):
                # fake entry with zero gradients:
                atom._parse_grad("0 0.0 0.0 0.0")
            else:
                atom._parse_grad(lines.next())
            #end if
        #end for
    except EOF:
        warn("gxread: gxfile incomplete")
    #end try

# # debug print:
# for atom in atoms:
#   print(atom)
    return atoms, loop, energy
#end def

def gxwrite(file,atoms,loop=1,energy=None):
    """Write gxfile in following format,
    (first and last lines are not parts of the file, rather for field width extimation):
    #12345 1234567890123456789012 1234567890123456789012 1234567890123456789012 123 123   123 123 123   123 123 123
    ^92.00         0.000000000000         0.000000000000         0.000000000000   1   1     0   0   0     0   0   0
    ^99.00         1.000000000000         0.000000000000         0.000000000000   0   2     1   0   0     0   0   0
    ^ 8.00         0.000000000000         0.000000000000         3.354241540616   2   3     1   2   0     1   0   0
    ^ 8.00         0.000000000000         0.000000000000        -3.354241540616   2   4     1   2   3     1   0   0
    ^ ...
    ^ 1.00        -5.184294677816         2.993153927796         0.427911418660   5  18     6   1   3     4   7   0
    ^ 1.00         0.000000000000        -5.986307855591         0.427911418660   5  19     7   1   3     4   7   0
    ^ 1.00         5.184294677816        -2.993153927796        -0.427911418660   5  20     8   1   4     4   7   0
    ^ 1.00        -5.184294677816        -2.993153927796        -0.427911418660   5  21     9   1   4     4   7   0
    ^ 1.00         0.000000000000         5.986307855591        -0.427911418660   5  22    10   1   4     4   7   0
    ^  -1.0         0.000000000000         0.000000000000         0.000000000000   0   0     0   0   0     0   0   0    0
    ^     -28546.827681881085     -28546.827681881085
    # 12345678901234567890123 12345678901234567890123
    ^    1        0.000000000000   0.000000000000   0.000000000000
    ^    2        0.000000000000   0.000000000000  -0.000493202624
    ^   ...
    ^   20       -0.000000596917  -0.000000344630   0.000002806462
    ^   21        0.000000000000   0.000000689261   0.000002806462
    #12345      1234567890123456 1234567890123456 1234567890123456
     """
    # special case for file=="-":
    if file == "-":
        write = stdout.write
    else:
        write = open(file,"w").write
    #end if

    #
    # First, the geometry section ...
    #
    for a in atoms:
        write( a._tostr_geom() )
    #end for

    #
    # ... closed by a line starting with negative number:
    #
    write( "%6.1f %22.12f %22.12f %22.12f   0   0     0   0   0     0   0   0    0\n" \
             % ( -loop, 0.0,0.0,0.0 ) )

    if energy == None: return
    #
    # Second, the energy and gradients (not needed if this is a new geometry):
    #
    write( " %23.12f %23.12f\n"  % (energy,energy) )

    for a in atoms:
        if is_dummy(a): continue
        if not hasattr(a,"grad"): break
        # FIXME: the case when only some of atoms dont have forces
        write( a._tostr_grad() )
    #end for
#end def

# You need to add "set modeline" and eventually "set modelines=5"
# to your ~/.vimrc for this to take effect.
# Dont (accidentally) delete these lines! Unless you do it intentionally ...
# Default options for vim:sw=4:expandtab:smarttab:autoindent:syntax
