CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE d2_energy(zvals, coordraw, tvec, ilist, imat, func,
     .                        n_atom, xyz, n_group,
     .                                   dftd2_energy)
C----------------------------------------------------------------------
      IMPLICIT none
C----------------------------------------------------------------------
      character*80                    :: func
      INTEGER                         :: n_atom, xyz, n_group
      INTEGER                         :: zvals(n_atom)
      INTEGER                         :: ilist(n_atom)
      LOGICAL                         :: imat(0:n_group,0:n_group)
      REAL*8                          :: coordraw(n_atom,xyz)
      REAL*8                          :: tvec(xyz)
      REAL*8                          :: dftd2_energy
C----------------------------------------------------------------------
Cf2py intent(in) zvals
Cf2py intent(in) coordraw
Cf2py intent(in) tvec
Cf2py intent(in) ilist
Cf2py intent(in) imat
Cf2py intent(in) func
Cf2py integer intent(hide),depend(coordraw) :: xyz=shape(coordraw,1)
Cf2py integer intent(hide),depend(coordraw) :: n_atom=shape(coordraw,0)
Cf2py integer intent(hide),depend(imat)     :: n_group=shape(imat,1)-1
Cf2py intent(out) dftd2_energy
C----------------------------------------------------------------------
C     Other Variables...
      REAL*8                          :: t_vec(xyz)
      REAL*8                          :: coords(xyz, n_atom)
      REAL*8                          :: cn(n_atom)
      REAL*8                          :: dcn(xyz, n_atom, n_atom)
C----------------------------------------------------------------------
C     Check input
      IF (xyz .ne. 3) THEN
          write(*,*) 'Something is wrong with coordinate input ...'
          stop
      END IF
C----------------------------------------------------------------------
C     get coords in a.u. and normal storage order
      t_vec  = tvec/0.52917726d0
      coords = transpose(coordraw)/0.52917726d0
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C
      cn  = 0.0d0
      dcn = 0.0d0
C
      call dftd(n_atom, n_group, coords, t_vec, zvals, ilist, imat,
     .           cn, dcn,
     .           func, 2,
     .                                    dftd2_energy)
C
C----------------------------------------------------------------------
      END SUBROUTINE d2_energy
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE d3_energy(zvals, coordraw, tvec, ilist, imat, func,
     .                        cn,
     .                        n_atom, xyz, n_group,
     .                                   dftd3_energy)
C----------------------------------------------------------------------
      IMPLICIT none
C----------------------------------------------------------------------
      character*80                    :: func
      INTEGER                         :: n_atom, xyz, n_group
      INTEGER                         :: zvals(n_atom)
      INTEGER                         :: ilist(n_atom)
      LOGICAL                         :: imat(0:n_group,0:n_group)
      REAL*8                          :: coordraw(n_atom,xyz)
      REAL*8                          :: tvec(xyz)
      REAL*8                          :: cn(n_atom)
      REAL*8                          :: dftd3_energy
C----------------------------------------------------------------------
Cf2py intent(in) zvals
Cf2py intent(in) coordraw
Cf2py intent(in) tvec
Cf2py intent(in) ilist
Cf2py intent(in) imat
Cf2py intent(in) func
Cf2py intent(in) cn
Cf2py integer intent(hide),depend(coordraw) :: xyz=shape(coordraw,1)
Cf2py integer intent(hide),depend(coordraw) :: n_atom=shape(coordraw,0)
Cf2py integer intent(hide),depend(imat)     :: n_group=shape(imat,1)-1
Cf2py intent(out) dftd3_energy
C----------------------------------------------------------------------
C     Other Variables...
      REAL*8                          :: t_vec(xyz)
      REAL*8                          :: coords(xyz, n_atom)
      REAL*8                          :: dcn(xyz, n_atom, n_atom)
C----------------------------------------------------------------------
C     Check input
      IF (xyz .ne. 3) THEN
          write(*,*) 'Something is wrong with coordinate input ...'
          stop
      END IF
C----------------------------------------------------------------------
C     get coords in a.u. and normal storage order
      t_vec  = tvec/0.52917726d0
      coords = transpose(coordraw)/0.52917726d0
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C
      dcn = 0.0d0
C
      call dftd(n_atom, n_group, coords, t_vec, zvals, ilist, imat,
     .           cn, dcn,
     .           func, 3,
     .                                    dftd3_energy)
C
C----------------------------------------------------------------------
      END SUBROUTINE d3_energy
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE d2_gradients(zvals, coordraw, tvec, ilist, imat, func,
     .                        n_atom, xyz, n_group,
     .                                   dftd2_energy, dftd2_gradients)
C----------------------------------------------------------------------
      IMPLICIT none
C----------------------------------------------------------------------
      character*80                    :: func
      INTEGER                         :: n_atom, xyz, n_group
      INTEGER                         :: zvals(n_atom)
      INTEGER                         :: ilist(n_atom)
      LOGICAL                         :: imat(0:n_group,0:n_group)
      REAL*8                          :: coordraw(n_atom,xyz)
      REAL*8                          :: tvec(xyz)
      REAL*8                          :: dftd2_energy
      REAL*8                          :: dftd2_gradients(n_atom,xyz)
C----------------------------------------------------------------------
Cf2py intent(in) zvals
Cf2py intent(in) coordraw
Cf2py intent(in) tvec
Cf2py intent(in) ilist
Cf2py intent(in) imat
Cf2py intent(in) func
Cf2py integer intent(hide),depend(coordraw) :: xyz=shape(coordraw,1)
Cf2py integer intent(hide),depend(coordraw) :: n_atom=shape(coordraw,0)
Cf2py integer intent(hide),depend(imat)     :: n_group=shape(imat,1)-1
Cf2py intent(out) dftd2_energy
Cf2py intent(out) dftd2_gradients
C----------------------------------------------------------------------
C     Other Variables...
      REAL*8                          :: t_vec(xyz)
      REAL*8                          :: coords(xyz, n_atom)
      REAL*8                          :: grads(xyz, n_atom)
      REAL*8                          :: cn(n_atom)
      REAL*8                          :: dcn(xyz, n_atom, n_atom)
C----------------------------------------------------------------------
C     Check input
      IF (xyz .ne. 3) THEN
          write(*,*) 'Something is wrong with coordinate input ...'
          stop
      END IF
C----------------------------------------------------------------------
C     get coords in a.u. and normal storage order
      t_vec  = tvec/0.52917726d0
      coords = transpose(coordraw)/0.52917726d0
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C
      cn  = 0.0d0
      dcn = 0.0d0
C
      call dftd(n_atom, n_group, coords, t_vec, zvals, ilist, imat,
     .           cn, dcn,
     .           func, 2,
     .                                    dftd2_energy, dispgrad=grads)
C
      dftd2_gradients = transpose(grads)
C----------------------------------------------------------------------
      END SUBROUTINE d2_gradients
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE d3_gradients(zvals, coordraw, tvec, ilist, imat, func,
     .                        cn, dcnraw,
     .                        n_atom, xyz, n_group,
     .                                   dftd3_energy, dftd3_gradients)
C----------------------------------------------------------------------
      IMPLICIT none
C----------------------------------------------------------------------
      character*80                    :: func
      INTEGER                         :: n_atom, xyz, n_group
      INTEGER                         :: zvals(n_atom)
      INTEGER                         :: ilist(n_atom)
      LOGICAL                         :: imat(0:n_group,0:n_group)
      REAL*8                          :: coordraw(n_atom,xyz)
      REAL*8                          :: tvec(xyz)
      REAL*8                          :: cn(n_atom)
      REAL*8                          :: dcnraw(n_atom,n_atom,xyz)
      REAL*8                          :: dftd3_energy
      REAL*8                          :: dftd3_gradients(n_atom,xyz)
C----------------------------------------------------------------------
Cf2py intent(in) zvals
Cf2py intent(in) coordraw
Cf2py intent(in) tvec
Cf2py intent(in) ilist
Cf2py intent(in) imat
Cf2py intent(in) func
Cf2py intent(in) cn
Cf2py intent(in) dcnraw
Cf2py integer intent(hide),depend(coordraw) :: xyz=shape(coordraw,1)
Cf2py integer intent(hide),depend(coordraw) :: n_atom=shape(coordraw,0)
Cf2py integer intent(hide),depend(imat)     :: n_group=shape(imat,1)-1
Cf2py intent(out) dftd3_energy
Cf2py intent(out) dftd3_gradients
C----------------------------------------------------------------------
C     Other Variables...
      REAL*8                          :: t_vec(xyz)
      REAL*8                          :: coords(xyz, n_atom)
      REAL*8                          :: dcn(xyz, n_atom, n_atom)
      REAL*8                          :: grads(xyz, n_atom)
C----------------------------------------------------------------------
C     Check input
      IF (xyz .ne. 3) THEN
          write(*,*) 'Something is wrong with coordinate input ...'
          stop
      END IF
C----------------------------------------------------------------------
C     get coords in a.u. and normal storage order
      t_vec  = tvec/0.52917726d0
      coords = transpose(coordraw)/0.52917726d0
      dcn(1,:,:)   = transpose(dcnraw(:,:,1))
      dcn(2,:,:)   = transpose(dcnraw(:,:,2))
      dcn(3,:,:)   = transpose(dcnraw(:,:,3))
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
      dftd3_energy = 0.0
      grads        = 0.0
C
      call dftd(n_atom, n_group, coords, t_vec, zvals, ilist, imat,
     .           cn, dcn,
     .           func, 3,
     .                                    dftd3_energy, dispgrad=grads)
C
      dftd3_gradients = transpose(grads)
C----------------------------------------------------------------------
      END SUBROUTINE d3_gradients
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
