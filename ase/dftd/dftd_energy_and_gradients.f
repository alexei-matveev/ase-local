CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Reimplementation of gradient module
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     subroutine g_disp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,
c    .                 s6,s18,rs6,rs8,alp6,alp8,noabc,num,
c    .                 version,echo,g,disp,gnorm,
c    .                 ngroup,ilist,imat,cn,dcn2,dcn3,tvec)
      subroutine g_disp(max_elem,maxc,n,ngroup,xyz,tvec,iz,ilist,imat,
     .                  cn, dcn,
     .                  mxc,c6ab,r2r4,r0ab,s6,s18,rs6,rs8,alp6,alp8,
     .                  version, numgrad,
     .                  gdsp, dispgrad)
      !
      implicit none
      !
      !------------------------------------------------------------------
      !
      logical, INTENT(in) :: numgrad
      logical, INTENT(in) :: imat(0:ngroup,0:ngroup)
      !
      integer, INTENT(in) :: max_elem, maxc, mxc(max_elem)
      integer, INTENT(in) :: n, ngroup
      integer, INTENT(in) :: version
      integer, INTENT(in) :: iz(n)
      integer, INTENT(in) :: ilist(n)
      !
      real*8, INTENT(in)  :: s6,s18,rs6,rs8,alp6,alp8
      real*8, INTENT(in)  :: tvec(3)
      real*8, INTENT(in)  :: xyz(3,n)
      real*8, INTENT(in)  :: cn(n), dcn(3,n,n)
      real*8, INTENT(in)  :: r0ab(max_elem,max_elem),r2r4(max_elem)
      real*8, INTENT(in)  :: c6ab(max_elem,max_elem,maxc,maxc,3)
      !
      !------------------------------------------------------------------
      !
      real*8, INTENT(out) :: gdsp
      real*8, INTENT(out) :: dispgrad(3,n)
      !
      !------------------------------------------------------------------
      !
      integer             :: iat, jat, kat, i123
      !
      real*8              :: dx, dy, dz
      real*8              :: R0, r, r2, r235
      real*8              :: c6
      real*8              :: damp1, damp6
      real*8              :: tmp1, tmp2
      real*8              :: gx1, gy1, gz1
      real*8              :: e6, e8, e10, e12, e6abc
      real*8              :: gdspr, gdspl
      real*8              :: shifted_xyz(3,n)
      !
      real*8, parameter   :: step = 2.d-5
      !
      !------------------------------------------------------------------
      !
      ! Initialize
      gdsp     = 0.0
      dispgrad = 0.0

      IF(numgrad) THEN
        !
        ! DFT-D2/3 NUMERICAL GRADIENTS
        !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !
        ! The initial step
        call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,
     .             rs6,rs8,alp6,alp8,version,.true.,
     .             e6,e8,e10,e12,e6abc,ngroup,ilist,imat,cn,tvec)
        gdsp = - s6 * e6 - s18 * e8
        !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !
        shifted_xyz = xyz
        !
        DO iat = 1, n
          DO i123 = 1, 3
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !
            shifted_xyz(i123,iat) = shifted_xyz(i123,iat) + step
            !
            call edisp(max_elem,maxc,n,shifted_xyz,iz,c6ab,mxc,
     .                 r2r4,r0ab,
     .                 rs6,rs8,alp6,alp8,version,.true.,
     .                 e6,e8,e10,e12,e6abc,ngroup,ilist,imat,cn,tvec)
            !
            gdspr = - s6 * e6 - s18 * e8
            !
            shifted_xyz(i123,iat) = shifted_xyz(i123,iat) - 2*step
            !
            call edisp(max_elem,maxc,n,shifted_xyz,iz,c6ab,mxc,
     .                 r2r4,r0ab,
     .                 rs6,rs8,alp6,alp8,version,.true.,
     .                 e6,e8,e10,e12,e6abc,ngroup,ilist,imat,cn,tvec)
            !
            gdspl = - s6 * e6 - s18 * e8
            !
            ! Combine to gradient
            dispgrad(i123,iat) = (gdspr - gdspl) / (step + step)
            if (iat.eq.2.and.i123.eq.1) Then
c           print*, 'num', shifted_xyz
            end if
            !
            shifted_xyz(i123,iat) = shifted_xyz(i123,iat) + step

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          END DO
c        print *, ',t', dispgrad(:,:)
        END DO
        !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !
      ELSE
      IF (version.eq.2) THEN
        !
        ! DFT-D2 ANALYTICAL GRADIENTS
        !
        DO iat = 1, n
          DO jat = 1, n
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !
            ! Screening interaction groups
            IF ( .not. imat(ilist(iat),ilist(jat))) cycle
            !
            dx = (xyz(1,iat) - xyz(1,jat) - tvec(1))
            dy = (xyz(2,iat) - xyz(2,jat) - tvec(2))
            dz = (xyz(3,iat) - xyz(3,jat) - tvec(3))
            !
            r2 = dx*dx + dy*dy + dz*dz
            !
            ! Screening for identical atoms
            IF (r2 .lt. 0.01d0) cycle
            !
            r235 = r2**3.5
            r    = sqrt(r2)
            !
            R0   = r0ab(iz(jat),iz(iat))*rs6
            !
            c6   = c6ab(iz(jat),iz(iat),1,1,1)*s6
            !
            ! Everything is set... start calculation
            damp6 = exp(-alp6*(r/R0-1.0d0))
            damp1 = 1.0d0 + damp6
            !
            tmp1  = damp6 / (damp1 * damp1 * r235 * R0)
            tmp2  = 6.0d0 / (damp1 * r * r235)
            !
            ! Now build up gradients
            gx1 = tmp1 *alp6 * dx - tmp2 * dx
            gy1 = tmp1 *alp6 * dy - tmp2 * dy
            gz1 = tmp1 *alp6 * dz - tmp2 * dz
            !
            dispgrad(1,iat) = dispgrad(1,iat) - gx1 * c6
            dispgrad(2,iat) = dispgrad(2,iat) - gy1 * c6
            dispgrad(3,iat) = dispgrad(3,iat) - gz1 * c6
            !
            ! Energy to control results
            gdsp = gdsp - c6 / (damp1 * r2**3) * 0.5d0
            !
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          END DO
        END DO
      ELSEIF (version.eq.3) THEN
        !
        ! DFT-D3 ANALYTICAL GRADIENTS
        !
      END IF
      END IF ! analytical / numerical gradients

      END SUBROUTINE g_disp
