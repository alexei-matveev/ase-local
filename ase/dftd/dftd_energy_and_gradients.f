CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Reimplementation of energy and gradient subroutines
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine e_disp(max_elem,maxc,n,ngroup,xyz,tvec,iz,ilist,imat,
     .                  cn,
     .                  mxc,c6ab,r2r4,r0ab,s6,s8,rs6,rs8,alp6,alp8,
     .                  version,
     .                  disp)
      !
      implicit none
      !
      !------------------------------------------------------------------
      !
      logical, INTENT(in) :: imat(0:ngroup,0:ngroup)
      !
      integer, INTENT(in) :: max_elem, maxc, mxc(max_elem)
      integer, INTENT(in) :: n, ngroup
      integer, INTENT(in) :: version
      integer, INTENT(in) :: iz(n)
      integer, INTENT(in) :: ilist(n)
      !
      real*8, INTENT(in)  :: tvec(3)
      real*8, INTENT(in)  :: xyz(3,n)
      real*8, INTENT(in)  :: s6,s8,rs6,rs8,alp6,alp8
      real*8, INTENT(in)  :: cn(n)
      real*8, INTENT(in)  :: r0ab(max_elem,max_elem),r2r4(max_elem)
      real*8, INTENT(in)  :: c6ab(max_elem,max_elem,maxc,maxc,3)
      !
      !------------------------------------------------------------------
      !
      real*8, INTENT(out) :: disp
      !
      !------------------------------------------------------------------
      !
      integer             :: iat, jat
      !
      real*8              :: dx, dy, dz
      real*8              :: r6, r8
      real*8              :: r, r2, R42, R0, rr
      real*8              :: C6, C8
      real*8              :: tmp1, tmp2
      real*8              :: damp6, damp8
      real*8              :: e6, e8
      !
      !------------------------------------------------------------------
      !
      disp = 0.0d0
      !
      IF (version.eq.2) THEN
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
            ! Screening for identical or to distant atoms
            IF (r2 .lt. 0.0001d0 .or. r2 .gt. 10000.0d0) cycle
            !
            r    = sqrt(r2)
            r6   = r2 * r2 * r2
            !
            R0   = r0ab(iz(jat),iz(iat))*rs6
            !
            C6   = c6ab(iz(jat),iz(iat),1,1,1)*s6
            !
            ! Everything is set... start calculation
            damp6 = 1.0d0
     .            + exp(-alp6 * (r / R0 - 1.0d0))
            !
            ! Energy
            disp = disp - c6 / (damp6 * r6) * 0.5d0
            !
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          END DO
        END DO
        !
      ELSEIF  (version.eq.3) THEN
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
            ! Screening for identical or to distant atoms
            IF (r2 .lt. 0.0001d0 .or. r2 .gt. 10000.0d0) cycle
            !
            ! get vdW parameters
            R0  = r0ab(iz(iat),iz(jat))
            R42 = r2r4(iz(iat)) * r2r4(iz(jat))
            !
            call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                 cn(iat),cn(jat),C6)
            C8 = 3.0d0*C6*R42
            !
            ! Calculate needed powers of r
            r   = dsqrt(r2)
            rr  = R0 / r
            r6  = r2 * r2 * r2
            r8  = r6 * r2
            !
            ! Calculate damping terms
            tmp1  = rs6 * rr
            tmp2  = rs8 * rr
            damp6 = 1.0d0 / (1.0d0 + 6.0d0 * tmp1**alp6)
            damp8 = 1.0d0 / (1.0d0 + 6.0d0 * tmp2**alp8)
            !
            ! Calculate e6 an e8
            e6 = c6*damp6/r6
            e8 = c8*damp8/r8
            !
            ! Add contribution of this pair
            disp = disp - 0.5d0 * (s6 * e6 + s8 * e8)
            !
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ENDDO
        ENDDO
      ELSE
        write(*,*) 'version', version,'does not exist'
        stop
      END IF ! version
      !
      end subroutine e_disp

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
      real*8              :: r6, r8, r10
      real*8              :: r, r2, r235, R42, R0, rr
      real*8              :: C6, C8
      real*8              :: t6, t8
      real*8              :: s8, s10
      real*8              :: damp1, damp6, damp8
      real*8              :: tmp1, tmp2, term
      real*8              :: gx1, gy1, gz1
      real*8              :: gdspr, gdspl
      real*8              :: dr2_dxyz(3)
      real*8              :: gC6(3)
      real*8              :: shifted_xyz(3,n)
      !
      real*8, parameter   :: step = 2.0d-5
      !
      !------------------------------------------------------------------
      !
      ! Initialize
      gdsp     = 0.0
      dispgrad = 0.0
      !
      s8 =s18
      s10=s18
      !
      IF(numgrad) THEN
        !
        ! DFT-D2/3 NUMERICAL GRADIENTS
        !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !
        ! The initial step
        call e_disp(max_elem,maxc,n,ngroup,xyz,tvec,iz,ilist,imat,cn,
     .              mxc,c6ab,r2r4,r0ab,s6,s8,rs6,rs8,alp6,alp8,version,
     .              gdsp)
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
            call e_disp(max_elem,maxc,n,ngroup,shifted_xyz,tvec,iz,
     .                  ilist,imat,cn,mxc,c6ab,r2r4,r0ab,s6,s8,rs6,rs8,
     .                  alp6,alp8,version,
     .                  gdspr)
            !
            shifted_xyz(i123,iat) = shifted_xyz(i123,iat) - 2*step
            !
            call e_disp(max_elem,maxc,n,ngroup,shifted_xyz,tvec,iz,
     .                  ilist,imat,cn,mxc,c6ab,r2r4,r0ab,s6,s8,rs6,rs8,
     .                  alp6,alp8,version,
     .                  gdspl)
            !
            ! Combine to gradient
            dispgrad(i123,iat) = (gdspr - gdspl) / (step + step)
            !
            shifted_xyz(i123,iat) = shifted_xyz(i123,iat) + step
            !
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          END DO
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
            ! Screening for identical or to distant atoms
            IF (r2 .lt. 0.0001d0 .or. r2 .gt. 10000.0d0) cycle
            !
            r235 = r2**3.5
            r    = sqrt(r2)
            r6   = r2 * r2 * r2
            !
            R0   = r0ab(iz(jat),iz(iat))*rs6
            !
            C6   = c6ab(iz(jat),iz(iat),1,1,1)*s6
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
        DO iat = 1, n
          DO jat = 1, n
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! LOOP FOR CALCULATION OF ALL DERIVATIVES OF THE FACTORS
            !  damp6, damp8, 1/r^6, 1/r^8
            !
            ! Screening interaction groups
            IF (.not. imat(ilist(iat),ilist(jat))) cycle
            !
            dx = (xyz(1,iat) - xyz(1,jat) - tvec(1))
            dy = (xyz(2,iat) - xyz(2,jat) - tvec(2))
            dz = (xyz(3,iat) - xyz(3,jat) - tvec(3))
            !
            r2 = dx*dx + dy*dy + dz*dz
            !
            ! Screening for identical or to distant atoms
            IF (r2 .lt. 0.0001d0 .or. r2 .gt. 10000.0d0) cycle
            !
            ! get vdW parameters
            R0  = r0ab(iz(iat),iz(jat))
            R42 = r2r4(iz(iat)) * r2r4(iz(jat))
            !
            ! C6 and C8 values for pair
            call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                 cn(iat),cn(jat),C6)
            C8 = 3.0d0*C6*R42
            !
            ! Calculate needed powers of r
            r   = dsqrt(r2)
            r6  = r2 * r2 * r2
            r8  = r6 * r2
            r10 = r8 * r2
            !
            ! Calculate damping terms
            t6    = (r / (rs6 * R0))**(-alp6)
            t8    = (r / (rs8 * R0))**(-alp8)
            damp6 = 1.d0 / (1.d0 + 6.d0 * t6)
            damp8 = 1.d0 / (1.d0 + 6.d0 * t8)
            !
            dr2_dxyz(1) = 2.0 * dx
            dr2_dxyz(2) = 2.0 * dy
            dr2_dxyz(3) = 2.0 * dz
            !
            DO i123 = 1, 3
              ! Calculate derivatives of damping, R^-6 and R^-8 terms
              term = (-3.D0*damp6*damp6*s6*C6/r8*t6*alp6
     .               +3.D0*damp6*s6*C6/r8
     .               -9.D0*damp8*damp8*s8*C6*R42/r10*t8*alp8
     .               +12.D0*damp8*s8*C6*R42/r10) * dr2_dxyz(i123)
              !
              ! Add to gradients
              dispgrad(i123,iat) = dispgrad(i123,iat) + term
              !
            END DO
            !
            ! Contribution of this pair to gdsp
            gdsp = gdsp
     .           - 0.5d0 * (s6 * damp6 * c6 / r6 + s8 * C8 * damp8 / r8)
            !
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          END DO ! jat
          !
          DO jat = 1, n
            DO kat = 1, n
              ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              ! LOOP FOR CALCULATION OF ALL DERIVATIVES OF THE FACTORS
              !  C6, C8
              !
              ! Screening interaction groups
              IF (.not. imat(ilist(jat),ilist(kat))) cycle
              !
              dx = (xyz(1,jat) - xyz(1,kat) - tvec(1))
              dy = (xyz(2,jat) - xyz(2,kat) - tvec(2))
              dz = (xyz(3,jat) - xyz(3,kat) - tvec(3))
              !
              r2 = dx*dx + dy*dy + dz*dz
              !
              ! Screening for identical or to distant atoms
              IF (r2 .lt. 0.0001d0 .or. r2 .gt. 10000.0d0) cycle
              !
              ! Calculate vdW parameters
              R0  = r0ab(iz(jat),iz(kat))
              R42 = r2r4(iz(jat)) * r2r4(iz(kat))
              !
              call c6_grads(max_elem,maxc,n,cn,dcn,iz,c6ab,mxc,
     .                      jat,kat,iat,gC6)
              !
              ! Calculate needed powers of r
              r   = dsqrt(r2)
              r6  = r2 * r2 * r2
              r8  = r6 * r2
              !
              ! Calculate damping terms
              rr = R0 / r
              t6 = (rr * rs6)**alp6
              t8 = (rr * rs8)**alp8
              damp6 =1.d0 / (1.d0 + 6.d0 * t6)
              damp8 =1.d0 / (1.d0 + 6.d0 * t8)
              !
              DO i123 = 1, 3
                ! Calculate derivatives of C6 and C8 terms
                term = -(damp6 * s6 / r6 + 3.D0 * damp8 * s8 * R42 / r8)
     .                 * gC6(i123)
                !
                ! Add to gradients
                dispgrad(i123,iat) = dispgrad(i123,iat) + 0.5d0 * term
              END DO
              ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            END DO ! kat
          END DO ! jat
        END DO ! iat
        !
      END IF
      END IF ! analytical / numerical gradients
      !
      END SUBROUTINE g_disp
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE getc6(maxc,max_elem,c6ab,mxc,iat,jat,nci,ncj,c6)
      implicit none
      integer maxc,max_elem
      integer iat,jat,i,j,mxc(max_elem)
      real*8  nci,ncj,c6,c6mem
      real*8  c6ab(max_elem,max_elem,maxc,maxc,3)
c the exponential is sensitive to numerics
c when nci or ncj is much larger than cn1/cn2
      real*8  cn1,cn2,r,rsum,csum,tmp1
      include 'param.inc'

      c6mem=-1.d+99
      rsum=0.0
      csum=0.0
      c6  =0.0
      do i=1,mxc(iat)
      do j=1,mxc(jat)
         c6=c6ab(iat,jat,i,j,1)
         if(c6.gt.0)then
            c6mem=c6
            cn1=c6ab(iat,jat,i,j,2)
            cn2=c6ab(iat,jat,i,j,3)
c distance
            r=(cn1-nci)**2+(cn2-ncj)**2
            tmp1=exp(k3*r)
            rsum=rsum+tmp1
            csum=csum+tmp1*c6
         endif
      enddo
      enddo

      if(rsum.gt.0)then
         c6=csum/rsum
      else
         c6=c6mem
      endif

      END SUBROUTINE getc6
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE limit(iat,jat,iadr,jadr)
      integer iat,jat
      iadr=1
      jadr=1
      i=100
 10   if(iat.gt.100) then
         iat=iat-100
         iadr=iadr+1
         goto 10
      endif

      i=100
 20   if(jat.gt.100) then
         jat=jat-100
         jadr=jadr+1
         goto 20
      endif

      END SUBROUTINE limit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE c6_grads(max_elem,maxc,n,cn,dcn,
     .           iz,c6ab,mxc,iat,jat,kat,anag)
      implicit   none
      include    'param.inc'
      integer    n,iz(*),max_elem,maxc,iat,jat,kat,mxc(max_elem)
      real*8     cn(n),anag(3)
      real*8     c6ab(max_elem,max_elem,maxc,maxc,3)
      real*8     term1,term2,term3
      real*8     dterm2(3),dterm3(3),dcn(3,n,n)
      real*8     zaehler,nenner,dzaehler(3),dnenner(3)
      integer    i,j,k
      !
      dterm2=dcn(:,kat,iat)
      dterm3=dcn(:,kat,jat)
      !
      zaehler=0.0d0
      nenner=0.0d0
      dzaehler=0.0d0
      dnenner=0.0d0
      do i=1,mxc(iz(iat))
        do j=1,mxc(iz(jat))
          term3=c6ab(iz(iat),iz(jat),i,j,3)-cn(jat)
          term2=c6ab(iz(iat),iz(jat),i,j,2)-cn(iat)
          term1=exp(k3*(term2*term2+term3*term3))
          zaehler=zaehler+c6ab(iz(iat),iz(jat),i,j,1)*term1
          nenner=nenner+term1
          do k=1,3
            dzaehler(k)=dzaehler(k)+c6ab(iz(iat),iz(jat),i,j,1)*term1*k3
     .                  *2.0d0*(term2*dterm2(k)+term3*dterm3(k))
             dnenner(k)=dnenner(k)+term1*k3
     .                  *2.0d0*(term2*dterm2(k)+term3*dterm3(k))
          enddo
        enddo
      enddo
      if (abs(nenner).gt.1.d-99) then
        do k=1,3
         anag(k)=(dzaehler(k)*nenner-(dnenner(k)*zaehler))
     .           /(nenner*nenner)
        enddo
      else
        anag=0.0d0
      endif
      END SUBROUTINE c6_grads
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

