CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C compute gradient
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     subroutine gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
c    .                 s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,num,
c    .                 version,echo,g,disp,gnorm,
c    .                 ngroup,ilist,imat)
      subroutine gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,
     .                 s6,s18,rs6,rs8,alp6,alp8,noabc,num,
     .                 version,echo,g,disp,gnorm,
     .                 ngroup,ilist,imat,cn,dcn2,dcn3,tvec)
      implicit none
      integer n,iz(*),max_elem,maxc,version,mxc(max_elem)
      real*8 xyz(3,*),r0ab(max_elem,max_elem),r2r4(max_elem)
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
c     real*8 g(3,*),s6,s18,rcov(max_elem)
      real*8 g(3,*),s6,s18
c     real*8 rs6,rs8,rs10,alp10,alp8,alp6
      real*8 rs6,rs8,alp8,alp6
      logical noabc,num,echo

      integer iat,jat,i,j,kat
c     real*8 R0,C6,alp,R42,disp,x1,y1,z1,x2,y2,z2,rr,e6abc
      real*8 R0,C6,C8,R42,disp,x1,y1,z1,x2,y2,z2,rr,e6abc
c     real*8 dx,dy,dz,r2,r,r4,r6,r8,r10,r12,t6,t8,t10,damp1
      real*8 dx,dy,dz,r2,r,r4,r6,r8,r10,t6,t8,damp1
c     real*8 damp6,damp8,damp10,e6,e8,e10,e12,gnorm,tmp1
      real*8 damp6,damp8,e6,e8,e10,e12,gnorm,tmp1
      real*8 s10,s8,gC6(3),term,step,dispr,displ,r235,tmp2
c     real*8 cn(n),gx1,gy1,gz1,gx2,gy2,gz2,rthr
      real*8 gx1,gy1,gz1,gx2,gy2,gz2,rthr
c     real*8 dcn2(3,n),dcn3(3,n,n)

c Additional variables for interaction groups
      integer, intent(in) :: ngroup
      integer, intent(in) :: ilist(n)
      logical, intent(in) :: imat(0:ngroup, 0:ngroup)
      real*8,  intent(in) :: cn(n), dcn2(3,n), dcn3(3,n,n)
c translation vector in au
      real*8, intent(in)  :: tvec(3)

c R^2 cut-off
      rthr=2500.

      if(echo)write(*,*)
      if(version.eq.2)then
      if(echo)write(*,*) 'doing analytical gradient O(N^2) ...'
      disp=0
      do iat=1,n-1
         do jat=iat+1,n
c           consider only if interaction defined by interaction groups
            IF (imat(ilist(iat),ilist(jat))) THEN
               R0=r0ab(iz(jat),iz(iat))*rs6
               dx=(xyz(1,iat)-xyz(1,jat)+tvec(1))
               dy=(xyz(2,iat)-xyz(2,jat)+tvec(2))
               dz=(xyz(3,iat)-xyz(3,jat)+tvec(3))
               r2  =dx*dx+dy*dy+dz*dz
c              if(r2.gt.rthr) cycle
               r235=r2**3.5
               r   =sqrt(r2)
               damp6=exp(-alp6*(r/R0-1.0d0))
               damp1=1.+damp6
               c6=c6ab(iz(jat),iz(iat),1,1,1)*s6
               tmp1=damp6/(damp1*damp1*r235*R0)
               tmp2=6./(damp1*r*r235)
               gx1=alp6* dx*tmp1-tmp2*dx
               gx2=alp6*(-dx)*tmp1+tmp2*dx
               gy1=alp6* dy*tmp1-tmp2*dy
               gy2=alp6*(-dy)*tmp1+tmp2*dy
               gz1=alp6* dz*tmp1-tmp2*dz
               gz2=alp6*(-dz)*tmp1+tmp2*dz
               g(1,iat)=g(1,iat)-gx1*c6
               g(2,iat)=g(2,iat)-gy1*c6
               g(3,iat)=g(3,iat)-gz1*c6
               g(1,jat)=g(1,jat)-gx2*c6
               g(2,jat)=g(2,jat)-gy2*c6
               g(3,jat)=g(3,jat)-gz2*c6
               disp=disp+c6*(1./damp1)/r2**3
            END IF
         enddo
      enddo
      disp=-disp
      goto 999
      endif

      if(num) then
      if(echo)write(*,*) 'doing numerical gradient O(N^3) ...'
      call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,!rcov,
     .     rs6,rs8,alp6,alp8,version,noabc,
     .     e6,e8,e10,e12,e6abc,ngroup,ilist,imat,cn,tvec)

      disp=-s6*e6-s18*e8-s6*e6abc

      step=2.d-5

      do i=1,n
      do j=1,3
      xyz(j,i)=xyz(j,i)+step
      call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,!rcov,
     .     rs6,rs8,alp6,alp8,version,noabc,
     .     e6,e8,e10,e12,e6abc,ngroup,ilist,imat,cn,tvec)
      dispr=-s6*e6-s18*e8-s6*e6abc
      xyz(j,i)=xyz(j,i)-2*step
      call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,!rcov,
     .     rs6,rs8,alp6,alp8,version,noabc,
     .     e6,e8,e10,e12,e6abc,ngroup,ilist,imat,cn,tvec)
      displ=-s6*e6-s18*e8-s6*e6abc
      g(j,i)=0.5*(dispr-displ)/step
      xyz(j,i)=xyz(j,i)+step
      enddo
      enddo

      else

      if(echo)write(*,*) 'doing analytical gradient O(N^3) ...'
c     call ncoord (n,rcov,iz,xyz,cn)
c precompute for analytical part
c Skipping this, as cn is provided from outside
c     call ncoorda(n,rcov,iz,xyz,cn,dcn2,dcn3)

      s8 =s18
      s10=s18

      disp=0

      do iat=1,n-1
         x1=xyz(1,iat)
         y1=xyz(2,iat)
         z1=xyz(3,iat)
         do jat=iat+1,n
            if(iat.eq.jat) cycle
c           consider only if interaction defined by interaction groups
            IF (.not. imat(ilist(iat),ilist(jat))) cycle
            x2 = xyz(1,jat) + tvec(1)
            y2 = xyz(2,jat) + tvec(2)
            z2 = xyz(3,jat) + tvec(3)
            R0=r0ab(iz(jat),iz(iat))
c stored as sqrt
            r42=r2r4(iz(iat))*r2r4(iz(jat))
            call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                    cn(iat),cn(jat),C6)
            C8 = 3.0d0*C6*r42

c dC6(iat,jat)/dxyz_iat
c     call grdc6iji(max_elem,maxc,n,xyz,cn,rcov,iz,c6ab,mxc,
c    .              iat,jat,iat,gC6)
c analytically
c     call anagrdc6(max_elem,maxc,n,xyz,cn,dcn2,dcn3,rcov,iz,c6ab,mxc,
c    .              iat,jat,iat,gC6)
      call anagrdc6(max_elem,maxc,n,cn,dcn2,dcn3,iz,c6ab,mxc,
     .              iat,jat,iat,gC6)


      dx = (x1-x2)**2
      dy = (y1-y2)**2
      dz = (z1-z2)**2
      r2 = dx+dy+dz
c     if(r2.gt.rthr) cycle
      r = dsqrt(r2)
      t6 = (r/(rs6*R0))**(-alp6)
      damp6 =1.d0/( 1.d0+6.d0*t6 )
      t8 = (r/(rs8*R0))**(-alp8)
      damp8 =1.d0/( 1.d0+6.d0*t8 )
c     t10 = (r/(rs10*R0))**(-alp10)
c     damp10=1.d0/( 1.d0+6.d0*t10 )

      r4 = r2**2
      r6 = r2*r4
      r8 = r4**2
      r10 = r4*r6
c     r12 = r6**2

      dx = 2.D0*x1-2.D0*x2

      term = -3.D0*damp6*damp6*s6*C6/r8*t6*alp6*dx
     &       -1.D0*damp6*s6*gC6(1)/r6
     &       +3.D0*damp6*s6*C6/r8*dx
     &       -9.D0*damp8*damp8*s8*C6*R42/r10*t8*alp8*dx
     &       -3.D0*damp8*s8*gC6(1)*R42/r8
     &       +12.D0*damp8*s8*C6*R42/r10*dx
c    &       -11.025D0*gC6(1)*R42**2*damp10*s10/r10
c    &       -33.075D0*C6*R42**2*damp10*damp10*s10/r12*t10*alp10*dx
c    &       +55.125D0*C6*R42**2*damp10*s10/r12*dx
      g(1,iat)=g(1,iat)+term
      g(1,jat)=g(1,jat)-term

      dy = 2.D0*y1-2.D0*y2

      term = -3.D0*damp6*damp6*s6*C6/r8*t6*alp6*dy
     &       -1.D0*damp6*s6*gC6(2)/r6
     &       +3.D0*damp6*s6*C6/r8*dy
     &       -9.D0*damp8*damp8*s8*C6*R42/r10*t8*alp8*dy
     &       -3.D0*damp8*s8*gC6(2)*R42/r8
     &       +12.D0*damp8*s8*C6*R42/r10*dy
c    &       -11.025D0*gC6(2)*R42**2*damp10*s10/r10
c    &       -33.075D0*C6*R42**2*damp10*damp10*s10/r12*t10*alp10*dy
c    &       +55.125D0*C6*R42**2*damp10*s10/r12*dy
      g(2,iat)=g(2,iat)+term
      g(2,jat)=g(2,jat)-term

      dz = 2.D0*z1-2.D0*z2

      term = -3.D0*damp6*damp6*s6*C6/r8*t6*alp6*dz
     &       -1.D0*damp6*s6*gC6(3)/r6
     &       +3.D0*damp6*s6*C6/r8*dz
     &       -9.D0*damp8*damp8*s8*C6*R42/r10*t8*alp8*dz
     &       -3.D0*damp8*s8*gC6(3)*R42/r8
     &       +12.D0*damp8*s8*C6*R42/r10*dz
c    &       -11.025D0*gC6(3)*R42**2*damp10*s10/r10
c    &       -33.075D0*C6*R42**2*damp10*damp10*s10/r12*t10*alp10*dz
c    &       +55.125D0*C6*R42**2*damp10*s10/r12*dz
      g(3,iat)=g(3,iat)+term
      g(3,jat)=g(3,jat)-term

c     term = -1.D0/(1.D0+6.D0*t6)*s6*C6/r6
c    &       -3.D0/(1.D0+6.D0*t8)*s8*C6*R42/r8
      term = -s8*C8*damp8/r8-s6*damp6*c6/r6

c    &       -11.025D0*C6*R42**2/(1.D0+6.D0*t10)*s10/r10
c     if(iat.gt.jat)then ! OMIT THIS! now doing distinct pairs
         disp=disp+term
c     endif

         enddo
         do jat=2,n
         if(iat.eq.jat) cycle
         x1=xyz(1,jat)+ tvec(1)
         y1=xyz(2,jat)+ tvec(2)
         z1=xyz(3,jat)+ tvec(3)
            do kat=1,jat-1
            if(iat.eq.kat) cycle
            IF (.not. imat(ilist(jat),ilist(kat))) cycle
            x2 = xyz(1,kat)
            y2 = xyz(2,kat)
            z2 = xyz(3,kat)
            R0=r0ab(iz(kat),iz(jat))
            R42=r2r4(iz(jat))*r2r4(iz(kat))

c dC6(jat,kat)/dxyz_iat
c     call grdc6ijk(max_elem,maxc,n,xyz,cn,rcov,iz,c6ab,mxc,
c    .              jat,kat,iat,gC6)
c analytically
c     call anagrdc6(max_elem,maxc,n,xyz,cn,dcn2,dcn3,rcov,iz,c6ab,mxc,
c    .              jat,kat,iat,gC6)
      call anagrdc6(max_elem,maxc,n,cn,dcn2,dcn3,iz,c6ab,mxc,
     .              jat,kat,iat,gC6)

      dx = (x1-x2)**2
      dy = (y1-y2)**2
      dz = (z1-z2)**2
      r2 = dx+dy+dz
      r = dsqrt(r2)
      rr=R0/r
      t6 = (rr*rs6)**alp6
      damp6 =1.d0/( 1.d0+6.d0*t6 )
      t8 = (rr*rs8)**alp8
      damp8 =1.d0/( 1.d0+6.d0*t8 )
c     t10 = (rr*rs10)**alp10
c     damp10=1.d0/( 1.d0+6.d0*t10 )

      r4 = r2**2
      r6 = r2*r4
      r8 = r4**2
c     r10 = r4*r6
c     r12 = r6**2

      term = -1.D0*damp6*s6*gC6(1)/r6
     &       -3.D0*damp8*s8*gC6(1)*R42/r8
c    &       -11.025D0*gC6(1)*R42**2*damp10*s10/r10
        g(1,iat)=g(1,iat)+term
        g(1,jat)=g(1,jat)-term

      term = -1.D0*damp6*s6*gC6(2)/r6
     &       -3.D0*damp8*s8*gC6(2)*R42/r8
c    &       -11.025D0*gC6(2)*R42**2*damp10*s10/r10
        g(2,iat)=g(2,iat)+term
        g(2,jat)=g(2,jat)-term

      term = -1.D0*damp6*s6*gC6(3)/r6
     &       -3.D0*damp8*s8*gC6(3)*R42/r8
c    &       -11.025D0*gC6(3)*R42**2*damp10*s10/r10
        g(3,iat)=g(3,iat)+term
        g(3,jat)=g(3,jat)-term

            enddo
         enddo

      enddo

      endif

 999  continue
c     do i=1,n
c        write(*,'(83F17.12)') g(1:3,i)
c     enddo
      gnorm=sum(abs(g(1:3,1:n)))
      if(echo)then
      write(*,*)
      write(*,*)'|G|=',gnorm
      endif

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c ncoord derivative
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine ncoorda(natoms,rcov,iz,xyz,cn,dcn2,dcn3)
      implicit none
      include 'param.inc'
c     integer iz(*),natoms,i,max_elem
      integer iz(*),natoms,i
      real*8 xyz(3,*),cn(*),dcn2(3,natoms),rcov(94)
      real*8 dcn3(3,natoms,natoms)
      integer iat
c     real*8 dx,dy,dz,r,damp,xn,rr,rrr,rco,tmp1,tmp2,tmp3
      real*8 dx,dy,dz,r,xn,rr,rrr,rco,tmp1,tmp2,tmp3

      dcn2=0.0d0
      dcn3=0.0d0
      do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            rco=rcov(iz(i))+rcov(iz(iat))
            rr=rco/r
            rrr=1.0d0/(r*r*r)

            tmp1=exp(-k1*(rr-1.0d0))
            tmp2=1.0d0/(tmp1+1.0d0)
            tmp3=tmp1*tmp2*tmp2*k1*rco*rrr

            xn=xn+tmp2
            dcn3(1,iat,i)=tmp3*dx
            dcn3(2,iat,i)=tmp3*dy
            dcn3(3,iat,i)=tmp3*dz
            dcn2(1,i)=dcn2(1,i)+tmp3*dx
            dcn2(2,i)=dcn2(2,i)+tmp3*dy
            dcn2(3,i)=dcn2(3,i)+tmp3*dz
         endif
      enddo
      cn(i)=xn
      enddo
      dcn2=-1.0d0*dcn2
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c calculate dC6/dr analytically
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     subroutine anagrdc6(max_elem,maxc,n,dum,cn,dcn2,dcn3,
c    .           rcov,iz,c6ab,mxc,iat,jat,kat,anag)
      subroutine anagrdc6(max_elem,maxc,n,cn,dcn2,dcn3,
     .           iz,c6ab,mxc,iat,jat,kat,anag)
      implicit   none
      include    'param.inc'
      integer    n,iz(*),max_elem,maxc,iat,jat,kat,mxc(max_elem)
c     real*8     cn(*),dcn2(3,n),dum(3,*),anag(3),rcov(max_elem)
      real*8     cn(*),dcn2(3,n),anag(3)
      real*8     c6ab(max_elem,max_elem,maxc,maxc,3)
      real*8     term1,term2,term3
      real*8     dterm2(3),dterm3(3),dcn3(3,n,n)
      real*8     zaehler,nenner,dzaehler(3),dnenner(3)
      integer    i,j,k
      if (iat.eq.kat) then
        dterm2=dcn2(:,iat)
        dterm3=dcn3(:,iat,jat)
      else
        dterm2=dcn3(:,kat,iat)
        dterm3=dcn3(:,kat,jat)
      endif
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
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C gradient of C6(iat,jat) wrt to xyz of kat or iat
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     subroutine grdc6iji(max_elem,maxc,n,dum,cn,
c    .           rcov,iz,c6ab,mxc,iat,jat,kat,g)
c     subroutine grdc6iji(max_elem,maxc,n,dum,
c    .           rcov,iz,c6ab,mxc,iat,jat,kat,g)
c     implicit none
c     integer n,iz(*),max_elem,maxc,iat,jat,kat,mxc(max_elem)
c     real*8  cn(*),dum(3,*),g(3),rcov(max_elem)
c     real*8  dum(3,*),g(3),rcov(max_elem)
c     real*8  c6ab(max_elem,max_elem,maxc,maxc,3)

c     real*8 xi,xj,c6r,c6l,st,xxx(2),yi,yj,x0i,x0j
c     real*8 xi,xj,c6r,c6l,st,xxx(2)
c     integer j,jjj(2)

c     jjj(1)=iat
c     jjj(2)=jat

c     st=1.d-5
c     do j=1,3
c       dum(j,kat)=dum(j,kat)+st

c       call ncoord12(n,rcov,iz,dum,jjj,xxx)
c       xi=xxx(1)
c       xj=xxx(2)
c       call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),xi,xj,c6r)

c       dum(j,kat)=dum(j,kat)-st*2.0d0

c       call ncoord12(n,rcov,iz,dum,jjj,xxx)
c       xi=xxx(1)
c       xj=xxx(2)
c       call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),xi,xj,c6l)
c       g(j)=0.5d0*(c6r-c6l)/st

c       dum(j,kat)=dum(j,kat)+st
c     enddo
c
c     end
c

c     subroutine grdc6ijk(max_elem,maxc,n,dum,cn,
c    .           rcov,iz,c6ab,mxc,iat,jat,kat,g)
c     subroutine grdc6ijk(max_elem,maxc,dum,cn,
c    .           rcov,iz,c6ab,mxc,iat,jat,kat,g)
c     implicit none
c     integer n,iz(*),max_elem,maxc,iat,jat,kat,mxc(max_elem)
c     integer iz(*),max_elem,maxc,iat,jat,kat,mxc(max_elem)
c     real*8  cn(*),dum(3,*),g(3),rcov(max_elem)
c     real*8  c6ab(max_elem,max_elem,maxc,maxc,3)

c     real*8 xi,xj,c6r,c6l,st,xxx(2),yi,yj,x0i,x0j
c     real*8 xi,xj,c6r,c6l,st,yi,yj,x0i,x0j
c     integer j

c x0 is the contribution of kat to the CN of iat and jat
c     call ncoord11(iat,kat,rcov,iz,dum,x0i)
c     call ncoord11(jat,kat,rcov,iz,dum,x0j)

c     st=1.d-5
c     do j=1,3
c       dum(j,kat)=dum(j,kat)+st

c       call ncoord11(iat,kat,rcov,iz,dum,yi)
c       call ncoord11(jat,kat,rcov,iz,dum,yj)
c       xi=cn(iat)-x0i+yi
c       xj=cn(jat)-x0j+yj

c       call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),xi,xj,c6r)

c       dum(j,kat)=dum(j,kat)-st*2.0d0

c       call ncoord11(iat,kat,rcov,iz,dum,yi)
c       call ncoord11(jat,kat,rcov,iz,dum,yj)
c       xi=cn(iat)-x0i+yi
c       xj=cn(jat)-x0j+yj

c       call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),xi,xj,c6l)
c       g(j)=0.5d0*(c6r-c6l)/st

c       dum(j,kat)=dum(j,kat)+st
c     enddo
c
c     end

