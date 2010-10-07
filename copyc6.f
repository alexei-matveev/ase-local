CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c copy from machine generated data statements inside pars.f
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine copyc6(maxc,max_elem,c6ab,maxci)
      implicit none
      integer maxc,max_elem,maxci(max_elem),nlines
      real*8  c6ab(max_elem,max_elem,maxc,maxc,3)
c     character*(*) fname
c     character*1  atmp
c     character*80 btmp
c     real*8  x,y,f,cn1,cn2,cmax,xx(10),kekse
c     integer iat,jat,i,n,l,j,k,il,iadr,jadr,nn,kk
      integer iat,jat,iadr,jadr,nn,kk
      include 'pars.inc'
      c6ab=-1
      maxci=0
c read file
      kk=1
      do nn=1,nlines
       iat=int(pars(kk+1))
       jat=int(pars(kk+2))
c       write(*,*)pars(kk+1),pars(kk+2)
c       write(*,*)iat,jat
       call limit(iat,jat,iadr,jadr)
       maxci(iat)=max(maxci(iat),iadr)
       maxci(jat)=max(maxci(jat),jadr)
c      write(*,*)pars(kk),pars(kk+1),pars(kk+2),pars(kk+3),pars(kk+4)
       c6ab(iat,jat,iadr,jadr,1)=pars(kk)  
       c6ab(iat,jat,iadr,jadr,2)=pars(kk+3)
       c6ab(iat,jat,iadr,jadr,3)=pars(kk+4)

       c6ab(jat,iat,jadr,iadr,1)=pars(kk) 
       c6ab(jat,iat,jadr,iadr,2)=pars(kk+4)
       c6ab(jat,iat,jadr,iadr,3)=pars(kk+3)
       kk=(nn*5)+1
      enddo
      end
