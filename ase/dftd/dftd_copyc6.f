CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c copy from machine generated data statements inside pars.f
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine copy_c6(maxc,max_elem,c6ab,maxci)
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c load old dftd2 parameters
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine d2_paramaters(max_elem,maxc,c6ab,r0ab)
      implicit none  
      integer max_elem,maxc
      real*8 r0ab(max_elem,max_elem)
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3)

      real*8 c6(86),r0(86)
      integer i,j

c the published radii in S.Grimme, J.Comput.Chem. 27, (2006), 1787-1799 (tab 1)
c refer to the following values multiplied by 1.1 (rs6 in this code)
c H, He
         r0(1:86) = (/ 0.91d0,0.92d0,
c Li-Ne
     .      0.75d0,1.28d0,1.35d0,1.32d0,1.27d0,1.22d0,1.17d0,1.13d0,
c Na-Ar
     .      1.04d0,1.24d0,1.49d0,1.56d0,1.55d0,1.53d0,1.49d0,1.45d0,
c K, Ca old
     .      1.35d0,1.34d0,
c Sc-Zn
     .      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,
     .      1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,
c Ga-Kr
     .      1.50d0,1.57d0,1.60d0,1.61d0,1.59d0,1.57d0,
c Rb, Sr
     .      1.48d0,1.46d0,
c Y-Cd
     .      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,
     .      1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,
c In, Sn, Sb, Te, I, Xe
     .      1.52d0,1.64d0,1.71d0,1.72d0,1.72d0,1.71d0,
c Cs,Ba,La,Ce-Lu
     .      1.638d0,1.602d0,1.564d0,1.594d0,1.594d0,1.594d0,1.594d0,
     .      1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,
     .      1.594d0,1.594d0,1.594d0,
c Hf, Ta-Au
     .      1.625d0,1.611d0,1.611d0,1.611d0,1.611d0,1.611d0,1.611d0,
     .      1.611d0,
c Hg,Tl,Pb,Bi,Po,At,Rn
     .      1.598d0,1.805d0,1.767d0,1.725d0,1.823d0,1.810d0,1.749d0/)

       c6(1:86) = (/0.14d0,0.08d0,
     .   1.61d0,1.61d0,3.13d0,1.75d0,1.23d0,0.70d0,0.75d0,0.63d0,
     .   5.71d0,5.71d0,10.79d0,9.23d0,7.84d0,5.57d0,5.07d0,4.61d0,
     .   10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,
     .   10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,16.99d0,
     .   17.10d0,16.37d0,12.64d0,12.47d0,12.01d0,24.67d0,24.67d0,
     .   24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,
     .   24.67d0,24.67d0,24.67d0,37.32d0,38.71d0,38.44d0,31.74d0,
     .   31.50d0,29.99d0,315.275d0,226.994d0,176.252d0,
     .  140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,
     .  140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,
     .  105.112d0,
     .  81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,
     .  57.364d0,57.254d0,63.162d0,63.540d0,55.283d0,57.171d0,56.64d0 /)

      c6ab = -1
      do i=1,86
         do j=1,i
            r0ab(i,j)=(r0(i)+r0(j))/0.52917726d0
            r0ab(j,i)=(r0(i)+r0(j))/0.52917726d0
            c6ab(i,j,1,1,1)=dsqrt(c6(i)*c6(j))
            c6ab(j,i,1,1,1)=dsqrt(c6(i)*c6(j))
         enddo
      enddo

      end
