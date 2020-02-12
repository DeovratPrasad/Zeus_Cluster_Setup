! a program to calculate differential emission measure
      PROGRAM stitch
      implicit NONE
! number of true grid points on each processor
      integer, parameter :: nx=64, ny=64, nz=1
      integer, parameter :: np1=8, np2=4, np3=1
      integer, parameter :: incr_min=0, incr_max=100
      real*4, parameter :: x1min=3.086e21,x1max=617.2e21 
      real*4, parameter :: x2min=0.0,x2max=3.1415926535897932
      real*4, parameter :: x3min=0.0,x3max=6.28318530717959
      logical, parameter :: xmhd=.false., xcosmic=.false.

      integer i, j, k, l, m, n, ig, jg, kg, incr
      character*15 :: binfile
      character*2 :: id

      real*4, dimension(nx,ny,nz) :: d, e, v1, v2, v3, b1, b2, b3, ecr

      real*4, dimension(nx*np1,ny*np2,nz*np3) :: dg, eg, v1g, v2g, v3g &
      , b1g, b2g, b3g, ecrg

      real*4 :: dx1, dx2, dx3, x1rat
      real*4 :: x1a(nx*np1+1), x1b(nx*np1+1), x2a(ny*np2+1), x2b(ny*np2+1)

      real*8 :: mu, mui, mue, TkeV, dvol, n_e, n_i, mp, pi

      integer, parameter :: nbins = 1000 
      real*8 :: log10TkeV(nbins), dem(nbins), dlog10TkeV, crate


      mu = 0.62; mue = 1.17; mui=1.0/(1./mu-1./mue); mp=1.67265e-24  
      pi = acos(-1.0)

!      x1rat = ( (x1max-x1min)/(sqrt(x1max*x1min)-x1min) - 1.0 )**(2.0/(nx*np1))

      x1rat = (x1max/x1min)**(1.0/(nx*np1))

      x1a(1) = x1min
      dx1 = (x1max-x1min)*(x1rat-1.0)/(x1rat**(nx*np1)-1.0) 

      do i = 1, np1*nx+1     
        if (i.ne.1) x1a(i) = x1a(i-1)+dx1 
        x1b(i) = x1a(i) + 0.5*dx1
        dx1 = x1rat*dx1
!        write(*,*) x1a(i), x1b(i)
!        write(18,*) i, x1a(i)
      enddo

      dx2 = (x2max-x2min)/(ny*np2) 
  
      do j = 1, np2*ny+1
        x2a(j) = x2min + (j-1)*dx2
        x2b(j) = x2a(j) + 0.5*dx2
      enddo

!      dx3 = (x3max-x3min)/(nz*np3)
 
      id = 'aa'

      do incr = incr_min, incr_max

      do n = 0, np3-1
      do m = 0, np2-1
      do l = 0, np1-1

        write(binfile,"(a3,a2,3i2.2,'.',i3.3)") 'bin',id,l,m,n,incr
        open(unit=42,file=binfile,status='unknown',form='binary' &
        ,convert='big_endian')

        read(42) (((d(i,j,k), i=1,nx), j=1,ny), k=1,nz)
        read(42) (((e(i,j,k), i=1,nx), j=1,ny), k=1,nz)

        dg(nx*l+1:nx*(l+1),ny*m+1:ny*(m+1),nz*n+1:nz*(n+1)) &
        = d(1:nx,1:ny,1:nz)

        eg(nx*l+1:nx*(l+1),ny*m+1:ny*(m+1),nz*n+1:nz*(n+1)) &
        = e(1:nx,1:ny,1:nz)

        close(42)

      enddo
      enddo
      enddo

      dlog10TkeV = 6./nbins
      do l = 1, nbins
        log10TkeV(l) = -4. + (l-0.5)*(2.+4.)/nbins
      enddo

      dem = 0.0

      do i = 1, np1*nx 
      do k = 1, np3*nz 
      do j = 1, np2*ny

        dvol = 2.0*pi*x1b(i)*x1b(i)*(x1a(i+1)-x1a(i))*sin(x2b(j))*(x2a(j+1)-x2a(j))/(np3*nz)

        n_e = dg(i,j,k)/mp/mue
        n_i = dg(i,j,k)/mp/mui
        TkeV = 2./3.*eg(i,j,k)/((n_e+n_i)*1.6022e-9)

        if (TkeV.gt.0.02) then
          crate = 1.0e-22*n_i*n_e*( 8.6e-3*TkeV**(-1.7) &
          + 5.8e-2*TkeV**0.5 + 6.3e-2)
        else if (TkeV.le.0.02.and.TkeV.ge.0.0017235) then
          crate = n_i*n_e*6.72e-22*(TkeV/0.02)**0.6
        else
          crate = n_i*n_e*1.544e-22*(TkeV/0.0017235)**6.0
        endif

        do l = 1, nbins-1
          if ( log10(TkeV).gt.log10TkeV(l) .and. log10(TkeV).le.log10TkeV(l+1) ) & 
             dem(l) = dem(l) + crate*dvol/dlog10TkeV 
        enddo

      enddo
      enddo
      enddo

      do l = 1, nbins-1 
        write(100*(incr+1),2001) l, log10TkeV(l), dem(l)
      enddo

      enddo !incr

2001  format(i7,2e20.7)

      END
