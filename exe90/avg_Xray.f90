      PROGRAM stitch
      implicit NONE
! number of true grid points on each processor
      integer, parameter :: nx=64, ny=64, nz=1
      integer, parameter :: np1=8, np2=4, np3=1
      real*4, parameter :: x1min=3.086e+21, x1max=6.172e+23
      real*4, parameter :: x2min=0.0,x2max=3.1415926535897932
      real*4, parameter :: x3min=0.0,x3max=6.28318530717959
      logical, parameter :: xmhd=.false., xcosmic=.false.

      integer i, j, k, l, m, n, ig, jg, kg, incr, iincr
      character*24 :: binfile, onefile
      character*6 :: twofile
      character*2 :: id

      real*4, dimension(nx,ny,nz) :: d, e, v1, v2, v3, b1, b2, b3, ecr

      real*4, dimension(nx*np1,ny*np2,nz*np3) :: dg, eg, v1g, v2g, v3g &
      , b1g, b2g, b3g, ecrg

      real*4 :: dx1, dx2, dx3, x1rat
      real*4 :: x1a(nx*np1), x1b(nx*np1), x2a(ny*np2), x2b(ny*np2)

      real*8 :: mu, mui, mue, tcool, crate, TkeV, vol, dvol, volx &
      , dvolx, n_e, n_i, mp, pi
      real*4 :: d1D(nx*np1), e1D(nx*np1), d1D_X(nx*np1), e1D_X(nx*np1)

      mu = 0.62; mue = 1.18; mui=1.0/(1./mu-1./mue); mp=1.67265e-24  
      pi = acos(-1.0)

!      x1rat = ( (x1max-x1min)/(sqrt(x1max*x1min)-x1min) - 1.0 )**(2.0/(nx*np1))

      x1rat = (x1max/x1min)**(1.0/(nx*np1))

      x1a(1) = x1min
      dx1 = (x1max-x1min)*(x1rat-1.0)/(x1rat**(nx*np1)-1.0)

      do i = 1, np1*nx     
        if (i.ne.1) x1a(i) = x1a(i-1)+dx1 
        x1b(i) = x1a(i) + 0.5*dx1
        dx1 = x1rat*dx1
!        write(*,*) x1a(i), x1b(i)
      enddo

      dx2 = (x2max-x2min)/(ny*np2) 
  
      do j = 1, np2*ny
        x2a(j) = x2min + (j-1)*dx2
        x2b(j) = x2a(j) + 0.5*dx2
      enddo

!      dx3 = (x3max-x3min)/(nz*np3)
 
      id = 'aa'

      do iincr = 0, 100

      eg = 0.0; dg = 0.0

      do incr = max(0,(iincr-1))*1, (iincr)*1

      do n = 0, np3-1
      do m = 0, np2-1
      do l = 0, np1-1

        write(binfile,"(a3,a2,3i2.2,'.',i3.3)") 'bin',id,l,m,n,incr
        open(unit=42,file=binfile,status='unknown',form='binary' &
        ,convert='big_endian')

        read(42) (((d(i,j,k), i=1,nx), j=1,ny), k=1,nz)
        read(42) (((e(i,j,k), i=1,nx), j=1,ny), k=1,nz)

        dg(nx*l+1:nx*(l+1),ny*m+1:ny*(m+1),nz*n+1:nz*(n+1)) &
        = dg(nx*l+1:nx*(l+1),ny*m+1:ny*(m+1),nz*n+1:nz*(n+1)) &
        + d(1:nx,1:ny,1:nz)

        eg(nx*l+1:nx*(l+1),ny*m+1:ny*(m+1),nz*n+1:nz*(n+1)) &
        = eg(nx*l+1:nx*(l+1),ny*m+1:ny*(m+1),nz*n+1:nz*(n+1)) &
        + e(1:nx,1:ny,1:nz)

        close(42)

      enddo
      enddo
      enddo

      enddo



      dg = dg/((iincr)*1-max(0,(iincr-1))*1+1)
      eg = eg/((iincr)*1-max(0,(iincr-1))*1+1) 


! 1D profiles

      k = 1
      write(onefile,"(a3,i3.3,a10)") 'one', iincr, '_X_C10.dat'
      open(unit=42,file=onefile,status='unknown')

      do i = 1, np1*nx
        d1D(i) = 0.0
        d1D_X(i) = 0.0
        vol = 0.0
        volx = 0.0
        do j = 1, np2*ny
          dvol = sin(x2b(j))
          vol = vol + dvol
          d1D(i) = d1D(i) + dg(i,j,k)*dvol
          e1D(i) = e1D(i) + eg(i,j,k)*dvol
 
          n_e = dg(i,j,k)/mp/mue
          n_i = dg(i,j,k)/mp/mui
          TkeV = 2./3.*eg(i,j,k)/((n_e+n_i)*1.6022e-9)

          if (TkeV.ge.0.1) then
            dvolx = sin(x2b(j))
            volx = volx + dvolx
            d1D_X(i) = d1D_X(i) + dvolx*dg(i,j,k)
            e1D_X(i) = e1D_X(i) + dvolx*eg(i,j,k) 
          endif 

        enddo
        d1D(i) = d1D(i)/vol 
        e1D(i) = e1D(i)/vol
        d1D_X(i) = d1D_X(i)/(volx+1.e-30)
        e1D_X(i) = e1D_X(i)/(volx+1.e-30)
        write(42,2001) x1b(i), d1D(i), d1D_X(i), e1D(i), e1D_X(i)  
      enddo

      close(42)

      enddo

2001  format(5e20.7)

      END
