      PROGRAM stitch
      implicit NONE
! number of true grid points on each processor
      integer, parameter :: nx=16, ny=8, nz=32
      integer, parameter :: np1=8, np2=8, np3=1
      real*8, parameter :: x1min=3.086e+21, x1max=3.086e24
      real*8, parameter :: x2min=0.0,x2max=3.1415926535897932
      real*8, parameter :: x3min=0.0,x3max=6.28318530717959
      logical, parameter :: xmhd=.false., xcosmic=.false.

      integer i, j, k, l, m, n, p, ig, jg, kg, incr, iincr
      character*24 :: binfile, onefile
      character*6 :: twofile
      character*2 :: id

      real*4, dimension(nx,ny,nz) :: d, e, v1, v2, v3, b1, b2, b3, ecr

      real*8, dimension(nx*np1,ny*np2,nz*np3) :: dg, eg, v1g, v2g, v3g &
      , b1g, b2g, b3g, ecrg

      real*8, dimension(nx*np1,ny*np2,nz*np3,0:100) :: eg_X, dg_X

      real*8 :: dx1, dx2, dx3, x1rat
      real*8 :: x1a(nx*np1), x1b(nx*np1), x2a(ny*np2), x2b(ny*np2)

      real*8 :: mu, mui, mue, tcool, crate, TkeV, vol, dvol, volx &
      , dvolx, n_e, n_i, mp, pi
      real*8 :: d1D(nx*np1), e1D(nx*np1), d1D_X(nx*np1), e1D_X(nx*np1)

      real*8 :: cooling_rate
      real*8 :: guniv, rs, mnot, gacc, tff

      mu = 0.61681; mue = 1.1765; mui=1.0/(1./mu-1./mue); mp=1.67265e-24


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

      do iincr = 0, 5

      eg = 0.0; dg = 0.0; eg_X = 0.0; dg_X = 0.0

      do incr = max(0,(iincr-1))*20, max(0,iincr)*20

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

        do k = 1, nz
        do j = 1, ny
        do i = 1, nx

          n_e = d(i,j,k)/mp/mue
          n_i = d(i,j,k)/mp/mui
          TkeV = 2./3.*e(i,j,k)/((n_e+n_i)*1.6022e-9)

          if (TkeV.ge.0.1) then
            dg_X(nx*l+i,ny*m+j,nz*n+k,incr) = d(i,j,k)
            eg_X(nx*l+i,ny*m+j,nz*n+k,incr) = e(i,j,k)
          endif

        enddo
        enddo
        enddo

        close(42)

      enddo
      enddo
      enddo

      enddo

      dg = dg/((iincr)*20-max(0,(iincr-1))*20+1)
      eg = eg/((iincr)*20-max(0,(iincr-1))*20+1) 

! 1D profiles

      write(onefile,"(a3,i3.3,a10)") 'one', iincr, '_Y_C10.dat'
      open(unit=42,file=onefile,status='unknown')

      mnot = 3.8d14*2.d33; guniv = 6.672e-8; rs=390.*3.086e21

      do i = 1, np1*nx
        d1D(i) = 0.0; e1D(i) = 0.0
        d1D_X(i) = 0.0; e1D_X(i) = 0.0
        vol = 0.0
        volx = 0.0
        do j = 1, np2*ny
        do k = 1, np3*nz

          dvol = sin(x2b(j))
          vol = vol + dvol  
 
          d1D(i) = d1D(i) + dg(i,j,k)*dvol
          e1D(i) = e1D(i) + eg(i,j,k)*dvol
 
          do incr = max(0,(iincr-1))*20, max(0,iincr)*20
            if ( dg_X(i,j,k,incr).gt.0.0 ) dvolx = sin(x2b(j))
            volx = volx + dvolx
            d1D_X(i) = d1D_X(i) + dvolx*dg_X(i,j,k,incr)
            e1D_X(i) = e1D_X(i) + dvolx*eg_X(i,j,k,incr) 
          enddo

        enddo
        enddo
        d1D(i) = d1D(i)/vol 
        e1D(i) = e1D(i)/vol
        d1D_X(i) = d1D_X(i)/(volx+1.e-30)
        e1D_X(i) = e1D_X(i)/(volx+1.e-30)


        TkeV = 0.6666666*e1D_X(i)*0.61681*1.67265e-24/(d1D_X(i)*1.6022d-9)
        n_e = d1D_X(i)/(1.67265e-24*1.1765)
        n_i = d1D_X(i)/(1.67265e-24*0.61681) - n_e

        if (TkeV.gt.0.02) then
          cooling_rate = 1.0d-22*n_i*n_e*( 8.6d-3*TkeV**(-1.7) &
         + 5.8d-2*TkeV**0.5 + 6.3d-2)
        else if (TkeV.le.0.02.and.TkeV.ge.0.0017235) then
          cooling_rate = n_i*n_e*6.72e-22*(TkeV/0.02)**0.6
        else
          cooling_rate = n_i*n_e*1.544e-22*(TkeV/0.0017235)**6.0
        endif

        tcool = e1D_X(i)/cooling_rate/pi/1.e16

        gacc = 2.*guniv*mnot*( log(1+x1b(i)/rs)-1./(1.+rs/x1b(i)))  &
        /x1b(i)**2

!        write(*,*) gacc, x1b(i), x1b(i)/rs, mnot, guniv

        tff = sqrt( 2.*x1b(i)/gacc )/pi/1.e7/1.e9

        write(42,2001) x1b(i)/3.086e21, d1D_X(i), e1D_X(i), (0.6666666 &
        *e1D_X(i)*mu/mue/1.6022e-9 &
        /(d1D_X(i)/mue/mp)**1.66666666), tcool, tff  
      enddo

      close(42)

      enddo

2001  format(6e20.7)

      END
