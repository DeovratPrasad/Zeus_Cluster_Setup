      PROGRAM stitch
      implicit NONE
! number of true grid points on each processor
      integer, parameter :: nx=64, ny=64, nz=1
      integer, parameter :: np1=8, np2=4, np3=1
      integer, parameter :: incr_min=81, incr_max=100
      real*8, parameter :: x1min=3.086e+21, x1max=6.172e+23
      real*8, parameter :: x2min=0.0,x2max=3.1415926535897932
      real*8, parameter :: x3min=0.0,x3max=6.28318530717959
      logical, parameter :: xmhd=.false., xcosmic=.false.

      integer i, l, m, n, incr
      character*15 :: onefile
      character*2 :: id

      real*8, dimension(nx*np1) :: dg, eg, davg, eavg

      real*8 :: x1b(nx*np1)
      real*8 :: fac_m
      
      id = 'aa'

!      davg = 0.0; eavg = 0.0

      do incr = incr_min, incr_max

      davg = 0.0; eavg = 0.0

!      write(*,*) incr

      do n = 0, np3-1
      do m = 0, np2-1
!
! the following only works if i have 4 processors in the theta direction
!
        if (m.eq.0.or.m.eq.3) then
          fac_m = 0.5*(1.-1./sqrt(2.))
        else
          fac_m = 0.5/sqrt(2.)
        endif
      do l = 0, np1-1

        write(onefile,"(a3,a2,3i2.2,'.',i3.3)") 'one',id,l,m,n,incr
!        write(*,*) onefile
        open(unit=42,file=onefile,status='old')

        read(42,*)
        read(42,*)

        do i = 1, nx

          read(42,2001) x1b(i+nx*l), dg(i+nx*l), eg(i+nx*l)

! this averaging is not quite right!

          davg(i+nx*l) = davg(i+nx*l) + fac_m*dg(i+nx*l)
          eavg(i+nx*l) = eavg(i+nx*l) + fac_m*eg(i+nx*l)

!          write(*,*) i

        enddo

        close(42)

      enddo
      enddo
      enddo

c      enddo !incr
 
c      davg = davg/(incr_max-incr_min+1); eavg = eavg/(incr_max-incr_min+1)
      

      write(onefile,"(a3,a2,i3.3)") 'one', id, incr
      open(unit=42,file=onefile,status='unknown')

      do i = 1, np1*nx

        write(42,2001) x1b(i), davg(i), eavg(i)

      enddo

      close(42) 

      enddo !incr

2001  format(3e20.7)

      END
