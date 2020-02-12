c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                  CRSOURCE                 \\\\\\\\\\
c
c                            Developed by
c                            Prateek Sharma 
c
c=======================================================================
c
       subroutine crsource
c
c implements a source term with a fixed mdot_in; eq. [19] of Guo & Oh
c
c......................................................................
      use real_prec
      use config
      use param
      use root
      use grid
      use field
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use cons
      use bndry
c
      implicit NONE
c
      integer  :: i, j, k
      real(rl) :: cr_rate
c mdot_in=6.36d25 is 1Mo/yr, eps=0.003 used in Guo
      real(rl), parameter :: mdot_in=6.36d25, nu=1.5, eps=0.00003
     &, rnot=6.16d22
      real(rl) :: theta0
c
      theta0=pi/6.0

      do k=ks, ke
      do j=js, je
      do i=is, ie

        if (x2b(j).lt.theta0.or.x2b(j).gt.pi-theta0) then

        cr_rate = nu*eps*mdot_in*clight*clight/(4.0*pi*rnot*rnot*rnot)
     &  *(rnot*x1bi(i))**(3.0+nu)*(1.0-exp(-x1b(i)*x1b(i)/(rnot*rnot)))

c        cr_rate = eps*mdot_in*clight*clight/(4.0*pi*rnot*rnot*rnot)
c     &  *exp(-x1b(i)/rnot)

        ecr(i,j,k) = ecr(i,j,k) + cr_rate*dt
 
        endif

      enddo
      enddo
      enddo
c
c making the boundary values out of date
c
      do i = 1,6
       bvstat(i,9) = 0      !  ecr
      enddo
c
      return
      end
