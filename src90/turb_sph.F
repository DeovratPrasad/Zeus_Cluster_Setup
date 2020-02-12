c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////               TURB_SPH                    \\\\\\\\\\
c
       subroutine turb_sph
c
c implement turbulent forcing in spherical geometry for the cluster problem
c using Ben Chandran's idea but applying in spherical geometry
c
c   dv/dt = F; thus forcing is proportional to density; i.e., roughly similar
c mach number is obtained everywhere.
c
c
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use root
      use field
      use grid
      use bndry
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use cons
c
      implicit NONE
      real(rl) :: lnot, alpha, cs0, amp
      integer :: i, j, k, l, m, n, iseed
      real(rl) :: rv1, rv2, rv3, kx, ky, kz, Fr, Ft, Fp, x, y, z, p1
      real(rl) :: Fx(in,jn,kn), Fy(in,jn,kn), Fz(in,jn,kn)
      real(rl) :: Ftx(in,jn,kn), Fty(in,jn,kn), Ftz(in,jn,kn)
c
      lnot = 20.0*3.086e21; alpha = 1.e-4; cs0 = 1.e7
      iseed = 1
      Ftx = 0.0; Fty = 0.0; Ftz = 0.0

      if (ldimen.eq.2) then

      do m = -4, 4 
      do l = -4, 4

      if (l.ne.0.and.m.ne.0) then

          kx = l*2.*pi/lnot
          ky = m*2.*pi/lnot
c
c choose amplitude of forcing s.t. F = alpha * c_s/t, t=20kpc/c_s; alpha ~ 0.01 say 
c
          amp = alpha*cs0*cs0/lnot/(l**2+m**2)**.5

          call ran2(iseed, rv1)
          call ran2(iseed, rv2) 

        Fx = 0.0; Fy = 0.0
        do k = ks, ke
        do j = js-1, je+1
        do i = is-1, ie+1
 
          x = x1b(i)*sin(x2b(j))*cos(x3b(k))
          y = x1b(i)*sin(x2b(j))*sin(x3b(k))

          p1 = kx*x + ky*y + 2.*pi*rv1

          Fx(i,j,k) = Fx(i,j,k) + amp*cos(p1) 

          p1 = kx*x + ky*y + 2.*pi*rv2

          Fy(i,j,k) = Fy(i,j,k) + amp*cos(p1)

        enddo
        enddo
        enddo 

        do k = ks, ke
        do j = js, je
        do i = is, ie
          Ftx(i,j,k) = Ftx(i,j,k) 
     &    + (Fx(i,j,k)*ky*ky - Fy(i,j,k)*kx*ky)/(kx*kx + ky*ky)
          Fty(i,j,k) = Fty(i,j,k)
     &    + (Fy(i,j,k)*kx*kx - Fx(i,j,k)*kx*ky)/(kx*kx + ky*ky) 
        enddo
        enddo
        enddo 

      endif

      enddo
      enddo
 
      do k = ks, ke
      do j = js, je
      do i = is, ie
c
c convert x,y,z force components into r, theta, phi
c        
        Fr = 0.5*(Ftx(i,j,k)+Ftx(i-1,j,k))*sin(x2b(j))*cos(x3b(k)) 
     &  + 0.5*(Fty(i,j,k)+Fty(i-1,j,k))*sin(x2b(j))*sin(x3b(k))

        Ft = 0.5*( Ftx(i,j,k)*cos(x2b(j))+Ftx(i,j-1,k)*cos(x2b(j-1)) )
     &  *cos(x3b(k))
     &  + 0.5*(Fty(i,j,k)*cos(x2b(j))+Fty(i,j-1,k)*cos(x2b(j-1)) )
     &  *sin(x3b(k))


        v1(i,j,k) = v1(i,j,k) + dt*Fr
        v2(i,j,k) = v2(i,j,k) + dt*Ft

      enddo
      enddo
      enddo

      endif

      if (ldimen.eq.3) then

      do n = -4, 4
      do m = -4, 4
      do l = -4, 4

      if (l.ne.0.and.m.ne.0) then

          kx = l*2.*pi/lnot
          ky = m*2.*pi/lnot
          kz = n*2.*pi/lnot
c
c choose amplitude of forcing s.t. F = alpha * c_s/t, t=20kpc/c_s; alpha ~ 0.01 say 
c
          amp = alpha*cs0*cs0/lnot/(l**2+m**2+n**2)**.5

          call ran2(iseed, rv1)
          call ran2(iseed, rv2)
          call ran2(iseed, rv3)

        Fx = 0.0; Fy = 0.0; Fz=0.0
        do k = ks, ke
        do j = js, je
        do i = is, ie

          x = x1b(i)*sin(x2b(j))*cos(x3b(k))
          y = x1b(i)*sin(x2b(j))*sin(x3b(k))

          p1 = kx*x + ky*y + 2.*pi*rv1

          Fx(i,j,k) = Fx(i,j,k) + amp*cos(p1)

          p1 = kx*x + ky*y + 2.*pi*rv2

          Fy(i,j,k) = Fy(i,j,k) + amp*cos(p1)

          p1 = kx*x + ky*y + 2.*pi*rv2

          Fz(i,j,k) = Fz(i,j,k) + amp*cos(p1)
 
        enddo
        enddo
        enddo

        do k = ks, ke
        do j = js, je
        do i = is, ie
          Ftx(i,j,k) = Ftx(i,j,k)
     &    + (Fx(i,j,k)*(ky*ky+kz*kz)-Fy(i,j,k)*kx*ky-Fz(i,j,k)*kx*kz)
     &    /(kx*kx + ky*ky + kz*kz)
          Fty(i,j,k) = Fty(i,j,k)
     &    + (Fy(i,j,k)*(kx*kx+kz*kz)-Fx(i,j,k)*kx*ky-Fz(i,j,k)*ky*kz)
     &    /(kx*kx + ky*ky + kz*kz)
          Ftz(i,j,k) = Ftz(i,j,k)
     &    + (Fz(i,j,k)*(kx*kx+ky*ky)-Fx(i,j,k)*kx*kz-Fy(i,j,k)*ky*kz)
     &    /(kx*kx + ky*ky + kz*kz)
        enddo
        enddo
        enddo

      endif

      enddo
      enddo
      enddo

      endif
c
      return
      end
c
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////               TURB_SPH                    \\\\\\\\\\
c
c=======================================================================
c
