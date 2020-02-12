c=======================================================================
c
      subroutine printd 
c
c  PURPOSE: dumps history file; conserved volume integrated quantities, etc.
c
c-----------------------------------------------------------------------
      use real_prec
      use param
      use config
      use field
      use grid
      use bndry
      use root
      use cons
      use gravmod
      use metallicity
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer, parameter :: nscal = 38
      integer :: i, j, k
      real(rl) :: scal(nscal)
      real(rl) :: dvol, darea, dm 
      real(rl) :: dmin, dav, vav, area, emin, emax
      real(rl) :: mui, n_e, n_i, TkeV, cooling_rate
      real(rl) :: dtcsmin, esnd, dsnd, drmin
      double precision :: buf_in2d(2,1), buf_out2d(2,1)   
      double precision :: buf_in2d1(2,1), buf_out2d1(2,1)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
      mui = 1.0/(1.0/mu-1.0/mue)
      scal = 0.0d0 
      scal(1) = time
      scal(2) = dt
      dmin = huge
      emin = huge; emax = -huge
      dtcsmin = huge

      do k = ks, ke
      do j = js, je
      do i = is, ie
        dvol = dvl1a(i)*dvl2a(j)*dvl3a(k)
        dm = dvol*d(i,j,k)
c total mass
        scal(3) = scal(3) + dm
c total internal energy
        scal(4) = scal(4) + dvol*e(i,j,k)
c total (r,t,p) KE
        scal(5) = scal(5) + dm*0.25*(v1(i,j,k)**2+v1(i+1,j,k)**2)
        scal(6) = scal(6) + dm*0.25*(v2(i,j,k)**2+v2(i,j+1,k)**2)
        scal(7) = scal(7) + dm*0.25*(v3(i,j,k)**2+v3(i,j,k+1)**2)
c total (r,t,p) magnetic energy
        if (xmhd) then
          scal(8) = scal(8) + dvol*0.25*(b1(i,j,k)**2+b1(i+1,j,k)**2)
          scal(9) = scal(9) + dvol*0.25*(b2(i,j,k)**2+b2(i,j+1,k)**2)
          scal(10) = scal(10) + dvol*0.25*(b3(i,j,k)**2+b3(i,j,k+1)**2)
        endif
c gravitational potential energy; for now, either xptmass or xdmpotl
        if (xptmass) then
          scal(11) = scal(11) - (guniv*ptmass/(x1b(i)-rg))*dm
        endif
        if (xdmpotl) then
          scal(11) = scal(11) + phi_dm(i)*dm
        endif
c minimum density
        dmin = min(dmin,d(i,j,k))
        scal(20) = dmin
        emin = min(emin,e(i,j,k))
        scal(24) = emin

        emax = max(emax,e(i,j,k))
        scal(25) = emax

        if (ldimen.eq.3) drmin = min(dx1a(i),g2b(i)*dx2a(j)
     &  ,g31b(i)*g32b(j)*dx3a(k))

        if (ldimen.eq.2)
     &   drmin = min(dx1a(i),g2b(i)*dx2a(j) )

        if (ldimen.eq.1)
     &   drmin = dx1a(i)


        if (dtcsmin.gt.drmin/sqrt(gamma*gamm1*e(i,j,k)/d(i,j,k)) ) then
          dtcsmin = drmin 
     &  /sqrt(gamma*gamm1*e(i,j,k)/d(i,j,k))
          scal(21) = dtcsmin 
          scal(22) = e(i,j,k)
          scal(23) = d(i,j,k)
        endif 
c angular momentum
        scal(14) = scal(14) + dm*0.5*(v3(i,j,k)+v3(i,j,k+1))
     &             *g31b(i)*g32b(j)
c cosmic ray energy
        if (xcosmic) then
          scal(15) = scal(15) + dvol*ecr(i,j,k)
        endif
c
c calculate amount of energy added to the box, which by defn. is amount of cooling for now.
c
        n_e = d(i,j,k)/(mp*mue)
        n_i = d(i,j,k)/(mp*mui)
        TkeV = gamm1*e(i,j,k)/((n_e+n_i)*1.6022d-9)

        if (TkeV.gt.0.02) then
          cooling_rate = 1.0d-22*n_i*n_e*(8.6d-3*TkeV**(-1.7)
     &    + 5.8d-2*TkeV**0.5 + 6.3d-2)
        else if (TkeV.le.0.02.and.TkeV.ge.0.0017235) then
          cooling_rate = n_i*n_e*6.72d-22*(TkeV/0.02)**0.6
        else
          cooling_rate = n_i*n_e*1.544d-22*(TkeV/0.0017235)**6.0
        endif
c
        if (TkeV.le.1.0d-4) cooling_rate = 0.0d0

        scal(16) = scal(16) + dvol*cooling_rate
c
        if (TkeV.gt.0.1) scal(17) = scal(17) + dm

        if (TkeV.le.0.1.and.TkeV.gt.0.01) scal(18) = scal(18) + dm

        if (TkeV.le.0.01) scal(19) = scal(19) + dm
c
      enddo
      enddo
      enddo
      scal(12) = eadd
      scal(13) = madd    
#ifdef MPI_USED
      buf_in(1:17) = scal(3:19)
      buf_in(18) = scal(20)
      buf_in(19) = scal(24)
      buf_in(20) = scal(25)


      buf_in2d(1,1) = scal(21)
      buf_in2d(2,1) = scal(22)

      buf_in2d1(1,1) = scal(21)
      buf_in2d1(2,1) = scal(23)

c
c MPI_2DOUBLE_PRECISION was not working fine with MPI_SUM
c
      call MPI_REDUCE(buf_in, buf_out, 17, MPI_FLOAT
     &, MPI_SUM, 0, comm3d, ierr)

      call MPI_REDUCE(buf_in(18), buf_out(18), 2, MPI_FLOAT
     &, MPI_MIN, 0, comm3d, ierr)

      call MPI_REDUCE(buf_in(20), buf_out(20), 1, MPI_FLOAT
     &, MPI_MAX, 0, comm3d, ierr)

      call MPI_REDUCE(buf_in2d, buf_out2d, 1,
     & MPI_2DOUBLE_PRECISION, MPI_MINLOC, 0, comm3d, ierr)

      call MPI_REDUCE(buf_in2d1, buf_out2d1, 1,
     & MPI_2DOUBLE_PRECISION, MPI_MINLOC, 0, comm3d, ierr)


      if (myid.eq.0) then
        scal(3:19) = buf_out(1:17)
        scal(20) = buf_out(18) 
        scal(24) = buf_out(19)
        scal(25) = buf_out(20)
        scal(21) = buf_out2d(1,1)
        scal(22) = buf_out2d(2,1)
        scal(23) = buf_out2d1(2,1)
#endif
      write(3,2001) (scal(i), i=1, 25)
#ifdef MPI_USED
      endif
#endif
2001  format(25e20.7)
c radial luminosities through inner boundary; integrated over theta and phi

      scal = 0.0d0
      scal(1) = time
      i = is
      dav = 0.0; vav = 0.0; area = 0.0
      do k = ks, ke
      do j = js, je
        darea = g2a(i)*g31a(i)*dvl2a(j)*dvl3a(k)
        dm = darea*d(i,j,k)*v1(i,j,k)
c mass
        scal(2) = scal(2) + dm
c internal energy
        scal(3) = scal(3) + darea*e(i,j,k)*v1(i,j,k)
c angular momentum
        scal(4) = scal(4) + dm*v3(i,j,k)*g31a(i)*g32b(j)
c kinetic energy
        scal(5) = scal(5) + dm*0.5*(v1(i,j,k)**2+
     &  0.5*(v2(i,j,k)**2+v2(i,j+1,k)**2)+v3(i,j,k)**2)
c Poynting E*B
        if (xmhd) then
          scal(6) = scal(6) +  darea*(
     &             v3(i,j,k)*b1(i,j,k)*b3(i,j,k)
     &             +0.25*(v2(i,j,k)+v2(i,j+1,k))*b1(i,j,k)
     &             *(b2(i,j,k)+b2(i,j+1,k))
     &             -v1(i,j,k)*(b3(i,j,k)**2
     &             +0.5*(b2(i,j,k)**2+b2(i,j+1,k)**2)))
        endif
        if (xcosmic) then
          scal(12) = scal(12) + darea*ecr(i,j,k)*v1(i,j,k)
        endif
        dav = dav + darea*d(i,j,k)
        vav = vav + darea*v1(i,j,k)
        area = area + darea
c
c diagnostics for inflow temperature structure; hot vs. cold feedback
c
        n_e = d(i,j,k)/(mp*mue)
        n_i = d(i,j,k)/(mp*mui)
        TkeV = gamm1*e(i,j,k)/((n_e+n_i)*1.6022d-9)

        if (TkeV.gt.0.1) then 
          scal(14) = scal(14) + dm
          scal(15) = scal(15) + darea*e(i,j,k)*v1(i,j,k)
          scal(16) = scal(16) + dm*v3(i,j,k)*g31a(i)*g32b(j)
          scal(17) = scal(17) + dm*0.5*(v1(i,j,k)**2+
     &  0.5*(v2(i,j,k)**2+v2(i,j+1,k)**2)+v3(i,j,k)**2)
        endif

        if (TkeV.le.0.1.and.TkeV.gt.0.01) then
          scal(18) = scal(18) + dm
          scal(19) = scal(19) + darea*e(i,j,k)*v1(i,j,k)
          scal(20) = scal(20) + dm*v3(i,j,k)*g31a(i)*g32b(j)
          scal(21) = scal(21) + dm*0.5*(v1(i,j,k)**2+
     &  0.5*(v2(i,j,k)**2+v2(i,j+1,k)**2)+v3(i,j,k)**2)
        endif

        if (TkeV.le.0.01) then
          scal(22) = scal(22) + dm
          scal(23) = scal(23) + darea*e(i,j,k)*v1(i,j,k)
          scal(24) = scal(24) + dm*v3(i,j,k)*g31a(i)*g32b(j)
          scal(25) = scal(25) + dm*0.5*(v1(i,j,k)**2+
     &  0.5*(v2(i,j,k)**2+v2(i,j+1,k)**2)+v3(i,j,k)**2)
        endif

      enddo
      enddo

      scal(32) = dav*vav/area

c radial luminosities through outer boundary; integrated over theta and phi

      i = ie+1
      dav = 0.0; vav = 0.0; area = 0.0
      do k = ks, ke
      do j = js, je
        darea = g2a(i)*g31a(i)*dvl2a(j)*dvl3a(k)
        dm = darea*d(i,j,k)*v1(i,j,k)
c mass
        scal(7) = scal(7) + dm
c internal energy
        scal(8) = scal(8) + darea*e(i,j,k)*v1(i,j,k)
c angular momentum
        scal(9) = scal(9) + dm*v3(i,j,k)*g31a(i)*g32b(j)
c kinetic energy
        scal(10) = scal(10) + dm*0.5*(v1(i,j,k)**2+
     &  0.5*(v2(i,j,k)**2+v2(i,j+1,k)**2)+v3(i,j,k)**2)
c Poynting E*B
        if (xmhd) then
          scal(11) = scal(11) +  darea*(
     &             v3(i,j,k)*b1(i,j,k)*b3(i,j,k)
     &             +0.25*(v2(i,j,k)+v2(i,j+1,k))*b1(i,j,k)
     &             *(b2(i,j,k)+b2(i,j+1,k))
     &             -v1(i,j,k)*(b3(i,j,k)**2
     &             +0.5*(b2(i,j,k)**2+b2(i,j+1,k)**2)))
        endif
        if (xcosmic) then
          scal(13) = scal(13) + darea*ecr(i,j,k)*v1(i,j,k)
        endif
        dav = dav + darea*d(i,j,k)
        vav = vav + darea*v1(i,j,k)
        area = area + darea
      enddo
      enddo

      scal(33) = dav*vav/area

c radial luminosities through r = sqrt(i_{min}i_{max}); integrated over theta and phi

      i = (is+ie)/2 
      dav = 0.0; vav = 0.0; area = 0.0
      do k = ks, ke
      do j = js, je
        darea = g2a(i)*g31a(i)*dvl2a(j)*dvl3a(k)
        dm = darea*d(i,j,k)*v1(i,j,k)
c mass
        scal(26) = scal(26) + dm
c internal energy
        scal(27) = scal(27) + darea*e(i,j,k)*v1(i,j,k)
c angular momentum
        scal(28) = scal(28) + dm*v3(i,j,k)*g31a(i)*g32b(j)
c kinetic energy
        scal(29) = scal(29) + dm*0.5*(v1(i,j,k)**2+
     &  0.5*(v2(i,j,k)**2+v2(i,j+1,k)**2)+v3(i,j,k)**2)
c Poynting E*B
        if (xmhd) then
          scal(30) = scal(30) +  darea*(
     &             v3(i,j,k)*b1(i,j,k)*b3(i,j,k)
     &             +0.25*(v2(i,j,k)+v2(i,j+1,k))*b1(i,j,k)
     &             *(b2(i,j,k)+b2(i,j+1,k))
     &             -v1(i,j,k)*(b3(i,j,k)**2
     &             +0.5*(b2(i,j,k)**2+b2(i,j+1,k)**2)))
        endif
        if (xcosmic) then
          scal(31) = scal(31) + darea*ecr(i,j,k)*v1(i,j,k)
        endif
        dav = dav + darea*d(i,j,k)
        vav = vav + darea*v1(i,j,k)
        area = area + darea
      enddo
      enddo
      scal(34) = dav*vav/area

      i = (3*is+ie)/4
      dav = 0.0; vav = 0.0; area = 0.0
      do k = ks, ke
      do j = js, je
        darea = g2a(i)*g31a(i)*dvl2a(j)*dvl3a(k)
        dm = darea*d(i,j,k)*v1(i,j,k)
c mass
        scal(35) = scal(35) + dm
        dav = dav + darea*d(i,j,k)
        vav = vav + darea*v1(i,j,k)
        area = area + darea
      enddo
      enddo
      scal(36) = dav*vav/area

      i = (is+3*ie)/4
      dav = 0.0; vav = 0.0; area = 0.0
      do k = ks, ke
      do j = js, je
        darea = g2a(i)*g31a(i)*dvl2a(j)*dvl3a(k)
        dm = darea*d(i,j,k)*v1(i,j,k)
c mass
        scal(37) = scal(37) + dm
        dav = dav + darea*d(i,j,k)
        vav = vav + darea*v1(i,j,k)
        area = area + darea
      enddo
      enddo
      scal(38) = dav*vav/area


      write(51,2002) (scal(i), i=1, 38)
2002  format(38e20.7)
c
      return
      end
