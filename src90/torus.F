c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 T O R U S                 \\\\\\\\\\
c
c=======================================================================
c
      subroutine torus 
c
c initialize Stone & Pringle 2001: P&P constant specific angular momentum 
c torus 
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use field
      use grid
      use cons
      use gravmod
      use root
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit none
c
      integer  :: i, j, k, l, iprob
      real(rl) :: eout(jn)
      real(rl) :: d0,acons,dist,w0,dens,fac,a3(in,jn),beta
      real(rl) :: w,z,eq29,cprime,en,bz0,rmin,rmax,amp,dmax,pmax 
      namelist  /pgen/ d0,acons,dist,w0,beta,rg,iprob,bz0
c
c-----------------------------------------------------------------------
c
c defining the properties of the central compact object; its located at origin 
c
       guniv = 1.0
       ptmass = 1.0
       x1ptm = 0.0
       x2ptm = 0.0
       x3ptm = 0.0
c
       amp = 1.0e-2
       d0 = 1.0e-4
       acons = 0.06667
       dist = 1.5
       w0 = 1.0
       beta = 200.0
       rg  = 0.5 ! r_S should be in code units
       iprob = 1
       bz0 = 0.0
      
       rmin = x1a(is)

       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
#ifdef MPI_USED
         buf_in(1) = d0
         buf_in(2) = acons
         buf_in(3) = dist
         buf_in(4) = w0
         buf_in(5) = beta
         buf_in(6) = rg
         buf_in(7) = bz0
         buf_in(8) = rmin
         ibuf_in(1) = iprob 
#endif
       endif
#ifdef MPI_USED
        call MPI_BCAST( buf_in, 8, MPI_FLOAT
     &                , 0, comm3d, ierr )
        call MPI_BCAST( ibuf_in, 1, MPI_INTEGER
     &                , 0, comm3d, ierr )
        if (myid .ne. 0) then
         d0 = buf_in(1)
         acons = buf_in(2)
         dist = buf_in(3)
         w0 = buf_in(4)
         beta = buf_in(5)
         rg = buf_in(6)
         bz0 = buf_in(7)
         rmin = buf_in(8)
         iprob = ibuf_in(1)
        endif
#endif
c
c  Setup background medium with constant density
c
       d = d0
       v1 = 0.0
       v2 = 0.0
       v3 = 0.0
       a3 = 0.0
       b3 = 0.0
c
c  Initialize thick torus using PP formulae
c  Note that:
c     dist = distortion of the torus, d=1 gives torus of small circular
c             cross section, d>1 gives fatter and fatter profiles.
c     w0 = radial position of density max
c  Set only density; pressure set using difference formula below
c    to achieve exact equilibrium
c  setup the vector potential A3 equal to p -- gives field lines
c    parallel to contours of p -- note how amplitude set
c
       cprime = 0.5/dist
       en = 1.0/gamm1

      do k=1,kn
      do j=1,jn
      do i=1,in
        w = x1b(i)*sin(x2b(j))
        z = x1b(i)*cos(x2b(j))
c recall rg is the Schwarzschild radius in code units
        eq29 = (guniv*ptmass)/(w0*(en + 1.))
     &       *(w0/(x1b(i)-rg)
     &     - 0.5*w0**3/((w0-rg)*w)**2 - cprime)
        fac = 0.0
        if (eq29 .gt. 0.0) then
          dens  = (eq29/acons)**en
          if (dens .gt. 100.*d(i,j,k)) fac = 1.0
          if (iprob .eq. 1 .and. dens .ge. 0.1)
     &        a3(i,j) = (dens-0.1)**2/(beta)
          if (iprob .eq. 3) then
            if (dens .ge. 0.1  .and. dens .lt. 0.54)
     &        a3(i,j) = (dens-0.1)**2/(beta)
            if (dens .ge. 0.54 .and. dens .lt. 0.73)
     &        a3(i,j) = (0.4-(dens-0.1)**2)/(beta)
            if (dens .ge. 0.73 .and. dens .lt. 0.87)
     &        a3(i,j) =(-0.4+(dens-0.1)**2)/(beta)
            if (dens .ge. 0.87)
     &        a3(i,j) = (0.8-(dens-0.1)**2)/(beta)
          endif
          v3(i,j,k) = sqrt(guniv*ptmass*w0**3/(w0-rg)**2)/w
        endif
        d(i,j,k) = (1.0-fac)* d(i,j,k) + fac*dens
      enddo
      enddo  
      enddo
c
c       setup the Magnetic field by differencing A
c
      if (iprob .eq. 1 .or. iprob .eq. 3) then
c
      do j=1,jn-1
      do i=1,in-1
        b1(i,j,:) = (g32a(j+1)*A3(i,j+1)-g32a(j)*A3(i,j))/(x2a(j+1)-
     &    x2a(j))/g31a(i)/g32b(j)
        b2(i,j,:) = -(g2a(i+1)*A3(i+1,j)-g2a(i)*A3(i,j))/(x1a(i+1)-
     &    x1a(i))/g2b(i)
      enddo
      enddo
c
c option for vertical field lines only through disk is not set yet 
c
      endif
c
c setting up exact equilibrium pressure
c
      e(in,:,:) = guniv*ptmass*d0/rmin

#ifdef MPI_USED
      do l = ntiles(1)-1, 0, -1
        if (l.eq.coords(1)) then
          if (l.ne.ntiles(1)-1) then
            nreq=nreq+1
            call MPI_RECV(eout, jn, MPI_FLOAT, n1p, 12345
     &      , comm3d, stat, ierr)
            call MPI_WAIT(stat, req(nreq), ierr)
            do k=1, kn
              e(ie+3,:,k) = eout(:)
            enddo
          endif
#endif
          do j = 1, jn
          do i = in-1, 1, -1
            e(i,j,:) = e(i+1,j,:) +
     &   (guniv*ptmass/(x1a(i+1)-rg)**2
     &    -0.5*(v3(i+1,j,:)**2+v3(i,j,:)**2)/g31a(i+1))
     &   *dx1b(i+1)*(d(i+1,j,:)+d(i,j,:))/(2.*gamm1)
          enddo
          enddo
          eout = e(is+2,:,(ks+ke)/2)
#ifdef MPI_USED
          if (l.ne.0) then
            nreq = nreq + 1
            call MPI_ISEND(eout, jn, MPI_FLOAT, n1m, 12345
     &      , comm3d, req(nreq), ierr)
          endif

        endif
      enddo 
#endif

      p = gamm1*e
c
      j = (1+jn)/2
      k = (1+kn)/2
      do i = 1, in
        write(10*(1+myid),1234) x1b(i), d(i,j,k), e(i,j,k), p(i,j,k)
     &, b1(i,j,k), b2(i,j,k)
      enddo
1234  format(6e20.7)
c
      return
      end
