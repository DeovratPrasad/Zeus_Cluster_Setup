c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 B O N D I                 \\\\\\\\\\
c
c=======================================================================
c
      subroutine bondi
c
c initialize an MHD spherical Bondi flow
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
      integer  :: i, j, k
      real(rl) :: ccode, a3(in+1,jn+1), rad, vr, temp
      real(rl) :: d_min, vr_min, r_min !density, vr at is0-2 
      namelist  /pgen/ rg 
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
       rg = 2.5d-6 ! r_S in code units 

       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
#ifdef MPI_USED
         buf_in( 1) = rg
#endif
       endif
#ifdef MPI_USED
        call MPI_BCAST( buf_in, 1, MPI_FLOAT
     &                , 0, comm3d, ierr )
        if (myid .ne. 0) then
         rg    = buf_in( 1)
        endif
#endif
c
c scaling parameters
c
       ccode = sqrt(2*guniv*ptmass/rg)

c
c      Set up initial vertical field 
c
       do j = 1, jn+1
       do i = 1, in+1
         a3(i,j) = 1e-6*x1a(i)*sin(x2a(j))
       enddo
       enddo

       do k = 1, kn
       do j = 1, jn
       do i = 1, in
         b1(i,j,k) = (g32a(j+1)*a3(i,j+1)-g32a(j)*a3(i,j))/(x2a(j+1)-
     &    x2a(j))/g31a(i)/g32b(j)
         b2(i,j,k) = -(g2a(i+1)*a3(i+1,j)-g2a(i)*a3(i,j))/(x1a(i+1)-
     &    x1a(i))/g2b(i)
         b3(i,j,k) = 0.0
         v2(i,j,k) = 0.0
         v3(i,j,k) = 0.0
       enddo
       enddo
       enddo 
c
c      Set up the flow; go over each process & read the file and initialize the flow 
c
      d_min = 1.0d0

#ifdef MPI_USED
      if (myid == 0) then
#endif
        rad = 0.0
        i = is-2
        k = is-2
        open(unit = 15)
        do while (rg*rad <= x1b(ie+3))
          read(15,*) rad, vr, temp
c line added for Eliot's file
          rad = 0.5*rad
          if (rad*rg >= x1a(k) .and. k <= ie+3) then
c remember Eliot's file has vr positive
            v1(k,:,:) = -vr*ccode
            k = k+1
          endif
          if (rad*rg >= x1b(i) .and. i <= ie+3) then
c density can be specified independently, specifies overall normalization, but
c MBH>>Mdot*time [doesn't matter if no self-gravity].
            d(i,:,:) = d_min*x1b(is-2)**2*v1(is-2,:,:)
     &               /(x1b(i)**2*v1(i,:,:))
            p(i,:,:) = d(i,:,:)*temp*ccode**2
            e(i,:,:) = p(i,:,:)/(gamma-1.0)
            i = i+1
          endif
        enddo
        close(15)
#ifdef MPI_USED
        r_min = x1b(is-2)
        d_min = d(is-2,(js+je)/2,(ks+ke)/2)
        vr_min = v1(is-2,(js+je)/2,(ks+ke)/2)
      endif
      buf_in(1) = r_min
      buf_in(2) = d_min
      buf_in(3) = vr_min
      call MPI_BCAST( buf_in, 3, MPI_FLOAT, 0, comm3d, ierr )
      if (myid .ne. 0) then
        r_min = buf_in(1)
        d_min = buf_in(2)
        vr_min = buf_in(3)
      endif

c open the file for each process, one at a time

      do j = 1, nprocs-1         
        if (j == myid) then
          rad = 0.0
          i = is-2
          k = is-2
          open(unit = 15)
          do while (rg*rad <= x1b(ie+3))
            read(15,*) rad, vr, temp
c line added for Eliot's file
            rad = 0.5*rad
            if (rad*rg >= x1a(k) .and. k <= ie+3) then
c remember Eliot's file has vr positive
              v1(k,:,:) = -vr*ccode
              k = k+1
            endif
            if (rad*rg >= x1b(i) .and. i <= ie+3) then 
c density can be specified independently, specifies overall normalization, but 
c MBH>>Mdot*time [doesn't matter if no self-gravity].
              d(i,:,:) = d_min*r_min**2*vr_min
     &               /(x1b(i)**2*v1(i,:,:))
              p(i,:,:) = d(i,:,:)*temp*ccode**2
              e(i,:,:) = p(i,:,:)/(gamma-1.0)
              i = i+1
            endif
          enddo
          close(15)
        endif
      enddo
#endif

c
c      j = (1+jn)/2
c      k = (1+kn)/2
c      do i = 1, in
c        write(10*(1+myid),1234) x1b(i), d(i,1,1), e(i,j,k), v1(i,j,k)
c     &, b1(i,j,k), b2(i,j,k)
c      enddo
c1234  format(6e20.7)
c
      return
      end
