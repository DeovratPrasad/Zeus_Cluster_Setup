c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                   CRDIFFUSE_D             \\\\\\\\\\
c
c                            Developed by Prateek Sharma
c
c
c driver for anisotropic cosmic ray diffusion. As in qupdate_d, we use 
c a dummy variable (w3de) to store old value of ecr and calculate cosmic 
c ray energy fluxes based on that. another dummy variable (w3dg) carries 
c the updated cosmic ray internal energy density. IMP: we assume 
c that each grid point is updated just once in one subcycle step.
c
c Since cosmic ray diffusion is poorly understood we use a phenomenological
c form for parallel diffusion of cosmic rays, Dcr = alpha*r*VA where alpha
c is a tuning parameter. The diffusion part of CR eq. is:
c
c d/dt(pcr) = div(D grad(pcr)) or in terms of ecr [our primary variable],
c d/dt(ecr) = div(D grad(ecr)), which looks the same.
c
c=======================================================================
c
       subroutine crdiffuse_d
c
c......................................................................
      use real_prec
      use param
      use cons
      use root
      use grid
      use field
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use cosmic 
      use config
      use bndry
      use metallicity
c
      implicit NONE
c
      integer  :: i, j, k, l
      integer  :: nsubcycle
      real(rl) :: dtdiff, dr2min
      integer  :: k1, k2 
      real(rl) :: valf
      real(rl) :: n_e, n_i, TkeV, mui 
c
c----------------------------------------------------------------------
c cosmic ray diffusion is subcycled;
c
c calculation is divided as in pdv_d, the standard way to overlap 
c communication and computation
c
c----------------------------------------------------------------------
c
c purpose of bvald is to update borders used in getting dtdiff for CRs
c
      nreq = 0
      nsub = nsub + 1
      call bvald (1,1,1,1,1,1,d)
#ifdef MPI_USED
         if(nreq .ne. 0)
     .      call MPI_WAITALL ( nreq, req, stat, ierr )
#endif

      mui = 1.0d0/(1.0d0/mu - 1.0d0/mue)

      do k=ks, ke
      do j=js, je  
      do i=is, ie
c        
c choosing Dcr [units are cm^2s^{-1}] from above prescription, Dcr=alpha*r*VA 
c is tricky. e.g., B1 is only defined from is:ie+1 and a zone centered Dcr
c will not be enough, we need Dcr(i-1,:,:), etc. Instead of doing like in
c case of qupdate we can define Dcr locally in crdiffuse.F and use it. Here
c for determining CR update timestep we use an upper-bound for Dcr.
c        
        valf = max( abs(b1(i,j,k)),abs(b1(i+1,j,k)),abs(b2(i,j,k))
     &  ,abs(b2(i,j+1,k)),abs(b3(i,j,k)),abs(b3(i,j,k+1)) )
     &  /sqrt( min( d(i,j,k),d(i-1,j,k),d(i+1,j,k),d(i,j-1,k),d(i,j+1,k)
     &  ,d(i,j,k-1),d(i,j,k+1) ) )
        Dcr(i,j,k) =  alpha_cr*x1b(i)*valf

c Fulai's form for Dcr

c        n_e = d(i,j,k)/(mp*mue)
c        n_i = d(i,j,k)/(mp*mui)
c        TkeV = gamm1*e(i,j,k)/((n_e+n_i)*1.6022d-9)

c        Dcr(i,j,k) = 3.0d28*sqrt(.12/n_e)*(TkeV/1.6)**0.1666666

      enddo
      enddo
      enddo 

      dtdiff = huge
      do k = ks, ke
      do j = js, je
      do i = is, ie
        dr2min = min(dx1a(i),g2b(i)*dx2a(j),g31b(i)*g32b(j)*dx3a(k))
        dr2min = dr2min*dr2min
c
c Courant stability in 3-D: dt<dx^2/6D
c
c        dtdiff = min(dtdiff, 0.166666*dr2min/Dcr(i,j,k))
        dtdiff = min(dtdiff, 0.25*dr2min/Dcr(i,j,k))
      enddo
      enddo
      enddo

#ifdef MPI_USED
      buf_in(1) = dtdiff
      call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
     &                      , MPI_2DOUBLE_PRECISION
     &                      , MPI_MINLOC, comm3d, ierr)
      dtdiff  =   buf_out(1)
#endif

      nsubcycle = int(dt/dtdiff) + 1
      dtdiff = dt/nsubcycle

      do l = 1, nsubcycle
c
c assigning previous ecr to the dummy variable w3de
c
        do k=ks,ke
        do j=js,je
        do i=is,ie
          w3de(i,j,k) = ecr(i,j,k)
        enddo
        enddo
        enddo
c
c-----------------------------------------------------------------------
        if(ldimen .eq. 2) go to 222
        if(ldimen .eq. 1) go to 111
c-----------------------------------------------------------------------
c
c Divide the computational volume into three equal pieces.  We must
c have at least 5 active zones in the 3-direction.
c
         k1 = int( real( ke - ks + 1 ) / 3.0 ) + ks
         k2 = int( real( ke - ks + 1 ) / 3.0 ) + k1
c
c i boundaries
c
c    1) Post sends and receives.
c
         nreq = 0
         nsub = nsub + 1
         call bvalecr (1,1,0,0,0,0,w3de)
c
c    2) Do first portion of the interior points.
c
         call crdiffuse (is+1,ie-1,js+1,je-1,ks+1,k1,w3de,w3dg,dtdiff)
c
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
         if(nreq .ne. 0)
     .      call MPI_WAITALL ( nreq, req, stat, ierr )
#endif      
c......................................................................
c
c j boundaries
c
c    1) Post sends and receives.
c
         nreq = 0
         nsub = nsub + 1
         call bvalecr (0,0,1,1,0,0,w3de)
c
c    2) Do middle 1/3 of the interior points, and some on borders.
c
         call crdiffuse (is,is,js+1,je-1,ks+1,k1,w3de,w3dg,dtdiff)
         call crdiffuse (ie,ie,js+1,je-1,ks+1,k1,w3de,w3dg,dtdiff)
         call crdiffuse (is,ie,js+1,je-1,k1+1,k2,w3de,w3dg,dtdiff)
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
         if(nreq .ne. 0)
     .      call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c k boundaries
c
c    1) Post sends and receives.
c
         nreq = 0
         nsub = nsub + 1
         call bvalecr (0,0,0,0,1,1,w3de)
c
c    2) Do last 1/3 of the interior points, and some on borders.
c
         call crdiffuse (is,ie,js,js,ks+1,k2,w3de,w3dg,dtdiff)
         call crdiffuse (is,ie,je,je,ks+1,k2,w3de,w3dg,dtdiff)
         call crdiffuse (is,ie,js,je,k2+1,ke-1,w3de,w3dg,dtdiff)
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
         if(nreq .ne. 0)
     .      call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c Finally, do the remaining border.
c
         call crdiffuse (is,ie,js,je,ks,ks,w3de,w3dg,dtdiff) 
         call crdiffuse (is,ie,js,je,ke,ke,w3de,w3dg,dtdiff)
c
        go to 999
c======================================================================
c     2D TRANSPORT
c======================================================================
c
222   continue
c
c i boundaries
c
c    1) Post sends and receives.
c
         nreq = 0
         nsub = nsub + 1
         call bvalecr (1,1,0,0,0,0,w3de)
c
c    2) Do first portion of the interior points.
c
         call crdiffuse (is+1,ie-1,js+1,je-1,ks,ks,w3de,w3dg,dtdiff)
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
         if(nreq .ne. 0)
     .      call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c j boundaries
c
c    1) Post sends and receives.
c
         nreq = 0
         nsub = nsub + 1
         call bvalecr (0,0,1,1,0,0,w3de)
c
         call crdiffuse (is,is,js+1,je-1,ks,ks,w3de,w3dg,dtdiff)
         call crdiffuse (ie,ie,js+1,je-1,ks,ks,w3de,w3dg,dtdiff)
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
       if(nreq .ne. 0)
     .    call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c
c Finally, do the remaining border.
c
         call crdiffuse (is,ie,js,js,ks,ks,w3de,w3dg,dtdiff)
         call crdiffuse (is,ie,je,je,ks,ks,w3de,w3dg,dtdiff)
c
      go to 999
c======================================================================
c     1D TRANSPORT
c======================================================================
c
111   continue
c
c i boundaries
c
c    1) Post sends and receives.
c
         nreq = 0
         nsub = nsub + 1
         call bvalecr (1,1,0,0,0,0,w3de)
c
c    2) Do first portion of the interior points.
c
         call crdiffuse (is+1,ie-1,js,js,ks,ks,w3de,w3dg,dtdiff)
c
c    3) Wait for communications to complete.
c
#ifdef MPI_USED
         if(nreq .ne. 0)
     .      call MPI_WAITALL ( nreq, req, stat, ierr )
#endif
c......................................................................
c
c Finally, do the remaining border.
c
         call crdiffuse (is,is,js,js,ks,ks,w3de,w3dg,dtdiff)
         call crdiffuse (ie,ie,js,js,ks,ks,w3de,w3dg,dtdiff)
c
c      Mark the boundary values out of date.
c
999    continue
       do i = 1,6
         bvstat(i,9) = 0      !  making border ecr out of date
       enddo
c
c update the internal energy
c
       do k=ks,ke
       do j=js,je
       do i=is,ie
           ecr(i,j,k) = w3dg(i,j,k)
       enddo
       enddo
       enddo
c
      enddo !subcycle
c
      return
      end
