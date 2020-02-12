c=======================================================================
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                  SOURCE_JET               \\\\\\\\\\
c
c                            Developed by
c                            Prateek Sharma 
c
c=======================================================================
c
c implements smooth jet source terms
c
       subroutine source_jet
c
c......................................................................
c
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
      use jet_vars
      use d_e_out
c
      implicit NONE
      real(rl) :: psi_jet, mdot_acc, mdot_jet, darea, dnew(in,jn,kn)
      integer :: i, j, k
c
       nreq = 0
       nsub = nsub + 1
       call bvald  (1,0,0,0,0,0,d)
       call bvalv1 (1,0,0,0,0,0,v1)
       dnew = d
c
c  apply mass source
c
      mdot_acc = 0.0d0
c   mass accretion through the inner radius
#ifdef MPI_USED
      if(coords(1).eq.0)then
#endif 
      do k=ks, ke
      do j=js, je
         darea = g2a(is)*g31a(is)*dvl2a(j)*dvl3a(k)
         mdot_acc=mdot_acc+darea*d(is,j,k)*v1(is,j,k)
      enddo
      enddo
#ifdef MPI_USED
      endif
      buf_in(1) = mdot_acc
      call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1, MPI_FLOAT
     &,                  MPI_SUM, comm3d, ierr)
      mdot_acc = buf_out(1)
#endif
c
c fix mdot so that it doesn't exceed the eddington valuye for 1e10 BH
c
c      mdot_acc = max(mdot_acc,-250.d0*msol/(pi*1.d7))
c
c positive mdot_jet
c
      mdot_jet = jet_eff*max(-mdot_acc,0.0d0)
c
      do k=ks, ke
      do j=js, je
      do i=is+1, ie

        dnew(i,j,k) = d(i,j,k) + dt*norm_jet*mdot_jet
     &  *psi_jet( x1b(i),x2b(j) )
c
        v1(i,j,k) = v1(i,j,k) + 2.0d0*dt*norm_jet*mdot_jet*v_jet
     &  *psi_jet( x1a(i),x2b(j) )/(d(i,j,k)+d(i-1,j,k))

        if ( v1(i,j,k).gt.clight ) then 
          v1(i,j,k)=clight
c          write(*,*) time, nhy, i, j, 'applied vfloor'
        endif
c
c        if (nspec.gt.1.and.x1b(i).lt.r_jet.and.x2b(j).lt.th_jet ) then
c          abun(i,j,k,1) = 1.0d15
c        endif
c        if(nspec.gt.1.and.x1b(i).lt.r_jet.and.x2b(j).gt.(pi-th_jet))then
c          abun(i,j,k,1) = 1.0d15
c        endif
c
      enddo
      enddo
      enddo
c
#ifdef MPI_USED
c
c    3) Wait for communications to complete.
c
       if (nreq .ne. 0) call MPI_WAITALL ( nreq, req, stat, ierr )
#endif

      i=is
      do k=ks, ke
      do j=js, je
c
        dnew(i,j,k) = d(i,j,k) + dt*norm_jet*mdot_jet
     &  *psi_jet( x1b(i),x2b(j) )
c
        v1(i,j,k) = v1(i,j,k) + 2.0d0*dt*norm_jet*mdot_jet*v_jet
     &  *psi_jet( x1a(i),x2b(j) )/(d(i,j,k)+d(i-1,j,k))

        if ( v1(i,j,k).gt.clight ) then 
           v1(i,j,k)=clight
c           write(*,*) time, nhy, i, j, 'applied vfloor'
        endif
c
c        if (nspec.gt.1.and.x1b(i).lt.r_jet.and.x2b(j).lt.th_jet ) then
c          abun(i,j,k,1) = 1.0d15
c        endif
c        if(nspec.gt.1.and.x1b(i).lt.r_jet.and.x2b(j).gt.(pi-th_jet))then
c          abun(i,j,k,1) = 1.0d15
c        endif
c
      enddo
      enddo

      d = dnew

#ifdef MPI_USED
      if (coords(1).eq.ntiles(1)-1) then
#endif
        d(ie,:,:) = d_outer
#ifdef MPI_USED
      endif
#endif
c
c making the boundary values out of date
c
      do i = 1,6
       bvstat(i,1) = 0      !  d
       bvstat(i,3) = 0      ! v1
      enddo

      return
      end

      function psi_jet( rad, theta)
      use real_prec
      use param
      use cons
      use jet_vars

      implicit NONE
      real(rl) :: psi_jet, rad, theta

      psi_jet = 0.25d0*(2.0d0+tanh((th_jet-theta)/sig_th)
     & +tanh((th_jet+theta-pi)/sig_th))*(1.0d0+tanh((r_jet-rad)/sig_r))

      return
      end
