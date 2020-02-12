c=======================================================================
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                  COOLING                  \\\\\\\\\\
c
c                            Developed by
c                            Prateek Sharma 
c
c=======================================================================
c
       subroutine cooling 
c
c Implements f-f (+metals) cooling; assuming that CGS units are used
c Density cannot simply be scaled out like in case when there is no cooling
c or self gravity. Tozzi & Norman 2001 cooling function; see Chandran 2007, 
c Guo & Peng Oh 2007
c
c......................................................................
      use real_prec
      use config
      use param
      use root
      use grid
      use field
      use d_e_out
      use brag, only : q_cond_tot
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use metallicity 
      use cons
      use bndry
c
      implicit NONE
c
      integer  :: i, j, k, n, nsubcycle
      real(rl) :: mui, n_e, n_i, TkeV, dvol, vol(in)
      real(rl) :: cooling_rate(in,jn,kn), const_heat(in)
      real(rl) :: darea, dtcool
      real(rl) :: buf_send(2*in), buf_recv(2*in)
c
      mui = 1.0d0/(1.0d0/mu - 1.0d0/mue)
c
      dtcool = huge

      do k=ks, ke
      do j=js, je
      do i=is, ie

        n_e = d(i,j,k)/(mp*mue)
        n_i = d(i,j,k)/(mp*mui)
        TkeV = gamm1*e(i,j,k)/((n_e+n_i)*1.6022d-9)
c
        if (TkeV.gt.0.02) then
          cooling_rate(i,j,k) = 1.0d-22*n_i*n_e*( 8.6d-3*TkeV**(-1.7)
     &    + 5.8d-2*TkeV**0.5 + 6.3d-2)
        else if (TkeV.le.0.02.and.TkeV.ge.0.0017235) then
          cooling_rate(i,j,k) = n_i*n_e*6.72e-22*(TkeV/0.02)**0.6 
        else 
          cooling_rate(i,j,k) = n_i*n_e*1.544e-22*(TkeV/0.0017235)**6.0
        endif 
c
        if (TkeV.le.1.0d-4) cooling_rate(i,j,k) = 0.0d0
c
        dtcool = min( dtcool
     &  , e(i,j,k)/(cooling_rate(i,j,k)+tiny) )

      enddo
      enddo
      enddo
c
#ifdef MPI_USED
      buf_in(1) = dtcool
      call MPI_ALLREDUCE( buf_in(1), buf_out(1), 1
     &                      , MPI_2DOUBLE_PRECISION
     &                      , MPI_MINLOC, comm3d, ierr)
      dtcool  =   buf_out(1) 
#endif
c
c factor of 4.0 chosen arbitrarily
c
      nsubcycle = min(int(4.0*dt/(dtcool+tiny) ) + 1, 100)

      dtcool = dt/float(nsubcycle)
c
      if (myid.eq.0.and.mod(nhy,100).eq.0) 
     & write(17,*) nsubcycle
c
      rad_cool=0.0
      feed_heat=0.0

      do n = 1, nsubcycle
c
      do i = is, ie 
        const_heat(i) = 0.0d0
        vol(i) = 0.0d0
      do k = ks, ke
      do j = js, je

        n_e = d(i,j,k)/(mp*mue)
        n_i = d(i,j,k)/(mp*mui)
        TkeV = gamm1*e(i,j,k)/((n_e+n_i)*1.6022d-9)
c
        if (TkeV.gt.0.02) then
          cooling_rate(i,j,k) = 1.0d-22*n_i*n_e*( 8.6d-3*TkeV**(-1.7)
     &    + 5.8d-2*TkeV**0.5 + 6.3d-2)
        else if (TkeV.le.0.02.and.TkeV.ge.0.0017235) then
          cooling_rate(i,j,k) = n_i*n_e*6.72e-22*(TkeV/0.02)**0.6
        else
          cooling_rate(i,j,k) = n_i*n_e*1.544e-22*(TkeV/0.0017235)**6.0
        endif

        if (TkeV.le.1.0d-4) cooling_rate(i,j,k) = 0.0d0

        dvol = dvl1a(i)*dvl2a(j)*dvl3a(k)
        vol(i) = vol(i) + dvol

        rad_cool(i) = rad_cool(i) + dvol*cooling_rate(i,j,k)
c
c
c turning off heating for now
c
c        const_heat(i) = const_heat(i) + dvol*cooling_rate(i,j,k)
c
         const_heat(i) = 0.0d0
      enddo
      enddo
c       
c FB htg. can in reality be cooling if conductive heating larger than
cooling
c
c        write(*,*) vol, (rad_cool(i)/float(n)), q_cond_tot(i)
c        heating_rate(i) = max((rad_cool(i)/float(n))-q_cond_tot(i)
c     &  ,0.0)/vol 

c      if (xbrag) const_heat(i) = const_heat(i) - q_cond_tot(i)  

      enddo

#ifdef MPI_USED
      buf_send(1:ie-is+1) = const_heat(is:ie)
      buf_send(ie-is+2:2*(ie-is)+2) = vol(is:ie) 

      call MPI_ALLREDUCE( buf_send, buf_recv, 2*(ie-is+1)
     &, MPI_DOUBLE_PRECISION, MPI_SUM, comm2d, ierr)

      const_heat(is:ie) = buf_recv(1:ie-is+1)
      vol(is:ie) = buf_recv(ie-is+2:2*(ie-is)+2)
#endif

      if (xbrag) const_heat(is:ie) = const_heat(is:ie) 
     & - q_cond_tot(is:ie)

      feed_heat(is:ie) = feed_heat(is:ie) + const_heat(is:ie)

      const_heat(is:ie) = const_heat(is:ie)/vol(is:ie)

      do k = ks, ke
      do j = js, je
      do i = is, ie
c
c combine heating rate w cooling rate
c
        cooling_rate(i,j,k) = cooling_rate(i,j,k) - const_heat(i)
c
c        if ( e(i,j,k).lt.0.0 ) write(*,*) 'cooling: eold<0', i, j, k
 
        if (cooling_rate(i,j,k).lt.0.0) then

        e(i,j,k) = e(i,j,k) - dtcool*cooling_rate(i,j,k)
c
c cooling done implicitly
c
        else

        e(i,j,k) = e(i,j,k)/( 1.0+dtcool*cooling_rate(i,j,k)/e(i,j,k) )

        endif

c        if ( e(i,j,k).lt.0.0 ) then
c          write(*,*) 'cooling: enew<0', i, j, k
c          stop
c        endif

      enddo 
      enddo
      enddo
c
#ifdef MPI_USED
      if (coords(1).eq.ntiles(1)-1) then
#endif
        e(ie,:,:) = e_outer
#ifdef MPI_USED
      endif
#endif

      enddo!subcycle

      rad_cool=rad_cool/float(nsubcycle)
      feed_heat=feed_heat/float(nsubcycle)

#ifdef MPI_USED
      buf_send(1:(ie-is)+1) = feed_heat(is:ie)
      buf_send((ie-is)+2:2*(ie-is)+2) = rad_cool(is:ie)

      call MPI_ALLREDUCE( buf_send(1), buf_recv(1), 2*(ie-is+1)
     &, MPI_DOUBLE_PRECISION, MPI_SUM, comm2d, ierr)

      feed_heat(is:ie) = buf_recv(1:(ie-is)+1)
      rad_cool(is:ie) = buf_recv((ie-is)+2:2*(ie-is)+2)
#endif

c
c making the boundary values out of date
c
      do i = 1,6
       bvstat(i,2) = 0      !  e
      enddo
c
      return
      end
