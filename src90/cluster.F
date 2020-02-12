c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 CLUSTER                   \\\\\\\\\\
c
c=======================================================================
c
      subroutine cluster
c
c initialize an MHD cluster in hydrodynamic equilibrium
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use field
      use grid
      use cons
      use metallicity
      use gravmod
      use root
      use domain
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
      real(rl) :: a3(in+1,jn+1)
      real(rl) :: rad, pres, gpot, den, pnot, pini, rmin
      real(rl) :: x, y, z, prefac
      integer :: iseed
      real(rl) :: rvalue
      namelist  /pgen/ x, y, z, prefac
c
c-----------------------------------------------------------------------
c
c defining the properties of the central compact object; its located at origin
c
c all quantities are in real cgs units, not in some other normalized form
c refer to Guo & Peng-Oh for the cluster parameters
c
       guniv = 6.673d-8
       x = 0.7
       y = 0.28
       z = 0.02
       prefac = 1.0d17

       v1=0.0d0; v2=0.0d0; v3=0.0d0
c
c density, pressure, phi_dm are read from a file
c
       rad=0.0d0
       i = is-2
c
#ifdef MPI_USED
       do j = 0, nprocs-1
        if (myid.eq.j) then
#endif
        open(unit=15)
        do while (rad<=x1b(ie+3))
          read(15,*) rad, den, pres, gpot
          if (rad>=x1b(i) .and. i<=ie+3) then
            d(i,:,:) = den
            p(i,:,:) = pres
            phi_dm(i) = gpot
            i=i+1
          endif
        enddo
        close(15)
#ifdef MPI_USED
        endif
       enddo
#endif
       e = p/gamm1
c
c
c set up a high-beta initial field; prefac is chosen in this way
c
c
       pnot = p(is,(js+je)/2,(ks+ke)/2)
       rmin = x1a(is)

       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
#ifdef MPI_USED
         buf_in(1) = x
         buf_in(2) = y
         buf_in(3) = z
         buf_in(4) = prefac
         buf_in(5) = p(is,(js+je)/2,(ks+ke)/2)
         buf_in(6) = x1a(is)
#endif
       endif
#ifdef MPI_USED
        call MPI_BCAST( buf_in, 6, MPI_FLOAT
     &                , 0, comm3d, ierr )
        if (myid .ne. 0) then
          x = buf_in(1)
          y = buf_in(2)
          z = buf_in(3)
          prefac = buf_in(4)
c
c pressure at r=rin
c
          pnot = buf_in(5)
          rmin = buf_in(6)
        endif
#endif

       mu = 1.0d0/(2.0*x+0.75*y+0.5625*z)
       mue = 2.0d0/(1.0d0+x)

      if (xmhd) then
c       do j = 1, jn+1
c       do i = 1, in+1
c         a3(i,j) = prefac*x1a(i)*sin(x2a(j))/(6.0d23) !roughly rout=200 kpc
c       enddo
c       enddo

       b2=0.0d0
       b3=0.0d0

       do k = 1, kn
       do j = 1, jn
       do i = 1, in
c         b1(i,j,k) = (g32a(j+1)*a3(i,j+1)-g32a(j)*a3(i,j))/(x2a(j+1)-
c     &    x2a(j))/g31a(i)/g32b(j)
c         b2(i,j,k) = -(g2a(i+1)*a3(i+1,j)-g2a(i)*a3(i,j))/(x1a(i+1)-
c     &    x1a(i))/g2b(i)
c         b3(i,j,k) = 0.0
         if (x2b(j).gt.0.5*pi) then
           b1(i,j,k) = prefac*1.d-12*x1ai(i)*x1ai(i)*(6.d23*6.d23)
         else
           b1(i,j,k) = -prefac*1.d-12*x1ai(i)*x1ai(i)*(6.d23*6.d23)
         endif
       enddo
       enddo
       enddo
      endif
c
c start with very small cosmic ray pressure
c    
      if (xcosmic) then
        j=(js+je)/2
        k=(ks+ke)/2 
        do i=1, in
c does not look steep enough
          pcr(i,:,:) = 0.005*pnot*(x1b(i)/rmin)**(-3.0)
c          pcr(i,:,:) = p(i,:,:)*0.2
          ecr(i,:,:) = pcr(i,:,:)/(gamma_cr-1.0)
          pini = p(i,j,k)
          p(i,:,:) = p(i,:,:) - pcr(i,:,:)
          e(i,:,:) = p(i,:,:)/(gamma-1.0)
          d(i,:,:) = d(i,:,:)*p(i,:,:)/pini
        enddo
      endif
c
c specifying abundance data in the innermost areas
c
      if (nspec.gt.1) then
      abun = 1.0d-15 
#ifdef MPI_USED
      if (coords(1).eq.0) then 
#endif
      do k=ks-2, ke+3
      do j=js-2, je+3
        abun(is-2:is,j,k,1) = 1.0d15
        abun(is-2:is+4,j,k,2) = 1.0d15
      enddo 
      enddo
#ifdef MPI_USED
      endif 
#endif
      endif 
c
c introduce large amplitude, random perturbations in density at small scales
c
      iseed = 1
      do k = ks, ke
      do j = js, je
      do i = is, ie
        call ran2(iseed,rvalue)
        d(i,j,k) = d(i,j,k)*(1.0 + 0.0*(rvalue-0.5) )
      enddo
      enddo
      enddo
c
c
c
c      j = (1+jn)/2
c      k = (1+kn)/2
c      do i = 1, in
c        write(10*(1+myid),1234) x1b(i), d(i,j,k), e(i,j,k), v1(i,j,k)
c     &, b1(i,j,k), b2(i,j,k), phi_dm(i)
c      enddo
c1234  format(6e20.7)
c
c      stop
c
      return
      end


c====================================================================
c  numerical recipes random number generator ran2
c    requires input seed value=iseed
c    returns real random number=rvalue
c    Also updates iseed for next call 
c
      subroutine ran2(iseed,rvalue)

      integer iseed
      real*8 rvalue
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,jj,kk,iv(NTAB),iy
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
c
      idum=iseed
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 jj=NTAB+8,1,-1
          kk=idum/IQ1
          idum=IA1*(idum-kk*IQ1)-kk*IR1
          if (idum.lt.0) idum=idum+IM1
          if (jj.le.NTAB) iv(jj)=idum
11      continue
        iy=iv(1)
      endif
      kk=idum/IQ1
      idum=IA1*(idum-kk*IQ1)-kk*IR1
      if (idum.lt.0) idum=idum+IM1
      kk=idum2/IQ2
      idum2=IA2*(idum2-kk*IQ2)-kk*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      jj=1+iy/NDIV
      iy=iv(jj)-idum2
      iv(jj)=idum
      if(iy.lt.1)iy=iy+IMM1
      rvalue=min(AM*iy,RNMX)
      iseed=idum
      return
      end
