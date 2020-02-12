c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 ICM_INITIAL               \\\\\\\\\\
c
c=======================================================================
c
      subroutine icm_initial
c
c initialize an MHD cluster in hydrodynamic equilibrium. Using the entropy 
c profile prescription of Cavagnolo et al.  K=T_{keV}/n_e^{2/3}
c
c K(r) = K_0 + K_100 (r/100kpc)^\alpha
c
c-----------------------------------------------------------------------
c
      use real_prec
      use config
      use param
      use root
      use grid
      use field
      use gravmod
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use cons
      use metallicity
      use domain
      use d_e_out
      use jet_vars
c
      implicit none
c
      integer  :: i, j, k, l, m, n, l1, m1, n1
      real(rl) :: x, y, z, prefac, pout_r, pout_s
      real(rl) :: pg, mnot, rs, rgv, rhonot, K0, K100, alpha, ne_out
      real(rl) :: fn, fnp, err, const, mnot_msun, rs_kpc, r100
     &, rmax
      real(rl) :: Kent(in), r_jet_kpc, sig_r_kpc
      integer :: iseed, mag_flag
      real(rl) :: rvalue, amp, aklmX, aklmY, aklmZ, xs, ys, zs, phiX
     &,phiY, phiZ, aklm, phi
      real(rl) :: kx, ky, kz, delrho(in,jn,kn)
      real(rl), dimension(in,jn,kn) :: Ax,Ay,Az,Ar,At,Ap
      namelist  /pgen/ amp, x, y, z, prefac, ne_out, K0, K100, alpha
     &, mnot_msun, rs_kpc, r100, rmax, cond_prefac, ncond_lim, mag_flag
     &, th_jet, r_jet_kpc, sig_th, sig_r_kpc, v_jet, jet_eff
c
c-----------------------------------------------------------------------
c
c all quantities are in real cgs units, not in some other normalized form
c refer to Guo & Peng-Oh for the cluster parameters
c
       guniv = 6.673d-8; mnot = 7.6d47; rs = 1.204d24; rgv = 8.733d21
       rhonot = 5.629d-23
c metallicity
       x = 0.7; y = 0.28; z = 0.02
c strength of B-field
       prefac = 1.0
c profile
       K0 = 30.; K100 = 150.; r100 = 100.; rmax = 200.; alpha = 1.2
       ne_out = 0.0015
c NFW potential parameters
       rs_kpc = 390.; mnot_msun = 3.8e14
c magnitude of density perturbations
       amp = 0.1
c conduction/subcycling-limiter parameters
       cond_prefac = 1.0; ncond_lim = 1000
c magnetic geometry
       mag_flag = 1
c jet parameters
       r_jet_kpc=2.5d0; th_jet=15.0d0*pi/180.; v_jet=5.d9
       sig_r_kpc=0.05d0; sig_th=0.05d0; jet_eff=6.d-3
       if (myid .eq. 0) then
         read (1, pgen)
         write (2, pgen)
#ifdef MPI_USED
         buf_in(1)=x; buf_in(2)=y; buf_in(3)=z; buf_in(4)=prefac
         buf_in(5)=K0; buf_in(6)=K100; buf_in(7)=alpha; buf_in(8)=ne_out
         buf_in(9)=rs_kpc; buf_in(10)=mnot_msun; buf_in(11)=amp
         buf_in(12)=r100; buf_in(13)=rmax; buf_in(14)=cond_prefac
         buf_in(15)=r_jet_kpc; buf_in(16)=th_jet; buf_in(17)=v_jet
         buf_in(18)=sig_r_kpc; buf_in(19)=sig_th; buf_in(20)=jet_eff
         ibuf_in(1)=ncond_lim; ibuf_in(2)=mag_flag
#endif
       endif
#ifdef MPI_USED
        call MPI_BCAST( buf_in, 20, MPI_FLOAT
     &                , 0, comm3d, ierr )
        call MPI_BCAST( ibuf_in, 2, MPI_INTEGER
     &                , 0, comm3d, ierr )
        if (myid .ne. 0) then
         x=buf_in(1); y=buf_in(2); z=buf_in(3); prefac=buf_in(4)
         K0=buf_in(5); K100=buf_in(6); alpha=buf_in(7); ne_out=buf_in(8)
         rs_kpc=buf_in(9); mnot_msun=buf_in(10); amp=buf_in(11) 
         r100=buf_in(12); rmax=buf_in(13); cond_prefac=buf_in(14)
         r_jet_kpc=buf_in(15); th_jet=buf_in(16); v_jet=buf_in(17)
         sig_r_kpc=buf_in(18); sig_th=buf_in(19); jet_eff=buf_in(20)
         ncond_lim=ibuf_in(1); mag_flag=ibuf_in(2)
        endif
#endif

       mnot = mnot_msun*msol; rs = rs_kpc*3.086e21;
       r_jet=r_jet_kpc*cmkpc; sig_r=sig_r_kpc*cmkpc
       norm_jet = 3.0d0/(2.0d0*pi*r_jet*r_jet*r_jet*(1.0d0-cos(th_jet)))

       phi_dm = -2.0*guniv*mnot*log(1.0+x1b/rs)/x1b
     & -0*4.0*pi*guniv*rhonot*rgv**3
     & *log(x1b/rgv+sqrt(1.0+(x1b/rgv)**2))/x1b

       v1=0.0d0; v2=0.0d0; v3=0.0d0
c
c entropy specified in terms of parameters K0, alpha, K100
c
       Kent = K0 + K100*(x1b/3.086e21/r100)**alpha
       mu = 1.0d0/(2.0*x+0.75*y+0.5625*z)
       mue = 2.0d0/(1.0d0+x)

c conversion of K: K = const.*p/d^gamma ! this is correct taking K=TkeV/n_e^{2/3}
       const = (mu*mp*(mue*mp)**gamm1)/(boltz*1.1604d7)
c
c now set up hydrostatic profile
c
#ifdef MPI_USED
      if (coords(1).eq.ntiles(1)-1) then
#endif
        d(in,:,:) = ne_out*mp*mue

        p(in,:,:) = Kent(in)*d(in,:,:)**gamma/const

#ifdef MPI_USED
      endif
      nreq = 0

      do l = ntiles(1)-1, 0, -1
        if (l.eq.coords(1)) then

          if (l.ne.ntiles(1)-1) then
            nreq=nreq+1

            call MPI_RECV(pout_r, 1, MPI_DOUBLE_PRECISION, n1p, 12345
     &      , comm3d, req(nreq), ierr)
   
            p(ie+3,:,:) = pout_r
          endif
#endif
c
c pressure solverd by Neton-Raphson. outer pressure value is the initial guess
c
          do i = in-1, 1, -1

            err = huge
            j = (1+jn)/2; k = (1+kn)/2
            pg = p(i+1,j,k)

            do while (err.ge.1.0e-10)

            fn = p(i+1,j,k) - pg + 0.5*(phi_dm(i+1)-phi_dm(i))
     & *((p(i+1,j,k)/Kent(i+1))**(1./gamma) + (pg/Kent(i))**(1./gamma))
     & *const**(1./gamma)

            fnp = -1. + 0.5*(phi_dm(i+1)-phi_dm(i))*(1./gamma)
     & *pg**(1./gamma-1.)*(const/Kent(i))**(1./gamma)


            pg = pg - fn/fnp             

            err = abs(fn/fnp)/pg

            enddo

          p(i,:,:) = pg
          d(i,:,:) = (const*p(i,:,:)/Kent(i))**(1./gamma)

        enddo
#ifdef MPI_USED
          pout_s = p(is+2,j,k)

          if (l.ne.0) then
            nreq = nreq + 1

            call MPI_ISEND(pout_s, 1, MPI_DOUBLE_PRECISION, n1m, 12345
     &      , comm3d, req(nreq), ierr)
          endif

        endif

c        call MPI_BARRIER(comm3d, ierr)

      enddo
#endif


      e = p/gamm1
c
      j = (1+jn)/2; k = (1+kn)/2
#ifdef MPI_USED
      if (coords(1).eq.ntiles(1)-1) then
#endif

        d_outer = d(ie,j,k) 
        e_outer = e(ie,j,k)

#ifdef MPI_USED
      endif
#endif


c
c introduce large amplitude, random perturbations in density at small
c scales
c
      iseed = 1
      delrho = 0.0
c
        do n = 4, 20
        do n1 = -n, n, 2*n
          kx = 2.0*pi*n1/(2.*rmax*3.086e21)
        do l = 4, 20
        do l1 = -l, l, 2*l
          ky = 2.0*pi*l1/(2.*rmax*3.086e21)
        do m = 4, 20
        do m1 = -m, m, 2*m
          kz = 2.0*pi*m1/(2.*rmax*3.086e21)
c
c          write(*,*) n1, l1, m1
c
          call ran2(iseed,rvalue)
          phi = 2.0*pi*rvalue

c          call ran2(iseed,rvalue) 

c equal power at all k's
          call ran2(iseed,rvalue)


          aklm = 0.15*amp*(0.5-rvalue)/sqrt(1.*n**2+l**2+m**2)
c
        do k = ks, ke
        do j = js, je
        do i = is, ie

          xs = x1b(i)*sin(x2b(j))*cos(x3b(k))
          ys = x1b(i)*sin(x2b(j))*sin(x3b(k))
          zs = x1b(i)*cos(x2b(j))

          delrho(i,j,k) = delrho(i,j,k) + aklm*d(i,j,k)
     &    *cos(phi + kx*xs + ky*ys + kz*zs)

        enddo
        enddo
        enddo
c        
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      d = d + delrho

c      write(*,*) 'OK'

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c magnetic field initialization
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (xmhd) then

c mag_flag.eq.1 for radial B

        if (mag_flag.eq.1) then

          b2=0.0d0
          b3=0.0d0

          do k = 1, kn
          do j = 1, jn
          do i = 1, in
            if (x2b(j).gt.0.5*pi) then
              b1(i,j,k) = prefac*x1ai(i)*x1ai(i)*(3.086d21*3.086d21)
     &        *(200./rmax)**2*1.e-6
            else
              b1(i,j,k) =-prefac*x1ai(i)*x1ai(i)*(3.086d21*3.086d21)
     &        *(200./rmax)**2*1.e-6
            endif
          enddo
          enddo
          enddo

        endif 

c mag_flag.eq.2 for azimuthal B

        if (mag_flag.eq.2) then

          b1=0.0d0
          b2=0.0d0

          do k = 1, kn
          do j = 1, jn
          do i = 1, in
            b3(i,j,k)=prefac*x1bi(i)*3.086e21*(200./rmax)*1.e-6
          enddo
          enddo
          enddo

        endif 

        if (mag_flag.eq.3) then

        if (ldimen.lt.3) then 
          write(*,*) 'tangled B-fields only allowed for 3-D'
          stop
        endif

        iseed = -1001
        Ax = 0.0
        Ay = 0.0
        Az = 0.0

        do n = 4, 10
        do n1 = -n, n, 2*n
          kx = 2.0*pi*n1/(2.*rmax*3.086e21)
        do l = 4, 10
        do l1 = -l, l, 2*l
          ky = 2.0*pi*l1/(2.*rmax*3.086e21)
        do m = 4, 10
        do m1 = -m, m, 2*m
          kz = 2.0*pi*m1/(2.*rmax*3.086e21)

          call ran2(iseed,rvalue)
          phiX = 2.0*pi*rvalue
          call ran2(iseed,rvalue)
          aklmX = (0.5-rvalue)/(1.*n**2+l**2+m**2)**1.33333
          call ran2(iseed,rvalue)
          phiY = 2.0*pi*rvalue
          call ran2(iseed,rvalue)
          aklmY = (0.5-rvalue)/(1.*n**2+l**2+m**2)**1.33333
          call ran2(iseed,rvalue)
          phiZ = 2.0*pi*rvalue
          call ran2(iseed,rvalue)
          aklmZ = (0.5-rvalue)/(1.*n**2+l**2+m**2)**1.33333

          do k = 1, kn
          do j = 1, jn
          do i = 1, in

            xs = x1b(i)*sin(x2b(j))*cos(x3b(k))
            ys = x1b(i)*sin(x2b(j))*sin(x3b(k))
            zs = x1b(i)*cos(x2b(j))
c
c prefac 1 should give beta~1000
c

            Ax(i,j,k) = Ax(i,j,k) + aklmX*cos(phiX + kx*xs +ky*ys+ 
     &      kz*zs)*1.e18*prefac

            Ay(i,j,k) = Ay(i,j,k) + aklmY*cos(phiY + kx*xs + ky*ys+
     &      kz*zs)*1.e18*prefac

            Az(i,j,k) = Az(i,j,k) + aklmZ*cos(phiZ + kx*xs + ky*ys+
     &      kz*zs)*1.e18*prefac
        
          enddo
          enddo
          enddo        

        enddo
        enddo
        enddo
        enddo
        enddo
        enddo

        do k = 1, kn
        do j = 1, jn
        do i = 1, in

          Ar(i,j,k) =Ax(i,j,k)*sin(x2a(j))*cos(x3a(k))+
     &       Ay(i,j,k)*sin(x2a(j))*sin(x3a(k))+Az(i,j,k)*cos(x2a(j)) 

          At(i,j,k) =Ax(i,j,k)*cos(x2b(j))*cos(x3a(k))+
     &       Ay(i,j,k)*cos(x2b(j))*sin(x3a(k))-Az(i,j,k)*sin(x2b(j))

          Ap(i,j,k) =  -Ax(i,j,k)*sin(x3b(k))+Ay(i,j,k)*cos(x3b(k))
       
        enddo
        enddo
        enddo
c
c reducing B at the poles
c

c        do j = 1, jn

c          Ar(:,j,:) = Ar(:,j,:)*exp( -abs(x2b(j)-0.5*pi)*1.5 )
c          At(:,j,:) = At(:,j,:)*exp( -abs(x2b(j)-0.5*pi)*1.5 )
c          Ap(:,j,:) = Ap(:,j,:)*exp( -abs(x2b(j)-0.5*pi)*1.5 )

c        enddo


c reduce theta and radial B

        do j = 1, jn

          Ar(:,j,:) = Ar(:,j,:)*exp( -.5*(x2b(j)-0.5*pi)**2 )
          At(:,j,:) = At(:,j,:)*exp( -.5*(x2b(j)-0.5*pi)**2 )
          Ap(:,j,:) = Ap(:,j,:)*exp( -.5*(x2b(j)-0.5*pi)**2 )

        enddo

        do i = 1, in

          Ar(i,:,:) = Ar(i,:,:)*exp( -2.3*(x1b(i)/rmax/3.086e21)**2 )
          At(i,:,:) = At(i,:,:)*exp( -2.3*(x1b(i)/rmax/3.086e21)**2 )
          Ap(i,:,:) = Ap(i,:,:)*exp( -2.3*(x1b(i)/rmax/3.086e21)**2 )

        enddo

c
c following is needed for safety at the pole; symmetry
c
        do k = 1, kn-1
          if (coords(2).eq.0) Ar(:,js,k+1) = Ar(:,js,k)
          if (coords(2).eq.ntiles(2)-1) Ar(:,je+1,k+1) = Ar(:,je+1,k)
        enddo

        do i = 1, in-1
          if (coords(2).eq.0) Ap(i+1,js,:)=g31a(i)*Ap(i,js,:)/g31a(i+1)
          if (coords(2).eq.ntiles(2)-1) Ap(i+1,je+1,:) 
     &    = g31a(i)*Ap(i,je+1,:)/g31a(i+1)
        enddo

        do k = ks-2, ke+2
        do j = js-2, je+2
        do i = is-2, ie+2

          b1(i,j,k)=( (g32a(j+1)*Ap(i,j+1,k)-g32a(j)*Ap(i,j,k))*dx2ai(j)
     &    -(At(i,j,k+1)-At(i,j,k))*dx3ai(k) )*g32bi(j)*g31ai(i)
                
          b2(i,j,k)=(Ar(i,j,k+1)-Ar(i,j,k))*dx3ai(k)*g32ai(j)*g31bi(i) 
     &    - (g31a(i+1)*Ap(i+1,j,k)-g31a(i)*Ap(i,j,k))*g31bi(i)*dx1ai(i)
                 

          b3(i,j,k) =((g2a(i+1)*At(i+1,j,k)-g2a(i)*At(i,j,k))*
     &    dx1ai(i) - (Ar(i,j+1,k)-Ar(i,j,k))*dx2ai(j) )*g2bi(i)

        enddo
        enddo
        enddo

        endif !mag_flag

    
      endif !xmhd
c
      return
      end
