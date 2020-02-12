c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                 CRDIFFUSE                 \\\\\\\\\\
c
c                            Developed by Prateek Sharma
c
c  called by driver crdiffuse_d. applying CR diffusion on selected grid 
c  points.
c
c=======================================================================
c
      subroutine crdiffuse(ibeg,iend,jbeg,jend,kbeg,kend
     &,ecrold,ecrnew,dtdiff)
c
c......................................................................
      use real_prec
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
      use config
      use cosmic 
c
      implicit NONE
c
      integer  :: i, j, k
      integer  :: ibeg,iend,jbeg,jend,kbeg,kend 
      real(rl) :: bh1(in,jn,kn), bh2(in,jn, kn), bh3(in,jn,kn)   
      real(rl) :: bh2cenx1(in,jn,kn), bh3cenx1(in,jn,kn)
      real(rl) :: bh1cenx2(in,jn,kn), bh3cenx2(in,jn,kn)
      real(rl) :: bh1cenx3(in,jn,kn), bh2cenx3(in,jn,kn) 
      real(rl) :: b1cenx2, b1cenx3, b2cenx1, b2cenx3, b3cenx1, b3cenx2
      real(rl) :: bctr, q1(in,jn,kn), q2(in,jn,kn), q3(in,jn,kn)
      real(rl) :: divq  
      real(rl) :: dtdiff
      real(rl) :: decrx1, decrx2, decrx3
      real(rl) :: lim 
      real(rl) :: ecrold(in,jn,kn), ecrnew(in,jn,kn)
c
c----------------------------------------------------------------------
c compute magnetic field unit vectors
c
c----------------------------------------------------------------------
c

       do k = kbeg, kend
       do j = jbeg, jend
       do i = ibeg, iend+1

         b2cenx1 = 0.25*(b2(i,j,k)+b2(i,j+1,k)
     &   +b2(i-1,j,k)+b2(i-1,j+1,k))
         b3cenx1 = 0.25*(b3(i,j,k)+b3(i,j,k+1)
     &   +b3(i-1,j,k)+b3(i-1,j,k+1))
         bctr = sqrt(b1(i,j,k)*b1(i,j,k)
     &   +b2cenx1*b2cenx1+b3cenx1*b3cenx1)
         bh1(i,j,k) = b1(i,j,k)/(bctr+tiny)
         bh2cenx1(i,j,k) = b2cenx1/(bctr+tiny)
         bh3cenx1(i,j,k) = b3cenx1/(bctr+tiny) 

         Dcr(i,j,k) = alpha_cr*x1a(i)*bctr 
     &   /sqrt(0.5*(d(i,j,k)+d(i-1,j,k)))

         decrx2 = lim( lim(ecrold(i,j+1,k)-ecrold(i,j,k)
     &            ,ecrold(i-1,j+1,k)-ecrold(i-1,j,k))
     &            , lim(ecrold(i,j,k)-ecrold(i,j-1,k)
     &            ,ecrold(i-1,j,k)-ecrold(i-1,j-1,k)) )
         decrx3 = lim( lim(ecrold(i,j,k+1)-ecrold(i,j,k)
     &            ,ecrold(i-1,j,k+1)-ecrold(i-1,j,k))
     &            , lim(ecrold(i,j,k)-ecrold(i,j,k-1)
     &            ,ecrold(i-1,j,k)-ecrold(i-1,j,k-1)) )
         decrx1 = (ecrold(i,j,k)-ecrold(i-1,j,k))*dx1bi(i)
         decrx2 = decrx2*g2ai(i)*dx2ai(j)
         decrx3 = decrx3*g31ai(i)*g32bi(j)*dx3ai(k)

         q1(i,j,k) = -Dcr(i,j,k)*bh1(i,j,k)
     &   *( bh1(i,j,k)*decrx1 + bh2cenx1(i,j,k)*decrx2
     &   + bh3cenx1(i,j,k)*decrx3 )

       enddo
       enddo
       enddo 

       do k = kbeg, kend
       do j = jbeg, jend+1
       do i = ibeg, iend

         b1cenx2 = 0.25*(b1(i,j,k)+b1(i+1,j,k)
     &   +b1(i,j-1,k)+b1(i+1,j-1,k))
         b3cenx2 = 0.25*(b3(i,j,k)+b3(i,j,k+1)
     &   +b3(i,j-1,k)+b3(i,j-1,k+1))
         bctr = sqrt(b1cenx2*b1cenx2+b2(i,j,k)*b2(i,j,k)
     &   +b3cenx2*b3cenx2)
         bh1cenx2(i,j,k) = b1cenx2/(bctr+tiny)
         bh2(i,j,k) = b2(i,j,k)/(bctr+tiny)
         bh3cenx2(i,j,k) = b3cenx2/(bctr+tiny)
         
         Dcr(i,j,k) = alpha_cr*x1b(i)*bctr
     &   /sqrt(0.5*(d(i,j,k)+d(i,j-1,k)))

         decrx1 = lim( lim(ecrold(i+1,j,k)-ecrold(i,j,k)
     &         ,ecrold(i+1,j-1,k)-ecrold(i,j-1,k))
     &         , lim(ecrold(i,j,k)-ecrold(i-1,j,k)
     &         ,ecrold(i,j-1,k)-ecrold(i-1,j-1,k)) )
         decrx3 = lim( lim(ecrold(i,j,k+1)-ecrold(i,j,k)
     &         ,ecrold(i,j-1,k+1)-ecrold(i,j-1,k))
     &         , lim(ecrold(i,j,k)-ecrold(i,j,k-1)
     &         ,ecrold(i,j-1,k)-ecrold(i,j-1,k-1)) )
         decrx1 = decrx1*dx1ai(i)
         decrx2 = (ecrold(i,j,k)-ecrold(i,j-1,k))*g2bi(i)*dx2bi(j)
         decrx3 = decrx3*g31bi(i)*g32ai(j)*dx3ai(k)
c
c keeping a spherical polar coordinate system in mind; avoid singularities at poles
c
          if (x2a(j).le.1d-10.and.x2a(j).ge.-1d-10) decrx3 = 0.0d0
          if (x2a(j).le.pi+1d-10.and.x2a(j).ge.pi-1d-10) decrx3 = 0.0d0
          q2(i,j,k) = -Dcr(i,j,k)*bh2(i,j,k)
     &   *( bh1cenx2(i,j,k)*decrx1 + bh2(i,j,k)*decrx2
     &   + bh3cenx2(i,j,k)*decrx3 )

       enddo
       enddo 
       enddo

       do k = kbeg, kend+1
       do j = jbeg, jend
       do i = ibeg, iend

         b1cenx3 = 0.25*(b1(i,j,k)+b1(i+1,j,k)
     &   +b1(i,j,k-1)+b1(i+1,j,k-1))
         b2cenx3 = 0.25*(b2(i,j,k)+b2(i,j+1,k)
     &   +b2(i,j,k-1)+b2(i,j+1,k-1))
         bctr = sqrt(b1cenx3*b1cenx3+b2cenx3*b2cenx3
     &   +b3(i,j,k)*b3(i,j,k))
         bh1cenx3(i,j,k) = b1cenx3/(bctr+tiny)
         bh2cenx3(i,j,k) = b2cenx3/(bctr+tiny)
         bh3(i,j,k) = b3(i,j,k)/(bctr+tiny)

         Dcr(i,j,k) = alpha_cr*x1b(i)*bctr
     &   /sqrt(0.5*(d(i,j,k)+d(i,j,k-1)))

         decrx1 = lim( lim(ecrold(i+1,j,k)-ecrold(i,j,k)
     &            ,ecrold(i+1,j,k-1)-ecrold(i,j,k-1))
     &            , lim(ecrold(i,j,k)-ecrold(i-1,j,k)
     &            ,ecrold(i,j,k-1)-ecrold(i-1,j,k-1)) )
         decrx2 = lim( lim(ecrold(i,j+1,k)-ecrold(i,j,k)
     &          ,ecrold(i,j+1,k-1)-ecrold(i,j,k-1))
     &          , lim(ecrold(i,j,k)-ecrold(i,j-1,k)
     &          ,ecrold(i,j,k-1)-ecrold(i,j-1,k-1)) )
         decrx1 = decrx1*dx1ai(i)
         decrx2 = decrx2*g2bi(i)*dx2ai(j)
         decrx3 = (ecrold(i,j,k)-ecrold(i,j,k-1))*g31bi(i)*g32bi(j)
     &  *dx3bi(k)
         q3(i,j,k) = -Dcr(i,j,k)*bh3(i,j,k)
     &  *( bh1cenx3(i,j,k)*decrx1 + bh2cenx3(i,j,k)*decrx2
     &  + bh3(i,j,k)*decrx3 )

        enddo
        enddo
        enddo
     
c
c update cosmic ray internal energy
c
       do k = kbeg, kend
       do j = jbeg, jend
       do i = ibeg, iend
         divq = (g2a(i+1)*g31a(i+1)*q1(i+1,j,k)
     &   -g2a(i)*g31a(i)*q1(i,j,k))*dvl1ai(i)
     &   + (g32a(j+1)*q2(i,j+1,k)-g32a(j)*q2(i,j,k))*dvl2ai(j)*g2bi(i)
     &   + (q3(i,j,k+1)-q3(i,j,k))*dvl3ai(k)*g31bi(i)*g32bi(j)
         ecrnew(i,j,k) = ecrold(i,j,k) - dtdiff*divq
       enddo
       enddo
       enddo
c
      return
      end
