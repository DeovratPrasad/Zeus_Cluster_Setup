c=======================================================================
c
      subroutine floor_d 
c
c  PURPOSE: applies density floor  
c
c-----------------------------------------------------------------------
      use config
      use field
      use grid
      use root
      use bndry
c
      implicit NONE
c
      integer :: i, j, k
      real(rl) dtsnd, newd, valf
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c

c      do k = ks-2, ke+3
c      do j = js-2, je+3
c      do i = is-2, ie+3

      do k = ks, ke
      do j = js, je
      do i = is, ie

#ifdef NODEF

        if (d(i,j,k).lt.dfloor) then
c
c madd is the mass added in real zones (not including the ghost zones)
c
          if (i.ge.is.and.i.le.ie.and.j.ge.js.and.j.le.je 
     &        .and.k.ge.ks.and.k.le.ke)
     &      madd = madd + (dfloor-d(i,j,k))*dvl1a(i)*dvl2a(j)*dvl3a(k) 
c
c but the floor is applied everywhere, including the ghosts
c
            d(i,j,k) = dfloor
          endif

#endif
c#ifdef NODEF
c a more well-controlled floor should use dtsnd 

         dtsnd = min( min( dx1a(i), g2b(i)*dx2a(j)
     &   , g31b(i)*g32b(j)*dx3a(k))
     &   /sqrt( gamma*gamm1*e(i,j,k)/d(i,j,k) ), huge )

         if (xmhd) then

           valf = sqrt( ( (b1(i,j,k)+b1(i+1,j,k))**2
     & + (b2(i,j,k)+b2(i,j+1,k))**2 + (b3(i,j,k)+b3(i,j,k+1))**2 )
     & /d(i,j,k) )


           dtsnd = min( dtsnd, min( dx1a(i), g2b(i)*dx2a(j),
     &     g31b(i)*g32b(j)*dx3a(k))/valf )
     &

         endif

c         write(*,*) dtmin, dfloor

         if ( dtsnd.lt. 10.*dtmin .or. d(i,j,k).lt.dfloor ) then

c       if ( dtsnd.lt. 3.*dtmin ) then

c         write(*,*) dtsnd, dx1a(i), g2b(i)*dx2a(j),
c     &   g31b(i)*g32b(j)*dx3a(k) 

         newd = max( d(i,j,k)*(10.*dtmin/dtsnd)**2, dfloor )
 
c          newd = d(i,j,k)*(3.*dtmin/dtsnd)**2

c         write(*,*) i, j, k, newd

          if (i.ge.is.and.i.le.ie.and.j.ge.js.and.j.le.je 
     &        .and.k.ge.ks.and.k.le.ke)
     &      madd = madd + (newd-d(i,j,k))*dvl1a(i)*dvl2a(j)*dvl3a(k)
c
c but the floor is applied everywhere, including the ghosts
c
            d(i,j,k) = newd 

         endif
c#endif

      enddo 
      enddo 
      enddo
c
      bvstat(:,2) = 0
c
      return
      end

c=======================================================================
c
      subroutine floor_e
c
c  PURPOSE: applies internal energy floor
c
c-----------------------------------------------------------------------
      use field
      use grid
      use root
c
      implicit NONE
c
      integer :: i, j, k
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c

      do k = ks-2, ke+3
      do j = js-2, je+3
      do i = is-2, ie+3

        if (e(i,j,k).lt.efloor) then
c
c madd is the mass added in real zones (not including the ghost zones)
c
          if (i.ge.is.and.i.le.ie.and.j.ge.js.and.j.le.je
     &        .and.k.ge.ks.and.k.le.ke)
     &      eadd = eadd + (efloor-e(i,j,k))*dvl1a(i)*dvl2a(j)*dvl3a(k)
c
c but the floor is applied everywhere, including the ghosts
c
          e(i,j,k) = efloor
        endif

      enddo
      enddo
      enddo
c
      return
      end
c
c=======================================================================
c
      subroutine ceil_e
c
c  PURPOSE: applies internal energy ceiling
c
c-----------------------------------------------------------------------
      use field
      use grid
      use root
      use bndry
c
      implicit NONE
c
      integer :: i, j, k
      real(rl) :: dtsnd, newe
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c

c      do k = ks-2, ke+3
c      do j = js-2, je+3
c      do i = is-2, ie+3

      do k = ks, ke
      do j = js, je
      do i = is, ie

#ifdef NODEF

        if (e(i,j,k).gt.eceil) then
c
c madd is the mass added in real zones (not including the ghost zones)
c
          if (i.ge.is.and.i.le.ie.and.j.ge.js.and.j.le.je
     &        .and.k.ge.ks.and.k.le.ke)
     &      eadd = eadd + (eceil-e(i,j,k))*dvl1a(i)*dvl2a(j)*dvl3a(k)
c
c but the floor is applied everywhere, including the ghosts
c
          e(i,j,k) = eceil
        endif

#endif

         dtsnd = min( min( dx1a(i), g2b(i)*dx2a(j) 
     &   , g31b(i)*g32b(j)*dx3a(k))
     &   /sqrt( gamma*gamm1*e(i,j,k)/d(i,j,k) ), huge )


         if ( dtsnd.lt. 10.*dtmin .or. e(i,j,k).gt.eceil ) then

c         if ( dtsnd.lt. 3.*dtmin ) then

          newe = min( e(i,j,k)/(10.*dtmin/dtsnd)**2, eceil )

c          newe = e(i,j,k)/(3.*dtmin/dtsnd)**2

          if (i.ge.is.and.i.le.ie.and.j.ge.js.and.j.le.je
     &        .and.k.ge.ks.and.k.le.ke)
     &      eadd = eadd + (newe-e(i,j,k))*dvl1a(i)*dvl2a(j)*dvl3a(k)

          e(i,j,k) = newe

         endif  

      enddo
      enddo
      enddo
c
      bvstat(:,3) = 0
c
      return
      end
