c=======================================================================
c/////////////////////////       OUTPUT2D       \\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine output2d 
c
c
c
c  FORMATTED 2D (phi-averaged) WRITES OF important field variables 
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
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
c
      implicit NONE
c
      integer :: i, j, k 
      real(rl) :: d2D(in,jn), e2D(in,jn), v12D(in,jn), v22D(in,jn)
      real(rl) :: v32D(in,jn), b12D(in,jn), b22D(in,jn), b32D(in,jn)
      real(rl) :: ecr2D(in,jn) 
      real(rl) :: dvol, vol
      real(rl) :: abun2D_1(in,jn), abun2D_2(in,jn) 
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
      do j = js, je
      do i = is, ie
      d2D(i,j)=0.0d0; e2D(i,j)=0.0d0; v12D(i,j)=0.0d0 
      v22D(i,j)=0.0d0; v32D(i,j)=0.0d0
      b12D(i,j)=0.0d0; b22D(i,j)=0.0d0; b32D(i,j)=0.0d0
      ecr2D(i,j)=0.0d0; abun2D_1(i,j)=0.0d0; abun2D_2(i,j)=0.0d0 
      vol=0.0d0
      do k = ks, ke
        dvol = dvl3a(k)
        vol = vol + dvol 
        d2D(i,j) = d2D(i,j) + dvol*d(i,j,k)
        e2D(i,j) = e2D(i,j) + dvol*e(i,j,k)
        v12D(i,j) = v12D(i,j) + dvol*0.5*(v1(i,j,k)+v1(i+1,j,k))
        v22D(i,j) = v22D(i,j) + dvol*0.5*(v2(i,j,k)+v2(i,j+1,k))
        v32D(i,j) = v32D(i,j) + dvol*0.5*(v3(i,j,k)+v3(i,j,k+1))
        if (xmhd) then
          b12D(i,j) = b12D(i,j) + dvol*0.5*(b1(i,j,k)+b1(i+1,j,k))
          b22D(i,j) = b22D(i,j) + dvol*0.5*(b2(i,j,k)+b2(i,j+1,k))
          b32D(i,j) = b32D(i,j) + dvol*0.5*(b3(i,j,k)+b3(i,j,k+1))
        endif
        if (xcosmic) then
          ecr2D(i,j) = ecr2D(i,j) + dvol*ecr(i,j,k)
        endif
        if (nspec.gt.1) then
          abun2D_1(i,j) = abun2D_1(i,j) + dvol*abun(i,j,k,1)
          abun2D_2(i,j) = abun2D_2(i,j) + dvol*abun(i,j,k,2)
        endif 
      enddo
      enddo
      enddo
      d2D=d2D/(vol+tiny); e2D=e2D/(vol+tiny); v12D=v12D/(vol+tiny) 
      v22D=v22D/(vol+tiny); v32D=v32D/(vol+tiny); b12D=b12D/(vol+tiny)
      b22D=b22D/(vol+tiny); b32D=b32D/(vol+tiny); ecr2D=ecr2D/(vol+tiny)
      abun2D_1=abun2D_1/(vol+tiny); abun2D_2=abun2D_2/(vol+tiny)
c      open(unit=42,file=twoDfile,status='unknown')
c      write(42,*) 'time = ', time, ' nhy =', nhy
c      write(42,*) 

      do j = js, je
      do i = is, ie
        write(42,2001) x1b(i), x2b(j), d2D(i,j), e2D(i,j), v12D(i,j)
     &  , v22D(i,j), v32D(i,j), b12D(i,j), b22D(i,j), b32D(i,j)
     &  , ecr2D(i,j), abun2D_1(i,j), abun2D_2(i,j)
      enddo
      enddo

      close(42) 
2001  format(13e20.7)
c
      return
      end
