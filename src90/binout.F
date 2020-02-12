c=======================================================================
c/////////////////////////       BINOUT       \\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine binout(filename)
c
c
c
c  UNFORMATTED 3D WRITES OF important field variables 
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
      character*15 :: filename
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
      open(unit=42,file=binfile,status='unknown'
     &,form='unformatted',convert='big_endian')

      write(42) (((real(d(i,j,k)), i=is,ie), j=js,je), k=ks,ke)
      write(42) (((real(e(i,j,k)), i=is,ie), j=js,je), k=ks,ke)
      write(42) (((real(v1(i,j,k)), i=is,ie), j=js,je), k=ks,ke)
      write(42) (((real(v2(i,j,k)), i=is,ie), j=js,je), k=ks,ke)
      write(42) (((real(v3(i,j,k)), i=is,ie), j=js,je), k=ks,ke)
      if (xmhd) write(42) 
     &(((real(b1(i,j,k)), i=is,ie), j=js,je), k=ks,ke)
      if (xmhd) write(42) 
     &(((real(b2(i,j,k)), i=is,ie), j=js,je), k=ks,ke)
      if (xmhd) write(42) 
     &(((real(b3(i,j,k)), i=is,ie), j=js,je), k=ks,ke)
      if (xcosmic) write(42) 
     &(((real(ecr(i,j,k)), i=is,ie), j=js,je), k=ks,ke)

      close(42) 
c
      return
      end
