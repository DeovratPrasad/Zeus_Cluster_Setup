! a program to take averages of differential emission measure
      PROGRAM stitch
      implicit NONE
      
      integer :: incr, iincr, l, p
      integer, parameter :: nbins=1000
      real*8 :: dem_avg(nbins), log10TkeV(nbins), dem(nbins)

      do iincr = 0, 5

      dem_avg = 0.0
      do incr = max(0,(iincr-1))*20, (iincr)*20
  
        write(*,*) incr       

        open( 100*(incr+1) ) 
        do l = 1, nbins-1 
          read(100*(incr+1),2001) p, log10TkeV(l), dem(l)
          dem_avg(l) = dem_avg(l) + dem(l)
        enddo

        close( 100*(incr+1) ) 

      enddo

      dem_avg = dem_avg/((iincr)*20-max(0,(iincr-1))*20+1)

      do l = 1, nbins-1
        write((iincr+1)*10,2001) l, log10TkeV(l), dem_avg(l)
      enddo

      enddo

2001  format(i7,2e20.7)

      END
