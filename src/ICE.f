      double precision function kernel(u)

c smoothing kernel

      implicit none
      double precision one,pi,two,ten,u,zero

      parameter(pi=3.1415927d0,one=1.0d0,two=2.0d0,ten=10.0d0,
     & zero=0.0d0)

c gaussian kernel:
      if (abs(u) .le. ten) then
       kernel = dexp(-u**2/two)/dsqrt(two*pi)
      else
       kernel=zero
      endif

      return
      end

      subroutine ickde(n,left,right,numgrid,gridpts,f0,h,niter,f1)

      implicit NONE
      INTEGER i, ii, j, k, m, n, numgrid, niter
      DOUBLE PRECISION left(n), right(n), gridpts(numgrid), f0(numgrid),
     & f1(numgrid), fnew, denom, h, meshwidth, numer, z
c function calls:
      DOUBLE PRECISION kernel

      m = numgrid

      meshwidth = gridpts(2)-gridpts(1)

      do 120 ii=1,niter

      do 100 i=1,m
        fnew = 0.0
        do 80 j=1,n
          denom = 0.0
          numer = 0.0
          if (abs(left(j)-right(j)) .ge. .00000001) then
           do 60 k=1,m
            if ((left(j) .lt. gridpts(k)) .AND.
     &       (right(j) .gt. gridpts(k))) then
              denom = denom + f0(k)
              z = (gridpts(i)-gridpts(k))/h
              numer = numer + kernel(z)*f0(k)/h
            end if
60         continue
          else
           z = (gridpts(i) - (left(j)+right(j))/2)/h
           numer = numer + kernel(z)/h
           denom = 1.0
          endif  
          if (abs(denom) .gt. .00000001) then
              fnew = fnew + numer/denom
          endif
80    continue
        f1(i) = fnew/n

100   continue
        do 110 i=1,m
          f0(i) = f1(i)
110   continue
120   continue

      return
      end

