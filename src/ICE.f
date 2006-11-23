      double precision function kernel(u, ker)

c smoothing kernel

      implicit none
      double precision one,pi,two,ten,u,zero
      integer ker

      parameter(pi=3.1415927d0,one=1.0d0,two=2.0d0,ten=10.0d0,
     & zero=0.0d0)

      if (ker .eq. 1) then 
c gaussian kernel:
          if (abs(u) .le. ten) then
              kernel = dexp(-u**2/two)/dsqrt(two*pi)
          else
              kernel=zero
          endif
      else
          if (ker .eq. 2) then
c epanechnikov kernel:
              if (abs(u) .le. one) then
                  kernel = (1-u**2)*3.0/4.0
              else
                  kernel=zero
              endif
          else
c biweight kernel
              if (abs(u) .le. one) then
                  kernel = (1-u**2)**2*15.0/16.0
              else
                  kernel=zero
              endif
          endif
      endif

      return
      end

      subroutine ickde(n,left,right,numgrid,gridpts,f0,h,niter,f1,ker)

      implicit NONE
      INTEGER i, ii, j, k, m, n, numgrid, niter, ker
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
          if (abs(left(j)-right(j)) .ge. meshwidth) then
           do 60 k=1,m
            if ((left(j) .lt. gridpts(k)) .AND.
     &       (right(j) .gt. gridpts(k))) then
              denom = denom + f0(k)
              z = (gridpts(i)-gridpts(k))/h
              numer = numer + kernel(z, ker)*f0(k)/h
            end if
60         continue
          else
           z = (gridpts(i) - (left(j)+right(j))/2)/h
           numer = numer + kernel(z, ker)/h
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

      subroutine icllde(n,left,right,numgrid,gridpts,f0,h,niter,f,ker)

      implicit NONE
      INTEGER              i, ii, j, k, m, n, numgrid, niter, ker
      DOUBLE PRECISION     left(n), right(n), gridpts(numgrid), 
     &                     f0(numgrid), f(numgrid), fder, fnew, fdenom, 
     &h, meshwidth, fnumer, kf,
     & fdernumer, z, fsum
      
      DOUBLE PRECISION kernel

      m = numgrid

      meshwidth = gridpts(2)-gridpts(1)

      do 120 ii=1,niter

      do 100 j=1,m
          fnew = 0.0
          fder = 0.0
          do 80 i=1,n
              fdenom = 0.0
              fnumer = 0.0
              fdernumer = 0.0
              do 60 k=1,m
                  if ((left(i) .lt. gridpts(k)) .AND.
     &            (right(i) .GT. gridpts(k))) then
                      fdenom = fdenom + f0(k)
                      z = (gridpts(k)-gridpts(j))/h
                      kf = kernel(z,ker)*f0(k)
                      fnumer = fnumer + kf/h
                      fdernumer = fdernumer + kf*z
                  end if
60            continue
              if (fdenom .gt. .00000001) then
                  fnew = fnew + fnumer/fdenom
                  fder = fder + fdernumer/fdenom
              endif
80        continue
          f(j) = fnew*EXP(-((fder/(h*fnew))**2)/2.0)/n

100   continue
      fsum=0.0
      do 110 j=1,m
          f0(j) = f(j)
          fsum = fsum+f(j)*meshwidth
110   continue
120   continue

      do 150 j=1,m
          f(j)=f(j)/fsum
150   continue

      return
      end

