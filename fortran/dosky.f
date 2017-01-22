      SUBROUTINE dosky
c
c determine sky for this frame
c mod 04-Aug-98 aah fix tek1k scattered light prob
c mod 02-Oct-99 aah 20K obj
c
      INTEGER MAXOBJ, NPTS, MAXXPTS, MAXYPTS
      PARAMETER (MAXOBJ=600000)
      PARAMETER ( NPTS=200000)
      PARAMETER ( MAXXPTS=4100)
      PARAMETER ( MAXYPTS=4100)
      REAL*8
     $  xra(MAXOBJ),xdec(MAXOBJ)
      REAL*4
     $  imbuf(MAXXPTS,MAXYPTS),
     $  thresh,rdnoise,skymn,skymed,
     $  coords(MAXOBJ,2),dx, pixel(NPTS),sky,sigsky,skew,
     $  apzero,gain,xfwhm(MAXOBJ),xmag(MAXOBJ),xerr(MAXOBJ)
      INTEGER
     $  nx,ny,nc, n, kk, jj,iap,iskyinr,iskyoutr,nap,np,
     $  ich
      CHARACTER cra*13,cdec*13
      COMMON /apblk/ apzero,rdnoise,gain,iap,iskyinr,iskyoutr,nap,ich
      COMMON /findblk/ coords,xra,xdec,xfwhm,xmag,xerr,nc,cra,cdec
      COMMON /skyblk/ thresh,sky,sigsky
      COMMON /imblk/ imbuf,nx,ny
c
c use 200,000 random points for determining sky background
c
      np = npts
      dx = float(nx*ny)/float(NPTS+1)
      do i=1,NPTS
        n = int(dx*float(i))
        n = min(n,nx*ny)
        ix = n/nx
        ix = n - ix*nx + 1
        iy = (n-ix-1)/nx + 1
        pixel(i) = imbuf(ix,iy)
      enddo
      call sort1 (pixel,np)
      call mmm (pixel,np,hibad,rdnoise,skymn,skymed,sky,
     $  sigsky,skew)
c correct error
      if (sigsky.lt.0.0) sigsky = 0.0
      return
      end

      FUNCTION NBLANK (string)
c
c find first blank character
c
      INTEGER
     *       nblank, i
      CHARACTER*1
     *       string(*)
c
      nblank = 0
      DO i=1,256
        nblank = nblank + 1
        IF (string(nblank).eq.' ') GOTO 100
      ENDDO
100   continue
      nblank = nblank - 1
      RETURN
      END

      SUBROUTINE SORT (n,x,index)
c
c heapsort of array x with corresponding index array index
c from numerical recipies
c
      DIMENSION x(n),index(n)
      k = n/2+1
      ir = n
10    continue
      IF (k.gt.1) THEN
        k = k-1
        xx = x(k)
        ii = index(k)
      ELSE
        xx = x(ir)
        ii = index(ir)
        x(ir) = x(1)
        index(ir) = index(1)
        ir = ir-1
        IF (ir.eq.1) THEN
          x(1) = xx
          index(1) = ii
          RETURN
        ENDIF
      ENDIF
      i = k
      j = k+k
20    IF (j.le.ir) THEN
        IF (j.lt.ir) THEN
          IF (x(j).lt.x(j+1)) j = j+1
        ENDIF
        IF (xx.lt.x(j)) THEN
          x(i) = x(j)
          index(i) = index(j)
          i = j
          j = j+j
        ELSE
          j = ir+1
        ENDIF
      GOTO 20
      ENDIF
      x(i) = xx
      index(i) = ii
      GOTO 10
      END
