      SUBROUTINE dophot
c
c perform concentric aperture photometry on image
c plagerized from Stetson
c 21-Dec-95 aah
c mod 04-Aug-98 aah add tek1k scattered light fix
c mod 02-Oct-99 aah 20K obj
c mod 02-Aug-2016 aah SEL/SExtractor format
c
      INTEGER MAXXPTS, MAXYPTS, MAXOBJ
      INTEGER MINSKY, MAXSKY, MAXAP
c note MAXXPTS/YPTS permit 4k4k detector maximum
      PARAMETER (MAXXPTS=4100)
      PARAMETER (MAXYPTS=4100)
c probably overkill; typical image has about 11K objects, but don't know about the plane
      PARAMETER (MAXOBJ=600000)
      PARAMETER (MAXSKY=1000 )
      PARAMETER (MINSKY=20)
      PARAMETER ( MAXAP=9)
      REAL*4
     $    imbuf(MAXXPTS,MAXYPTS),xc,yc,par(MAXAP),
     $    ut,exptime, coords(MAXOBJ,2),apzero,rdnoise,
     $    sky,fwhm,lobad,hibad,gain,skysig,skymod,skyvar,
     $    fwhmx(MAXOBJ), fwhmy(MAXOBJ), wscale, thresh,sigsky,
     $    apsky(MAXSKY), apmax, err1, err2, err3, dum,
     $    skymn,skymed,skew,epoch,sharp(MAXOBJ),
     $    ratio(1024,1024),xmag(MAXOBJ),xerr(MAXOBJ),
     $    xfwhm(MAXOBJ)
      REAL*8
     $    apmag(MAXAP), area(MAXAP), aperr(MAXAP),
     $    xra(MAXOBJ),xdec(MAXOBJ)
      INTEGER
     $    i, j, k, n, nx, ny, nc, ixdir, iydir,iap(MAXAP),nap,
     $    iskyinr,iskyoutr, izdir,ich, ix, iy, nbias, nlin, nscat,
     $    nn
      CHARACTER
     $    date*10,filter*6,object*80,ra*20,dec*20,cra*13,cdec*13
      COMMON /keyword/ ut,exptime,epoch,date,filter,object,ra,dec
      COMMON /findblk/ coords,xra,xdec,xfwhm,xmag,xerr,nc,cra,cdec
      COMMON /fwhmblk/ fwhmx,fwhmy,sharp,fwhm
      COMMON /imblk/ imbuf,nx,ny
      COMMON /apblk/ apzero,rdnoise,gain,iap,iskyinr,iskyoutr,nap,ich
      COMMON /skyblk/ thresh,sky,sigsky
      COMMON /initblk/ wscale,ixdir,iydir,izdir,lobad,hibad
c
      apmxsq = -1.
      do i=1,nap
        par(i) = iap(i)
        apmxsq = amax1(apmxsq,(par(i)+0.5)**2)
      enddo
c     print *,'par',apmxsq,par
c
c loop over number of objects
c
      if (nc.le.0) return
      nn = 0
c      print *,'nc ',nc,coords(3,1),coords(3,2),nx,ny,iskyinr,iskyoutr
      DO n=1,nc
        xc = coords(n,1)
        yc = coords(n,2)
        if ((xc.ge.1.0).and.(xc.le.float(nx))
     $    .and.(yc.ge.1.0).and.(yc.le.float(ny))) then
        rin = iskyinr
        rinsq = rin*rin
c set outer radius to either input iskyoutr or else the radius
c that includes MAXSKY pixels, whichever is smaller
        routsq = float(MAXSKY)/3.142 + rinsq
        dum = float(iskyoutr)**2
        if (dum.lt.routsq) routsq = dum
        rout = sqrt(dum)
c limits of submatrix to be searched for this object
        lx = max0(1,int(xc-rout)+1)
        mx = min0(nx,int(xc+rout))
        ly = max0(1,int(yc-rout)+1)
        my = min0(ny,int(yc+rout))
        edge =amin1(xc-0.5,(nx+0.5)-xc,yc-0.5,(ny+0.5)-yc)
        apmax = -1.1E38
c       if (n.eq.3) print *,'lx ',lx,mx,ly,my,edge,rinsq,routsq,rout
c if star aperture extends outside array, its magnitude is bad
        DO j=1,nap
          apmag(j) = 0.0
          if (edge.lt.par(j)) apmag(j) = -1.1D38
          area(j) = 0.0
        ENDDO
        nsky = 0
c
c loop through submatrix for this object
c
        DO j=ly,my
          dysq =(j-yc)**2
            DO k=lx,mx
              rsq = dysq + (k-xc)**2
              datum = imbuf(k,j)
c is pixel in sky annulus?
              if ((rsq.lt.rinsq).or.(rsq.gt.routsq).or.
     $           (nsky.gt.maxsky).or.(datum.lt.lobad).or.
     $           (datum.gt.hibad)) goto 100
                 nsky = nsky+1
                 apsky(nsky) = datum
                 goto 200
100           continue
              if (rsq.gt.apmxsq) goto 200
c pixel is in at least one of the apertures.
              r = sqrt(rsq)-0.5
              DO m=1,nap
c if pixel completely outside aperture, skip
                if (r.gt.par(m)) goto 150
c determine fraction of pixel inside aperture
                fractn=amax1(0.0,amin1(1.0,par(m)-r))
                if ((datum.ge.lobad).and.(datum.le.hibad).and.
     $            (apmag(m).gt.-1.0D38)) then
                  apmag(m) = apmag(m)+dble(fractn*datum)
                  area(m) = area(m) + dble(fractn)
                  apmax = max(apmax,datum)
                else
                  apmag(m) = -1.1D38
                endif
150           continue
              ENDDO
200           continue
         ENDDO
      ENDDO
c determine sky for this object
      if (nsky.lt.MINSKY) then
c      if (n.eq.3) print *,' Too few sky pixels',nsky,MINSKY
        do i=1,nap
          apmag(i) = 99.999D0
          aperr(i) = 9.999
        enddo
        goto 500
      endif
c     call calcsky(apsky,nsky,skymod,skysig)
      call sort1(apsky,nsky)
      call mmm(apsky,nsky,hibad,rdnoise,skymn,skymed,skymod,
     $    skysig,skew)
c     if (n.eq.3) print *,'nsky ',nsky,MINSKY,xc,yc,apmag(1),skymod
c kludge to handle error
      if (skysig.lt.0.0) skysig = 0.0
      skyvar = skysig**2
      sigsq = skyvar/float(nsky)
c subtract sky from each aperture
        DO i=1,nap
          apmag(i) = apmag(i) - dble(skymod)*area(i)
c if star+sky fainter than sky, or earlier badap occurred...
          if (apmag(i).le.0.0D0) goto 300
          err1 = area(i)*skyvar
          err2 = apmag(i)/gain
          err3 = sigsq*(area(i)**2)
          aperr(i) = amin1(9.999,1.0857*sqrt(err1+err2+err3)/
     $        sngl(apmag(i)))
c derive error and instrumental magnitude
            apmag(i) = apzero - 2.5D0*dlog10(apmag(i))
            if (apmag(i).gt.99.999D0) goto 300
            goto 400
300       continue
          apmag(i) = 99.999D0
          aperr(i) = 9.999
400       continue
        ENDDO
500     continue
        if (apmax.lt.-1.0E38.or.apmax.gt.99999.99) apmax=99999.99
c        if (n.eq.3) print *,'apmag ',apmag(1),fwhmx(1),fwhmy(1)
        if (apmag(1).lt.99.0) then
          if (fwhmx(n).lt.90.0.and.fwhmy(n).lt.90.0) then
            nn =min((nn + 1),99999)
            write (2,901) xra(n),xdec(n),xc,yc,fwhmx(n),fwhmy(n),
     $       sharp(n),apmax,
     $       skymod,(apmag(i),i=1,nap),xmag(n),
     $       (aperr(i),i=1,nap),xerr(n)
901         format (2f12.7,2f8.2,3f6.2,2f9.2,20f7.3)
c            write (2,900) (apmag(i),i=1,nap),xmag(n)
c            write (2,900) (aperr(i),i=1,nap),xerr(n)
c900         format (10f7.3)
          endif
        endif
        endif
      ENDDO
      RETURN
      END

      SUBROUTINE calcsky(sky,nsky,val,sig)
c
c determine sky from input array
c
      REAL*8 sum
      REAL*4 sky(*),val,sig
      INTEGER nsky
c
      call sort1 (sky,nsky)
      val = sky(nsky/2+1)
      sum=0.
      do i=1,nsky
        sum = sum + (sky(i)-val)**2
      enddo
      sig = dsqrt(sum)/dfloat(nsky-1)
      RETURN
      END

      SUBROUTINE SORT1 (x,n)
c
c heapsort of array x 
c from numerical recipies
c
      REAL*4  x(*)
      INTEGER n
      if (n.eq.1) return
      k = n/2+1
      ir = n
10    continue
      IF (k.gt.1) THEN
        k = k-1
        xx = x(k)
      ELSE
        xx = x(ir)
        x(ir) = x(1)
        ir = ir-1
        IF (ir.eq.1) THEN
          x(1) = xx
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
          i = j
          j = j+j
        ELSE
          j = ir+1
        ENDIF
      GOTO 20
      ENDIF
      x(i) = xx
      GOTO 10
      END
