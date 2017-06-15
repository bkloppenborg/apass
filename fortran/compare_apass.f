      program compare_apass
c
c compare APASS DR7 calibration fields with new APASS epp files
c to determine the residual offsets for each filter
c hardcoded for 4 calibration fields
c note: only currently works for BVgri
c uses central 600x600 region to set zeropoint
c written 2014-11-21 aah
c mod 2017-06-10 aah new file formats
c
      integer MP,NFIELD,NFILT
      parameter (MP=500000)
      parameter (NFIELD=4)
      parameter (NFILT=12)
      real*8 xra,xdec,cra(MP),cdec(MP),errsq,errmax
      real*8 ra,dec,hjd,xval,dra
      real*4 cmag(MP,NFILT),cerr(MP,NFILT),
     $      xdif(MP),radif(MP),xx(MP),yy(MP),
     $      decdif(MP),xerr(MP),xmag(MP)
      real*4 amass,ccdx,ccdy,fwhmx,fwhmy,amax,sky,apmag,aperr,
     $    psfmag,psferr,xd(MP)
      real*4 xbeg,xend,ybeg,yend,rad,rad2,dm,xsum,xoff
      real*4 apm(7),ape(7),cm(5),ce(5)
      integer nstd,i,j,k,n,ifilt
      integer iap,photflag,propid,propfield,isystem,jd0,kset,
     $    kser,istar,nmatch,indx(MP),id(MP)
      character file1*80,file2*80,file3*80,file4*15
      character nfile*4,object*25
      character*5 ifield(NFIELD)
      data xbeg,xend,ybeg,yend/1700.,2299.,1700.,2299./
      data ifield/'10006','10040','10068','10106'/
c
c errmax is error radius for match.  assume decent astrometry here
c
      errmax = (1.3/3600.)**2
      dra = 1.5/3600.
      rad = 0.
      rad2 = 0.
c
c read the four calibration fields
c
      i=1
      do j=1,NFIELD
      file4 = 'field_'//ifield(j)//'.txt'
      print *,'Reading file ',file4
      open (unit=1,file=file4,status='old')
10    continue
        read (1,900,end=20) cra(i),cdec(i),
     $   cm,ce

900     format (10x,f11.6,7x,f11.6,17x,f7.3,
     $    7x,9f7.3)
c need to put calib mags into epp order
        cmag(i,2) = cm(2)
        cmag(i,3) = cm(1)
        cmag(i,8) = cm(3)
        cmag(i,9) = cm(4)
        cmag(i,10) = cm(5)
        cerr(i,2) = ce(2)
        cerr(i,3) = ce(1)
        cerr(i,8) = ce(3)
        cerr(i,9) = ce(4)
        cerr(i,10) = ce(5)
        if (cmag(i,2).lt.16.0) i=i+1
        goto 10
20    continue
      close(1)
      enddo
      nstd = i-1
      print *,'nstd = ',nstd
c
      print *,'Enter epp file list: '
      read (5,'(a)') file1
      print *,'Choose filter(1-10): '
      read (5,*) ifilt
      print *,'Choose aperture (1-7; usually 2): '
      read (5,*) iap
      print *,'Enter output file name (raster_SR.txt): '
      read (5,'(a)') file2
      open (unit=1,file=file1,status='old')
      open (unit=2,file=file2,status='new')
c loop over multiple epp files to combine results
100   continue
        read (1,'(a)',end=400) file3
        print *,'Processing file ',file3
        open (unit=3,file=file3,status='old')
c read header
          read (3,*)
          read (3,*)
c loop over stars in current epp file
          n=1
200       continue
          read (3,9003,end=300) hjd,ra,dec,ifil,amass,
     $     ccdx,ccdy,fwhmx,fwhmy,amax,sky,
     $     (apm(i),ape(i),i=1,7),psfmag,psferr,photflag,
     $     object,
     $     jd0,nfile,kset,kser,istar
9003      format (f12.5,f12.7,f12.7,i5,f7.3,
     $     f9.3,f9.3,f7.3,f7.3,
     $     f8.1,f7.1,16f8.4,
     $     i5,a25,i6,1x,
     $     a4,i5,i5,i7)
           apmag = apm(iap)
           aperr = ape(iap)
           if ((dec.gt.-1.7).and.(dec.lt.1.7).and.
     $        (amass.lt.1.7).and.
     $        (ifil.eq.ifilt).and.abs(aperr).lt.0.1) then
             xval = ra-dra
             call locate (cra,xval,nstd,k)
             if(k.ne.0.and.k.ne.nstd) then
             do i=k,nstd
               if ((cra(i)-ra).gt.dra) goto 200
               errsq = (cra(i)-ra)**2 + (cdec(i)-dec)**2
               if (errsq.lt.errmax.and.abs(cerr(i,ifilt)).lt.0.1) then
                 xdif(n) = cmag(i,ifilt)-apmag
                 xmag(n) = apmag
                 xerr(n) = aperr
                 radif(n) = (cra(i)-ra)*3600.
                 decdif(n) = (cdec(i)-dec)*3600.
                 xx(n) = ccdx
                 yy(n) = ccdy
                 indx(n) = i
                 n = n+1
                 goto 200
               endif
             enddo
             endif
           endif
           goto 200
300       continue
          nmatch = n-1
          print *,'nmatch ',nmatch
c
c have all matches for current ifield, epp file. now generate
c residuals and write to file
c
c first, find mean offset at center of field
c replace this with SORT to get median
c
          xsum = 0.
          nx = 1
          do i=1,nmatch
            if (xx(i).ge.xbeg.and.xx(i).le.xend.and.
     $        yy(i).ge.ybeg.and.yy(i).le.yend) then
              xd(nx) = xdif(i)
              id(nx) = nx
              nx = nx+1
            endif
          enddo
          nx = nx-1
          if (nx.ge.4) then
            call sort (nx,xd,id)
            k = nx/2
            xoff = xd(k)
            print *,'xoff,nx= ',xoff,nx
            do i=1,nx
               print *,i,xd(i)
            enddo
          else
            goto 350
          endif
c
c now write all residuals after correcting for offset
         do i=1,nmatch
           k = indx(i)
           dm = xdif(i)-xoff
           if (dm.lt.0.5) then
             write (2,903) cra(k),cdec(k),xx(i),yy(i),
     $        cmag(k,ifilt),cerr(k,ifilt),xmag(i),xerr(i),
     $        radif(i),decdif(i),dm
903          format (2f12.6,2f10.3,4f8.3,f8.3,f8.3,f8.3)
           endif
         enddo
350      continue
        close(3)
c process next fred file
        goto 100
400   continue
      close(1)
      close(2)
      stop
      end

      SUBROUTINE SORT (n,x,indx)
c
c heapsort of array x with corresponding index array index
c from numerical recipies
c
      DIMENSION x(n),indx(n)
      k = n/2+1
      ir = n
10    continue
      IF (k.gt.1) THEN
        k = k-1
        xx = x(k)
        ii = indx(k)
      ELSE
        xx = x(ir)
        ii = indx(ir)
        x(ir) = x(1)
        indx(ir) = indx(1)
        ir = ir-1
        IF (ir.eq.1) THEN
          x(1) = xx
          indx(1) = ii
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
          indx(i) = indx(j)
          i = j
          j = j+j
        ELSE
          j = ir+1
        ENDIF
      GOTO 20
      ENDIF
      x(i) = xx
      indx(i) = ii
      GOTO 10
      END

      subroutine locate(x,xval,n,indx)
c
c locate xval inside of ra using bisection method
c indx is such that x(indx) <= xval <= x(indx+1)
c indx=0 or indx=n indicates out of range
c limitations: x must be real*8 and ascending order
c
      real*8 xval,x(n),rn,r1
      integer n,jl,ju,jm,indx
c
      jl = 0
      ju = n+1
      rn = x(n)
      r1 = x(1)
100   if ((ju-jl).gt.1) then
        jm = (ju+jl)/2
        if ((rn.gt.r1).eqv.
     $    (xval.gt.x(jm))) then
           jl = jm
        else
          ju = jm
        endif
        goto 100
      endif
      indx = jl
      return
      end
