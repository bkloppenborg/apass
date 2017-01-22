      PROGRAM HEPP 
c
c match astrometry files into .epp format
c maximum of 5 filters, 20K stars per frame
c written 14-Jan-1997 aah
c this version does not use slalib...fix later
c algorithm is to use the first observation of an object as a
c 'seed'.  When looking at a new data point, see which seed
c is the closest to that object, and if the distance error
c is less than 2arcsec, add it to that seed.
c mod 14-June-1998 aah
c mod 09-Sept-1998 aah check proper number of filters
c mod 08-Mar-1999 aah improved neighbor match
c mod 27-May-1999 aah curve of growth calcs
c mod 04-Jan-2000 aah y2k
c mod 17-Aug-2000 aah user input error radius
c mod 17-Oct-2002 aah g77
c mod 17-Feb-2003 aah add hjd offsets
c mod 24-Feb-2003 aah add IR
c mod 26-Aug-2007 aah add init file read for lat/lon
c mod 06-May-2009 aah add x,y
c mod 30-Aug-2009 aah add filenum
c mod 30-Oct-2009 aah epoch photometry file format
c mod 14-Dec-2009 aah add "system" to output file
c mod 29-Jan-2009 aah add apass convenience
c mod 01-Apr-2010 aah sel's format
c mod 18-Nov_2011 aah add stacking
c mod 07-Mar-2013 aah add ZS,Y support
c mod 10-Aug-2016 aah new ast file format
c NOTE: THIS VERSION ONLY FOR APASS
c input on command line: S/N, base filename
c
      INTEGER MAXFILT,MAXPRO,MAXAP
      PARAMETER (MAXAP = 8)
      PARAMETER (MAXFILT = 15)
      PARAMETER (MAXPRO = 65000)
      REAL*4 latitude,longitude,errrad,apsize(MAXAP),wscale
      INTEGER i,j,firstfile,lastfile,iap,ierr,kindx,kset,
     $   iapcor,ilatlon,iapset(7),propid(MAXPRO),prid
      INTEGER photflag,prfield(MAXPRO),prfd,jj,apflag
      INTEGER istack(MAXPRO)
      CHARACTER
     *      base*60,fname*60,fnum*4,
     $      cstring*80, fil*8,chr*80,hemisphere
      character*25 prnam(MAXPRO),prname,object,obj
      LOGICAL ok
      common /blkx/ latitude,longitude,errrad,fnum,iap,iapcor,apsize
      common /blky/ npr,istack,prnam,prfield
      DATA kset /1/
c
c     print *,' HEPP version 3.0 10-Aug-2016'
c
c get hemisphere and directory from command line and do checks
c hemisphere is single character, either n or s
c filename is image directory name, format yymmdd
c
      narg = iargc()
      call getarg(1,hemisphere)
      call getarg(2,base)
c
c set lat/lon
c
      isystem = 40
      if (hemisphere.eq.'s') then
       latitude = -30.165278
       longitude =  70.815
      else
       latitude = 32.80333
       longitude = 105.4747
      endif
      iapset(1) = 2
      iapset(2) = 3
      iapset(3) = 4
      iapset(4) = 5
      iapset(5) = 7
      iapset(6) = 9
      iapset(7) = 11
      iapset(8) = 11
      wscale = 2.566
      do i=1,8
        apsize(i) = wscale*float(iapset(i))
      enddo
c
c note: for the purposes of APASS, we're not going to use a proposal file.
c all fields will be named acoording to the "object"
c
c get all file names and open appropriately
c
      nbase = nblank(base)
c
c kludge for loop.  Never had more than 2000 files in a night
c
      firstfile = 101
      lastfile = 2000
c assume 16arcsec aperture
      iap = 2
      iapcor = 2
      errrad = 3.
      photflag = 0
      apflag = 1
      fname = base(1:nbase)//'.epp'
      open (unit=3,file=fname,status='new')
c write header
      write (3,9001)
9001  format ('#INSTRUMENTAL MAGNITUDES, EPOCHPHOT 3.0'/
     $  '#',4x,'HJD',7x,'RA',9x,'DEC',9x,'FILT',1x,'AMASS',
     $   2x,'CCDX',5x,'CCDY',5x,'fwmx',3x,'fwmy',4x,'Peak',
     $   3x,'Sky',3x,'apr',2x,'apmag',3x,'aperr',
     $   1x,'pflg',
     $   2x,'proid',1x,'profd',2x,'syst',1x,'night',1x,'file',
     $   1x,'kset',1x,'kser',3x,'star')
c
c main loop over number of frames
c
c ierr = 0 and object = blank so first file will be processed correctly
c
      ierr = 0
      object = ' '
      fil = ' '
      j = firstfile - 1
100   CONTINUE
        j = j + 1
        if (j.gt.lastfile) goto 200
c
c try "b" file first
c
        write (fnum,9010) j
9010    format (i4.4)
        fname = base(1:nbase)//'b.'//fnum
c       write (6,9011) fname
9011    format (' Processing file ',a60)
c       get all filters for this set
c       one set:  defined by having the same objname, and no filter dups
        inquire (file=fname,exist=ok)
        IF (ok) THEN
          open (unit=1,file=fname,status='old')
          call readstars(ierr,object,obj,fil)
          close(1)
          if (ierr.eq.3) goto 200
        ELSE
c
c then try "a" file if "b" file does not exist
c
          write (fnum,9012) j
9012      format (i4.4)
          fname = base(1:nbase)//'a.'//fnum
          inquire (file=fname,exist=ok)
          IF (ok) THEN
            open (unit=1,file=fname,status='old')
            call readstars(ierr,object,obj,fil)
            close(1)
            if (ierr.eq.3) goto 200
          ENDIF
        ENDIF
        goto 100
200     CONTINUE
c
c ierr = 0 for first entry into readstars
c ierr = 1 for first file to process
c ierr = 2 for current matchset, not finished
c ierr = 3 for finished matchset
c ierr = 4 for file error, abort
c
c continue until all frames of this set have been gathered
c then match objects and write to file
c
        call matchset(kindx)
        prid = 99999
        prfd = 88888
        call writeset(kindx,kset,iapset,prid,prfd,photflag,
     $    isystem,object)
        kset = kset + 1
        ierr = 1
        if (j.gt.lastfile) goto 300
c
c note: have read one too many files so need to back up
c
        j = j - 1
        object = obj
        goto 100
300   CONTINUE
      close(3)
      STOP
      END

      SUBROUTINE READSTARS (ierr,object,obj,fil)
c
c read astrometry file on channel 1
c mod 27-May-1999 aah to calculate curve-of-growth
c mod 22-Nov-2011 aah stacking options
c ierr = 0 for first entry into readstars
c ierr = 1 for first file to process
c ierr = 2 for current matchset, not finished
c ierr = 3 for finished matchset
c ierr = 4 for file error, abort
c
      INTEGER MAXSTAR,MAXFILT,MAXAP,MAXPRO
      PARAMETER (MAXSTAR = 600000)
      PARAMETER (MAXFILT = 15)
      PARAMETER (MAXAP = 8)
      PARAMETER (MAXPRO = 65000)
      REAL*8 r(MAXSTAR,MAXFILT),d(MAXSTAR,MAXFILT),
     $  jd(MAXSTAR,MAXFILT),jdx,rx,dx
      REAL*4 fx(MAXSTAR,MAXFILT),fy(MAXSTAR,MAXFILT),
     $  sh(MAXSTAR,MAXFILT),apm(MAXSTAR,MAXFILT),
     $  apmag(MAXSTAR,MAXFILT,MAXAP),
     $  aperr(MAXSTAR,MAXFILT,MAXAP),
     $  rerr(MAXFILT),xmag(MAXAP),xerr(MAXAP),
     $  ut(MAXFILT),utd,raerr,ymag(MAXSTAR,MAXAP),
     $  yerr(MAXSTAR,MAXAP),diff(MAXSTAR),apcor(MAXAP),
     $  apsize(MAXAP),xx,yy,fxx,fyy,x1(MAXSTAR,MAXFILT),
     $  y1(MAXSTAR,MAXFILT),skydn(MAXSTAR,MAXFILT)
      REAL*4 latitude,longitude,xx1,xy1,xfx,xfy,
     $  apmx,skydnx,errrad,errmax,xxerr,yyerr
      INTEGER indx(MAXSTAR,MAXFILT),
     $  npts(MAXSTAR,MAXFILT),ii,kk,
     $  i,j,k,ierr,ifilts(26),jfilt,jx,iapcor,
     $  nstar(MAXFILT),iap,firstfile,indx1(MAXSTAR)
      INTEGER iyear,imonth,iday,jd0,iplate
      INTEGER npr,istack(MAXPRO),prfield(MAXPRO)
      CHARACTER*8 filts(26),fil
      CHARACTER*4 filen(MAXSTAR,MAXFILT),fnum
      CHARACTER object*25,obj*25,ch*80,filter*8
      character*25 prnam(MAXPRO),prname
      COMMON /blka/ r,d,fx,fy,sh,apm,apmag,aperr,skydn,
     $  indx,npts
      COMMON /blkb/ jd,ut,rerr,nstar,x1,y1,jd0,filen
      common /blkx/ latitude,longitude,errrad,fnum,iap,iapcor,apsize
      common /blky/ npr,istack,prnam,prfield
      DATA
     *       filts /'U       ','B       ','V       ','R       ',
     *              'F675W   ','I       ','HALC    ','F814W   ',
     *              'OPEN    ','Z       ','J       ','H       ',
     $              'K       ','KP      ','SU      ','SG      ',
     $              'SR      ','SI      ','SZ      ','HA      ',
     $              'TB      ','TG      ','TR      ','OIII    ',
     $              'ZS      ','Y       '/
      DATA ifilts /1,2,3,4,4,5,4,5,6,3,1,2,3,3,7,8,9,10,11,
     $   12,24,25,26,16,13,14/
      DATA nstar /MAXFILT*0/
c
      picon = datan(1.D0)/45.0
      errmax = (errrad/3600.)**2  ! user input radius
      iplate = 0
c for first entry, initialize arrays
      if (ierr.eq.0) then
        do i=1,MAXSTAR
          do j=1,MAXFILT
            npts(i,j) = 0
            r(i,j) = 0.0
            d(i,j) = 0.0
            jd(i,j) = 0.0
            fx(i,j) = 0.0
            fy(i,j) = 0.0
            do k=1,MAXAP
              apmag(i,j,k) = 0.0
              aperr(i,j,k) = 0.0
            enddo
            apm(i,j) = 0.0
            skydn(i,j) = 0.0
            x1(i,j) = 0.0
            y1(i,j) = 0.0
            filen(i,j) = "   "
c           sh(i,j) = 0.0
          enddo
        enddo
      endif
c
c read header (all comments starting with #)
c
10    continue
      read (1,'(a)') ch
      if (ch(1:1).ne.'#') then
        backspace(1)
        goto 20
      endif
      if (ch(2:7).eq.'Filter') then
        read(ch,9001) filter
9001    format (9x,a8)
      elseif (ch(2:3).eq.'UT') then
        read (ch,9002) uth
9002    format (5x,f7.4)
      elseif (ch(2:8).eq.'Exptime') then
        read (ch,9003) exptime
9003    format (10x,f7.2)
      elseif (ch(2:5).eq.'Date') then
        if (ch(12:12).eq.'-') then
          read (ch,9009) iyear,imonth,iday
9009      format (9x,i2,1x,i2,1x,i2)
        else
          read (ch,9004) iday,imonth,iyear
9004      format (7x,i2,1x,i2,1x,i2)
        endif
c y2k mod
        if (iyear.lt.50) then
          iyr = iyear + 2000
        else
          iyr = iyear + 1900
        endif
c get the mean RA error for the frame
      elseif (ch(2:4).eq.'No.') then
        read (ch,9005) raerr
9005    format(37x,f8.3)
      elseif (ch(2:7).eq.'Object') then
        obj = " "
        read (ch,9022) obj
9022    format (9x,a)
       elseif (ch(2:9).eq.'PlateFit') then
        iplate = 1
      endif
      goto 10
20    continue
c
c header has been read...now check to see if this is a new set.
c note: new object is always new set.
c if new set, don't read the star data since we will be rereading
c entire file anyway
c
      if (ierr.eq.1) fil = filter
      if (obj.ne.object.and.ierr.eq.2) then
        ierr = 3
        fil = filter
        return
      endif
c
c now get stacking option
c istk = 0 no stacking
c istk = 1 stacking all frames
c istk = 2 multifilter time series, stack until new filter seen
c
        istk = 0
c
c get filter index
c note: 26 is hardcoded, kludge
c
      jfilt = 1
      do i=1,26
        if (filter.eq.filts(i)) then
          jfilt = ifilts(i)
          goto 30
        endif
      enddo
30    continue
c if no stacking and repeat of filter, then end of this set
      if (nstar(jfilt).ne.0.and.istk.eq.0) then
        ierr = 3
        fil = filter
        return
      endif
c if stacking a time series and repeat of filter, end of set
      if (fil.ne.filter.and.nstar(jfilt).ne.0.and.istk.eq.2) then
        ierr = 3
        fil = filter
        return
      endif
      object = obj
c utx is the midpoint exposure time for this frame
c jdx is the midpoint exposure in JD
      utx = uth + exptime/(2.*3600.)
      exptime = 2.5*alog10(exptime)
      rerr(jfilt) = raerr
      ut(jfilt) = utx
      call jday (iday,imonth,iyr,utx,jdx)
cXXX this sets jd0 to first JD in file,
cXXX even if UTD rolls over later
      if (ierr.eq.0) then
        jd0 = int(jdx)
      endif
c
c read all stars for this frame
c note: hardcoded to 8 "apertures" (7 real + Sextractor)
c
      i = 0
100   continue
        i = min((i+1),MAXSTAR)
        if (iplate.eq.0) then
c       read (1,9024,end=200) r(i,jfilt),d(i,jfilt),
c    $    x1(i),y1(i),
c    $    fx(i,jfilt),fy(i,jfilt),sh(i,jfilt),apm(i,jfilt),
c    $    skydn(i,jfilt)
        read (1,9024,end=200) rx,dx,xx1,xy1,xfx,xfy,shx,
     $    apmx,skydnx,(xmag(j),j=1,8),(xerr(j),j=1,8)
c     $    apmx,skydnx,(xmag(j),j=1,7),(xerr(j),j=1,7)
c9024    format (2f12.8,2f8.2,3f6.2,f9.2,f9.2/7f7.3/7f7.3)
9024    format (2f12.8,2f8.2,3f6.2,f9.2,f9.2,20f7.3)
c       read (1,9025) (ymag(i,j),j=1,7)
9025    format (9f7.3)
c       read (1,9025) (yerr(i,j),j=1,7)
        else
c         read (1,9026,end=200) r(i,jfilt),d(i,jfilt),
c    $    x1(i),y1(i),
c    $    fx(i,jfilt),fy(i,jfilt),sh(i,jfilt),apm(i,jfilt),
c    $    skydn(i,jfilt),
c    $    (ymag(i,j),j=1,7),(yerr(i,j),j=1,7)
        read (1,9026,end=200) rx,dx,xx1,xy1,xfx,xfy,shx,
     $    apmx,skydnx,(xmag(j),j=1,8),(xerr(j),j=1,8)
c     $    apmx,skydnx,(xmag(j),j=1,7),(xerr(j),j=1,7)
c9026    format (2f12.8,2f8.2,3f6.2,f9.2,f9.2,
c     $    46x,18f7.3)
9026    format (2f12.8,2f8.2,3f6.2,f9.2,f9.2,
     $    46x,20f7.3)
        endif
c
c crude matching to older objects; needs sorting
c if xmag = 99, then data bad and skip object
c
        if (xmag(iap).gt.60.00) goto 100
c
c if already stars in list, and stacking enabled, look to
c see if we have a match with an existing object, and if so,
c sum the data
c
        ii = 0
        if (nstar(jfilt).gt.0.and.istk.ne.0) then
          kk = 0
          xxerr = 1.E30
          cosdjk = cos(dx*picon)**2
          do k=1,nstar(jfilt)
            xnpts = dble(npts(k,jfilt))
            yyerr = ((rx-r(k,jfilt)/xnpts)**2)*cosdjk +
     $             (dx-d(k,jfilt)/xnpts)**2
            if (yyerr.lt.xxerr) then
              xxerr = yyerr
              kk = k
            endif
          enddo
          if (xxerr.lt.errmax) then
            ii = kk
          endif
        endif
        xmag(iap) = xmag(iap) + exptime
        if (ii.eq.0) then
          nstar(jfilt) = nstar(jfilt) + 1
          ii = nstar(jfilt)
          npts(ii,jfilt) = 1
          r(ii,jfilt) = rx
          d(ii,jfilt) = dx
          jd(ii,jfilt) = jdx
          x1(ii,jfilt) = xx1
          y1(ii,jfilt) = xy1
          fx(ii,jfilt) = xfx
          fy(ii,jfilt) = xfy
          sh(ii,jfilt) = shx
          apm(ii,jfilt) = apmx
          do k=1,MAXAP
            apmag(ii,jfilt,k) = xmag(k)
            aperr(ii,jfilt,k) = xerr(k)
          enddo
          skydn(ii,jfilt) = skydnx
          filen(ii,jfilt) = fnum
          do j=1,8
            ymag(ii,j) = xmag(j)
            yerr(ii,j) = xerr(j)
          enddo
        else
          npts(ii,jfilt) = npts(ii,jfilt) + 1
          r(ii,jfilt) = r(ii,jfilt) + rx
          d(ii,jfilt) = d(ii,jfilt) + dx
          jd(ii,jfilt) = jd(ii,jfilt) + jdx
          x1(ii,jfilt) = x1(ii,jfilt) + xx1
          y1(ii,jfilt) = y1(ii,jfilt) + xy1
          fx(ii,jfilt) = fx(ii,jfilt) + xfx
          fy(ii,jfilt) = fy(ii,jfilt) + xfy
          sh(ii,jfilt) = sh(ii,jfilt) + shx
          apm(ii,jfilt) = apm(ii,jfilt) + apmx
          do k=1,MAXAP
            apmag(ii,jfilt,k) = apmag(ii,jfilt,k) + xmag(iap)
            aperr(ii,jfilt,k) = aperr(ii,jfilt,k) + xerr(iap)
          enddo
          skydn(ii,jfilt) = skydn(ii,jfilt) + skydnx
          filen(ii,jfilt) = fnum
          do j=1,8
            ymag(ii,j) = xmag(j)
            yerr(ii,j) = xerr(j)
          enddo
        endif
        goto 100
200   continue
c
c now calculate median aperture differences
c
      do j=1,8
        np = 1
        do i=1,nstar(jfilt)
          kflag = 0
          xx = 1.E30
          fxx = x1(i,j)
          fyy = y1(i,j)
          do k=1,nstar(jfilt)
            yy = (x1(k,j)-fxx)**2 +
     $               (y1(k,j)-fyy)**2
            if (yy.lt.xx.and.k.ne.i) then
              xx = yy
            endif
          enddo
          xx = sqrt(xx)
          if (xx.lt.apsize(iap)) kflag = 1
          if (kflag.eq.0.and.ymag(i,j).lt.90.0.and.
     $      ymag(i,iap).lt.90.0
     $      .and.yerr(i,j).lt.0.02.and.yerr(i,iap).lt.0.02
     $      .and.yerr(i,j).gt.0.001.and.yerr(i,iap)
     $      .gt.0.001) then
            diff(np) = ymag(i,j) - ymag(i,iap)
            indx1(np) = np
            np = np+1
          endif
        enddo
        np = np-1
c        if (j.eq.2) then
c          open (unit=33,file='gmtest.txt',status='new')
c          write (33,9098) (diff(kk),kk=1,np)
c9098      format (f7.3)
c          close(33)
c        endif
        if (np.gt.2) then
          call sort(np,diff,indx1)
          apcor(j) = diff(np/2)
        else
          apcor(j) = 0.0
        endif
      enddo
c
c now correct magnitudes.  If no neighbors, reasonable s/n,
c then mag = mag(iap).  If companion within iap*2 pixels,
c then mag = mag(1) - apcor.  Likewise, if err > 0.1,
c then mag = mag(1) - apcor.
c
      fil = filter
      ierr = 2
      return
      do i=1,nstar(jfilt)
        xmag(iap) = ymag(i,iap)
        xerr(iap) = yerr(i,iap)
c check for s/n
        if (ymag(i,iap).lt.17.0) then
           xmag(iap) = ymag(i,iapcor) - apcor(iapcor)
           xerr(iap) = yerr(i,iapcor)
        else
c check for neighbors
          xx = 1.E30
          fxx = x1(i,jfilt)
          fyy = y1(i,jfilt)
          do k=1,nstar(jfilt)
            yy = (x1(k,jfilt)-fxx)**2 +
     $               (y1(k,jfilt)-fyy)**2
            if (yy.lt.xx.and.k.ne.i) then
              xx = yy
            endif
          enddo
          xx = sqrt(xx)
          if (xx.lt.apsize(iap)) then
             xmag(iap) = ymag(i,iapcor) - apcor(iapcor)
             xerr(iap) = yerr(i,iapcor)
          endif
        endif
        do k=1,MAXAP
          apmag(i,jfil,kt) = xmag(k) + exptime
          aperr(i,jfilt,k) = xerr(k)
          if (xmag(k).gt.90.0) apmag(i,jfilt,k) = 99.999
        enddo
      enddo
      RETURN
      END

      SUBROUTINE MATCHSET (kindx)
c
c match all stars for this set
c no good way, so here is the first crack:
c  (1) if V present, use its starlist as master.
c  (2) if V not present, use reddest filter as master.
c  (3) compare all other filters with master, find closest star
c       to master star, write its index into indx.
c
      INTEGER MAXSTAR,MAXFILT,MAXAP
      PARAMETER (MAXSTAR = 600000)
      PARAMETER (MAXFILT = 15)
      PARAMETER (MAXAP = 8)
      REAL*8 r(MAXSTAR,MAXFILT),d(MAXSTAR,MAXFILT),
     $  jd(MAXSTAR,MAXFILT),picon,cosdjk
      REAL*4 fx(MAXSTAR,MAXFILT),fy(MAXSTAR,MAXFILT),
     $  sh(MAXSTAR,MAXFILT),apm(MAXSTAR,MAXFILT),
     $  apmag(MAXSTAR,MAXFILT,MAXAP),
     $  aperr(MAXSTAR,MAXFILT,MAXAP),
     $  rerr(MAXFILT),ut(MAXFILT),errmax,skydn(MAXSTAR,MAXFILT),
     $  yerr,err,x1(MAXSTAR,MAXFILT),y1(MAXSTAR,MAXFILT),
     $  apsize(MAXAP)
      REAL*4 latitude,longitude,errrad,xnpts,ynpts
      real*8 xx1,xx2,xx3
      INTEGER indx(MAXSTAR,MAXFILT),jd0,
     $  i,j,k,ierr,ifilts(8),jfilt,npts(MAXSTAR,MAXFILT),
     $  nstar(MAXFILT),indxv
      INTEGER iap,iapcor
      CHARACTER*4 filen(MAXSTAR,MAXFILT),fnum
      COMMON /blka/ r,d,fx,fy,sh,apm,apmag,aperr,skydn,
     $  indx,npts
      COMMON /blkb/ jd,ut,rerr,nstar,x1,y1,jd0,filen
      common /blkx/ latitude,longitude,errrad,fnum,iap,iapcor,apsize
c
      indxv = 3
      errmax = (errrad/3600.)**2  ! user input radius
      picon = datan(1.D0)/45.0
c match with respect to V frame, or reddest frame if
c V not present
      if (nstar(indxv).ne.0) then
        kindx = indxv
      else
        do j=MAXFILT,1,-1
          if (nstar(j).ne.0) then
            kindx = j
            goto 10
          endif
        enddo
10      continue
      endif
      do i=1,nstar(kindx)
        indx(i,kindx) = i
      enddo
c
c main loop
c slow, kludge, must check every kindx star for new nearest
c
c loop over number of filters
c
      DO i=1,MAXFILT
        if (nstar(i).ne.0.and.i.ne.kindx) then
          do j=1,nstar(kindx)
            indx(j,i) = 0
            xerr = 1.E30
c loop over stars within each filter to find match w/V star
c note: assumption is small frame, so that can just subtract
c RA's without doing true spherical trig
c note: for stacking, ynpts and xnpts refer to number of stacked points
c for nonstacking, they are =1
            ynpts = dble(npts(j,kindx))
            cosdjk = cos(d(j,kindx)*picon/ynpts)**2
c
c find closest star for other filter to the current V star
c
            do k=1,nstar(i)
              xnpts = dble(npts(k,i))
              yerr = ((r(j,kindx)/ynpts-r(k,i)/xnpts)**2)*cosdjk +
     $               (d(j,kindx)/ynpts-d(k,i)/xnpts)**2
              if (yerr.lt.xerr) then
                xerr = yerr
                kk = k
              endif
            enddo
c if this is within errmax of the V star, assume it is same and update
            if (xerr.lt.errmax) then
              indx(j,i) = kk
            endif
          enddo
        endif
      ENDDO
      RETURN
      END

      SUBROUTINE WRITESET (kindx,kset,iapset,prid,prfd,photflag,
     $  isystem,object)
c
c write this set of data to disk file
c
      INTEGER MAXSTAR,MAXFILT,MAXAP
      PARAMETER (MAXSTAR = 600000)
      PARAMETER (MAXFILT = 15)
      PARAMETER (MAXAP = 8)
      REAL*8 r(MAXSTAR,MAXFILT),d(MAXSTAR,MAXFILT),
     $  jd(MAXSTAR,MAXFILT),jdelt(MAXFILT),dt
      REAL*8 rasum,rasq,decsum,decsq,xnr,fxsum,fysum,
     $  shsum,rdsum,jdsum,jdx,rax,decx,hjd
      real*8 xx1,xx2,xx3,ray
      REAL*4 fx(MAXSTAR,MAXFILT),fy(MAXSTAR,MAXFILT),
     $  sh(MAXSTAR,MAXFILT),apm(MAXSTAR,MAXFILT),
     $  apmag(MAXSTAR,MAXFILT,MAXAP),
     $  aperr(MAXSTAR,MAXFILT,MAXAP),
     $  rerr(MAXFILT),smag(MAXFILT),serr(MAXFILT),
     $  smax(MAXFILT),ut(MAXFILT),decy,amass,
     $  x1(MAXSTAR,MAXFILT),y1(MAXSTAR,MAXFILT),
     $  xs(MAXFILT),ys(MAXFILT),psfmag,psferr,
     $  skydn(MAXSTAR,MAXFILT),apsize(MAXAP)
      REAL*4 latitude,longitude,errrad,xnpts,utx
      INTEGER indx(MAXSTAR,MAXFILT),ifil,iapv,
     $  i,j,k,kk,ierr,ifilts(8),jfilt,jd0,
     $  nstar(MAXFILT),kindx,kset,kser,photflag,
     $  npts(MAXSTAR,MAXFILT),nmin(MAXFILT)
      INTEGER iap,iapcor,iapset(7),prid,prfd,isystem
      CHARACTER*8 filts(8)
      CHARACTER*4 filen(MAXSTAR,MAXFILT),fils(MAXFILT)
      CHARACTER object*25,obj*25,ch*80,filt*8,fnum*4
      DATA filts/'       V','       B','       U',
     $   '       R','       I','       Z',
     $   '       X','       Y'/
      COMMON /blka/ r,d,fx,fy,sh,apm,apmag,aperr,skydn,
     $  indx,npts
      COMMON /blkb/ jd,ut,rerr,nstar,x1,y1,jd0,filen
      common /blkx/ latitude,longitude,errrad,fnum,iap,iapcor,apsize
c
c calculate means
c
      do j=1,MAXFILT
         nmin(j) = 0
      enddo
      do j=1,MAXFILT
        if (nstar(j).gt.0) then
        do i=1,nstar(j)
          nmin(j) = max(nmin(j),npts(i,j))
          xnpts = float (npts(i,j))
          if (npts(i,j).lt.1) xnpts = 1.0
          r(i,j) = r(i,j)/xnpts
          d(i,j) = d(i,j)/xnpts
          jd(i,j) = jd(i,j)/xnpts
          fx(i,j) = fx(i,j)/xnpts
          fy(i,j) = fy(i,j)/xnpts
          do k=1,MAXAP
            apmag(i,j,k) = apmag(i,j,k)/xnpts
            aperr(i,j,k) = aperr(i,j,k)/xnpts
          enddo
          apm(i,j) = apm(i,j)/xnpts
          skydn(i,j) = skydn(i,j)/xnpts
          x1(i,j) = x1(i,j)/xnpts
          y1(i,j) = y1(i,j)/xnpts
c         sh(i,j) = sh(i,j)/xnpts
        enddo
        endif
      enddo
      do j=1,MAXFILT
        if (nstar(j).gt.0) nmin(j) = max ((nmin(j)-1),1)
      enddo
c
c loop over stars in this field
c
      psfmag = 0.0
      psferr = 0.0
      kser = kset  ! KLUDGE
      iapv = iapset(iap)
      utx = 0.0
c
      ifil = iapset(iap)
      DO i=1,nstar(kindx)
        do j=1,MAXFILT
c        if (npts(i,j).ge.nmin(j)) then
          k = indx(i,j)
          if (k.ne.0) then
            indx(i,j) = 0
            call airmass(latitude,longitude,r(k,j),d(k,j),utx,jd(k,j),
     $           amass)
            ray = r(k,j)/15.
            call helcor (jd(k,j),ray,d(k,j),dt)
            hjd = jd(k,j) + dt
            do kk=1,MAXAP
              if (apmag(k,j,kk).gt.99.999) apmag(k,j,kk) = 99.999
              if (aperr(k,j,kk).gt.9.999) aperr(k,j,kk) = 9.999
            enddo
            write(3,9005) hjd,r(k,j),d(k,j),
     $       j,amass,
     $       x1(k,j),y1(k,j),fx(k,j),fy(k,j),
     $       apm(k,j),skydn(k,j),iapv,
     $       (apmag(k,j,kk),aperr(k,j,kk),kk=1,MAXAP),
     $       photflag,object,
     $       jd0,filen(k,j),kset,kser,i
9005         format (f12.5,f12.7,f12.7,
     $         i5,f7.3,
     $         f9.3,f9.3,f7.3,f7.3,
     $         f8.1,f7.1,i5,
     $         16f8.4,
     $         i5,a25,
     $         i6,1x,a4,i5,i5,i7)
       endif
c      endif
      enddo
      enddo
      j = 0
      DO i=1,MAXFILT
        if (nstar(i).ne.0) j = j+1
      ENDDO

      do j=1,MAXFILT
        if (nstar(j).ne.0) then
          do i=1,nstar(j)
            npts(i,j) = 0
            r(i,j) = 0.0
            d(i,j) = 0.0
            jd(i,j) = 0.0
            fx(i,j) = 0.0
            fy(i,j) = 0.0
            do k=1,MAXAP
              apmag(i,j,k) = 0.0
              aperr(i,j,k) = 0.0
            enddo
            apm(i,j) = 0.0
            skydn(i,j) = 0.0
            x1(i,j) = 0.0
            y1(i,j) = 0.0
            filen(i,j) = "   "
c           sh(i,j) = 0.0
          enddo
          nstar(j) = 0
        endif
      enddo
      RETURN
      END

      SUBROUTINE AIRMASS (latitude,longitude,ra,dec,ut,jdx,xm)
c
c calculate airmass for given position
c ra,dec in degrees, ut in hrs
c
      real*8 hjd0,ra,dec,picon,jdx
      real*4 ut,st,secz,xm,alon,alat,latitude,longitude,ha,d
      real*4 errrad
      real*8 xi
      integer ist,iap,iapcor
      integer ix
      character fnum*4
c
      picon = datan(1.D0)/45.
      alon = longitude/15.
      ix = int(jdx)
      xi = ix + 0.5
      ut = (jdx - xi)*24.0
      hjd0 = jdx-ut/24.0
      call stime (hjd0,ut,alon,st,ist)
      ha = st-ra/15.
      if (ha.gt.12.0) ha=ha-24.0
      if (ha.lt.-12.0) ha=ha+24.0
      ha = ha*15.*picon
      d = dec*picon
      alat = latitude*picon
      alon = longitude*picon
      secz = 1./(sin(alat)*sin(d)+cos(alat)*cos(d)*cos(ha))
      xm = secz - 0.0018167*(secz-1.)-0.002875*(secz-1.)**2
     $  - 0.0008083*(secz-1.)**3
      RETURN
      END

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

      subroutine stime (date,ut,alon,st,ist)
c
c	*** calculate sidereal time ***
c
c	inputs -
c		date	julian date - 2400000. (at 0h ut)
c		ut	ut in decimal hours
c		alon	observer longitude in hours
c			  ( longitude in deg / 15)
c	outputs -
c		st	sidereal time in decimal hours
c		ist	s. t. in format hhmm  (integer)
c
c	*** programmed by arne a. henden 1977 ***
c
      REAL*8 date
      REAL*4 sthr,ut,alon,st
      INTEGER ist
c
c	sthr is # sidereal hrs since jan 1, 1900
      sthr=6.6460556+2400.0512617*(date-15020.)/36525.+
     $   ut*1.00273-alon
c	get s. t. from sthr - # whole s. t. days since jan 1, 1900
      st=sthr-int(sthr/24.0)*24.0
c	now combine to get ist
      ist=int(st)
      ist=ist*100+int((st-float(ist))*60.)
      return
      end

      SUBROUTINE SORT (n,x,indx)
c
c heapsort of array x with corresponding index array index
c from numerical recipies
c
      REAL*4 x(n),xx
      INTEGER indx(n),n,k,ir
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

        subroutine jday(id,im,iyr,uthr,date)
c this sub calculates julian date from ut date
c inputs are id,im,iy, uthr--day,month,year,and ut hour
c output is date, julian date-2,400,000
c Henden, 1977
c
	double precision date,jdo,iiy
        dimension mo(12)
        integer id,im,iyr,iy
        data mo /0,31,59,90,120,151,181,212,243,273,304,334/
c
c leap is number of leap days since 1900
c
        iy = iyr - 1900
        leap=iy/4
        if((4*leap-iy) .eq. 0 .and. im .lt. 3) leap=leap-1
        iiy=365.0*float(iy)
        jdo=15020.0+iiy+float(leap)+float(mo(im))+float(id)
        date=jdo+uthr/24.0-0.5
        return
        end
	SUBROUTINE HELCOR (JD,ALPHA,DELTA,DT)
C
C  THIS SUBROUTINE CALCULATES THE HELIOCENTRIC CORRECTION
C   JD IS THE JULIAN DATE (MINUS 2400000) AT THE UT OF THE OBSERVATION
C   ALPHA IS THE RA IN HOURS
C   DEC IS THE DECLINATION IN DEGREES
C   DT IS THE HELCOR IN DAYS
C
C   FORTRAN 77 COMPATIBLE
C
C   LAST REVISED 21-APR-89
C   
	REAL*8 JD,ALPHA,DELTA,DT
	PICON=6.28318531/360.0
	TT=(JD-15020.)/36525.
	P=(1.396041+0.000308*(TT+0.5))*(TT-0.499998)
	RL=279.696678+36000.76892*TT+0.000303*(TT**2)-P
	G=358.475833+35999.04975*TT-0.00015*(TT**2)
8	IF(G .LT. 360.) GO TO 11
	G=G-360.
	GO TO 8
11	IF(RL .LT. 360.) GO TO 10
	RL=RL-360.
	GO TO 11
10	RL=RL*PICON
	G=G*PICON
	GLM=G-RL
	GLP=G+RL
	TGLP=2.0*G+RL
	TGLM=2.0*G-RL
	X=.99986*COS(RL)-0.025127*COS(GLM)+0.008374*COS(GLP)+.000105*
     1  COS(TGLP)+.000063*TT*COS(GLM)+.000035*COS(TGLM)
	Y=.917308*SIN(RL)+.023053*SIN(GLM)+.007683*SIN(GLP)+.000097*
     1  SIN(TGLP)-.000057*TT*SIN(GLM)-.000032*SIN(TGLM)
	A=ALPHA*15.0*PICON
	DEL=DELTA*PICON
	DT=-0.0057755*((COS(DEL)*COS(A))*X+(.4337751*SIN(DEL)+
     1  COS(DEL)*SIN(A))*Y)
	RETURN
	END
