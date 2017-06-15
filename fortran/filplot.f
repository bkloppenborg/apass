      PROGRAM filplot
c
c read APASS means file, plot scatter diagrams for each field
c written 28-July-1998 aah
c mod 16-Oct-1998 aah closest match + ubvri output
c mod 11-Jan-1999 aah nobs and b-v additions
c mod 15-Oct-1999 aah abs(xerr)
c mod 20-Oct-2002 aah g77
c mod 30-Mar-2009 aah add zplot features
c mod 11-Aug-2010 aah apass format
c mod 25-Aug-2010 aah hierr mods
c mod 22-Jul-2015 aah limit points
c
      INTEGER MP,MAXFILT
      PARAMETER (MP=80000)
      PARAMETER (MAXFILT=6)
c     REAL*8 ymn,ymx,xmn,xmx,xx,yy
      REAL*4 ymn,ymx,xmn,xmx,xx,yy
      REAL*4 ymin,ymax,xmin,xmax,aspect,zscale1
      REAL*4 zscale2,err,dist,xscale,yscale,verno
      REAL*4 x(MP),y(MP),z(mp),dif,difi
      REAL*4 z4(mp),z1(mp,maxfilt),z2(mp,maxfilt)
      REAL*4 x1(MP),x2(MP),x3(MP),y1(MP),y2(MP),y3(MP)
      REAL*4 raerr,decerr,xmag(MAXFILT),xerr(MAXFILT)
      REAL*8 ra(mp),dec(mp),picon,conv,rax,decx
      REAL*8 ra1,dec1,ra2,dec2,xi,eta,hierr,errlim
      INTEGER npts,k,nr,mr,iflag,eflag,zflag,nfield
      INTEGER i,j,nobs(MP),icolmag,icolerr,ilen,nerr
      INTEGER limflag,indx(MP),n1,n2,n3
      CHARACTER sname*10,iname*10,file1*80,flag*1,psfile*40,
     $  encapstr*52,psname*10,oldname*10,xname*10,file2*80,
     $  file3*80
      CHARACTER*10 fnames(4000),newname
      DATA aspect /1.35/  ! 23/17
      DATA psname /'postencap '/
c
      verno = 2.1
      verno = 1.9
      picon = datan(1.D0)/45.
      eflag = 0
10    continue
      print *,'FILPLOT Version 2.0 22-Jul-2015'
      print *,' Enter master file name (/uw3/vp/photometry.res): '
      read (5,'(a)') file1
      open (unit=1,file=file1,status='old')
      print *,' Enter 0(display), 1(printer) 2(file): '
      read (5,*) iflag
      print *,' Enter output file name: '
      read (5,'(a)') file2
      open (unit=22,file=file2,status='new')
      print *,'Normal (0) or use fieldtable(1): '
      read (5,*) zflag
      if (zflag.eq.1) then
        print *,'Enter fieldtable name: '
        read (5,'(a)') file3
        open (unit=33,file=file3,status='old')
        oldname = ' '
        i=0
30      continue
          read (33,9077,end=40) newname
9077      format (a10)
          if (newname.ne.oldname) then
            i = i+1
            fnames(i) = newname
            oldname = newname
          endif
          goto 30
40      continue
        close(33)
        nfield = i
        print *,'nfields= ',nfield
      endif
50    continue
      print *,' all stars (0) or 800 stars (1): '
      read (5,*) limflag
      print *,' Enter filter number for mag (1=V): '
      read (5,*) icolmag
      icolerr = icolmag
      print *,'Enter mean error, below which no plot: '
      read (5,*) errlim
c cols are V, B-V, B, g', r', i'
      print *,' Enter star name, xxx=all: '
      read (5,'(a)') sname
c read header plus first record to get starting name
      read (1,*)
      read (1,904,end=200) iname
      oldname = iname
      backspace (1)
c
c loop over records
c
60    continue
      npts = 0
      xmax = -1.e38
      xmin = 1.e38
      ymax = -1.e38
      ymin = 1.e38
c read one field
100   continue
        read (1,904,end=200) iname,rax,raerr,decx,decerr,nr,
     $    mr,(xmag(j),j=1,MAXFILT),(xerr(j),j=1,MAXFILT)
904     format (a10,f11.6,f7.3,f11.6,f7.3,2i5,12f7.3)
        xname = oldname
        if (iname.eq.oldname) then
          if (abs(xerr(icolerr)).gt.9.0) goto 100
          npts = min(mp,(npts + 1))
          indx(npts) = npts
          ra(npts) = rax
          dec(npts) = decx
          x(npts) = xmag(icolmag)
          y(npts) = abs(xerr(icolerr))
          if (y(npts).gt.1.0.and.y(npts).lt.9.0)
     $      y(npts) = 1.0
          nobs(npts) = nr
          do i=1,MAXFILT
            z1(npts,i) = xmag(i)
            z2(npts,i) = abs(xerr(i))
          enddo
          ymax = amax1(y(npts),ymax)
          ymin = amin1(y(npts),ymin)
          xmax = amax1(xmag(icolmag),xmax)
          xmin = amin1(xmag(icolmag),xmin)
        else
          oldname = iname
          backspace (1)
          goto 150
        endif
        goto 100
c
c finished reading one field. now check to see if it is the right one
c
150     continue
        j = 0
        if (zflag.eq.1) then
          do i=1,nfield
            if (xname.eq.fnames(i)) j = 1
          enddo
          if (j.ne.1) goto 210
        else
          if (sname.eq.'xxx'.or.xname.eq.sname) goto 210
        endif
        goto 60
c
c yup, so plot it
c
200   continue
      eflag = 1
      print *,'Field not found!'
210   continue
      print *,'field: ',xname,npts,xmin,xmax,ymin,ymax,iflag
c trap to prevent plotting for basically blank fields
      if (npts.le.3) goto 320
      nerr = 0
      hierr = 0.0
      do i=1,npts
        if (x(i).gt.11.0.and.x(i).lt.14.0) then
          hierr = hierr + y(i)
          nerr = nerr + 1
        endif
      enddo
      print *,'hierr ',npts,nerr,hierr
      if (nerr.gt.0) then
         hierr = hierr / dble(nerr)
         if (hierr.lt.errlim) goto 60
      endif
      xscale = abs(xmax - xmin)
      yscale = abs(ymax - ymin)
      IF (iflag.eq.0) THEN
        i = sm_device('X11 -g 800x800')
      ELSEIF (iflag.eq.1) THEN
        i = sm_device('postscript')
      ELSE
        print *,'Enter output filename: '
        read (5,'(a40)') psfile
        encapstr = psname//psfile
        i = sm_device(encapstr)
      ENDIF
      call sm_graphics
      xmn = xmin
      xmx = xmax
      ymn = ymin
      ymx = ymax
c     print *,xmin,xmax,ymin,ymax
c     call sm_limits (xmn,xmx,ymn,ymx)
      call sm_limits (xmin,xmax,ymin,ymax)
      call sm_ctype ('default')
      call sm_grelocate (11000,32000)
      call sm_label (xname)
      call sm_box (1,2,0,0)
      call sm_xlabel('magnitude')
      call sm_ylabel('err')
c     call sm_ptype (2.43d2,1)
      call sm_ptype (2.43e2,1)
c     call sm_expand(1.0d0)
      call sm_expand(1.0e0)
c
c kludge to limit number of plotted points
c
      call sort (npts,x,indx)
      if (limflag.eq.1) npts = min(npts,800)
c create three writes to speed up x11
      n1 = 0
      n2 = 0
      n3 = 0
      do j=1,npts
        if (nobs(indx(j)).eq.1) then
          n1 = n1 + 1
          x1(n1) = x(j)
          y1(n1) = y(indx(j))
        else if (nobs(indx(j)).eq.2) then
          n2 = n2 + 1
          x2(n2) = x(j)
          y2(n2) = y(indx(j))
        else
          n3 = n3 + 1
          x3(n3) = x(j)
          y3(n3) = y(indx(j))
        endif
      enddo
      if (n1.ne.0) then
        call sm_ctype ('red')
        do j=1,n1
          xx = x1(j)
          yy = y2(j)
          call sm_points(xx,yy,1)
        enddo
      endif
      if (n2.ne.0) then
        call sm_ctype ('green')
        do j=1,n2
          xx = x2(j)
          yy = y2(j)
          call sm_points(xx,yy,1)
        enddo
      endif
      if (n3.ne.0) then
        call sm_ctype ('default')
        do j=1,n3
          xx = x3(j)
          yy = y(j)
          call sm_points(xx,yy,1)
        enddo
      endif
      call sm_ctype ('default')
c     call sm_expand(1.0d0)
      call sm_expand(1.0e0)
      call sm_gflush
250   continue
      n = 0
      print *,'Enter s=star, n=next, x=exit: '
300   continue
      if (iflag.eq.0) then
        call sm_curs(xx,yy,k)
        flag=char(k)
      else
        flag = 'x'
        call sm_hardcopy
      endif
      if (flag.eq.'s') then
        n = n + 1
        dif = 1.e38
        do i=1,npts
          difi = ((xx-x(i))/xscale)**2 + ((yy-y(indx(i)))/yscale)**2
          if (difi.le.dif) then
            k = i
            dif = difi
          endif
        enddo
310     continue
        n1 = k
        k = indx(k)
        print *,' Closest: ',ra(k),dec(k),x(n1),y(k)
        write (22,9003) xname,n,ra(k),dec(k),(z1(k,j),j=1,maxfilt),
     $   (z2(k,j),j=1,maxfilt)
9003    format (a10,i5,1x,2f13.8,12f7.3)
        goto 300
      elseif (flag.eq.'n') then
c       if (eflag.eq.1) goto 810
        call sm_erase
320     continue
        if (sname.ne.'xxx') then
          rewind(1)
          goto 50
        else
          goto 60
        endif
      endif
810   continue
      close(1)
      stop
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

      SUBROUTINE SORT (n,x,indx)
c
c heapsort of array x with corresponding indx array indx
c from numerical recipies
c
      REAL*4 x(n)
      INTEGER indx(n)
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

