	program FILSTAT
c
c  This program calculates mean magnitudes/colors/coordinates
c version for newphot 29-Jan-98 aah
c mod 20-Apr-98 aah improved match
c mod 14-June-1998 aah errrad=2arcsec
c mod 29-June-1998 aah remove points with flags=11111
c mod 21-Dec-1998 aah sorted + binary search
c mod 09-Apr-1999 aah coord match uses cos(dec)
c mod 15-Oct-1999 aah use initial errs for single points
c mod 17-Aug-2000 aah user input matching radius
c mod 03-Feb-2006 aah v-i format
c mod 12-Jun-2006 aah fix header
c mod 17-Jun-2007 aah add m n
c mod 19-Feb-2009 aah intermediate storage; 'add' option
c mod 04-Jul-2010 aah filter format
c mod 25-Aug-2010 aah declination ranges
c note: subtract mode does not yet work
c
c space = 94*MAXSTAR Bytes
c
c ***********************************************************************
c
c note: 5,000,000 stars permitted for now
      PARAMETER (MAXSTAR=5000000)
      PARAMETER (MAXFILT=11)
      REAL*8 ra(MAXSTAR),dec(MAXSTAR),ra2(MAXSTAR),dec2(MAXSTAR),
     $  rax,decx,dist,err,raz,err2,cosd,picon,declo,dechi
      REAL*8 rsq,rsqmax
      REAL*4 smag(MAXFILT,MAXSTAR),
     $  xmag(MAXFILT),xerr(MAXFILT),serr(MAXFILT,MAXSTAR),
     $  smag2(MAXFILT,MAXSTAR),raerr,decerr,errrad
      REAL*4 tmag(7),terr(7),xmag1,xerr1,ccdx,ccdy
      INTEGER npts,i,j,flags(11),kset,iflag,jflag,indx(MAXSTAR),
     $  ix
c     INTEGER*2 xnm(MAXFILT,MAXSTAR),xnr(MAXSTAR),ferr(MAXSTAR)
      INTEGER xnm(MAXFILT,MAXSTAR),xnr(MAXSTAR),ferr(MAXSTAR)
      INTEGER xmr(MAXSTAR), optype, csv, jcsdss, flag1, flag2
      CHARACTER file1*60,file2*60,file3*60,file4*60,star*10
      CHARACTER*10 sname(MAXSTAR),names(100),oldname
      character*80 ch
      LOGICAL fflag, nflag(MAXSTAR)
c
c ***********************************************************************
c
      print *,'    Program FILSTAT version 3.2   05-Jun-2012'
      print *
      print *
c
c ************************************************************************
c
      picon = datan(1.D0)/45.0
      rsqmax = 2.*(2048.-100.)**2 ! note for 4096x4096 only
      raerr = 0.0
      decerr = 0.0
      print *,'Enter type of operation (0=rebuild,1=add,2=del,',
     $   '3=summary): '
      read (5,*) optype
      print *,'Lower declination value: '
      read (5,*) declo
      print *,'Upper declination value: '
      read (5,*) dechi
      if (optype.ne.0) then
        print *,'Enter intermediate file: '
        read (5,'(a)') file4
        open (unit=4,file=file4,status='old')
      endif
      if (optype.ne.3) then
        print *,'Enter file list: '
        read (5,'(a)') file1
        open (unit=1,file=file1,status='old')
      endif
      print *,'Enter list of star names, xxx=all: '
      read (5,'(a)') file2
      if (file2.eq.'xxx') then
        jflag = 1
      else
        jflag = 0
        open (unit=2,file=file2,status='old')
        nst = 1
10      continue
        read (2,'(a)',end=20) names(nst)
        nst = min((nst+1),100)
        goto 10
20      continue
        nst = nst - 1
        close(2)
      endif
      print *,'Enter output file name: '
      read (5,'(a)') file3
      open (unit=3,file=file3,status='new')
      if (optype.gt.2) then
        print *,'J/C (0), SDSS (1), APASS (2), 7color (3): '
        read (5,*) jcsdss
        print *,'column delimited (0) or csv(1): '
        read (5,*) csv
        write (3,908)
908     format ('#',2x,'Name',4x,'RA(J2000)',3x,'raerr',2x,
     $    'DEC(J2000)',1x,'decerr',1x,'nobs',2x,'mobs',7x,
     $    'filt  mag  err')
        print *,'Enter min number of nights: '
        read (5,*) minobs
      endif
      if (optype.ne.3) then
        print *,'Enter matching error radius (arcsec): '
        read (5,*) errrad
c
        err = (errrad/3600.)  ! user input position err
        err2 = err*err
      endif
      npts = 0
      fflag = .true.
c
c read in intermediate file if already exists
c
      j=1
      if (optype.gt.0) then
60      continue
          read (4,9015,end=70) sname(j),ra(j),ra2(j),dec(j),
     $     dec2(j),xnr(j),xmr(j),
     $     (smag(i,j),smag2(i,j),serr(i,j),xnm(i,j),i=1,MAXFILT)
9015     format (a10,4d24.16,2i8,11(3e16.8,i10))
         indx(j) = j
         j=j+1
         goto 60
70      continue
        npts = j-1
        print *,'npts= ',npts
        iflag = 1
        close(4)
      endif
      if (optype.eq.3) goto 700
c
c loop over number of files
c
100   continue
      read (1,'(a)',end=700) file2
      write (6,905) file2
905   format (' Processing file: ',a50)
      open (unit=2,file=file2,status='old')
c skip over header
      read (2,*)
      read (2,*)
      read (2,*)
c
c for each new file, set the night flag to true
c this enables update to xnr for the first obs of a field on a night
c
      if (npts.ne.0) then
        do i=1,npts
         nflag(i) = .true.
        enddo
      endif
c loop over stars
200   continue
      read (2,901,end=600) rax,decx,ccdx,ccdy,flag1,flag2,
     $   hjd,avexx,kset,igroup,
     $  star,ifil,xmag1,xerr1
901   format(2f12.7,2f10.3,1x,i1.1,1x,i1.1,1x,f12.6,
     $  1x,f5.3,i5,1x,i10,1x,a10,1x,i5,1x,2f8.4)
      if (decx.lt.declo.or.decx.gt.dechi) goto 200
      if (flag1.ne.0) goto 200
c kludge - if outside corrected image circle, ignore
      rsq = (2048.-ccdx)**2 + (2048.-ccdy)**2
c     if (rsq.gt.rsqmax) goto 200
c
      cosd = dcos(decx*picon)
      if (jflag.eq.0) then
        do i=1,nst
          if (star.eq.names(i)) goto 220
        enddo
        goto 200
      endif
220   continue
c special check for first point in first file
      if (npts.eq.0) then
        npts = 1
        indx(npts) = npts
        sname(npts) = star
        ra(npts) = rax
        ra2(npts) = rax*rax
        dec(npts) = decx
        dec2(npts) = decx*decx
        xnr(npts) = 1
        xmr(npts) = 1
        do i=1,MAXFILT
          smag(i,npts) = 0.
          smag2(i,npts) = 0.
          xnm(i,npts) = 0
          serr(i,npts) = 0.
        enddo
        nflag(npts) = .false.
        if (xmag1.lt.90.0.and.flag1.eq.0) then
          smag(ifil,npts) = xmag1
          smag2(ifil,npts) = xmag1*xmag1
          xnm(ifil,npts) = 1
          serr(ifil,npts) = xerr1
        else
          smag(ifil,npts) = 0.
          smag2(ifil,npts) = 0.
          xnm(ifil,npts) = 0
          serr(ifil,npts) = xerr1
        endif
c
c now more general check for subsequent objects
c
      else
        iflag = 0
c first, find appropriate region to check
        raz = rax - 2.*err/cosd
        call locate(raz,ra,xmr,indx,npts,ix)
        if (ix.eq.0) ix = 1
        raz = rax + 2.*err/cosd
        i = ix
c then loop to see if object already in array
230     continue
        j = indx(i)
c if gt raz, beyond err box and we have new object
        if (ra(j)/dble(xmr(j)).gt.raz) goto 250
          dist = ((ra(j)/dble(xmr(j))-rax)*cosd)**2
     $       + (dec(j)/dble(xmr(j))-decx)**2
          if (dist.lt.err2) then
c object already in array, so add it to sums
            iflag = 1
            sname(j) = star     ! updates names
            xmr(j) = xmr(j) + 1
            if (nflag(j)) then
              xnr(j) = xnr(j) + 1
              nflag(j) = .false.
            endif
            ra(j) = ra(j) + rax
            ra2(j) = ra2(j) + rax*rax
            dec(j) = dec(j) + decx
            dec2(j) = dec2(j) + decx*decx
c     NOTE: This introduces a bug, a measurement will be recorded in the RA/DEC
c     averages, but the magnitude may not be recorded.
            if (xmag1.lt.90.0.and.flag1.eq.0) then
              smag(ifil,j) = smag(ifil,j) + xmag1
              smag2(ifil,j) = smag2(ifil,j) + xmag1*xmag1
c     NOTE: xnm appears to be a counter for measurements in a given filter.
              xnm(ifil,j) = xnm(ifil,j) + 1
            endif
            goto 250
          endif
          i = i+1
          if (i.le.npts) goto 230
c we have updated ra. need to sort locally since ordering may
c have changed.
c
250     continue
        if (iflag.eq.0) then
c new object, so create sums for it
          npts = min((npts+1),MAXSTAR)
          sname(npts) = star
          ra(npts) = rax
          ra2(npts) = rax*rax
          dec(npts) = decx
          dec2(npts) = decx*decx
          xnr(npts) = 1
          xmr(npts) = 1
          nflag(npts) = .false.
          do i=1,MAXFILT
            smag(i,npts) = 0.
            smag2(i,npts) = 0.
            xnm(i,npts) = 0
            serr(i,npts) = 0.
          enddo 
          if (xmag1.lt.90.0.and.flag1.eq.0) then
              smag(ifil,npts) = xmag1
              smag2(ifil,npts) = xmag1*xmag1
              xnm(ifil,npts) = 1
              serr(ifil,npts) = xerr1
          endif
          ferr(npts) = 0
c now insert object into ra-ordering
          i = ix
400       continue
          if (rax.lt.ra(indx(i))/dble(xmr(indx(i)))) then
            do j = npts,i+1,-1
              indx(j) = indx(j-1)
            enddo
            indx(i) = npts
            goto 500
          endif
          i = i+1
          if (i.ge.npts) then
             indx(npts) = npts
             goto 500
          endif
          goto 400
500       continue
        endif
      endif
      goto 200
600   continue
      write (6,906) npts
906   format ('   npts now: ',i7)
      goto 100
700   continue
c
c when get here, have all stars entered.
c now form means and standard deviations (s.e., not m.e.)
c note: add write to intermediate file here
c field 10 ra 8 dec 8 rasq 8 decsq 8 mag 4*5 magsq 4*5
c ran 4 decn 4 magn 4*5 ferr 2; total=112 bytes/record
c or 23 longwords
c make last record have ra=-xxx so that we can stop
c the initial read easily.
c
      if (optype.eq.0.or.optype.eq.1.or.optype.eq.2) then
        do j=1,npts
          k = indx(j)
          write (3,9015) sname(k),ra(k),ra2(k),dec(k),dec2(k),
     $     xnr(k),xmr(k),(smag(i,k),smag2(i,k),
     $     serr(i,k),xnm(i,k),i=1,MAXFILT)
        enddo
        close(3)
        goto 850
      endif
c optype =3, summary
      do j=1,npts
        i = indx(j)
c       print *,' pts: ',j,i,xnr(i),xmr(i),ra(i),dec(i),smag(3,i),
c    $   sname(i)
        if (xnr(i).lt.minobs) goto 800
        rax = ra(i)/dble(xmr(i))
        decx = dec(i)/dble(xmr(i))
        if (xmr(i).lt.2) then
          raerr = 0.0
          decerr = 0.0
        else
          raerr = dsqrt(abs(dble(xmr(i))*ra2(i)-ra(i)*ra(i))/
     $   (dble(xmr(i))*(dble(xmr(i))-1.)))
          decerr = dsqrt(abs(dble(xmr(i))*dec2(i)-dec(i)*dec(i))/
     $   (dble(xmr(i))*(dble(xmr(i))-1.)))
         raerr = raerr*3600.*dcos(decx*picon)
         decerr = decerr*3600.
         if (raerr.gt.100.0) raerr = 99.999
         if (decerr.gt.100.0) decerr = 99.999
        endif
        do k=1,MAXFILT
          if (xnm(k,i).ne.0) then
            xmag(k) = smag(k,i)/xnm(k,i)
            if (xnm(k,i).lt.2) then
              xerr(k) = -serr(k,i)
            else
             xerr(k) = sqrt(abs(float(xnm(k,i))*smag2(k,i)-
     $        smag(k,i)*smag(k,i))/(float(xnm(k,i))*
     $        (float(xnm(k,i))-1.)))
c              if (k.eq.3.and.sname(i).eq.'0020110035') then
c              print *,rax,decx,xmag(3),xerr(3),xnm(3,i),smag(3,i),
c     $           smag2(3,i)
c              endif
            endif
          else
           xmag(k) = 99.999
           xerr(k) = 9.999
          endif
         enddo
         if (csv.eq.0) then
           if (jcsdss.eq.0) then
c Johnson/Cousins output. DOES NOT WORK WITH FILCON FILES
             if (xmag(3).gt.90.0) goto 800
             tmag(1) = xmag(3)
             if (xmag(2).lt.90.0.and.xmag(3).lt.90.0) then
               tmag(2) = xmag(2) - xmag(3)
               terr(2) = sqrt (xerr(2)**2 + xerr(3)**2)
             else
               tmag(2) = 99.999
               terr(2) = 9.999
             endif
             if (xmag(2).lt.90.0.and.xmag(3).lt.90.0) then
               tmag(3) = xmag(1) - xmag(2)
               terr(3) = sqrt (xerr(1)**2 + xerr(2)**2)
             else
               tmag(3) = 99.999
               terr(3) = 9.999
             endif
             if (xmag(2).lt.90.0.and.xmag(3).lt.90.0) then
               tmag(4) = xmag(3) - xmag(4)
               terr(4) = sqrt (xerr(3)**2 + xerr(4)**2)
             else
               tmag(4) = 99.999
               terr(4) = 9.999
             endif
             if (xmag(2).lt.90.0.and.xmag(3).lt.90.0) then
               tmag(5) = xmag(4) - xmag(5)
               terr(5) = sqrt (xerr(4)**2 + xerr(5)**2)
             else
               tmag(5) = 99.999
               terr(5) = 9.999
             endif
             if (xmag(2).lt.90.0.and.xmag(3).lt.90.0) then
               tmag(6) = xmag(3) - xmag(5)
               terr(6) = sqrt (xerr(3)**2 + xerr(5)**2)
             else
               tmag(6) = 99.999
               terr(6) = 9.999
             endif
             write (3,903) sname(i),rax,raerr,decx,decerr,xnr(i),
     $        (tmag(k),k=1,6),(terr(k),k=1,6),xmr(i)
903           format (a10,f11.6,f7.3,f11.6,f7.3,i5,12f7.3,i5)
           elseif (jcsdss.eq.1) then
c SDSS output
           elseif (jcsdss.eq.2) then
c APASS output
c            if (xmag(3).gt.90.0) goto 800
             tmag(1) = xmag(3)
             terr(1) = xerr(3)
             if (xmag(2).lt.90.0.and.xmag(3).lt.90.0) then
               tmag(2) = xmag(2) - xmag(3)
               terr(2) = sqrt (xerr(2)**2 + xerr(3)**2)
             else
               tmag(2) = 99.999
               terr(2) = 9.999
             endif
             tmag(3) = xmag(2)
             terr(3) = xerr(2)
             tmag(4) = xmag(8)
             terr(4) = xerr(8)
             tmag(5) = xmag(9)
             terr(5) = xerr(9)
             tmag(6) = xmag(10)
             terr(6) = xerr(10)
             write (3,904) sname(i),rax,raerr,decx,decerr,xnr(i),
     $         xmr(i),(tmag(k),k=1,6),(terr(k),k=1,6)
904          format (a10,f11.6,f7.3,f11.6,f7.3,2i5,12f7.3)
           elseif (jcsdss.eq.3) then
c 7color
             tmag(1) = xmag(2)
             terr(1) = xerr(2)
             tmag(2) = xmag(3)
             terr(2) = xerr(3)
             tmag(3) = xmag(7)
             terr(3) = xerr(7)
             tmag(4) = xmag(8)
             terr(4) = xerr(8)
             tmag(5) = xmag(9)
             terr(5) = xerr(9)
             tmag(6) = xmag(10)
             terr(6) = xerr(10)
             tmag(7) = xmag(11)
             terr(7) = xerr(11)
             write (3,907) sname(i),rax,raerr,decx,decerr,xnr(i),
     $         xmr(i),(tmag(k),k=1,7),(terr(k),k=1,7)
907          format (a10,f11.6,f7.3,f11.6,f7.3,2i5,14f7.3)
          endif
        elseif (csv.eq.1) then
          if (jcsdss.eq.2) then
             if (xmag(3).gt.90.0) goto 800
             tmag(1) = xmag(3)
             terr(1) = xerr(3)
             if (xmag(2).lt.90.0.and.xmag(3).lt.90.0) then
               tmag(2) = xmag(2) - xmag(3)
               terr(2) = sqrt (xerr(2)**2 + xerr(3)**2)
             else
               tmag(2) = 99.999
               terr(2) = 9.999
             endif
             tmag(3) = 99.999
             terr(3) = 9.999
c transform from Jester et al. (2005)
c     V-R    =    1.09*(r-i) + 0.22        0.03
c     Rc-Ic  =    1.00*(r-i) + 0.21        0.01
             if (xmag(9).lt.90.0.and.xmag(10).lt.90.0) then
               tmag(4) = 1.09 * (xmag(9)-xmag(10)) + 0.22
               terr(4) = sqrt (xerr(9)**2 + xerr(10)**2)
             else
               tmag(4) = 99.999
               terr(4) = 9.999
             endif
             if (xmag(9).lt.90.0.and.xmag(10).lt.90.0) then
               tmag(5) = 1.00 * (xmag(9)-xmag(10)) + 0.21
               terr(5) = sqrt (xerr(9)**2 + xerr(10)**2)
             else
               tmag(5) = 99.999
               terr(5) = 9.999
             endif
             if (tmag(4).lt.90.00.and.tmag(5).lt.90.0) then
               tmag(6) = tmag(4) + tmag(5)
               terr(6) = sqrt (terr(4)**2 + terr(5)**2)
             else
               tmag(6) = 99.999
               terr(6) = 9.999
             endif
           endif
           write (3,9013) sname(i),rax,raerr,decx,decerr,xnr(i),
     $      (tmag(k),k=1,6),(terr(k),k=1,6),xmr(i)
9013       format (a10,',',f11.6,',',f7.3,',',f11.6,',',f7.3,',',
     $       i5,',',12(f7.3,','),i5,',')
        endif
800    continue
      enddo
850   continue
      close (3)
      stop
888   continue
      backspace(2)
      read (2,'(a)') ch
      print *,' file read err: ',ch
      stop
      end

      subroutine locate(rax,ra,xmr,indx,n,ira)
c locate rax inside of ra using bisection method
c ira is such that ra(ira) <= rax <= ra(ira+1)
c ira=0 or ira=n indicates out of range
      real*8 rax,ra(n),rn,r1
c     integer*2 xmr(n)
      integer xmr(n)
      integer n,ira,jl,ju,jm,indx(n)
c
      jl = 0
      ju = n+1
      rn = ra(indx(n))/dble(xmr(indx(n)))
      r1 = ra(indx(1))/dble(xmr(indx(1)))
100   if ((ju-jl).gt.1) then
        jm = (ju+jl)/2
        if ((rn.gt.r1).eqv.
     $    (rax.gt.ra(indx(jm))/dble(xmr(indx(jm))))) then
           jl = jm
        else
          ju = jm
        endif
        goto 100
      endif
      ira = jl
      return
      end
