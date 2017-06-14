      program filtran
c
c   Written by A. Henden
c     from an original program by R. Kaitchuck
c   Fortran 77 compatible
c
c   External subroutines: Solve2, Lin, Prompt, Bell, ucase, rdnames, wrtnames,
c                         extrd, disply
c
c   Revised 04-May-98 aah
c  this version works with newphot .ins files
c  revised to include (V-I) 17-Feb-99 aah
c revised for g77 17-October-2002 aah
c major revision 22 August 2003 aah
c minor fixes 26 November 2006 aah
c mod 07-Oct-2009 aah new file format 
c mod 01-Nov-2009 aah to work with epp files
c mod 18-Jun-2010 aah to use magnitudes instead of mag+multicolor
c mod 25-Jan-2014 aah to use apass standard file or landolt standard file
c mod 18-Dec-2014 aah to add zs,Y
c mod 14-Jun-2017 aah for tera, 64bit SM
c
      integer MAXSTARS,MAXC,MAXFILT,MAXFC
      PARAMETER (MAXSTARS = 70000)
      PARAMETER (MAXFILT = 14)
      PARAMETER (MAXFC = 2)
      PARAMETER (MAXC = 1)
      real*8 sumdev(MAXFC),sumdev2(MAXFC),sqerr
      real*8 rax,decx,hjd,sra(MAXSTARS),sdec(MAXSTARS)
      real*4 smag(MAXSTARS,MAXFILT),
     $   sext(MAXFC),fext(MAXFC),serr(MAXSTARS,MAXFILT),
     $   coef(MAXFC),zero(MAXFC),x1(MAXSTARS),y1(MAXSTARS),
     $   umag(MAXSTARS,MAXFILT),dmag(MAXFC),
     $   c(MAXSTARS,MAXC),sm(MAXFILT),se(MAXFILT),
     $   xx(MAXSTARS,MAXFILT),delc(MAXFC),fmag(MAXSTARS),
     $   uerr(MAXSTARS,MAXFILT),x(MAXSTARS,MAXFC),
     $   xmag(MAXSTARS,MAXFC),ymag(MAXSTARS,MAXFC)
      real*4 tmag(5),terr(5)
      real*4 xk,yk,amag,amagerr,am,xbeg,ybeg,xlen,ylen
      real*4 slope,b,zstar
      integer icolor(MAXC,MAXFC),korder(MAXFC),icolor1,icolor2
      INTEGER filno,iset,isetnew,istar,istarnew,nstd
      INTEGER iiflag,jjflag,jj,magflag
      integer indx(MAXSTARS),ka,kb,nstar,nfirst,nindx
      integer sindx(MAXSTARS,2),nfilt1,nfilt2,nfilt3
      integer iflagapass
      character*30 xlabel,ylabel
      character*50 fileid,raw,input
      character*50 file0,file1,file2,file3
      character*6 filt(MAXFILT),fil1,fil2,ichr
      character*11 star(MAXSTARS),star1,sname,xstar(MAXSTARS)
      character*13 clabel(5),fil3
      character*3 iflg(6),char1,char2,char3
      character doit,skip,lab(90),zpt,warn,adummy,ifl
      character tname*15
      logical connect,mflag,erflg
      common /names/ raw, file1, file2
      common /plotit/ x1,y1,nfit,xbeg,ybeg,xlen,ylen,xlabel,
     $   ylabel,icolor1,connect,slope,b,icolor2,mflag
      common /stds/ sra,sdec,smag,serr,nstd,star
      data xdif /0.0/
c
      sqerr = (2./3600.)**2
      erflg=.false.
      nc=1
      jjflag = 0
c ***********************************************************************
      print *,'                      Program FILTRAN version 2.1'
      print *,'                          December 17, 2014'
      print *
c ************************************************************************
c
      print *,'This program computes the transformation coefficients',
     $          ' to the standard'
      print *,'photometric system, and optionally computes the first',
     $          ' order extinction'
      print *,'coefficients and/or the zero point terms alone.  ',
     $          'The input consists of the' 
      print *,'ouput file from HEPP and the star list file which ',
     $        'contains the magnitude and '
      print *,'colors of the standard stars.  The results are ',
     $  'written to a FILTRAN file. '
      print *
      print *
      print *,' use bright standards (0=no,1=yes): '
      read *,magflag
c **********************************************************************
c  open the input files and check the file ID strings
c **********************************************************************
      call rdnames
100   if(file1 .eq. ' ') then
        print *,'  Enter the instrumental magnitude file name --- '
        read (5,'(a)') file1
      else
        print *,'  Enter the instrumental magnitude file [',
     $          file1(1:nblank(file1)),']: '
        read(5,'(a)') input
        if(input .ne. ' ') file1=input
      endif
      iflag=1
c
      open(unit=1,file=file1,status='old',err=200)
      read (1,'(a)') fileid
      if(fileid(1:24) .ne. '#INSTRUMENTAL MAGNITUDES') then
        call bell
        print *,'  This is the wrong kind of file ... try again'
        close (1)
        file1 =' '
        go to 100
      endif
c
103   if(file2 .eq. ' ') then
        print *,'  Enter the standard star file --- '
        read (5,'(a)') file2
      else
        print *,'  Enter the standard star file [',
     $          file2(1:nblank(file2)),']: '
        read(5,'(a)') input
        if(input .ne. ' ') file2=input
      endif
      iflag=2
c
      open(unit=2,file=file2,status='old',err=200)
      read (2,'(a)') fileid
      if(fileid(1:15) .ne. '#STAR LIST FILE') then
        call bell
        print *,'  This is the wrong kind of file ... try again: '
        close (2)
        if(iflag .eq. 2) file2=' '
        if(iflag .eq. 1) file1=' '
        go to 103
      endif
      read (2,*)
c
      call wrtnames
c
      print *,'  Enter the FILTRAN file name --- '
      read (5,'(a)') file3
c
      i=1
      go to 290
200   call bell
      print *,'  File cannot be opened, try again --- '
      if(iflag .eq. 1) then
        file1=' '
        go to 100
      endif
      if(iflag .eq. 2) then
        file2=' '
        go to 103
      endif
c
c *************************************************************************
c  Get filter info and bring in the data
c *************************************************************************
c
c first, read in all standard stars, convert to magnitude format,
c  and store
c
290   continue
      iflagapass = 0
      read (2,'(a)') ichr
      if (ichr.eq.'000201') iflagapass = 1
      backspace (2)
300   continue
        if (iflagapass.eq.0) then
          read (2,900,err=350,end=400) star(i),sra(i),sdec(i),
     $      (sm(j),j=1,6),(se(j),j=1,6)
900       format(a10,2x,f12.8,1x,f12.8,28x,12f7.3)
c else is apass stds
        else
          read (2,9080,err=350,end=400) star(i),sra(i),sdec(i),
     $     smag(i,3),smag(i,2),smag(i,8),smag(i,9),smag(i,10),
     $     serr(i,3),serr(i,2),serr(i,8),serr(i,9),serr(i,10)
9080      format (a10,1x,f10.6,8x,f10.6,17x,f7.3,7x,4f7.3,f7.3,7x,4f7.3)
          smag(i,1) = 99.999
          serr(i,1) = 9.999
          do k=11,MAXFILT
            smag(i,k) = 99.999
            serr(i,k) = 9.999
          enddo
          do k=4,7
            smag(i,k) = 99.999
            serr(i,k) = 9.999
          enddo
          i = i+1
          goto 300
        endif
        if (magflag.eq.0) then
           if (sm(1).lt.10.0) goto 300
        endif
        if (sm(1).gt.60.0) then
           smag(i,3) = 99.999
           serr(i,3) = 9.999
        else
           smag(i,3) = sm(1)
           serr(i,3) = se(1)
        endif
        if (sm(2).gt.60.0.or.sm(1).gt.60.0) then
           smag(i,2) = 99.999
           serr(i,2) = 9.999
        else
           smag(i,2) = sm(2) + sm(1)
           serr(i,2) = sqrt (se(2)**2 + se(1)**2)
        endif
        if (sm(3).gt.60.0.or.smag(i,2).gt.60.0) then
           smag(i,1) = 99.999
           serr(i,1) = 9.999
        else
           smag(i,1) = sm(3) + smag(i,2)
           serr(i,1) = sqrt(serr(i,2)**2 + se(3)**2)
        endif
        if (sm(1).gt.60.0.or.sm(4).gt.60.0) then
           smag(i,4) = 99.999
           serr(i,4) = 9.999
        else
           smag(i,4) = sm(1) - sm(4)
           serr(i,4) = sqrt(se(1)**2 + se(4)**2)
        endif
        if(smag(i,4).gt.60.0.or.sm(5).gt.60.0) then
           smag(i,5) = 99.999
           serr(i,5) = 9.999
        else
           smag(i,5) = smag(i,4) - sm(5)
           serr(i,5) = sqrt (smag(i,4)**2 + se(5)**2)
        endif
        if (sm(6).lt.90.0.and.sm(1).lt.60.0) then
            smag(i,5) = sm(1) - sm(6)
            serr(i,5) = sqrt(se(1)**2 + se(6)**2)
        endif
        do j=6,MAXFILT
          smag(i,j) = 99.999
          serr(i,j) = 9.999
        enddo
        i=i+1
        goto 300
350   continue
        print *,"read error in standard file",i
        close(2)
        stop
400   continue
      close(2)
      nstd = i-1
c
c now read in SDSS standards
c
      open (unit=2,file="sdss_stds.txt",status='old')
      read (2,*)
      read (2,*)
405   continue
        read (2,935,end=410) tname,rax,decx,(tmag(j),terr(j),j=1,5)
935     format (a15,2x,f11.5,f11.4,5(2x,2f9.3,5x))
c       print *,tname,rax,decx,tmag(3),terr(5)
        if (magflag.eq.0) then
          if (tmag(3).lt.10.0) goto 405
        endif
        do j=1,5
          if (tmag(j).lt.-50.0) tmag(j) = 99.999
          if (terr(j).lt.-50.0) terr(j) = 9.999
        enddo
        do j=1,nstd
          sq = (rax - sra(j))**2 + (decx-sdec(j))**2
          if (sq.lt.sqerr) then
          k = j
          goto 408
          endif
         enddo
c we have a new star
      nstd = nstd+1
      k = nstd
      do j=1,6
       smag(k,j) = 99.999
       serr(k,j) = 9.999
       enddo
       star(k) = tname(1:8)//tname(13:15)
       sra(k) = rax
       sdec(k) = decx
c we have an existing star
408   continue
        do j=1,5
        jj = j+6
        smag(k,jj) = tmag(j)
        serr(k,jj) = terr(j)
        enddo
409   continue
      goto 405
410   continue
      close(2)
c
c read in special extinction stars
c
      open (unit=2,file='sdssext.txt',status='old')
415   continue
         nstd=nstd+1
         read (2,978,end=417) sra(nstd),sdec(nstd),smag(nstd,2),
     $      smag(nstd,3),smag(nstd,8),smag(nstd,9),smag(nstd,10)
978      format (2f14.6,5f8.3)
            star(nstd) = 'SDSSEXT'
         goto 415
417   continue
      close(2)
      nstd = nstd-1
      print *,"number of standards= ",nstd
c     do i=1,nstd
c       print *,star(i),sra(i),sdec(i),smag(i,2),smag(i,3)
c     enddo
c
c need to read through file, find all standard star matches.
c
      warn = 'N'
      n=1
      jj = 0
      nfirst = 1
c skip over header
      read(1,*)
      read (1,901,end=500) hjd,ray,decy,filno,am,amag,amagerr,
     $  iset,istar
      backspace(1)
c
c  start loop over stars (n) and their filters (k)
c  
450   continue
      do i=1,MAXFILT
        umag(n,i) = 99.999
        uerr(n,i) = 9.999
        xx(n,i) = 9.999
      enddo
c read in one star
460   continue
      read (1,901,end=500) hjd,ray,decy,filno,am,amag,amagerr,
     $  isetnew,istarnew
901   format (f12.5,f12.7,f12.7,i5,f7.3,
     $  9x,9x,7x,7x,8x,7x,5x,f8.4,f8.4,
     $  8x,8x,5x,2x,25x,6x,5x,i5,5x,i7)
c    $  8x,8x,5x,6x,6x,6x,6x,5x,i5,5x,i7)
      if (istarnew.eq.istar.and.iset.eq.isetnew) then
        if(filno.le.MAXFILT) then
        umag(n,filno) = amag
        uerr(n,filno) = amagerr
        xx(n,filno) = am
        rax=ray
        decx=decy
        endif
        goto 460
      else
       istar = istarnew
       iset = isetnew
       backspace(1)
      endif
c
c  get the standard star data from the star list file
c
      call search(warn,rax,decx,sname,dmag,nindx,nfirst)
      if (sname.eq.'NEXT') goto 450
      sindx(n,1) = nindx
      sindx(n,2) = n
      n=min((n+1),MAXSTARS)
      go to 450
500    nstar=n-1
      print *,'nstar for fits: ',nstar
c
c check to see what kind of transformation the observer wants
c
      print *
510   print *,'  Do you want to find the zero points ONLY (y/n)?'
      print *,'  (The transformation coefficients must already be ',
     $          'known.) --- '
      call prompt(zpt,flag)
      if(flag .eq. 0) go to 510 
c
520   if(zpt .eq. 'Y') then
        open(unit=3,file=file3,status='old',err=9999)
        read (3,'(a)') fileid
        if(fileid .ne. 'TRANSFORMATION COEFFICIENTS') then
          call bell
          print *,'  The TRAN file has an improper file type.'
          print *,'  Enter the TRAN file name --- '
          read *,file3
          close (3)
          go to 520 
        endif
        call backup(file3)
c       call zcheck(filt,imag,nc,icolor,fext,sext,coef)
      else
        call fcheck(file3)
      endif
4030  print *,'Do you want standards written to file? (Y/N):'
      call prompt(ifl,flag)
c
      go to 600
530   call bell
      print *,'  TRAN file cannot be opened.'
      print *,'  Enter the TRAN file name --- '
      read *,file3
      go to 520 

c
c now do fit for each magnitude
c
600   continue
        print *,'which filter to fit (0 quits): '
        read (5,*) nfilt1
        if (nfilt1.eq.0) goto 9999
        print *,'1st filter for color: '
        read (5,*) nfilt2
        print *,'2nd filter for color: '
        read (5,*) nfilt3
        write (char1,902) nfilt1
        write (char2,902) nfilt2
        write (char3,902) nfilt3
902     format (I3)
c
c now check and see if we have enough observations for fit using these
c 3 filters
c
        n = 1
        do i=1,nstar
           j = sindx(i,1)
           k = sindx(i,2)
           if (smag(j,nfilt1).lt.90.0.and.smag(j,nfilt2).lt.90.0.and.
     $        smag(j,nfilt3).lt.90.0.and.umag(k,nfilt1).lt.90.0.and.
     $        umag(k,nfilt2).lt.90.0.and.umag(k,nfilt3).lt.90.0) then
             indx(n) = 1
             xstar(n) = star(j)
             xmag(n,1) = smag(j,nfilt1)
             xmag(n,2) = smag(j,nfilt2) - smag(j,nfilt3)
             ymag(n,1) = umag(k,nfilt1)
             ymag(n,2) = umag(k,nfilt2) - umag(k,nfilt3)
             x(n,1) = xx(k,nfilt1)
             x(n,2) = (xx(k,nfilt2) + xx(k,nfilt3))/2.
             n = n + 1
           endif
         enddo
         nfit = n - 1
      print *,'nfit ',nfit
c     do i=1,nfit
c     print *,i,xstar(i),xmag(i,1),ymag(i,1),x(i,1)
c     enddo
c
c  make sure that there are enough observations to do a legitimate solution.
c   The number of equations per star is nc+1 with 3(nc+1) unknowns.  So at
c   least 3 observations required (only 2 if you don't solve for first order
c   extinction.
c
      if(zpt .ne. 'Y') then
        if(nfit .lt. 3) then
          call bell
          print *,'   Insufficient number of observations, at least ',
     $            'three stars required.'
          go to 9999
        endif
      endif
      if(flag .eq. 0) go to 4030
c
      if (ifl.eq.'Y') then
        open (unit=66,file='temp.yy',status='new')
        do ka=1,nfit
          write (66,9070) xstar(ka)
9070      format (a11)
          write (66,9071) xmag(ka,1),ymag(ka,1),x(ka,1),
     $       xmag(ka,2),ymag(ka,2),x(ka,2)
9071      format(5x,5f7.3)
        enddo
        close(66)
      endif
c
c main iteration loop
c
      iiflag = 0
75    continue
c
c  check range of colors
c
      IF(zpt .eq. 'N') then
        colmax=ymag(1,2)
        colmin=ymag(1,2)
        do m=2,nfit
          if(ymag(m,2) .gt. colmax) colmax=ymag(m,2)
          if(ymag(m,2) .lt. colmin) colmin=ymag(m,2)
        enddo
        if(colmax-colmin .lt. 0.5) then
          print *,'  The first color index only spans a range of ',
     $            colmax-colmin
          print *,'  This may lead to a poor determination of the ',
     $            'transformation coefficients.'
8010      print *,'  Procede anyway (y/n)? --- '
          call prompt(doit,flag)
          if(flag .eq. 0) go to 8010
          if(doit .eq. 'N') go to 999
        endif
      ENDIF
c
c  get the 2nd order extinction coeff
c
      if(zpt .eq. 'N') then
        print *,'  Enter the second order extinction coefficients: '
8001    print *,' k"('//char1//')= '
        read(5,*,err=8001) sext(1)
8002    print *,' k"('//char2//'-'//char3//')= '
        read(5,*,err=8002) sext(2)
      endif
c    
c  find out what kind of solution to do
c
379   print *,'  Do you wish to solve for 1st order extinction',
     $          ' (y/n)? --- '
      call prompt(doit,flag)
      if(flag .eq. 0) go to 379
c
      IF(doit .eq. 'Y') then
c
c  check air mass range
c
        xmax=x(1,1)
        xmin=x(1,1)
        do m=1,nfit
          if(x(m,1) .gt. xmax) xmax=x(m,1)
          if(x(m,1) .lt. xmin) xmin=x(m,1)
        enddo
        xdif=xmax-xmin
        if(xdif .lt. .7) then
          print *,'  Your air masses only span a range of ',xdif
          print *,'  This may lead to poor determination of ',
     $              'the extinction coefficients.'
161       print *,'  Do you wish to proceed anyway (y/n)? --- '
          call prompt(doit,flag)
          if(flag .eq. 0) go to 161
          if(doit .eq. 'N') then
1611        print *,'  Do you wish to enter the extinction ',
     $                  'coefficent (y/n)? --- '
            call prompt(doit,flag)
            if(flag .eq. 0) go to 1611
            if(doit .eq. 'N') go to 999
            doit = "Z"
          endif
        endif
        if (doit.ne."Z") then
c *********
c  do the solution with extinction, or only find zero pts + extinction
c *********
         nn=2
         if(zpt .eq. 'N') then
           do i=1,nfit
              print *,i,x(i,1),xmag(i,1),xmag(i,2),ymag(i,1),ymag(i,2)
           enddo
           call solve3(nfit,nn,x,xmag,ymag,sext,fext,coef,zero)
         else
           call zextinc(1,nfit,nn,x,xmag,ymag,sext,fext,coef,zero,
     $                 filt,icolor,imag)
         endif
        endif
      ELSE
        doit = "Z"
      ENDIF
      print *,'sext',sext
      print *,'fext',fext
      print *,'coef',coef
      print *,'zero',zero
c     goto 30080
c *********
c  Solve with extinction coefficients fixed
c *********
c  Enter the extinction coefficients
c
        if (doit.eq."Z") then
8999    continue
8003      print *,'  Enter  k('//char1//'): '
          read(5,*,err=8003) fext(1)
8004      print *,'  Enter  k('//char2//'-'//char3//'): '
          read(5,*,err=8004) fext(2)
c
c  Solve for coefficients + zero points
c
7500  continue
        if(zpt .eq. 'N') then
          do i=1,nfit
c
c    solve for (smag-umag)o vs (color)
c
            y1(i)=xmag(i,1)-ymag(i,1)+fext(1)*x(i,1)+sext(1)*
     $                ymag(i,2)*x(i,2)
            x1(i)=xmag(i,2)
          enddo
          call lin (x1,y1,nfit,coef(1),zero(1))
c
c  plot results and overlay the fit
c
        ylabel = char1//'-Inst. Mag.'
        fil3 = char2//'-'//char3
        write(xlabel,'(a13)') fil3
c
        xbeg=2.0
        ybeg=3.0
        xlen=8.0
        ylen=5.0
        icolor1=15
        icolor2=12
        slope=coef(1)
        b=zero(1)
        connect=.false.
        mflag=.false.
        iiflag = 0
        call disply2 (indx,iiflag)
c
c   solve for (stand color-obs color)o vs stand color
c
                do 3002 i=1,nfit
                y1(i)=xmag(i,2)-ymag(i,2)+fext(2)*x(i,2)+sext(2)*
     $                ymag(i,2)*x(i,2)
3002            x1(i)=xmag(i,2)
                call lin (x1,y1,nfit,a,b)
c
c  plot results and overlay the fit
c
      fil3 = char2//'-'//char3
      ylabel = '('//fil3(1:nblank(fil3))//')-Int Col'
      write(xlabel,'(a13)') fil3
c
                slope=a
                print *,'slope ',slope,b
        call disply2 (indx,iiflag)
c
                coef(2)=1./(1.-a)
                zero(2)=coef(2)*b
57              continue
c
c now remove bad stars for next iteration
c
      if (iiflag.ne.0) then
        i = 1
9000    continue
          if (indx(i).eq.0) then
            print *,'Deleting star: ',xstar(i)
            if (i.ne.nfit) then
              do j=i,nfit-1
                indx(j) = indx(j+1)
                do k=1,nc+1
                  x(j,k) = x(j+1,k)
                  ymag(j,k) = ymag(j+1,k)
                  xmag(j,k) = xmag(j+1,k)
                  xstar(j) = xstar(j+1)
                enddo
              enddo
            endif
            nfit = nfit - 1
            i = i - 1
          endif
          i = i + 1
          if (i.le.nfit) goto 9000
        iiflag = 0
        goto 7500
      endif
	else
c
c  Solve for zero points only
c
	    call zextinc(0,nfit,nn,x,xmag,ymag,sext,fext,coef,zero,filt,
     $                 icolor,imag)
	endif
	endif
c **********************************************************************
c look at errors
c **********************************************************************
30080 continue
      print *,'sext',sext
      print *,'fext',fext
      print *,'coef',coef
      print *,'zero',zero
	print *,'                          Calculated Values and Errors'
        write (6,157) nfilt1,(char2//' -'//char3)
157	format(/,19x,i6,2x,a10)
	line=1
        do k=1,MAXFC
         sumdev(k) = 0.0
         sumdev2(k) = 0.0
        enddo
	do 7000 k=1,nfit
		c(k,1)=coef(2)*(ymag(k,2)*(1.-sext(2)*x(k,2))-
     $               fext(2)*x(k,2))+zero(2)
		fmag(k)=ymag(k,1)-fext(1)*x(k,1)-sext(1)*ymag(k,2)*
     $               x(k,1)+coef(1)*c(k,1)+zero(1)
c
		delmag= xmag(k,1)-fmag(k)
                sumdev(1) = sumdev(1) + delmag
                sumdev2(1) = sumdev2(1) + delmag*delmag
                delc(1)=xmag(k,2)-c(k,1)
                sumdev(2) = sumdev(2) + delc(1)
                sumdev2(2) = sumdev2(2) + delc(1)*delc(1)
		write(6,7005) xstar(k),fmag(k),c(k,1)
7005		format(1x,a14,1x,6(f7.3,11x))
                iflg(1) = '   '
                if (abs(delmag).gt.0.05) iflg(1)='***'
                do m=1,nc
                iflg(m+1) = '   '
                if (abs(delc(m)).gt.0.05) iflg(m+1)='***'
                enddo
		write(6,7007) delmag,iflg(1),(delc(m),iflg(m+1),m=1,nc)
7007		format(1x,'errors = ',6x,6(f7.3,a3,8x))
c 	line=line+1
	print *
c       if(line .ge. 6) then
c               print *,'Press any key to continue:'
c               read (5,'(a)') adummy
c		line=1
c	endif
7000	continue
c
c now calculate standard deviations
c
        zstar = float(nfit)
        DO m=1,nc+1
         sumdev(m) = dsqrt(abs(zstar*sumdev2(m) - sumdev(m)*sumdev(m))/
     $       (zstar*(zstar-1.)))
        ENDDO
        write (6,7008) (sumdev(m),m=1,nc+1)
7008    format (' St deviations: ',6(f7.3,11x))
               print *,'Press any key to continue:'
               read (5,'(a)') adummy
c  *********************************************************************
c  Make plots of Calc Mag vs True Mag and Calc color vs True Color
c  *********************************************************************
c  do mag first ...
c
	write(xlabel,"('True ',a3)") char1
	write(ylabel,"('Your ',a3)") char1
c
	xbeg=2.0
	ybeg=3.0
	xlen=8.0
	ylen=5.0
	icolor1=15
	icolor2=12
	slope=1				! overlay a 45 deg line for reference
	b=0
	connect=.false.
	mflag=.false.
		do k=1,nfit
		y1(k)=fmag(k)
		x1(k)=xmag(k,1)
		enddo
c
	call disply2(indx,iiflag)
c
c  now do colors
c
	do m=1,nc
c   write(xlabel,9066) char2,char3
9066     format ('True ('//a3//' - '//a3//')')
c   write(ylabel,9067) char2,char3
9067     format ('Your ('//a3//' - '//a3//')')
		do k=1,nfit
		y1(k)=c(k,m)
		x1(k)=xmag(k,m+1)
		enddo
	call disply2(indx,iiflag)
c
	enddo	
c
c now remove bad stars for next iteration
c
      if (iiflag.ne.0) then
        i = 1
90      continue
          if (indx(i).eq.0) then
            print *,'Deleting star: ',xstar(i)
            if (i.ne.nfit) then
              do j=i,nfit-1
                indx(j) = indx(j+1)
                do k=1,nc+1
                  x(j,k) = x(j+1,k)
                  ymag(j,k) = ymag(j+1,k)
                  xmag(j,k) = xmag(j+1,k)
                  xstar(j) = xstar(j+1)
                enddo
              enddo
            endif
            nfit = nfit - 1
            i = i - 1
          endif
          i = i + 1
          if (i.le.nfit) goto 90
        iiflag = 0
        goto 75
      endif
5000  continue
c
c  *********************************************************************
c  write the transformation coeff to screen and to TRAN file
c  *********************************************************************
      write (6,3000) nfilt1,nfilt2,nfilt3,fext,sext,coef,
     $    zero,xdif,nfit,zpt,sumdev
3000   format (3i5,2f9.4,2f9.4,2f9.4,2f9.4,f9.4,i6,1x,a1,
     $   2f9.4,2x,a40)
c
c
c
5060	if(zpt .eq. 'Y') rewind (3)
c
          if (jjflag.eq.0) then
            jjflag = 1
	  write(3,2999) file1
2999	  format('TRANSFORMATION COEFFICIENTS'/
     $    'filename: ',a40/' fil1',' fil2',
     $    ' fil3','   fext1','    fext2','    sext1','    sext2',
     $    '   coef1','   coef2','   zero1','   zero2','    xdif',
     $    '   nfit', '   zp ','   std1 ','   std2 ')
          endif
      write (3,3000) nfilt1,nfilt2,nfilt3,fext,sext,coef,
     $    zero,xdif,nfit,zpt,sumdev
c
               print *,'Press any key to continue:'
               read (5,'(a)') adummy
        goto 600
9999    continue
999	close(1)
	close(3)
c
	end
c
c
      subroutine info(numfil,filt,imag,nc,icolor,korder)
c
c  this subroutine prompts the user for information about which filters
c  are to be used for mag and colors
c
c  input:   numfil=number of filters, 
c           filt= char array of filters in user data set.
c  output:  imag = filter number which corresponds to the magnitude
c           nc = number of colors to form
c           icolor = an integer array which contains the numbers of the 
c                    filters to difference to form the colors.  icolor(1,1)
c                    and icolor(1,2) correspond to the filters used to make
c                    the color index which contains the mag filter (e.g.
c                    B-V and V).  This is required for the transformation eq.
c	    korder = an array containing the order the mag/colors are to
c                     be used from the standard star list.
c
      integer MAXFC,MAXC,MAXFILT
      PARAMETER (MAXFC = 2)
      PARAMETER (MAXC = 1)
      PARAMETER (MAXFILT = 11)
      integer  icolor(MAXC,MAXFC),korder(MAXFC)
      integer numfil,nc
      character*6 filt(MAXFILT)
      character text(40),c1,c2,c3,ok,adummy
      character*3 ihead(MAXFILT)
      character*5 epoch,ra,dec
      data ihead/MAXFILT*'   '/
c
c  read in the filter labels from the star list file
c
      do 400 m=1,MAXFC
400   korder(m)=0
      rewind (2)
      read (2,*)
      read (2,1000) text
1000  format(68x,40a1)
      rewind (2)
c parse, assuming 3 character names
      kl=1
      j=1
2500  if(text(j) .eq. ' ') go to 2000
        c1=text(j)
        c2=text(j+1)
        c3=text(j+2)
        ihead(kl)=c1//c2//c3
        j=j+2
        kl=kl+1
        if (kl.gt.6) goto 2510
2000  continue
        j=j+1
        if(j .le. 40) go to 2500
2510  continue
c  kstand is the mag + number of colors in the star list file (standard stars)
      kstand=kl-1
c
c  Test to make sure that the standard star file contains at least one
c   mag and one color
c
      if(kstand .lt. 2) then
        call bell
        print *,'  The standard star file must contain at least 1 ',
     $            'magnitude and 1 color'
        print *,'Press any key to continue:'
        read (5,'(a)') adummy
        stop
      endif
c
      print *,'  The star list file contains the following '
     $  ,'magnitude and colors...'
      print 3000,(m,ihead(m),m=1,kstand)
3000  format(5x,5(5x,i1,'...',a3))
      print *
      print *,'  These are the only allowed magnitude and colors you '
     $  ,'can form from '
      print *,'  your observations.  The minimum configuration '
     $  ,'you can solve for '
      print *,'  is one magnitude and a color which contains that'
     $  ,' magnitude.'
      print *,'  For example, V and B-V.'
c
c  ID magnitude filter
c    this relys on the fact that filter designations from the star list
c    file will contain 3 characters, except for the magnitude
c
      imag=0
      do i=1,numfil
        do k=1,kstand
          if(filt(i) .eq. ihead(k)) then
            korder(1)=k
            imag=i
          endif
        enddo
      enddo
c
c  check to make sure a magnitude filter was found
c
      if(imag .eq. 0) then
        call bell
        print *,'  The magnitude filter in the star-list file',
     $                  ' is not found in your data file.'
        print *,'Press any key to continue:'
        read (5,'(a)') adummy
        stop
      endif
c
c  ID the color which contains the mag filter and make its filters the
c   first elements of icolor
c
      print *
781   print *,'  Which item number above corresponds to the COLOR',
     $          ' index that contains '
        print 2001,kstand
2001  format('   the magnitude filter? (1-',i1,'): ')
      read *,item
      if((item.lt.1).or.(item.gt.kstand)) then
        call bell
        go to 781
      endif
c
c		******	does the selected color contain the mag filter?
c
      if((ihead(item)(1:1).ne.filt(imag)).and.(ihead(item)(3:3).ne.
     $      filt(imag))) then
        call bell
        print *,'  This color does not contain the magnitude ',
     $                  'filter'
        go to 781
      endif
c
c    		******	are both required filters present in users data?
c
      nc=0
      do i=1,numfil
        if(ihead(item)(1:1).eq.filt(i)) then
          do i2=1,numfil
            if(ihead(item)(3:3).eq.filt(i2)) then
              icolor(1,1)=i
              icolor(1,2)=i2
              nc=1
            endif
          enddo
        endif
      enddo
c
      if(nc .eq. 0) then
        call bell
        print *,'  Your data file does not contain the needed'
     $                 ,' filters'
        go to 781
      endif
c
c  ID the color filters
c
      korder(2)=item
      print *
15    nc=1
      do 101 k=1,kstand
        if(k .eq. korder(1)) go to 101
        if(k .eq. korder(2)) go to 101
        do i=1,numfil
          if(ihead(k)(1:1) .eq. filt(i)) then
            do i2=1,numfil
              if(ihead(k)(3:3) .eq. filt(i2)) then
505             print *,'  Do you wish to solve for (',ihead(k),
     $                             ') (y/n) ? --- '
                call prompt(ok,flag)
                if(flag .eq. 0) go to 505
                if(ok .eq. 'N') go to 506
                nc=nc+1
                icolor(nc,1)=i
                icolor(nc,2)=i2
                korder(nc+1)=k
506             continue
              endif
            enddo
          endif
        enddo
101   continue
c
      if(nc .lt. 1) then
        call bell
        print *,'  You MUST solve for at least one color'
        go to 15
      endif
c
      print *
      print *,'  The magnitude and color(s) to be processed from ',
     $          'your data file...'
      print 3001,(ihead(korder(m)),m=1,nc+1)
3001  format(10x,5(5x,a3))
c
      return
      end

      subroutine search(warn,rax,decx,sname,dmag,nindx,nfirst)
c
c  This routine searches the standard star list for object identified
c  by the character array STAR and returns the magnitude and colors
c
c  Input: star name (star), number of filters (nf), warning flag (warn)
c  Output: mag and colors (dmag)
c
      INTEGER MAXSTD,MAXFC,MAXFILT
      PARAMETER (MAXSTD = 70000)
      PARAMETER (MAXFC = 2)
      PARAMETER (MAXFILT = 14)
      real*8 sra(MAXSTD),sdec(MAXSTD),err,errmax
      real*8 rax,decx
      real*4  smag(MAXSTD,MAXFILT), dmag(MAXFC),
     $    serr(MAXSTD,MAXFILT)
      integer nlast,nfirst,nstd,i,nindx
      character*11 star(MAXSTD),sname,star2
      character warn,ds
      common /stds/ sra,sdec,smag,serr,nstd,star
c
      sname = " "
      errmax = 2.*(5./3600.)**2
300   continue
      iflag = 0
      DO i=1,nstd
        err = (rax-sra(i))**2 + (decx-sdec(i))**2
        IF (err.lt.errmax) THEN
          DO j=1,MAXFILT
           dmag(j) = smag(i,j)
          ENDDO
          sname = star(i)
          iflag = 1
          nindx = i
          goto 320
        ENDIF
      ENDDO
320   continue
      if (iflag.eq.1) RETURN
c
c  not found
c
999   if(warn .eq. 'Y') then
        call bell
        print *,'  ',sname,'not found in standard star list'
        print *,'  These stars were found ...'
c       call display(items,sname)
c       call display(1,sname)
        print *
        print *,'  Enter a new choice (type NEXT to skip) --- '
        read (5,'(a14)') sname
        star2=sname
        call ucase(star2,14)
        if(star2 .eq. 'NEXT') return
        go to 300
      else
        sname='NEXT'
        return
      endif
      end

	subroutine display(nstar,star1)
c
c  this displays nstar names on the screen, six across and pauses when
c   the screen is full
c
      integer MAXNC
      PARAMETER (MAXNC = 10000)
      integer nstar
	character*14 star1(maxnc)
      character adummy
c
c  irows = number of complete rows of 6 stars each
c  iover = number of stars in last partial row
c  iscrns = number of full screens of 20 lines
c  ilast = number of complete lines on last partial screen
c
	irows=nstar/6
	iover=nstar-irows*6
	iscrns=irows/20
	ilast=irows-iscrns*20
c
	if(irows .lt. 20) then
		print 200, (star1(k),k=1,nstar)
200		format(6(3x,a10))
		return
	else
		k=1
		iend=iscrns
		ilines=20
107		do 105 ir=1,iend
		do 100 il=1,ilines
		print 200, (star1(ik),ik=k,k+5)
100		k=k+6
105            print *,'Press any key to continue:'
               read (5,'(a)') adummy
		iscrns=iscrns-1
		if(iscrns .eq. 0) then
		iend=1
		ilines=ilast
		if(ilast .ge. 1) go to 107
		endif
		if(iover .gt. 0) print 200,(star1(ik),ik=k,k+iover-1)
	endif
	return
	end
	

      subroutine zextinc(iext,nstar,nc,x,smag,urmag,sext,fext,coef,zero,
     $              filt,icolor,imag)
c
c  this subroutine determines the zero point values with or without
c   an extinction determination
c
c   iext=0 extinction fixed
c   iext=1 determine the extinction
c NOTE: dimensions sext,fext,coef,zero to match other routines
c then used index (1) for all uses of these parameters
c
      integer MAXSTARS,MAXFC,MAXC,MAXFILT
      PARAMETER (MAXSTARS = 70000)
      PARAMETER (MAXC=1)
      PARAMETER (MAXFILT=11)
      PARAMETER (MAXFC=2)
      real*4 x(MAXSTARS,MAXFC),smag(MAXSTARS,MAXFC),
     $    urmag(MAXSTARS,MAXFC),
     $       sext(MAXFC),fext(MAXFC),coef(MAXFC),zero(MAXFC),
     $       x1(MAXSTARS),y1(MAXSTARS)
      integer icolor(MAXC,MAXFC),iext,nstar,nc
      character*6 filt(MAXFILT)
      character*30 xlabel,ylabel
      character*13 fil3
      logical connect,mflag
      common /plotit/ x1,y1,npts,xbeg,ybeg,xlen,ylen,xlabel,ylabel,
     $                    icolor1,connect,slope,b,icolor2,mflag
c
c  Do extinction case first
c
      if(iext .eq. 1) then
c
c  do mag first
c
      do 100 i=1,nstar
        x1(i)=x(i,1)
100   y1(i)=smag(i,1)-urmag(i,1)+sext(1)*urmag(i,2)*x(i,1)
     $        -coef(1)*smag(i,2)
      call lin(x1,y1,nstar,fext,zero)
      write (6,980) nstar,fext,zero
980   format ('linfit: ',i2,2f12.6)
      write (6,981) (i,x1(i),y1(i),i=1,nstar)
981   format (i5,2f12.6)
c
c  plot results and overlay the fit
c
      write(ylabel,"('Stand-Inst ',a6)") filt(imag)
      xlabel='Air Mass'
c
      npts=nstar
      xbeg=2.0
      ybeg=3.0
      xlen=8.0
      ylen=5.0
      icolor1=15
      icolor2=12
      slope=fext(1)
      b=zero(1)
      connect=.false.
      mflag=.false.
c
      call disply
      fext=-fext
c
c  Now do the non-extinction case
c
      ELSE
c
c  do mag first
c
      tot=0.
      do 500 i=1,nstar
500   tot=tot+smag(i,1)-urmag(i,1)+sext(1)*urmag(i,2)*x(i,1)
     $        -coef(1)*smag(i,2)+fext(1)*x(i,1)
      zero=tot/nstar
c
c  do colors
c
      do 600 km=1,nc
        tot=0.0
        do 700 kl=1,nstar
          x1(kl)=x(kl,km+1)
700     tot=tot+smag(kl,km+1)-coef(1)*urmag(kl,km+1)
     $             *(1.-sext(1)*x1(kl))
     $             +coef(1)*fext(1)*x1(kl)
        zero=tot/nstar
600   continue
      ENDIF
      return
      end
c
	subroutine fcheck(file3)
c
c  This routine checks for the existence of a file, asks the user about
c   overwriting and making a backup file
c
	character*30 file3
	character doit
	logical found
c
700	inquire(file=file3,exist=found)
c
	if(found) then
720	  	print *,'  This TRAN file already exists, do you wish ',
     $                    'to Overwrite (y/n)? --- '
	  	call prompt(doit,flag)
	        if(flag .eq. 0) go to 720
	  	if(doit .eq. 'Y') then 
722		    print *,'  Are you sure (y/n)? --- '
		    call prompt(doit,flag)
		    if(flag .eq. 0) go to 722
	  	    call backup(file3)
		else
		print *,'  Enter TRAN file name --- '
		read *,file3
		go to 700
		endif
	  endif
		open(unit=3,file=file3,status='unknown')
c
	return
	end
c
c
	subroutine backup(file1)
c
c  this routine creates a backup copy of a file before overwriting it
c
	character*30 file1,fx,input
	character*80 str
c
	do j=1,30
	fx(j:j)= ' '
	enddo
c
c  locate start of file extension
c
	i=index(file1,'.')
	if(i .ne. 0) then
		do k=1,i
		fx(k:k)=file1(k:k)
		enddo
		fx(i+1:i+3)='BAK'
	else
		do k=1,30
		if(file1(k:k) .eq. ' ') go to 100
		fx(k:k)=file1(k:k)
		enddo
100		fx(k:k+3)='.BAK'
	endif
c
c
	print *,'  Enter the backup file name [',fx(1:nblank(fx)),']: '
	read(5,'(a)') input
	if(input .ne. ' ') fx=input
c
      str='copy '//file1//' '//fx
      return
      end

       subroutine disply2 (indx,iiflag)
c
c  Common block paramaters (a common is used to prevent the system stack
c   from exceeding 64KB.)
c
c       x        the x data array              
c       y        the y data array              
c       npts      number of data points         
c       xbeg      x starting position in inches 
c       ybeg      y starting position in inches 
c       xlen      length of x axis in inches    
c       ylen      length of y axis in inches    
c       xlab      x axis label                  
c       ylab      y axis label                  
c       icolor1   color index to use for data
c       connect   logical true to connect data points
c       slope     slope of fit
c       b         intercept of fit
c       icolor2   color index of fit
c       mflag     logical variable indicating if this is a magnitude plot
c
c
      integer MAXSTARS
      PARAMETER (MAXSTARS=70000)
      real*4  x(MAXSTARS),y(MAXSTARS),rarray(7)
      real*8 dist,dmin,x1(2),y1(2)
      real*8 xmx,xmn,ymn,ymx,xz(MAXSTARS),yz(MAXSTARS)
      real*8 xx,yy
      integer indx(MAXSTARS),iiflag
      integer*2 iarray(9)
      character*30 xlab,ylab
      character*10 str
      character*40 string
      character key*1
      logical connect,mflag
      common /plotit/ x,y,npts,xbeg,ybeg,xlen,ylen,xlab,ylab,
     $  icolor1,connect,slope,b,icolor2,mflag
      data iflag /0/
c
c  init the graphics package, enter graphics mode
c
      IF (iflag.eq.0) THEN
        i = sm_device('X11 -g 800x800')
        iflag = 1
        call sm_graphics
      ENDIF
      call sm_erase
      xmin = 5.e20
      ymin = 5.e20
      xmax = -5.e20
      ymax = -5.e20
      DO i=1,npts
        xmin = min(xmin,x(i))
        xmax = max(xmax,x(i))
        ymin = min(ymin,y(i))
        ymax = max(ymax,y(i))
        xz(i) = x(i)
        yz(i) = y(i)
      ENDDO
      xmx = xmax
      xmn = xmin
      ymx = ymax
      ymn = ymin
      call sm_limits(xmn,xmx,ymn,ymx)
c     call sm_limits(xmin,xmax,ymin,ymax)
      call sm_ctype ('white')
      call sm_box (1,2,0,0)
      call sm_xlabel(xlab)
      call sm_ylabel(ylab)
      call sm_gflush
      call sm_ctype ('green')
      IF (connect) THEN
        call sm_conn(xz,yz,npts)
c       call sm_conn(x,y,npts)
      ELSE
        call sm_ptype(2.43d2,1)
c       call sm_ptype(2.43e2,1)
        call sm_expand(1.1d0)
c       call sm_expand(1.1e0)
        call sm_points(xz,yz,npts)
c       call sm_points(x,y,npts)
        call sm_ptype(1.1d1,1)
c       call sm_ptype(1.1e1,1)
        call sm_expand(1.0d0)
c       call sm_expand(1.0e0)
      ENDIF
      call sm_gflush
      x1(1) = xmin
      y1(1) = b + slope*x1(1)
      x1(2) = xmax
      y1(2) = b + slope*x1(2)
      call sm_ctype('red')
      call sm_lweight(2.d0)
c     call sm_lweight(2.e0)
      call sm_conn(x1,y1,2)
c     call sm_conn(x1,y1,2)
c
      call sm_gflush
      x1(1) = xmin
      y1(1) = b + slope*x1(1)
      x1(2) = xmax
      y1(2) = b + slope*x1(2)
      call sm_ctype('red')
      call sm_conn(x1,y1,2)
c     call sm_conn(x1,y1,2)
      call sm_gflush
c
c loop to remove discrepant points
c
800   continue
      print *,'d=delete,u=undelete,e=end:'
      call sm_curs(xk,yk,k)
      key=char(k)
      IF (key.eq.'d') THEN
        dmin=1.E38
        indx2 = 1
        iiflag = 1
        DO i=1,npts
          dist=(xk-x(i))**2+(yk-y(i))**2
          IF (dist.le.dmin) THEN
            dmin=dist
            indx2=i
          ENDIF
        ENDDO
        indx(indx2) = 0
        print 9090,x(indx2),y(indx2)
9090     format ('Closest point is: ',2f9.2)
      ELSEIF (key.eq.'u') THEN
        dmin=30000.
        indx2 = 1
        DO i=1,npts
          dist=(xk-x(i))**2+(yk-y(i))**2
          IF (dist.le.dmin) THEN
            dmin=dist
            indx2=i
          ENDIF
        ENDDO
        indx(indx2) = 1
        print 9090,x(indx2),y(indx2)
      ELSEIF (key.eq.'e') THEN
        RETURN
      ENDIF
      goto 800
      end

	subroutine solve3 (n,nfilt,x,smag,urmag,sext,fext,coef,zero)
c
c	a modified version of SOLVE by A. Henden 1980
c       last modified June 13, 1989 by Kaitchuck
c
c  	least squares solution of the UBVRI transformation equations
c   input -
c	smag(n,5)   standard magnitudes and colors in the form:
c		smag(n,1) = mag
c		smag(n,2) = color 1
c		smag(n,3) = color 2
c		smag(n,4) = color 3
c		smag(n,5) = color 4
c	urmag(n,5)  your corresponding instrumental magnitudes & colors
c	x(n,5)      air mass values
c	n           number of observations to reduce
c	nfilt       mag + number of colors = (1-5)
c	sext(5)     second order extinction
c   output -
c	fext(5)     first order extinction
c	coef(5)     transformation color coefficients
c	zero(5)     zero points
c
c	to have fext coef and zero to all be valid, you
c       must use standards with both a wide range of colors
c	and a wide range of airmases.  with a cluster, you
c	will get valid coef parameters, but fext and zero may
c	be invalid...so beware!
c
c
      INTEGER MAXSTARS,MAXFC
      PARAMETER (MAXSTARS = 70000)
      PARAMETER (MAXFC = 2)
	real*4 smag(maxstars,MAXFC),urmag(maxstars,MAXFC),
     $  x(maxstars,MAXFC),
     $  sext(MAXFC),fext(MAXFC),coef(MAXFC),zero(MAXFC)
c   loop over colors
	do 20 k=1,nfilt
c   zero matrix elements
	a1=0.
	a2=0.
	a3=0.
	a4=0.
	a5=0.
	a6=0.
	a7=0.
	a8=0.
	a9=0.
	a10=0.
	a11=0.
	a12=0.
c  loop over number of standard stars
	do 10 i=1,n
c  calculate matrix elements sums
	temp=urmag(i,k)*(1.-sext(k)*x(i,k))
	std=smag(i,k)
	if(k .ne. 1) go to 5
c  note: the mag equation has slightly different form
	temp=smag(i,2)
	std=std-urmag(i,k)
5	continue
c
c  the matrix looks like:
c	( a1  a2   a3   . a4  )
c       ( a5  a6   a7   . a8  )
c       ( a9  a10  a11  . a12 )
c  note: a2=a5, a3=a9, and a7=a10.  but explicit here for clarity
c
	a1=a1+temp*temp
	a2=a2+temp*x(i,k)
	a3=a3+temp
	a4=a4+std*temp
	a5=a5+temp*x(i,k)
	a6=a6+x(i,k)*x(i,k)
	a7=a7+x(i,k)
	a8=a8+std*x(i,k)
	a9=a9+temp
	a10=a10+x(i,k)
	a11=a11+1
	a12=a12+std
10	continue
c   calculate minors
	aa=a7*a10-a6*a11
	bb=a5*a11-a7*a9
	cc=a6*a9-a5*a10
	dd=a8*a11-a7*a12
	ee=a6*a12-a8*a10
	ff=a5*a12-a8*a9
c  calculate determinant
	det=a1*aa+a2*bb+a3*cc
c  solve for wanted values
	coef(k)=(a4*aa+a2*dd+a3*ee)/det
	zero(k)=(a2*ff-a1*ee+a4*cc)/det
	fext(k)=(a1*dd-a4*bb+a3*ff)/det
	if(k .ne. 1) fext(k)=fext(k)/coef(k)
20	continue
	return
	end
