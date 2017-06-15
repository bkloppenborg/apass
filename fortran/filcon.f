      program filcon
c
c  This program reads an instrumental data file, reads the TRAN file
c   and converts the observations to a standard photometic system
c
c   Written by R. Kaitchuck
c   modified 6-Apr-90 AAH
c   modified 7-May-93 AAH for nofs format
c   modified 29-Jul-93 AAH to be more flexible
c   modified 9-Nov-93 AAH to add errs
c   modified 18-Nov-93 AAH for arbitrary number of filters
c   mod 08-Sep-98 first filter fix for new format
c   mod 17-Feb-2003 aah for filter vs. color format
c   mod 07-Oct-2009 aah new file format
c   mod 03-Jul-2010 aah magnitude format
c   mod 07-Mar-2013 aah add ZS,Y
c   mod 14-Jun-2017 aah new pepp file format
c
c   Fortran 77 compatible
c
c   External subroutines: Prompt, Bell, openit, getfilt, match, ucase
c
      INTEGER MAXTRAN,MAXFILT
      PARAMETER (MAXFILT=14)
      PARAMETER (MAXTRAN=20)
      REAL*8 hjd(MAXFILT),ra(MAXFILT),dec(MAXFILT),rax,decx,hjdx
      REAL*4 coef(MAXTRAN,2),zero(MAXTRAN,2),fext(MAXTRAN,2),
     $   sext(MAXTRAN,2),umag(MAXFILT),urmag(MAXFILT),
     $   uerr(MAXFILT),urerr(MAXFILT),xx(MAXFILT)
      REAL*4 ut,xn,raerr,decerr,amver
      REAL*4 xccdx,xccdy,ccdx(MAXFILT),ccdy(MAXFILT)
      REAL*4 amass,amax,apmag,aperr,dmagx,dmag(MAXFILT)
      INTEGER photflag,propid,propfield,kset,ksetnew,istar,istarnew
      INTEGER satflags(MAXFILT),flags(MAXFILT),fflag,tflag
      INTEGER i,j,k,igroup,isys,isystem(MAXFILT)
      INTEGER iparm(20,3),ngt,night(MAXFILT)
      CHARACTER file3*50,file2*50,fileid*80
      CHARACTER object*25,objectnew*25
      character*30 type
      CHARACTER*10 star,skip,akcess
c
c ***********************************************************************
      print *,'                      Program FILCON version 3.2'
      print *,'                         14-Jun-2017 '        
      print *
c ************************************************************************
c
      print *,'This program converts an instrumental data file to ',
     $          'the standard photometric'
      print *,'system.  The input file comes from program AMATCH',
     $          ' and the transformation '
      print *,'coefficients come from the TRAN.DAT file.'
      print *,'  This version for magnitudes only.'
      print *
      print *
c
c  open and read the TRAN file
c
      open(unit=1,file='tran.dat',status='old',err=10)
      read (1,'(a)') fileid
      IF (fileid .ne. 'TRANSFORMATION COEFFICIENTS') THEN
        call bell
        print *,'  The TRAN file has an improper file ID.'
        close (1)
        GOTO 999
      ENDIF
      GOTO 20
10    call bell
      print *,'  TRAN file cannot be opened.'
      go to 999
20    continue
c
c  read all coefs
c
      read (1,*)
      read (1,*)
      i=1
25    continue
        read (1,3000,end=27) iparm(i,1),iparm(i,2),iparm(i,3),
     $  fext(i,1),fext(i,2),sext(i,1),sext(i,2),coef(i,1),
     $  coef(i,2),zero(i,1),zero(i,2)
3000    format (3i5,2f9.4,2f9.4,2f9.4,2f9.4)
        i = min (i+1,MAXTRAN)
        goto 25
27    continue
      close(1)
c note: can only have maxtran-1 transformations; kludge
      ntran = i - 1
      print *,'number of trans equations: ',ntran
c
c  open the input data file
c
30    print *,'  Enter the input data file name --- '
      read (5,'(a)') file3
      open(unit=3,file=file3,status='old',err=35)
      read (3,'(a)') fileid
      IF (fileid(1:24) .ne. '#INSTRUMENTAL MAGNITUDES') THEN
        call bell
        print *,'  This is the wrong file type, try again.'
        close (3)
        GOTO 30
      ENDIF
      go to 40
35    continue
      call bell
      print *,'  Unable to open this file, try again'
      GOTO 30
c
c  get output file name
c
40    print *,'  Enter the output file name --- '
      read (5,'(a)') file2
      IF(file2 .eq. file3) THEN
      	call bell
      	print *,'  Input and output files can not be the same.'
      	close (3)
      	GOTO 30
      ENDIF
c
      print *,'Saturation limit (fp): '
      read (5,*) satur
      type = 'STANDARD MAGNITUDES AND COLORS'
c     call openit(2,file2,type,.true.,akcess,num,filter)
      open(unit=2,file=file2,status='new')
c
c
c  write output header id
c
      IF(akcess .ne. 'append') THEN
          write(2,9008)
9008      format('# STANDARD MAGNITUDES ONLY')
          write(2,9018)
9018      format('# FILCON ver 3.2'/
     $   ,'# RA (J2000)',4x,'DEC',4x,'CCDX',6x,'CCDY'5x,
     $   'Flags',6x,'HJD',3x,
     $   'Airmass',2x,'Set',5x,'Group',4x,'Field',
     $   6x,'Filt',3x,'Mag',5x,'Error',4x,'dmag',
     $   4x,'sys',1x,'night')
      ENDIF
c
c read first line and backspace for initialization
c
        read (3,*)
        read (3,9003,end=999) hjdx,rax,decx,j,amass,
     $   xccdx,xccdy,amax,apmag,aperr,
     $   photflag,object,ngt,kset,istar
9003   format (f12.5,f12.7,f12.7,i5,f7.3,2f9.3,14x,
     $   f8.1,12x,f8.4,f8.4,16x,i5,2x,a25,
     $   i6,5x,i5,5x,i7,f8.4)
        backspace(3)
c
c  start loop over stars and their filters (k)
c note this version misses last star of file
c  
      igroup = 1
100   continue
c
c read in one star
c
      iflag = 0
        do i=1,MAXFILT
          umag(i) = 99.999
          uerr(i) = 9.999
          xx(i) = 1.
          flags(i) = photflag
        enddo
c loop until reach end of this star's data
110     continue
         read (3,9003,end=999) hjdx,rax,decx,j,amass,
     $   xccdx,xccdy,amax,apmag,aperr,photflag,
     $   objectnew,ngt,ksetnew,istarnew,dmagx
         if (objectnew.eq.object.and.kset.eq.ksetnew) then
           umag(j) = apmag
           uerr(j) = aperr
           xx(j) = amass
           ra(j) = rax
           dec(j) = decx
           hjd(j) = hjdx
           ccdx(j) = xccdx
           ccdy(j) = xccdy
           dmag(j) = dmagx
           isystem(j) = isys
           night(j) = ngt
           satflags(j) = 0
           if (amax.ge.satur) satflags(j) = 1
           flags(j) = photflag
           goto 110
         else
           object = objectnew
           kset = ksetnew
           backspace(3)
         endif
c
c at this point, have all magnitudes for this kset.  Now transform.
c
      do i=1,MAXFILT
        urmag(i) = 99.999
        urerr(i) = 9.999
      enddo
      do i=1,ntran
        i1 = iparm(i,1)
        i2 = iparm(i,2)
        i3 = iparm(i,3)
         if (umag(i1).lt.90.0.and.umag(i2).lt.90.0
     $     .and.umag(i3).lt.90.0.and.urmag(i1).gt.90.0) then
c found transform match
c calculate color
        urc = umag(i2) - umag(i3)
        urx = (xx(i2) + xx(i3))/2.
        stdc = coef(i,2)*(urc*(1.-sext(i,2)*urx) -
     $         fext(i,2)*urx) + zero(i,2)
c calculate magnitude
        urmag(i1) = umag(i1) - fext(i,1)*xx(i1) -
     $    sext(i,1)*stdc*xx(i1) + coef(i,1)*stdc +
     $    zero(i,1)
        urerr(i1) = uerr(i1)
        endif
      enddo
c all possible magnitudes are transformed.  Write results
      do i=1,MAXFILT
        if (urmag(i).lt.90.0) then
      write(2,9005) ra(i),dec(i),
     $  ccdx(i),ccdy(i),flags(i),satflags(i),hjd(i),xx(i),
     $  kset,igroup,object,i,urmag(i),urerr(i),dmag(i),
     $  isystem(i),night(i)
9005  format(2f12.7,2f10.3,1x,i1.1,1x,i1.1,1x,f12.6,
     $  1x,f5.3,i5,1x,i10,1x,a25,1x,i5,1x,3f8.4,2i6)
        endif
      enddo
      igroup = igroup + 1
      GOTO 100
c
999   continue
      close (2)
      stop
      end
