      program filepp
c
c  This program reads an instrumental epp file
c  and the 7 correction matrices.  It then steps through the
c  epp file, applies the photometric correction, and writes
c  a line to the equivalent pepp file.
c  written 11-Jul-2010 aah
c  mod 07-Mar-2013 aah include ZS,Y
c  mod 04-Nov-2016 aah use filelist
c  mod 06-Nov-2016 aah use new epp file format
c  note: hardcoded to ap=2 (rad=3) for now
c
c   Fortran 77 compatible
c
c   External subroutines: Prompt, Bell, openit, getfilt, match, ucase
c
      INTEGER MAXTRAN,MAXFILT
      PARAMETER (MAXFILT=14)
      REAL*8 hjd,ra,dec
      REAL*4 ut,xn,raerr,decerr,amver
      REAL*4 ccdx,ccdy,t,u,xmag,dmag
      REAL*4 amass,amax,apmag,aperr
      REAL*4  cormat(43,43,14),rowcol(43)
      REAL*4 fwhmx,fwhmy,sky,psfmag,psferr
      INTEGER photflag,propid,propfield,kset,ksetnew,istar,istarnew
      INTEGER satflags(MAXFILT),flags(MAXFILT),fflag,tflag
      INTEGER i,j,k,igroup,iap,isystem,jd0,nfile,kser
      CHARACTER file3*50,file2*50,fileid*80,object*25
      CHARACTER residpath*15,eppfile*14,rfile*80,pfile*80
      character*30 type
      CHARACTER*10 star,skip,akcess
c
c ***********************************************************************
      print *,'                      Program FILEPP version 2.0'
      print *,'                         04-Nov-2016     '     
      print *
c ************************************************************************
c
      print *,'This program corrects the photometry of an',
     $          ' instrumental data file.'
      print *,'The input file comes from program HEPP',
     $          ' along with the 7 correction matrices.'
      print *,'  This version for APASS only.'
      print *
      print *
      print *,'Enter filelist to process: '
      read (5,'(a)') pfile
c format pfile:  eppfile resid_path
      open (unit=22,file=pfile,status='old')
5     continue
        read (22,9000,end=999) eppfile,residpath
9000    format(a14,2x,a15)
c
c  open and read the correction matrices
c
      do i=1,43
         do j=1,43
           do k=1,14
              cormat(i,j,k) = 0.0
           enddo
         enddo
      enddo
      rfile = residpath//'/resid_sm_B.txt'
      open (unit=1,file=rfile,status='old')
      do i=1,41
        do j=1,41
          read (1,*) cormat(i+1,j+1,2)
        enddo
      enddo
      close(1)
      rfile = residpath//'/resid_sm_V.txt'
      open (unit=1,file=rfile,status='old')
      do i=1,41
        do j=1,41
          read (1,*) cormat(i+1,j+1,3)
        enddo
      enddo
      close(1)
      rfile = residpath//'/resid_sm_su.txt'
      open (unit=1,file=rfile,status='old')
      do i=1,41
        do j=1,41
          read (1,*) cormat(i+1,j+1,7)
        enddo
      enddo
      close(1)
      rfile = residpath//'/resid_sm_sg.txt'
      open (unit=1,file=rfile,status='old')
      do i=1,41
        do j=1,41
          read (1,*) cormat(i+1,j+1,8)
        enddo
      enddo
      close(1)
      rfile = residpath//'/resid_sm_sr.txt'
      open (unit=1,file=rfile,status='old')
      do i=1,41
        do j=1,41
          read (1,*) cormat(i+1,j+1,9)
        enddo
      enddo
      close(1)
      rfile = residpath//'/resid_sm_si.txt'
      open (unit=1,file=rfile,status='old')
      do i=1,41
        do j=1,41
          read (1,*) cormat(i+1,j+1,10)
        enddo
      enddo
      close(1)
      rfile = residpath//'/resid_sm_sz.txt'
      open (unit=1,file=rfile,status='old')
      do i=1,41
        do j=1,41
          read (1,*) cormat(i+1,j+1,11)
        enddo
      enddo
      close(1)
      rfile = residpath//'/resid_sm_zs.txt'
      open (unit=1,file=rfile,status='old')
      do i=1,41
        do j=1,41
          read (1,*) cormat(i+1,j+1,13)
        enddo
      enddo
      close(1)
      rfile = residpath//'/resid_sm_y.txt'
      open (unit=1,file=rfile,status='old')
      do i=1,41
        do j=1,41
          read (1,*) cormat(i+1,j+1,14)
        enddo
      enddo
      close(1)
      do i=1,41
        do k=1,14
          cormat(i,1,k) = cormat(i,2,k)
          cormat(i,43,k) = cormat(i,42,k)
        enddo
      enddo
      do j=1,41
        do k=1,14
          cormat(1,j,k) = cormat(2,j,k)
          cormat(43,j,k) = cormat(42,j,k)
        enddo
      enddo
      do i=1,43
        rowcol(i) = float(i)*100. - 150.
       enddo
c
c open input file
c
      file2 = eppfile//".epp"
      open(unit=2,file=file2,status='old')
      read (2,'(a)') fileid
      read (2,*)
c
c  get output file name
c
      file3 = eppfile//".pepp"
      open(unit=3,file=file3,status='new')
c
c
c  write output header id
c
      write (3,9001)
9001  format ('#INSTRUMENTAL MAGNITUDES, FILEPP 2.0'/
     $  '#',4x,'HJD',7x,'RA',11x,'DEC',7x,'FILT',1x,'AMASS',
     $   2x,'CCDX',5x,'CCDY',5x,'fwmx',3x,'fwmy',4x,'Peak',
     $   3x,'Sky',3x,'apr',2x,'apmag',3x,'aperr',
     $   3x,'psfmag',2x,'psferr',1x,'pflg',
     $   4x,'object',19x,'jd0',3x,'file',
     $   1x,'kset',1x,'kser',3x,'star',2x,'magcor')
c
c loop over records in file
c
100   continue
        read (2,9003,end=999) hjd,ra,dec,ifil,amass,
     $   ccdx,ccdy,fwhmx,fwhmy,amax,sky,iap,
     $   apmag,aperr,psfmag,psferr,photflag,
     $   object,jd0,nfile,kset,kser,istar
9003         format (f12.5,f12.7,f12.7,i5,f7.3,
     $         f9.3,f9.3,f7.3,f7.3,
     $         f8.1,f7.1,i5,16x,f8.4,f8.4,
     $         f8.4,f8.4,80x,i5,a25,i6,1x,
     $         a4,i5,i5,i7)
        do i=1,42
           if (ccdx.gt.rowcol(i).and.ccdx.lt.rowcol(i+1)) then
              i1 = i
              i2 = i+1
            endif
           if (ccdy.gt.rowcol(i).and.ccdy.lt.rowcol(i+1)) then
              j1 = i
              j2 = i+1
            endif
          enddo
c do bilinear interpolation from Numerical Recipies, page 96
        t = (ccdx-rowcol(i1))/(rowcol(i2)-rowcol(i1))
        u = (ccdy-rowcol(j1))/(rowcol(j2)-rowcol(j1))
        dmag = (1.-t)*(1.-u)*cormat(i1,j1,ifil) +
     $         t*(1.-u)*cormat(i2,j1,ifil) +
     $         t*u*cormat(i2,j2,ifil) +
     $         (1.-t)*u*cormat(i1,j2,ifil)
        xmag = apmag + dmag
        write (3,9004) hjd,ra,dec,ifil,amass,
     $   ccdx,ccdy,fwhmx,fwhmy,amax,sky,iap,
     $   xmag,aperr,psfmag,psferr,photflag,
     $   object,jd0,nfile,kset,kser,istar,dmag
9004         format (f12.5,f12.7,f12.7,i5,f7.3,
     $         f9.3,f9.3,f7.3,f7.3,
     $         f8.1,f7.1,i5,f8.4,f8.4,
     $         f8.4,f8.4,i5,2x,a25,i6,1x,
     $         a4,i5,i5,i7,f8.4)
        goto 100
c
999   continue
      close (2)
      close (3)
      stop
      end
