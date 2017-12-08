	program FILEXTR
c
c This program extracts individual star data from files
c version for newphot 29-Jan-98 aah
c mod 20-Apr-98 aah improved match
c mod 14-June-1998 aah errrad=2arcsec
c mod 29-June-1998 aah remove points with flags=11111
c mod 03-Mar-2007 aah 6 or 5 colors
c mod 11-Aug-2010 aah filcon format
c mod 08-Dec-2017 aah 25-char names
c
c space = 94*MAXSTAR Bytes
c
c ***********************************************************************
c
      PARAMETER (MAXSTAR=3000)
      PARAMETER (MAXFILT=6)
      REAL*8 ra(MAXSTAR),dec(MAXSTAR),
     $  rax,decx,dist,err,hjd,cosd
      REAL*4 
     $  xmag,xerr,dmag,ccdx,ccdy,
     $  raerr,decerr
      INTEGER npts,i,j,flags(11),kset,iflag,jflag,names(MAXSTAR)
      INTEGER pflag,sflag,kgroup,filt,isys,night,oldname
      INTEGER field1,field2
      CHARACTER file1*60,file2*60,file3*60,star*25,ch*90
      LOGICAL fflag
c
c ***********************************************************************
c
	print *,'    Program FILEXTR version 1.1   08-Dec-2017'
	print *
	print *
c
c ************************************************************************
c
      picon = datan(1.D0)/45.0
      print *,'Enter file list: '
      read (5,'(a)') file1
      open (unit=1,file=file1,status='old')
      print *,'Enter list of star names and positions: '
      read (5,'(a)') file2
      open (unit=2,file=file2,status='old')
      nst = 1
10    continue
        read (2,9000,end=20) names(nst),ra(nst),dec(nst)
9000    format (i10,6x,2f13.8)
        nst = min((nst+1),1000)
        goto 10
20    continue
      nst = nst - 1
      close(2)
      print *,'Number of stars to search: ',nst
      print *,'Enter output file name: '
      read (5,'(a)') file3
      open (unit=3,file=file3,status='new')
      write (3,908)
908   format ('#',2x,'Name',19x,'RA(J2000)',3x,'raerr',2x,
     $  'DEC(J2000)',2x,'decerr',3x,'Flags',3x,'V',
     $  5x,'B-V',5x,'U-B',5x,'V-R',5x,'R-I',5x,'Errors')
c
c loop over number of files
c
      err = (2.5/3600.)**2  ! position error = 2.5arcsec
      nm1 = MAXFILT-1
      npts = 0
100   continue
      read (1,'(a)',end=400) file2
      write (6,905) file2
905   format (' Processing file: ',a50)
      open (unit=2,file=file2,status='old')
c skip over header
      read (2,*)
      read (2,*)
      read (2,*)
      xmag = 99.999
      xerr = 9.999
c loop over stars
200   continue
        read (2,901,end=300) rax,decx,ccdx,ccdy,pflag,sflag,
     $    hjd,avexx,kset,kgroup,star,filt,xmag,xerr,dmag,
     $    isys,night
901   format (2f12.7,2f10.3,2i2,f13.6,f6.3,i5,i11,1x,a25,i6,
     $    1x,3f8.4,i6,i6)
        do j=1,nst
          cosd = dcos(dec(j)*picon)
          dist = ((ra(j)-rax)*cosd)**2 + (dec(j)-decx)**2
          if (dist.lt.err) then
            write (3,901) rax,decx,ccdx,ccdy,pflag,sflag,
     $       hjd,avexx,kset,kgroup,field1,field2,filt,xmag,xerr,dmag,
     $       isys,night
          endif
        enddo
        goto 200
300   continue
      goto 100
400   continue
      close (3)
      stop
      end
