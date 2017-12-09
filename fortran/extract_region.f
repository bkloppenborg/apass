      program extract_region
c
c extract ra/dec region from "sum" file
c written 09-Dec-2017 aah
c
      INTEGER MAXFILT
      PARAMETER (MAXFILT = 6)
      real*8 rax,decx,ralo,rahi,declo,dechi
      real*4 ccdx,ccdy,avexx,raerr,decerr
      real*4 xmag(MAXFILT),xerr(MAXFILT)
      integer nr,mr
      character file1*80,file2*80,file3*80,namex*25
      character ch1*139
c
      print *,'extract_region ver 1.0'
      print *,'note: does not cross 0/360 boundary'
      print *,'Name of sum file: '
      read (5,'(a)') file1
      print *,'Name of output file: '
      read (5,'(a)') file2
      print *,'Lower RA limit (degrees): '
      read (5,*) ralo
      print *,'Upper RA limit (degrees): '
      read (5,*) rahi
      print *,'Lower Dec limit: '
      read (5,*) declo
      print *,'Upper Dec limit: '
      read (5,*) dechi 
c
      open (unit=1,file=file1,status='old')
      open (unit=2,file=file2,status='new')
c read header and copy to output file
      read (1,'(a)') ch1
      write(2,'(a)') ch1
c now loop over records in sum file
300   continue
        read (1,904,end=400) namex,rax,raerr,decx,decerr,nr,
     $    mr,(xmag(j),j=1,MAXFILT),(xerr(j),j=1,MAXFILT)
904     format (a10,f11.6,f7.3,f11.6,f7.3,2i5,12f7.3)
        if (rax.lt.ralo.or.rax.gt.rahi.or.
     $    decx.lt.declo.or.decx.gt.dechi) goto 300
        write (2,904) namex,rax,raerr,decx,decerr,nr,
     $    mr,(xmag(j),j=1,MAXFILT),(xerr(j),j=1,MAXFILT)
         goto 300
400   continue
      close(1)
      close(2)
      stop
      end
