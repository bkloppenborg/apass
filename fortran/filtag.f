      program filtag
c
c tag observations as photometric or nonphotometric
c written 11-Aug-2010 aah
c mod 09-Dec-2017 aah new file format
c
      INTEGER MAXST,MAXJD
      PARAMETER (MAXST = 100)
      PARAMETER (MAXJD = 100)
      real*8 crax,cdecx,cra(MAXST,MAXJD),cdec(MAXST,MAXJD)
      real*8 hjd
      real*4 ccdx,ccdy,avexx,xmag1,xerr1,dmag
      integer isys,ifil,flag1,flag2,igroup,i,j,iflag,njd
      integer nst(MAXST),iprop(MAXST,MAXJD),jd(MAXJD)
      integer istar(MAXST,MAXJD),kset(MAXST,MAXJD)
      character file1*80,file2*80,file3*80
      character ch1*139,ch2*139,ch3*139
      character*25 name(MAXST,MAXJD),namex
c
      print *,'Input list of fred files: '
      read (5,'(a)') file1
      print *,'Input list of fields to tag: '
      read (5,'(a)') file2
      print *,'Mark fields photometric(0) or nonphotometric(1): '
      read (5,*) itag
c
      open (unit=1,file=file2,status='old')
      read (1,900) cra(1,1),cdec(1,1),ccdx,ccdy,flag1,flag2,
     $    hjd,avexx,kset(1,1),igroup,name(1,1),ifil,
     $    xmag1,xerr1,dmag,isys,jd(1)
      njd = 1
      nst(1) = 1
100   continue
        read (1,900,end=200) crax,cdecx,ccdx,ccdy,flag1,flag2,
     $    hjd,avexx,ksetx,igroup,namex,ifil,
     $    xmag1,xerr1,dmag,isys,jdx
900     format(2f12.7,2f10.3,1x,i1.1,1x,i1.1,1x,f12.6,
     $  1x,f5.3,i5,1x,i10,1x,a25,1x,i5,1x,3f8.4,2i6)
        iflag = 0
        do i=1,njd
          if (jd(i).eq.jdx) then
             nst(i) = nst(i) + 1
             iflag = i
           endif
         enddo
         if (iflag.eq.0) then
c new star
             njd = njd + 1
             nst(njd) = 1
             iflag = njd
             jd(njd) = jdx
          endif
          j = nst(iflag)
          cra(j,iflag) = crax
          cdec(j,iflag) = cdecx
          name(j,iflag) = namex
          kset(j,iflag) = ksetx
        goto 100
200   continue
      close(1)
      print *,'njds to search: ',njd
      do i=1,njd
           do j=1,nst(i)
            print *,i,j,kset(j,i),jd(i)
           enddo
       enddo
c
c loop over fred files
c
      open (unit=1,file=file1,status='old')
300   continue
        read (1,'(a)',end=600) file3
        print *,'Processing file ',file3
        open (unit=3,file=file3,status='old')
        read (3,'(a)') ch1
        read (3,'(a)') ch2
        read (3,'(a)') ch3
        read (3,900) crax,cdecx,ccdx,ccdy,flag1,flag2,
     $    hjd,avexx,ksetx,igroup,namex,ifil,
     $    xmag1,xerr1,dmag,isys,jdx
        iflag=0
        do j=1,njd
          if (jd(j).eq.jdx) iflag = 1
        enddo
        if (iflag.eq.0) then
           close(3)
           goto 300
        endif
c need to correct this file
        backspace(3)
        file1 = file3(1:nblank(file3))//'x'
        open (unit=2,file=file1,status='new')
        write (2,'(a)') ch1
        write (2,'(a)') ch2
        write (2,'(a)') ch3
400     continue
        read (3,900,end=500) crax,cdecx,ccdx,ccdy,flag1,flag2,
     $    hjd,avexx,ksetx,igroup,namex,ifil,
     $    xmag1,xerr1,dmag,isys,jdx
        do j=1,njd
          if (jd(j).eq.jdx) then
             do i=1,nst(j)
c               if(kset(i,j).eq.ksetx.and.name(i,j).eq.
c     $           namex) then
                if (kset(i,j).eq.ksetx) then
                 flag1 = itag
                endif
              enddo
           endif
        enddo
        write (2,900) crax,cdecx,ccdx,ccdy,flag1,flag2,
     $    hjd,avexx,ksetx,igroup,namex,ifil,
     $    xmag1,xerr1,dmag,isys,jdx
        goto 400
500    continue
       close(2)
       close(3)
       goto 300
600   continue
      close(1)
      close(2)
      close(3)
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
