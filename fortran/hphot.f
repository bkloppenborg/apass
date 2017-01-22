      PROGRAM hphot
c
c command line input:  s 130901.130
c where s = s or n
c 130901.130 is fits image to process
c  photometric extraction program for fits images
c  read fits headers and extract information
c  then read image into memory
c  determine sky median and rms
c  find all objects in frame
c  do aperture photometry on the objects
c  create .raw file
c
c  written 19-Jun-95 aah fileread
c  mod 21-dec-95 aah photometry
c mod 29-Jun-98 aah npix in marginal sum 21 -> 11
c mod 04-Aug-98 aah add tek1k scattered light fix
c mod 25-Nov-98 aah limit apmag to 99.9
c mod 30-Nov-98 aah tek1k linearity and scaling probs
c mod 21-Jan-2003 aah add astrocam, filelist
c mod 09-May-2003 aah improved input
c mod 07-Nov-2003 aah 60K objects
c mod 11-Dec-2007 aah objra, objdec, date-obs mods
c mod 01-Nov-2009 aah new formats
c mod 18-Sep-2014 aah logic to prevent photometry with blank images
c mod 02-Aug-2016 aah SEL file structure
c mod 09-Aug-2016 aah single command line file input
c
      PARAMETER (MAXXPTS=4100)
      PARAMETER ( MAXYPTS=4100)
      INTEGER
     $    i,j,first,last,nbase1,nbase2,curfile,srcflag,
     $    ierr,extflag,narg,lenchr,n
      CHARACTER
     $    fname1*80,base1*80,base2*80,fnum*4,fname2*80,
     $    fname3*80,chr*80,infile*80,astfile*80,outfile*80,
     $    hemisphere*80
      LOGICAL
     $    ok
c
c      print *,'HPHOT Version 3.1  10-Aug-2016'
c
c get filename from command line and do checks
c filename is image file name, format yymmdd.xxxx
c
      narg = iargc()
c      print *,'narg ',narg
      call getarg(1,hemisphere)
      call getarg(2,chr)
c     print *,'narg ',narg,chr
c
c create the three required file names
c
      lenchr = nblank(chr)
      infile = chr(1:lenchr)
      astfile = hemisphere(1:1)//chr(1:lenchr)//'.asu'
      n = lenchr - 5
      outfile = chr(1:n)//'a'//chr((n+1):lenchr)
c      print *,'files ',lenchr,infile,astfile,outfile
      ok = .true.
      inquire(file=infile,exist=ok)
      if (ok) then
        inquire (file=astfile,exist=ok)
        if (ok) then
          goto 10
        endif
      endif
      goto 800
10    continue
      srcflag = 2
      call doinit
c      write (6,9011) infile
9011  format ('  Processing file ',a60)
      call doread(infile,ierr)
      if (ierr.eq.1) goto 700
      if (ierr.eq.2) goto 800
      call dosky
c added ierr here to prevent photometry of blank image
      call dofind(astfile,ierr)
      if (ierr.eq.1) goto 700
      call doparm
      open (unit=2,file=outfile,status='new',err=800)
      call writehead(infile)
      call dophot
      close(2)
700   continue
      STOP
800   continue
      print *,'Problem with files: ',infile,astfile
      STOP
      END

      subroutine doread(infile,jerr)
c MODIFIED FOR KPNO
c
c...Load a FITS file
c  modified from one by Blaise Canzian
c    by A. Henden 19-June-1995
c  mod 12-Sept-1998 aah to remove bias columns
c  mod 30-Nov-1998 aah tek1k linearity fix + sc- fix
c  mod 06-Dec-1998 aah add ti800
c  mod 21-Jan-2003 aah add astrocam
c  mod 28-Nov-2005 aah SRO, mac osx
c  mod 21-Aug-2007 aah  to read nmsu headers
c
c  inputs:
c       infile (char*80)  input file name
c  outputs
c       twodimage (maxxpts,maxypts)  real*4  the data
c       nx  number of columns
c       ny  number of rows
c       jerr  integer flag;  =0 if ok, =1 if bias/flat; =2err
c	plus keywords
c
      integer MAXXPTS, MAXYPTS, MAXAP
      parameter( MAXXPTS=4100)
      parameter( MAXYPTS=4100)
      parameter( MAXAP=9 )
c
      real*8 rval, bzero, bscale
      real*4 amap(MAXXPTS*MAXYPTS), imbuf(MAXXPTS,MAXYPTS),
     $     ut,exptime,epoch,apzero,rdnoise,gain,wscale,buffr(720)
      integer ierr, intval, bitpix, bytpix, ipix, pixrec, items,
     $ npix, npixleft, iaxnum, iline, bch, ech, jerr, nbias,
     $  nlin,nscat,nap,iap(MAXAP),iskyinr,iskyoutr,ich
      integer*2 buffw(1440),buffx
      integer*4 mapdm(2), buffl(720), naxis, recnum, nx, ny
      character header*2880, keyw*8, aline*80, chval*1, infile*80,
     $      date*10,filter*6,object*80,ra*20,dec*20
      character*6 ff(5),nmsufilt(5)
      logical bitp, simple, theend
      integer*1 buffer(2880)
      equivalence( buffer, buffw, buffl, buffr )
      COMMON /apblk/ apzero,rdnoise,gain,iap,iskyinr,iskyoutr,nap,ich
      COMMON /imblk/ imbuf,nx,ny
      COMMON /keyword/ ut,exptime,epoch,date,filter,object,ra,dec
      COMMON /flags/ nbias,nscat,nlin
      COMMON /dum1/ amap
      DATA ff /'U     ','B     ','V     ','R     ','I     '/
      DATA nmsufilt /'U     ','B     ','V     ','R     ','I     '/
c
  1   format(a)
  2   format(t4,(a),$)
  3   format(t4,(a))
  4   format(t4,(a),i6)
c
      jerr = 0
      bitp = .false.
      simple = .false.
      theend = .false.
      bscale = 1.0d0
      bzero = 0.0d0
      epoch = 2000.0
      filter = 'V'
c
c...Open The FITS File
c
      open( unit=11, file=infile, form='unformatted', recl=2880,
     & access='direct', status='old', err=99 )  
c
      recnum = 1
      do while( .not. theend )
         read( unit=11, rec=recnum, iostat=ierr ) header
         if( ierr .ne. 0 ) then
            write(6,4) 'Read error #', ierr
            write(6,3) 'In header:  not correct format for FITS file.'
            jerr = 2
            go to 99                    
         end if
c
c Header format is 36 lines of 80 characters per line ('card').
c
         iline = 0
         do while( iline.le.35 .and. .not.theend )
            bch = 80*iline+1
            ech = 80*iline+80
            aline(1:80) = header(bch:ech)
c
c Extract keyword from positions 1-8 in the line.
c
            keyw = aline(1:8)
            if( keyw .eq. 'END' ) then
               theend = .true.
            else if( keyw .eq. 'SIMPLE' ) then
               read( aline(30:30), 1 ) chval
               if( chval .eq. 'T' ) simple = .true.
            else if( keyw .eq. 'NAXIS' ) then
               read( aline(11:30), * ) intval
               naxis = intval
               if( naxis .gt. 2 ) then
                  write(6,25) 'Number of axes NAXIS=', naxis,
     $             'too large.'
 25               format( t4, (a), i2, 1x, (a) )
                  jerr = 2
                  go to 99
               end if
            else if( keyw .eq. 'BITPIX' ) then
               read( aline(11:30), * ) intval
c
c Number of bits per pixel.  This program handles 16 or 32.
c
               bitpix = intval
               bytpix = abs(bitpix/8)
               bitp = .true.
            else if( keyw(1:5) .eq. 'NAXIS' ) then
               read( keyw(6:8), * ) iaxnum
               read( aline(11:30), * ) intval
               mapdm(iaxnum) = intval
            else if( keyw .eq. 'BSCALE' ) then
               read( aline(11:30), * ) rval
               bscale = rval
            else if( keyw .eq. 'BZERO' ) then
               read( aline(11:30), * ) rval
               bzero = rval
c ich=50 is nmsu 1.0m
            else if (keyw .eq. 'UT') then
              if (ich.eq.90.or.ich.eq.5.or.ich.eq.8.or.ich.eq.16
     $           .or.ich.eq.50) then
                read(aline(11:30),900) iuth,iutm,uts
900             format (1x,i2,1x,i2,1x,f2.0)
              elseif (ich.eq.40) then
                read(aline(11:30),910) iuth,iutm,uts
910             format (1x,i2,1x,i2,1x,f2.0)
              elseif (ich.eq.4) then
                read(aline(11:30),921) iuth,iutm,uts
921             format (1x,i2,1x,i2,1x,f5.2)
              else
                read(aline(11:30),920) iuth,iutm,uts
920             format (2x,i2,1x,i2,1x,f5.2)
              endif
              ut = float(iuth) + float(iutm)/60. + uts/3600.
            else if (keyw .eq. 'UTSHUT') then
              read(aline(11:30),930) iuth,iutm,uts
930           format (1x,i1,1x,i2,1x,f4.1)
c901           format (1x,i2,1x,i2,1x,f4.1)
              ut = float(iuth) + float(iutm)/60. + uts/3600.
            else if (keyw .eq. 'RA') then
              call getstr(aline(11:30),20,ra,20)
c             ra = ra(1:(nblank(ra(2:20))+1))//'.00'
            else if (keyw .eq. 'OBJCTRA') then
              call getstr(aline(11:30),20,ra,20)
              ra(3:3) = ':'
              ra(6:6) = ':'
              ra(9:11) = '.00'
            else if (keyw .eq. 'OBJCTDEC') then
              call getstr(aline(11:30),20,dec,20)
              dec(4:4) = ':'
              dec(7:7) = ':'
              dec(10:12) = '.00'
            else if (keyw .eq. 'DEC') then
              call getstr(aline(11:30),20,dec,20)
c             dec = dec(1:(nblank(dec(2:20))+1))//'.00'
            else if (keyw .eq. 'FILTER') then
               if (ich.eq.50) then
                 read( aline(11:30), * ) intval
                 filter = nmsufilt(intval)
               else
                 call getstr(aline(11:30),20,filter,6)
               endif
            else if (keyw .eq. 'FILTERS') then
              call getstr(aline(11:30),20,filter,6)
c             read (aline(12:13),*) i
c             filter = ff(i)
c             filter = aline(19:19)//'     '
            else if (keyw .eq. 'FILTER_1') then
              call getstr(aline(11:30),20,filter,6)
            else if (keyw .eq. 'OBJECT') then
              if (aline(11:15).eq.'    ') then
                 object = 'NGC6811'
              else
                call getstr(aline(11:30),20,object,80)
              endif
              if (object(1:4).eq.'FLAT'.or.object(1:4).eq.
     $           'BIAS') then
                jerr = 1
                return
              endif
              call makename(object)
            else if (keyw .eq. 'FEXPTIME' .or.
     $         keyw.eq.'EXPTIME'.or.keyw.eq.'EXPOSURE') then
              read (aline(11:30),*) exptime
            else if (keyw .eq. 'EPOCH   ') then
              read (aline(11:30),*) epoch
c             epoch = 1950.0
            else if (keyw .eq. 'DATE-OBS') then
              call getstr(aline(11:30),20,date,10)
              if (aline(22:22).eq.'T') then
                read(aline(22:30),910) iuth,iutm,uts
              ut = float(iuth) + float(iutm)/60. + uts/3600.
              endif
            else if (keyw .eq. 'DATE_OBS') then
              call getstr(aline(11:30),20,date,10)
            end if
            iline = iline + 1
         end do
         recnum = recnum + 1
      end do
c
      nx = mapdm(1)
      ny = mapdm(2)
c
      if( .not. simple ) then
         write(6,3) 'Not simple FITS format.'
         jerr = 2
         go to 99
      end if
c
c Compute number of pixels in map.
c
      npix = 1
      do i=1,naxis
         npix = npix*mapdm(i)
      end do
      npixleft = npix
      if( npix .gt. MAXXPTS*MAXYPTS ) then
         write(6,3) 'Map too big.'
         jerr = 2
         go to 99
      end if
c
      if( .not. bitp ) then
         write(6,3) 'Number of bits per pixel missing from header.'
         jerr = 2
         go to 99
      end if
      ipix = 0
      pixrec = 2880/bytpix
c
c Read until there are no more pixels left.  Note that the pixels may not fill
c the entire record due to incommensurability of lengths.
c
      do while( npixleft .gt. 0 )
         read( unit=11, rec=recnum, iostat=ierr ) buffer
         if( ierr .ne. 0 ) then
            write(6,4) 'Read error #', ierr
            write(6,3) 'In buffer:  wrong format for FITS file.'
            jerr = 2
            go to 99
         end if
c
c Byte swapping for VMS
c
         if( bitpix .eq. 16 ) then
            call fits_bswap(buffer,2880)
         else if( bitpix .eq. 32 ) then
            call fits_lbswap(buffer,2880)
         else if( bitpix .eq.-32 ) then
            call fits_lbswap(buffer,2880)
         else
            write(6,3) 'Bits per pixel not 16 or 32.'
            jerr = 2
            goto 99
         end if
c
c Number of pixels containing data in this line.
c
         items = min( npixleft, pixrec )
         if( bitpix .eq. 16 ) then
            do i=1,items
               ipix = ipix + 1
               buffx = 256*buffer(i*2) + buffer(i*2-1)
               amap(ipix) = sngl(bscale*buffw(i) + bzero)
c              amap(ipix) = sngl(bscale*buffx + bzero)
            end do
         else if( bitpix .eq. 32 ) then
            do i=1,items
               ipix = ipix + 1
               amap(ipix) = sngl(bscale*buffl(i) + bzero)
            end do
         else if( bitpix .eq.-32 ) then
            do i=1,items
               ipix = ipix + 1
               amap(ipix) = bscale*buffr(i) + bzero
            end do
         else
            write(6,3) 'Bits per pixel not 16 or 32.'
c
c It is possible to handle 8 bits per pixel, but in that case there is concern
c for how UNIX interprets BYTE variables versus 8-bit integers.  On the VAX,
c BYTE variables are SIGNED.
c
            jerr = 2
            go to 99
         end if
         npixleft = npixleft - pixrec
         recnum = recnum + 1
      end do
      nbias=0
      do j=1,mapdm(2)
        do i=1,mapdm(1)
            imbuf(i,j) = amap(i + (j-1)*mapdm(1))
         end do
      end do
      jerr = 0
c
 99   close(11)
      return
      end

c
c Subroutines to byte swap.  Modified from code by Jeff Pier.
c
      subroutine fits_bswap (inpt,num_bytes)
c
      integer num_bytes,indx1,indx2
      integer*1 inpt(*), tempo
c
c swap pairs of bytes:
c
      do i=1,num_bytes/2
          indx1 = 2*i-1
          indx2 = 2*i
          tempo = inpt(indx1)
          inpt(indx1) = inpt(indx2)
          inpt(indx2) = tempo
      enddo
c
      return
      end
c
      subroutine fits_lbswap (inpt,num_bytes)
c
        integer num_bytes,indx1,indx2,indx3,indx4
      integer*1 inpt(*), tempo
c
c swap quartets of bytes:
c
      do i=1,num_bytes/4
          indx1 = 4*(i-1)+1
          indx2 = 4*(i-1)+2
          indx3 = 4*(i-1)+3
          indx4 = 4*(i-1)+4
          tempo = inpt(indx4)
          inpt(indx4) = inpt(indx1)
          inpt(indx1) = tempo
          tempo = inpt(indx2)
          inpt(indx2) = inpt(indx3)
          inpt(indx3) = tempo
      enddo
c
      return
      end

        SUBROUTINE getstr (str1,n1,str2,n2)
c
c copy str1 into str2 delimited by ticks
c
      integer n1,n2,i,j,k
      character*(*) str1,str2
c
      do i=1,n1
        if (str1(i:i).eq.'''') then
          do j=i+1,n1
            if (str1(j:j).eq.'''') goto 100
            k=j-i
            if (k.le.n2) then
              str2(k:k) = str1(j:j)
            endif
          enddo
         endif
      enddo
100   continue
      if (k.lt.n2) then
       do i=k+1,n2
        str2(i:i) = ' '
       enddo
      endif
      return
      end

      subroutine makename (object)
      integer i,j
      character object*80
      do i=1,20
        j = i + 1
        if (object(i:i).eq.' '.and.object(j:j).ne.' ') then
           object(i:i) = "_"
        endif
        if (object(i:i).ge.'a'.and.object(i:i).le.'z') then
          j = ichar(object(i:i)) - 32
          object(i:i) = char(j)
        endif
      enddo
      return
      end

      SUBROUTINE dofind(file3,ierr)
c
c read ast file, create starlist for aperture photometry
c written 02-Aug-2016 aah
c
c outputs:
c    coords(x+y,n)  coordinates of located objects
c
      INTEGER MAXOBJ
      PARAMETER (MAXOBJ=600000)
      real*8 xra(MAXOBJ),xdec(MAXOBJ)
      real*4 xmag(MAXOBJ),xerr(MAXOBJ)
      real*4 coords(MAXOBJ,2),xfwhm(MAXOBJ)
      integer nc,n
      character file3*80,cra*13,cdec*13,chr*230
      common /findblk/coords,xra,xdec,xfwhm,xmag,xerr,nc,cra,cdec
c
      ierr = 0
c
c read header
c
      open (unit=33,file=file3,status='old')
      cra = ' 00:00:00.000'
      cdec ='-00:00:00.000'
10    continue
        read(33,'(a)',end=30) chr
        if (chr(1:6).eq.'#RANEW') then
          read (chr,900) cra
900       format (9x,a13)
        elseif (chr(1:7).eq.'#DECNEW') then
          read (chr,901) cdec
901       format (10x,a13)
        elseif (chr(1:3).eq.'#1N') then
          goto 30
        endif
        goto 10
30    continue
c
c read starlist
c
      n = 1
40    continue
        read (33,902,end=50) coords(n,1),coords(n,2),xmag(n),xerr(n),
     $    xfwhm(n),xra(n),xdec(n)
902     format (f14.3,f11.3,7x,f10.4,f10.4,f7.2,20x,f14.7,f13.7)
        if (xfwhm(n).gt.99.999) xfwhm(n) = 99.99
        if (xfwhm(n).lt.0) xfwhm(n) = 99.99
        n = n+1
        goto 40
50    continue
      nc = n-1
c      print *,'nc ',nc
      close(33)
      return
      end

      SUBROUTINE doparm
c
c calculate image parameters (centroid, fwhm)
c written 19 Dec 1997 aah
c mod 20K objects 02-Oct-99 aah
c
      INTEGER MAXXPTS, MAXYPTS, MAXOBJ, MAXAP
      PARAMETER (MAXXPTS=4100) 
      PARAMETER (MAXYPTS=4100)
      PARAMETER (MAXOBJ=600000)
      PARAMETER ( MAXAP=9)
      REAL*8 sumt,xra(MAXOBJ),xdec(MAXOBJ)
      REAL*4 fwhmave(MAXOBJ), imbuf(MAXXPTS,MAXYPTS),
     $  coords(MAXOBJ,2),sky,val(100),x,y,err,
     $  sigsky,apzero,rdnoise,gain
      REAL*4 xmag(MAXOBJ),xerr(MAXOBJ),xfwhm(MAXOBJ)
      REAL*4 fwhmx(MAXOBJ),fwhmy(MAXOBJ),sharp(MAXOBJ),fwhm
      INTEGER npix,ixmin,ixmax,iymin,iymax,npx,npy,nx,ny,nc,
     $  iap(MAXAP),iskyinr,iskyoutr,nap,ixc,iyc,i,j,k,n
      INTEGER ixdir,iydir,izdir,ich
      CHARACTER cra*13,cdec*13
      COMMON /imblk/ imbuf,nx,ny
      COMMON /findblk/ coords,xra,xdec,xfwhm,xmag,xerr,nc,cra,cdec
      COMMON /fwhmblk/ fwhmx,fwhmy,sharp,fwhm
      COMMON /apblk/ apzero,rdnoise,gain,iap,iskyinr,iskyoutr,nap,ich
      COMMON /skyblk/ thresh,sky,sigsky
c
c     npix = iap(nap)*2
      npix = iap(nap-1)
      n = 0
      DO i=1,nc
        fwhmx(i) = xfwhm(i)
        fwhmy(i) = xfwhm(i)
        fwhmave(i) = xfwhm(i)
      ENDDO
c
c now sort to find mean fwhm for image
c
      if (n.eq.0) n=1
      call sort1 (fwhmave,n)
      if (n.eq.1) then
        fwhm = fwhmave(1)
      else
        fwhm = fwhmave(n/2+1)
      endif
      RETURN
      END

      SUBROUTINE dophot
c
c perform concentric aperture photometry on image
c plagerized from Stetson
c 21-Dec-95 aah
c mod 04-Aug-98 aah add tek1k scattered light fix
c mod 02-Oct-99 aah 20K obj
c mod 02-Aug-2016 aah SEL/SExtractor format
c
      INTEGER MAXXPTS, MAXYPTS, MAXOBJ
      INTEGER MINSKY, MAXSKY, MAXAP
c note MAXXPTS/YPTS permit 4k4k detector maximum
      PARAMETER (MAXXPTS=4100)
      PARAMETER (MAXYPTS=4100)
c probably overkill; typical image has about 11K objects, but don't know about the plane
      PARAMETER (MAXOBJ=600000)
      PARAMETER (MAXSKY=1000 )
      PARAMETER (MINSKY=20)
      PARAMETER ( MAXAP=9)
      REAL*4
     $    imbuf(MAXXPTS,MAXYPTS),xc,yc,par(MAXAP),
     $    ut,exptime, coords(MAXOBJ,2),apzero,rdnoise,
     $    sky,lobad,hibad,gain,skysig,skymod,skyvar,
     $    wscale, thresh,sigsky,
     $    apsky(MAXSKY), apmax, err1, err2, err3, dum,
     $    skymn,skymed,skew,epoch,edge,
     $    ratio(1024,1024),xmag(MAXOBJ),xerr(MAXOBJ),
     $    xfwhm(MAXOBJ)
      REAL*4 fwhmx(MAXOBJ),fwhmy(MAXOBJ),sharp(MAXOBJ),fwhm
      REAL*8
     $    apmag(MAXAP), area(MAXAP), aperr(MAXAP),
     $    xra(MAXOBJ),xdec(MAXOBJ)
      INTEGER
     $    i, j, k, n, nx, ny, nc, ixdir, iydir,iap(MAXAP),nap,
     $    iskyinr,iskyoutr, izdir,ich, ix, iy, nbias, nlin, nscat,
     $    nn
      CHARACTER
     $    date*10,filter*6,object*80,ra*20,dec*20,cra*13,cdec*13
      COMMON /keyword/ ut,exptime,epoch,date,filter,object,ra,dec
      COMMON /findblk/ coords,xra,xdec,xfwhm,xmag,xerr,nc,cra,cdec
      COMMON /fwhmblk/ fwhmx,fwhmy,sharp,fwhm
      COMMON /imblk/ imbuf,nx,ny
      COMMON /apblk/ apzero,rdnoise,gain,iap,iskyinr,iskyoutr,nap,ich
      COMMON /skyblk/ thresh,sky,sigsky
      COMMON /initblk/ wscale,ixdir,iydir,izdir,lobad,hibad
c
      apmxsq = -1.
      do i=1,nap
        par(i) = iap(i)
        apmxsq = amax1(apmxsq,(par(i)+0.5)**2)
      enddo
c      print *,'par',apmxsq,par
c
c loop over number of objects
c
      if (nc.le.0) return
      nn = 0
c      print *,'nc ',nc,coords(3,1),coords(3,2),nx,ny,iskyinr,iskyoutr
      DO n=1,nc
        xc = coords(n,1)
        yc = coords(n,2)
        if ((xc.ge.1.0).and.(xc.le.float(nx))
     $    .and.(yc.ge.1.0).and.(yc.le.float(ny))) then
        rin = iskyinr
        rinsq = rin*rin
c set outer radius to either input iskyoutr or else the radius
c that includes MAXSKY pixels, whichever is smaller
        routsq = float(MAXSKY)/3.142 + rinsq
        dum = float(iskyoutr)**2
        if (dum.lt.routsq) routsq = dum
        rout = sqrt(dum)
c limits of submatrix to be searched for this object
        lx = max0(1,int(xc-rout)+1)
        mx = min0(nx,int(xc+rout))
        ly = max0(1,int(yc-rout)+1)
        my = min0(ny,int(yc+rout))
        edge =amin1(xc-0.5,(nx+0.5)-xc,yc-0.5,(ny+0.5)-yc)
        apmax = -1.1E38
c         if (n.eq.3) print *,'lx ',lx,mx,ly,my,edge,rinsq,routsq,rout
c if star aperture extends outside array, its magnitude is bad
        DO j=1,nap
          apmag(j) = 0.0
          if (edge.lt.par(j)) apmag(j) = -1.1D38
          area(j) = 0.0
        ENDDO
        nsky = 0
c
c loop through submatrix for this object
c
        DO j=ly,my
          dysq =(j-yc)**2
            DO k=lx,mx
              rsq = dysq + (k-xc)**2
              datum = imbuf(k,j)
c is pixel in sky annulus?
              if ((rsq.lt.rinsq).or.(rsq.gt.routsq).or.
     $           (nsky.gt.maxsky).or.(datum.lt.lobad).or.
     $           (datum.gt.hibad)) goto 100
                 nsky = nsky+1
                 apsky(nsky) = datum
                 goto 200
100           continue
              if (rsq.gt.apmxsq) goto 200
c pixel is in at least one of the apertures.
              r = sqrt(rsq)-0.5
              DO m=1,nap
c if pixel completely outside aperture, skip
                if (r.gt.par(m)) goto 150
c determine fraction of pixel inside aperture
                fractn=amax1(0.0,amin1(1.0,par(m)-r))
                if ((datum.ge.lobad).and.(datum.le.hibad).and.
     $            (apmag(m).gt.-1.0D37)) then
                  apmag(m) = apmag(m)+dble(fractn*datum)
                  area(m) = area(m) + dble(fractn)
                  apmax = max(apmax,datum)
                else
                  apmag(m) = -1.1D38
                endif
150           continue
              ENDDO
200           continue
         ENDDO
      ENDDO
c determine sky for this object
      if (nsky.lt.MINSKY) then
c      if (n.eq.3) print *,' Too few sky pixels',nsky,MINSKY
        do i=1,nap
          apmag(i) = 99.999D0
          aperr(i) = 9.999
        enddo
        goto 500
      endif
c     call calcsky(apsky,nsky,skymod,skysig)
      call sort1(apsky,nsky)
      call mmm(apsky,nsky,hibad,rdnoise,skymn,skymed,skymod,
     $    skysig,skew)
c      if (n.eq.3) print *,'nsky ',nsky,MINSKY,xc,yc,apmag(1),skymod,
c     $   lobad,hibad
c kludge to handle error
      if (skysig.lt.0.0) skysig = 0.0
      skyvar = skysig**2
      sigsq = skyvar/float(nsky)
c subtract sky from each aperture
        DO i=1,nap
          apmag(i) = apmag(i) - dble(skymod)*area(i)
c if star+sky fainter than sky, or earlier badap occurred...
          if (apmag(i).le.0.0D0) goto 300
          err1 = area(i)*skyvar
          err2 = apmag(i)/gain
          err3 = sigsq*(area(i)**2)
          aperr(i) = amin1(9.999,1.0857*sqrt(err1+err2+err3)/
     $        sngl(apmag(i)))
c derive error and instrumental magnitude
            apmag(i) = apzero - 2.5D0*dlog10(apmag(i))
            if (apmag(i).gt.99.999D0) goto 300
            goto 400
300       continue
          apmag(i) = 99.999D0
          aperr(i) = 9.999
400       continue
        ENDDO
500     continue
        if (apmax.lt.-1.0E38.or.apmax.gt.99999.99) apmax=99999.99
c        if (n.eq.3) print *,'apmag ',apmag(1),fwhmx(1),fwhmy(1)
        if (apmag(1).lt.99.0) then
          if (fwhmx(n).lt.90.0.and.fwhmy(n).lt.90.0) then
            nn =min((nn + 1),99999)
            write (2,901) xra(n),xdec(n),xc,yc,fwhmx(n),fwhmy(n),
     $       sharp(n),apmax,
     $       skymod,(apmag(i),i=1,nap),xmag(n),
     $       (aperr(i),i=1,nap),xerr(n)
901         format (2f12.7,2f8.2,3f6.2,2f9.2,20f7.3)
c            write (2,900) (apmag(i),i=1,nap),xmag(n)
c            write (2,900) (aperr(i),i=1,nap),xerr(n)
c900         format (10f7.3)
          endif
        endif
        endif
      ENDDO
      RETURN
      END

      SUBROUTINE calcsky(sky,nsky,val,sig)
c
c determine sky from input array
c
      REAL*8 sum
      REAL*4 sky(*),val,sig
      INTEGER nsky
c
      call sort1 (sky,nsky)
      val = sky(nsky/2+1)
      sum=0.
      do i=1,nsky
        sum = sum + (sky(i)-val)**2
      enddo
      sig = dsqrt(sum)/dfloat(nsky-1)
      RETURN
      END

      SUBROUTINE SORT1 (x,n)
c
c heapsort of array x 
c from numerical recipies
c
      REAL*4  x(*)
      INTEGER n
      if (n.eq.1) return
      k = n/2+1
      ir = n
10    continue
      IF (k.gt.1) THEN
        k = k-1
        xx = x(k)
      ELSE
        xx = x(ir)
        x(ir) = x(1)
        ir = ir-1
        IF (ir.eq.1) THEN
          x(1) = xx
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
          i = j
          j = j+j
        ELSE
          j = ir+1
        ENDIF
      GOTO 20
      ENDIF
      x(i) = xx
      GOTO 10
      END

      SUBROUTINE doinit
c
c do all initialization stuff
c mod 04-Aug-98 aah add chip keyword
c mod 02-Oct-99 aah 20K obj
c
      INTEGER MAXOBJ, MAXAP
      PARAMETER (MAXOBJ=600000)
      PARAMETER ( MAXAP=9)
      REAL*8 xra(MAXOBJ),xdec(MAXOBJ)
      REAL*4
     $    coords(MAXOBJ,2),apzero,rdnoise,gain,
     $    sky,sigsky,
     $    lobad,hibad,wscale,thresh,
     $    xfwhm(MAXOBJ),xmag(MAXOBJ),xerr(MAXOBJ)
      REAL*4 fwhmx(MAXOBJ),fwhmy(MAXOBJ),sharp(MAXOBJ),fwhm
      INTEGER
     $   ixdir,iydir,izdir,i,nc,nap,iap(MAXAP),iskyinr,iskyoutr,ich,
     $   nbias,nscat,nlin
      CHARACTER txt*80,cra*13,cdec*13
      COMMON /apblk/ apzero,rdnoise,gain,iap,iskyinr,iskyoutr,nap,ich
      COMMON /findblk/ coords,xra,xdec,xfwhm,xmag,xerr,nc,cra,cdec
      COMMON /fwhmblk/ fwhmx,fwhmy,sharp,fwhm
      COMMON /skyblk/ thresh,sky,sigsky
      COMMON /flags/ nbias,nscat,nlin
      COMMON /initblk/ wscale,ixdir,iydir,izdir,lobad,hibad
c
      iap(1) = 2
      iap(2) = 3
      iap(3) = 4
      iap(4) = 5
      iap(5) = 7
      iap(6) = 9
      iap(7) = 11
      iap(8) = 13
      iap(9) = 15
      nap = 7
      nscat = 0
      iskyinr = 12
      iskyoutr = 21
      apzero = 17.0
      ich = 40
      ixdir = 1
      iydir = 0
      izdir = 1
      wscale = 2.5660
      rdnoise = 9.210
      gain = 2.6900
      lobad = -200.0
      hibad = 65000.
      thresh = 3.000
      fwhm = 3.000
      RETURN
      END

      SUBROUTINE dosky
c
c determine sky for this frame
c mod 04-Aug-98 aah fix tek1k scattered light prob
c mod 02-Oct-99 aah 20K obj
c
      INTEGER MAXOBJ, NPTS, MAXXPTS, MAXYPTS, MAXAP
      PARAMETER (MAXOBJ=600000)
      PARAMETER ( NPTS=200000)
      PARAMETER ( MAXXPTS=4100)
      PARAMETER ( MAXYPTS=4100)
      PARAMETER (MAXAP=9)
      REAL*8
     $  xra(MAXOBJ),xdec(MAXOBJ)
      REAL*4
     $  imbuf(MAXXPTS,MAXYPTS),
     $  thresh,rdnoise,skymn,skymed,
     $  coords(MAXOBJ,2),dx, pixel(NPTS),sky,sigsky,skew,
     $  apzero,gain,xfwhm(MAXOBJ),xmag(MAXOBJ),xerr(MAXOBJ)
      INTEGER
     $  nx,ny,nc, n, kk, jj,iskyinr,iskyoutr,nap,np,
     $  ich,iap(MAXAP)
      CHARACTER cra*13,cdec*13
      COMMON /apblk/ apzero,rdnoise,gain,iap,iskyinr,iskyoutr,nap,ich
      COMMON /findblk/ coords,xra,xdec,xfwhm,xmag,xerr,nc,cra,cdec
      COMMON /skyblk/ thresh,sky,sigsky
      COMMON /imblk/ imbuf,nx,ny
c
c use 200,000 random points for determining sky background
c
      np = npts
      dx = float(nx*ny)/float(NPTS+1)
      do i=1,NPTS
        n = int(dx*float(i))
        n = min(n,nx*ny)
        ix = n/nx
        ix = n - ix*nx + 1
        iy = (n-ix-1)/nx + 1
        pixel(i) = imbuf(ix,iy)
      enddo
      call sort1 (pixel,np)
      call mmm (pixel,np,hibad,rdnoise,skymn,skymed,sky,
     $  sigsky,skew)
c correct error
      if (sigsky.lt.0.0) sigsky = 0.0
      return
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

      SUBROUTINE SORT (n,x,index)
c
c heapsort of array x with corresponding index array index
c from numerical recipies
c
      DIMENSION x(n),index(n)
      k = n/2+1
      ir = n
10    continue
      IF (k.gt.1) THEN
        k = k-1
        xx = x(k)
        ii = index(k)
      ELSE
        xx = x(ir)
        ii = index(ir)
        x(ir) = x(1)
        index(ir) = index(1)
        ir = ir-1
        IF (ir.eq.1) THEN
          x(1) = xx
          index(1) = ii
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
          index(i) = index(j)
          i = j
          j = j+j
        ELSE
          j = ir+1
        ENDIF
      GOTO 20
      ENDIF
      x(i) = xx
      index(i) = ii
      GOTO 10
      END

      SUBROUTINE writehead(fname)
c MODIFIED FOR KPNO
c
c write file header
c mod 04-aug-98 aah fix tek1k scatttered light prob
c mod 09-sep-98 aah add version number on write
c mod 02-Oct-99 aah 20K obj
c
      INTEGER MAXOBJ,MAXAP,MAXXPTS,MAXYPTS
      PARAMETER (MAXOBJ=600000)
      PARAMETER ( MAXAP=9)
      PARAMETER (MAXXPTS=4100)
      PARAMETER ( MAXYPTS=4100)
      REAL*8 xra(MAXOBJ),xdec(MAXOBJ)
      REAL*4 
     $  ut,exptime,coords(MAXOBJ,2),apzero,rdnoise,gain,
     $  sky,thresh,sigsky,wscale,
     $  imbuf(MAXXPTS,MAXYPTS),epoch,
     $  xmag(MAXOBJ),xerr(MAXOBJ),
     $  lobad,h,ibad,xfwhm(MAXOBJ)
      REAL*4 fwhmx(MAXOBJ),fwhmy(MAXOBJ),sharp(MAXOBJ),fwhm
      INTEGER
     $  nc,nap,iap(MAXAP),iskyinr,iskyoutr,nx,ny,
     $  ich, nbias, nlin, nscat, ixdir, iydir, izdir
      CHARACTER fname*80,date*10,filter*6,object*80,ra*20,dec*20,
     $  cra*13,cdec*13
      COMMON /keyword/ ut,exptime,epoch,date,filter,object,ra,dec
      COMMON /findblk/ coords,xra,xdec,xfwhm,xmag,xerr,nc,cra,cdec
      COMMON /apblk/ apzero,rdnoise,gain,iap,iskyinr,iskyoutr,nap,ich
      COMMON /skyblk/ thresh,sky,sigsky
      COMMON /flags/ nbias,nscat,nlin
      COMMON /fwhmblk/ fwhmx,fwhmy,sharp,fwhm
      COMMON /imblk/ imbuf,nx,ny
      COMMON /initblk/ wscale,ixdir,iydir,izdir,lobad,hibad
c
      write (2,915)
915   format ('#HPHOT Version 1.0')
      write (2,900) fname
900   format ('#File= ',a60)
      write (2,914) nx,ny
914   format ('#nx,ny= ',2i5)
      write (2,901) ut
901   format ('#UT= ',f7.4)
      write (2,902) date
902   format ('#Date= ',a10)
c902   format ('#Date= ',a10)
      write (2,903) filter
903   format ('#Filter= ',a6)
      write (2,904) object
904   format ('#Object= ',a60)
c     write (2,905) ra(1:20)
      write (2,905) cra(1:12)
c905   format ('#RA= ',a19)
905   format ('#RA= ',a12,'       ')
c     write (2,906) dec
      write (2,906) cdec(1:13)
c906   format ('#DEC= ',a20)
906   format ('#DEC= ',a13,'       ')
      write (2,912) epoch
912   format ('#EPOCH=',f9.3)
      write (2,907) exptime
907   format ('#Exptime= ',f7.2)
      wscale = 2.67
      write (2,908) wscale
908   format ('#Scale= ',f7.4)
      write (2,909) ixdir,iydir,izdir
909   format ('#Orient= ',i2,i2,i2)
      write (2,910) ich
910   format ('#CHIP= ',i3)
      write (2,911) fwhm
911   format ('#FWHM= ',f7.2)
      write (2,913) sky,sigsky,thresh
913   format ('#sky= ',3f7.2)
      if (nbias.ne.0) write (2,923)
923   format ('#PHOT removed bias columns b4 proc')
      if (nscat.ne.0) write (2,924)
924   format ('#PHOT fixed scattered light problem')
      if (nlin.ne.0) write (2,925)
925   format ('#PHOT fixed linearity problem')
      if (ich.ge.1.and.ich.le.3) write (2,926)
926   format ('#PHOT fixed linearity and scaling problems')
      RETURN
      END

      SUBROUTINE  MMM (SKY, NSKY, HIBAD, READNS, SKYMN, SKYMED, 
     .     SKYMOD, SIGMA, SKEW)
C
C=======================================================================
C
C               Official DAO version:  1988 July 1
C
C This version of MMM (modified by PBS 1984.IV.10ff) assumes that
C the sky brightnesses in the one-dimensional array SKY are already
C sorted on entering this routine, and that pixels outside the "bad"
C limits have already been eliminated.
C
C This particular version of MMM also takes cognizance of the fact that,
C pixels falling below the LOBAD threshold already having been 
C eliminated, the contaminated sky pixels values overwhelmingly display
C POSITIVE departures from the true value.
C
C If for some reason it is impossible to obtain the mode of the sky
C distribution, this will be flagged by setting SIGMA = -1.0.
C
C Arguments
C
C     SKY (INPUT) is a real vector containing actual sorted sky values.
C    NSKY (INPUT) is the number of defined elements in SKY.
C  SKYMOD (OUTPUT) is the estimated mode of the sky values.
C   SIGMA (OUTPUT) is the computed standard deviation of the peak in
C         the sky histogram.
C    SKEW (OUTPUT) is the computed skewness of the peak in the sky
C         histogram.
C
C=======================================================================
C
      IMPLICIT NONE
      INTEGER NSKY
      REAL SKY(NSKY)
C
      DOUBLE PRECISION DSQRT, DBLE
      REAL ALOG10, AMIN1, AMAX1
C
      DOUBLE PRECISION SUM,SUMSQ
      REAL CUT, CUT1, CUT2, DELTA, SKYMID, SKYMED, SKYMN, SKYMOD
      REAL SIGMA, SKEW, R, SIGN, HIBAD, CENTER, SIDE, READNS
      REAL DMOD, OLD, CLAMP
      INTEGER I, J, K, L, M
      INTEGER MINIMM, MAXIMM, NITER, ISTEP, MAXIT, MINSKY, JSTEP
      LOGICAL REDO
      DATA MAXIT / 30 /, MINSKY / 20 /
C
C-----------------------------------------------------------------------
C
C SECTION 1
C
      IF (NSKY .LE. 0) THEN
         GO TO 9900
      END IF
      SKYMID=0.5*(SKY((NSKY+1)/2)+SKY(NSKY/2+1))
C
C SKYMID is the median value for the whole ensemble of sky pixels.
C Notice that if NSKY is odd, then (NSKY+1)/2 and (NSKY/2)+1 will be the
C same number, whereas if NSKY is even, they will be different numbers.
C This same trick will be used again later.
C
C Initialize the variables for accumulating the mean and standard
C deviation, and initialize the rejection limits.
C
      SUM=0.D0
      SUMSQ=0.D0
      CUT1=AMIN1(SKYMID-SKY(1), SKY(NSKY)-SKYMID, HIBAD-SKYMID)
C
C For the first pass we will consider only pixels in a symmetric 
C interval of brightness values about the median value.  This exploits
C the assumption that all the bad pixels are already rejected from the
C lower end of the brightness range.
C
      CUT2=SKYMID + CUT1
      CUT1=SKYMID - CUT1
C
      MINIMM=0
      DO 1010 I=1,NSKY
         IF (SKY(I) .LT. CUT1) THEN
            MINIMM=I
            GO TO 1010
         END IF
         IF (SKY(I) .GT. CUT2) GO TO 1020
         DELTA=SKY(I)-SKYMID
         SUM=SUM+DELTA
         SUMSQ=SUMSQ+DELTA**2
         MAXIMM=I
 1010 CONTINUE
C
C Henceforth in this subroutine, MINIMM will point to the highest value
C rejected at the lower end of the vector, and MAXIMM will point to the
C highest value accepted at the upper end of the vector.
C MAXIMM-MINIMM is the number of pixels within the acceptance range.
C
C Compute mean and sigma (from the first pass).
C
 1020 CONTINUE
      SKYMED=0.5*(SKY((MINIMM+MAXIMM+1)/2)+SKY((MINIMM+MAXIMM)/2+1))
      SKYMN=SUM/DBLE(MAXIMM-MINIMM)
      SIGMA=DSQRT(SUMSQ/DBLE(MAXIMM-MINIMM)-SKYMN**2)
      SKYMN=SKYMN+SKYMID
C
C The middle sky value, SKYMID, was subtracted off up above and added 
C back in down here to reduce the truncation error in the computation 
C of SIGMA.
C Note that this definition of SIGMA is incorrect by a factor of
C SQRT [NSKY/(NSKY-1.)], but for all but pathological cases (where none
C of this can be trusted anyway), it's close enough.
C
      SKYMOD=SKYMN
      IF (SKYMED .LT. SKYMN) SKYMOD=3.*SKYMED-2.*SKYMN
C
C If the mean is less than the mode, that means the contamination is
C slight, and the mean value is what we really want.  Note that this
C introduces a slight bias toward underestimating the sky when
C the scatter in the sky is caused by random fluctuations rather than
C by contamination, but I think this bias is negligible compared to the
C problem of contamination.
C
C-----------------------------------------------------------------------
C
C SECTION 2
C
C Rejection and recomputation loop:
C
      NITER=0
      OLD = 0.
      CLAMP = 1.
 2000 NITER=NITER+1
      IF ((NITER .GT. MAXIT) .OR. (MAXIMM-MINIMM .LT. MINSKY)) THEN
         GO TO 9900
      END IF
C
C Compute Chauvenet rejection criterion.
C         
      R=ALOG10(FLOAT(MAXIMM-MINIMM))
      R=AMAX1(2., (-.1042*R+1.1695)*R+.8895)
C
C Compute rejection limits (symmetric about the current mode).
C
      CUT=R*SIGMA+0.5*ABS(SKYMN-SKYMOD)
      CUT=AMAX1(1.5,CUT)
      CUT1=SKYMOD-CUT
      CUT2=SKYMOD+CUT
C
C Recompute mean and sigma by adding and/or subtracting sky values
C at both ends of the interval of acceptable values.
C
C At each end of the interval, ISTEP will show the direction we have to 
C step through the vector to go from the old partition to the new one.
C Pixels are added or subtracted depending upon whether the limit is 
C moving toward or away from the mode.
C
      REDO=.FALSE.
C
C Is CUT1 above or below the minimum currently-accepted value?
C
      ISTEP=INT(SIGN(1.0001, CUT1-SKY(MINIMM+1)))
      JSTEP=(ISTEP+1)/2
C
C If ISTEP = +1, JSTEP = 1.  If ISTEP = -1, JSTEP=0.  If ISTEP = +1, 
C then we know that at least one pixel must be deleted at the low end.
C
      IF (ISTEP .GT. 0) GO TO 2120
 2100 IF ((ISTEP .LT. 0) .AND. (MINIMM .LE. 0)) GO TO 2150
C
C Quit when SKY(MINIMM) < CUT1 <= SKY(MINIMM+1)
C
      IF ((SKY(MINIMM) .LE. CUT1) .AND. (SKY(MINIMM+1) .GE. CUT1))
     .     GO TO 2150
C
C If ISTEP is positive, subtract out the sky value at MINIMM+1; if 
C ISTEP is negative, add in the sky value at MINIMM.
C
 2120 CONTINUE
      DELTA=SKY(MINIMM+JSTEP)-SKYMID
      SUM=SUM-REAL(ISTEP)*DELTA
      SUMSQ=SUMSQ-REAL(ISTEP)*DELTA**2
      MINIMM=MINIMM+ISTEP
      REDO=.TRUE.                                 ! A change has occured
      GO TO 2100
C
 2150 CONTINUE
C
C Is CUT2 above or below the current maximum?
C
      ISTEP=INT(SIGN(1.0001, CUT2-SKY(MAXIMM)))
      JSTEP=(ISTEP+1)/2
C
C If ISTEP = +1, JSTEP = 1.  If ISTEP = -1, JSTEP=0.  If ISTEP = -1, 
C then we know that we must subtract at least one pixel from the high 
C end.
C
      IF (ISTEP .LT. 0) GO TO 2220
 2200 IF ((ISTEP .GT. 0) .AND. (MAXIMM .GE. NSKY)) GO TO 2250
C
C Quit when SKY(MAXIMM) <= CUT2 < SKY(MAXIMM+1)
C
      IF ((SKY(MAXIMM) .LE. CUT2) .AND. (SKY(MAXIMM+1) .GE. CUT2))
     .     GO TO 2250
C
C If ISTEP is positive, add in the sky value at MAXIMM+1; if ISTEP is 
C negative, subtract off the sky value at MAXIMM.
C
 2220 DELTA=SKY(MAXIMM+JSTEP)-SKYMID
      SUM=SUM+REAL(ISTEP)*DELTA
      SUMSQ=SUMSQ+REAL(ISTEP)*DELTA**2
      MAXIMM=MAXIMM+ISTEP
      REDO=.TRUE.                                 ! A change has occured
      GO TO 2200
C
 2250 CONTINUE
C
C Compute mean and sigma (from this pass).
C
      SKYMN=SUM/DBLE(MAXIMM-MINIMM)
      SIGMA=DSQRT(SUMSQ/DBLE(MAXIMM-MINIMM)-SKYMN**2)
      SKYMN=SKYMN+SKYMID
C
C Obtain the median.  To first approximation, the median would be the
C value of the sky in the middle pixel in the sorted data (if the
C total number is odd) or the mean of the two pixels straddling
C the middle (if the total number of pixels is even).
C
C     SKYMED=0.5*(SKY((MINIMM+MAXIMM+1)/2)+SKY((MINIMM+MAXIMM)/2+1))
C
C However, this is not good enough.  If you look at the estimator for
C the mode, you will note that a tiny change in the list of sky pixels,
C just sufficient to alter the median value of the sky brightness by
C one unit, will change the estimator of the mode by three units.  We
C really want something more robust than this.  As a first attempt
C at a more robust median estimator, I estimated the median
C of the distribution by the mean of the central twenty percent of sky
C values.  This involved considerable care to make sure you get
C a perfectly symmetric sample of pixels about the median, whether
C there is an even or an odd number of pixels within the acceptance
C interval.  However, even this is not good enough when you have 
C quantized intensities and detectors with low enough read noise that 
C a large fraction of the sky pixels have the same (quantized) intensity:  
C a given fixed percentage on either side of the median does not 
C necessarily give a good estimator of how the values differ to the two 
C sides. So I add another criterion:  if either limiting value is not
C sufficiently different from the median, expand the range by one pixel 
C on either side until we reach a significantly different intensity value.
C In a normal distribution, the central 20% corresponds to +/- 0.25 sigma,
C so expand the range until both limiting values differ from the central
C value by at least 0.25 times the read noise.
C
      SKYMED=0.0
      CENTER = REAL(MINIMM+1 + MAXIMM)/2.
      SIDE = REAL(NINT(0.2*REAL(MAXIMM-MINIMM)))/2. + 0.25
      J = NINT(CENTER-SIDE)
      K = NINT(CENTER+SIDE)
      L = NINT(CENTER-0.25)
      M = NINT(CENTER+0.25)
      R = 0.25*READNS
 2305 CONTINUE
      IF ((J .GT. 1) .AND. (K .LT. NSKY) .AND. (
     .     (SKY(L)-SKY(J) .LT. R) .OR. (SKY(K)-SKY(M) .LT. R) )) THEN
         J = J-1
         K = K+1
         GO TO 2305
      END IF
C
      DO 2310 I=J,K
 2310 SKYMED=SKYMED+SKY(I)
C
      SKYMED=SKYMED/REAL(K-J+1)
      IF (SKYMED .LT. SKYMN) THEN
         DMOD=3.*SKYMED-2.*SKYMN - SKYMOD
      ELSE
         DMOD=SKYMN - SKYMOD
      END IF
C
C If the mean is less than the mode, that means the contamination is
C slight, and the mean value is what we really want.  Note that this
C introduces a slight bias toward underestimating the sky when
C the scatter in the sky is caused by random fluctuations rather than
C by contamination, but I think this bias is negligible compared to the
C problem of contamination.
C
C If the limits have not yet stopped moving, try again.
C
      IF (DMOD*OLD .LT. 0.) CLAMP = 0.5*CLAMP
      SKYMOD = SKYMOD + CLAMP*DMOD
      OLD = DMOD
      IF (REDO) GO TO 2000
C
C-----------------------------------------------------------------------
C
C Normal return.
C
      SKEW=(SKYMN-SKYMOD)/AMAX1(1., SIGMA)
      NSKY=MAXIMM-MINIMM
      RETURN
C
C-----------------------------------------------------------------------
C
C An error condition has been detected.
C
 9900 SIGMA=-1.0
      SKEW=0.0
      RETURN
C
      END!
