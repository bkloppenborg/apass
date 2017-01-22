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
c
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
