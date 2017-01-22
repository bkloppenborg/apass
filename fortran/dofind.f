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
      close(33)
      return
      end
