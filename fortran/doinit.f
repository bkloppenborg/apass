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
     $    sky,fwhmx(MAXOBJ),fwhmy(MAXOBJ),sigsky
     $    fwhm,lobad,hibad,wscale,thresh,sharp(MAXOBJ),
     $    xfwhm(MAXOBJ),xmag(MAXOBJ),xerr(MAXOBJ)
      INTEGER
     $   ixdir,iydir,i,nc,nap,iap(MAXAP),iskyinr,iskyoutr,ich,
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
