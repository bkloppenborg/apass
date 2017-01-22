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
      REAL fwhmave(MAXOBJ), imbuf(MAXXPTS,MAXYPTS),
     $  coords(MAXOBJ,2),sky,lobad,hibad,val(100),x,y,err,
     $  fwhmx(MAXOBJ),fwhmy(MAXOBJ),fwhm,sigsky,sharp(MAXOBJ)
      REAL*4 xmag(MAXOBJ),xerr(MAXOBJ),xfwhm(MAXOBJ)
      INTEGER npix,ixmin,ixmax,iymin,iymax,npx,npy,nx,ny,nc,
     $  iap(MAXAP),iskyinr,iskyoutr,nap,ixc,iyc,i,j,k,n
      INTEGER ixdir,iydir,izdir,ich
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
