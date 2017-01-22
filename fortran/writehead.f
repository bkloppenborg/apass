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
     $  sky,thresh,sigsky,fwhmx(MAXOBJ),wscale
     $  fwhmy(MAXOBJ),fwhm,imbuf(MAXXPTS,MAXYPTS),epoch,
     $  sharp(MAXOBJ),xfwhm(MAXOBJ),xmag(MAXOBJ),xerr(MAXOBJ),
     $  lobad,hibad
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
