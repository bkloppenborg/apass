       subroutine disply
c
c  Common block parameters (a common is used to prevent the system stack
c   from exceeding 64KB.)
c
c       x        the x data array              
c       y        the y data array              
c       npts      number of data points         
c       xbeg      x starting position in inches 
c       ybeg      y starting position in inches 
c       xlen      length of x axis in inches    
c       ylen      length of y axis in inches    
c       xlab      x axis label                  
c       ylab      y axis label                  
c       icolor1   color index to use for data
c       connect   logical true to connect data points
c       slope     slope of fit
c       b         intercept of fit
c       icolor2   color index of fit
c       mflag     logical variable indicating if this is a magnitude plot
c
c
      common /plotit/ x,y,npts,xbeg,ybeg,xlen,ylen,xlab,ylab,
     $  icolor1,connect,slope,b,icolor2,mflag
      real*4  x(3000),y(3000),rarray(7)
      real*4 xbeg,ybeg,xlen,ylen,slope,b
      real*4 xx(3000),yy(3000),xmin,xmax,ymin,ymax,x1(2),y1(2)
      integer icolor1,icolor2,npts
      integer*2 iarray(9)
      character*30 xlab,ylab
      character*10 str
      character*40 string
      logical connect,mflag
      data iflag /0/
c
c  init the graphics package, enter graphics mode
c
      IF (iflag.eq.0) THEN
        i = sm_device('X11')
        iflag = 1
        call sm_graphics
      ENDIF
      call sm_erase
      xmin = 5.e20
      ymin = 5.e20
      xmax = -5.e20
      ymax = -5.e20
      DO i=1,npts
        xx(i) = x(i)
        yy(i) = y(i)
        xmin = min(xmin,xx(i))
        xmax = max(xmax,xx(i))
        ymin = min(ymin,yy(i))
        ymax = max(ymax,yy(i))
      ENDDO
      call sm_limits(xmin,xmax,ymin,ymax)
      call sm_ctype ('white')
      call sm_box (1,2,0,0)
      call sm_xlabel(xlab)
      call sm_ylabel(ylab)
      call sm_gflush
      call sm_ctype ('green')
      IF (connect.eqv.(.true.)) THEN
        call sm_conn(xx,yy,npts)
        call sm_conn(xx,yy,npts)
      ELSE
        call sm_ptype(2.43e2,1)
        call sm_expand(1.2e0)
        call sm_points(xx,yy,npts)
        call sm_expand(1.0e0)
        call sm_ptype(1.1e1,1)
      ENDIF
      call sm_gflush
      x1(1) = xmin
      y1(1) = b + slope*x1(1)
      x1(2) = xmax
      y1(2) = b + slope*x1(2)
      call sm_ctype('green')
      call sm_conn(x1,y1,2)
      call sm_conn(x1,y1,2)
      call sm_gflush
      call sm_alpha
c     write (6,900)
900   format (' Press any key to continue: ')
c      call sm_redraw(0)
c      read (5,*) i
      return
      end
