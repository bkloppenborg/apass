	subroutine solve2 (n,nfilt,x,smag,urmag,sext,fext,coef,zero)
c
c	a modified version of SOLVE by A. Henden 1980
c       last modified June 13, 1989 by Kaitchuck
c
c  	least squares solution of the UBVRI transformation equations
c   input -
c	smag(n,5)   standard magnitudes and colors in the form:
c		smag(n,1) = mag
c		smag(n,2) = color 1
c		smag(n,3) = color 2
c		smag(n,4) = color 3
c		smag(n,5) = color 4
c	urmag(n,5)  your corresponding instrumental magnitudes & colors
c	x(n,5)      air mass values
c	n           number of observations to reduce
c	nfilt       mag + number of colors = (1-5)
c	sext(5)     second order extinction
c   output -
c	fext(5)     first order extinction
c	coef(5)     transformation color coefficients
c	zero(5)     zero points
c
c	to have fext coef and zero to all be valid, you
c       must use standards with both a wide range of colors
c	and a wide range of airmases.  with a cluster, you
c	will get valid coef parameters, but fext and zero may
c	be invalid...so beware!
c
c
      PARAMETER (maxstars = 3000)
	dimension smag(maxstars,5),urmag(maxstars,5),x(maxstars,5),
     $  sext(5),fext(5),coef(5),zero(5)
c   loop over colors
	do 20 k=1,nfilt
c   zero matrix elements
	a1=0.
	a2=0.
	a3=0.
	a4=0.
	a5=0.
	a6=0.
	a7=0.
	a8=0.
	a9=0.
	a10=0.
	a11=0.
	a12=0.
c  loop over number of standard stars
	do 10 i=1,n
c  calculate matrix elements sums
	temp=urmag(i,k)*(1.-sext(k)*x(i,k))
	std=smag(i,k)
	if(k .ne. 1) go to 5
c  note: the mag equation has slightly different form
	temp=smag(i,2)
	std=std-urmag(i,k)
5	continue
c
c  the matrix looks like:
c	( a1  a2   a3   . a4  )
c       ( a5  a6   a7   . a8  )
c       ( a9  a10  a11  . a12 )
c  note: a2=a5, a3=a9, and a7=a10.  but explicit here for clarity
c
	a1=a1+temp*temp
	a2=a2+temp*x(i,k)
	a3=a3+temp
	a4=a4+std*temp
	a5=a5+temp*x(i,k)
	a6=a6+x(i,k)*x(i,k)
	a7=a7+x(i,k)
	a8=a8+std*x(i,k)
	a9=a9+temp
	a10=a10+x(i,k)
	a11=a11+1
	a12=a12+std
10	continue
c   calculate minors
	aa=a7*a10-a6*a11
	bb=a5*a11-a7*a9
	cc=a6*a9-a5*a10
	dd=a8*a11-a7*a12
	ee=a6*a12-a8*a10
	ff=a5*a12-a8*a9
c  calculate determinant
	det=a1*aa+a2*bb+a3*cc
c  solve for wanted values
	coef(k)=(a4*aa+a2*dd+a3*ee)/det
	zero(k)=(a2*ff-a1*ee+a4*cc)/det
	fext(k)=(a1*dd-a4*bb+a3*ff)/det
	if(k .ne. 1) fext(k)=fext(k)/coef(k)
20	continue
	return
	end
