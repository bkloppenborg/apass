	subroutine lin (x,y,m,a,b)
c
c  This subroutine is call SMOOTH in appendix I
c
c  *** Linear least squares routine from Nielson ***
c  Simple solution with no weighting factors
c
c  Inputs:
c		x,y   data arrays
c		m     number of points in arrays
c  Outputs:
c		a,b   where y=ax+b
c
c  *** written by A. Henden 1973 ***
c
	dimension x(m),y(m)
	double precision a1,a2,a3,c1,c2,det
c  initialize summing parameters
	a2=0.
	a3=0.
	c1=0.
	c2=0.
	a1=m
c  loop to set up matrix coefficients
	do 10 i=1,m
	a2=a2+x(i)
	a3=a3+x(i)*x(i)
	c1=c1+y(i)
	c2=c2+y(i)*x(i)
10	continue
c  solve matrix - simple since only 2x2
	det=1./(a1*a3-a2*a2)
	b=-(a2*c2-c1*a3)*det
	a=(a1*c2-c1*a2)*det
	return
	end
