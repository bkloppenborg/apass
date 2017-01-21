	subroutine wrtnames
c
c  this routine writes the input/output file names used by the next program
c  
	common /names/ raw, int, starlis
	character*50 raw,int,starlis
	open(unit=66,file='names.dat',status='old',err=100)
	write(66,150) raw,int,starlis
150	format(3(a50,/))
	close (66)
	return
100	open(unit=66,file='names.dat',status='new',err=200)
	write(66,150) raw,int,starlis
	close(66)
200	return
	end
