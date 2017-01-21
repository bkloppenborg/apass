	subroutine rdnames
c
c  this routine reads the input/output file names used by previous programs
c  
	common /names/ raw, int, starlis
	character*50 raw,int,starlis
	open(unit=66,file='names.dat',status='old',err=100)
	read(66,150) raw,int,starlis
150	format(3(a50,/))
	close (66)
	return
100	raw=' '
	int=' '
	starlis=' '
	return
	end
