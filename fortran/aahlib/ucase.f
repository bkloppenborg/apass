	subroutine ucase(name,len)
c
c  this function converts a character string of length LEN to upper case
c   LEN can not exceed 80
c
c   Written by R. Kaitchuck
c   Fortran 77 compatible
c
c   External subroutines: none
c
c   Last revised 25-Apr-89
c
	character*80 name
	dimension ihold(80)
c
	do 10 i=1,len
	ihold(i)=ichar(name(i:i))
	if((ihold(i).ge.97).and.(ihold(i).le.122)) ihold(i)=ihold(i)-32
10	name(i:i)=char(ihold(i))
	return
	end
