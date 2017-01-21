	subroutine prompt(change,flag)
	character change
	flag=0
	read (5,'(a1)') change
	if((change .eq. 'n').or.(change .eq. 'N')) then
		flag = 1
		change = 'N'
	endif
	if((change .eq. 'y').or.(change .eq. 'Y')) then
		flag = 1
		change = 'Y'
	endif
	if(flag .eq. 0) then
		call bell
		write(6,10)
10		format(4x,'Invalid response --- try again.'/)
	endif
	return
	end
