	subroutine extrd(file,numfil,filt,fext,erflg)
c
c  This routine reads an extinction coefficient file and returns the
c   coefficients.
c
c	input:  file   -- file name
c		filt   -- filter array
c		numfil -- number of elements in filter array (5 max)
c	output:	fext   -- ext. coeff
c		erflg  -- true is read successful
c
	character c1,c2,c3,text(80),ok
	character*3 filt(5),xfilt(5)
	character*30 file,id
	dimension fext(5),xext(5)
	logical erflg
	erflg=.true.
c
c  open file
c
8	open(unit=97,file=file,status='old',err=10)
	go to 15
10	call bell
	close (97)
	print *,'  Error opening file, try again, file name ... '
	read *, file
	go to 8
c
c  Check file type
c
15	read(97,'(a30)') id
	if(id .ne. 'FIRST ORDER EXTINCTION') then
		call bell
	        close (97)
		print *,'  Wrong file type, try again, file name ... '
		read *,file
		go to 8
	endif
c
c  Read the filter designations from the file
c
	read (97,1000) text
1000	format(80a1)
	kl=1
	j=1
2500	if(text(j) .eq. ' ') go to 2000
		c1=text(j)
		c2=text(j+1)
		c3=text(j+2)
		xfilt(kl)=c1//c2//c3
		j=j+2
		kl=kl+1
2000	j=j+1
	if(j .le. 80) go to 2500
	kfilt=kl-1
c
c  kfilt is the mag + number of colors in the extinction file 
c
c  Read in the extinction coefficients
c
	read(97,*) (xext(m),m=1,kfilt)
	close (97)
c
c  See if the file contains any values for colors
c
	nc=0
	do i=1,5
		do j=1,3
	   	if(xfilt(i)(j:j) .eq. '-') nc=nc+1
		enddo
	enddo
c
c  Does the filter designations passed to this sub contain any colors?
c
	kc=0
	do i=1,5
		do j=1,3
		if(filt(i)(j:j) .eq. '-') kc=kc+1
		enddo
	enddo
c
c  Make the the file contains enough colors
c
	if((kc .gt. nc).and.(nc .ne. 0)) then
		call bell
		print *,'  Not enough colors in the extinction file'
		erflg=.false.
		return
	endif
c
c  See if colors need to be built up from filters in file
c
	if((kc .gt. 0).and.(nc .eq. 0)) go to 800
c
c  Match filter designations in file with those passed to sub
c
	i=1
30	ifnd=0
	do j=1,kfilt
		if(filt(i).eq.xfilt(j)) then
			fext(i)=xext(j)
			ifnd=1
		endif
	enddo
	if(ifnd .eq. 0) then
	   call bell
	   print *,'  The extinction coefficient for ',filt(i),
     $             ' was not found'
32	   print *,'  Please enter value --- '
	   read(6,*,err=32) fext(i)
	endif
	i=i+1
	if(i .gt. numfil) go to 852
	go to 30
c
c  Build the required color coefficients
c      First find the magnitude filter (its the only single char filter)
c
800	i=1
801	do j=1,kfilt
		if(filt(i).eq.xfilt(j)) then
			fext(i)=xext(j)
			imag=i
			go to 810
		endif
	enddo
	i=i+1
	if(i .le. numfil) go to 801
   	call bell
	Print *,'  The magnitude filter was not found in the file'
	erflg=.false.
	return
c
c
810	do i=1,numfil
	if(i .eq. imag) go to 825			! skip mag
		do j=1,kfilt
		if(filt(i)(1:1) .eq. xfilt(j)) then	      ! match 1st char
			do jj=1,kfilt
			if(filt(i)(3:3) .eq. xfilt(jj)) then      ! match 3rd
			fext(i)=xext(j)-xext(jj)                  ! calc coeff
			endif
			enddo
		endif
		enddo
825	continue
	enddo
c
c  Display results
c
852	call system('cls')
	do kk=1,numfil
	  write(5,855) filt(kk),fext(kk)
855	  format(3x,'k(',a3,') = ',f6.3)
	enddo
c
	print *
830	print *,'  Are these values from the extinction file ',
     $          'correct (y/n)? '
	call prompt(ok,flag)
	if(flag .eq. 0) go to 830
	if(ok .eq. 'N') erflg=.false.
	return
	end
