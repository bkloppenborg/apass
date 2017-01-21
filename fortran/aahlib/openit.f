	subroutine openit(iunit,file1,type,wrt,akcess,numfil,filt)
c
c  This routines opens the file and
c  optionally appends data to an existing file.
c
c  input: iunit = unit number
c         file1 =  file name
c         wrt   = true if writing to the file, false if read only
c         type = expected file id code
c
c  output: akcess = 'append' if user wishes to append to the file
c          numfil,filt = in the case of appending, will return the number
c                        of filters and the filter designations
c
c  external subroutines: ucase,getfilt
c
c  Written by R. Kaitchuck
c  Fortran 77 compatible
c  **** Last Revised 6-Aug-90
c
	character*30 file1,type
	character*10 ftype,akcess
	character*3 filt(5)
	character doit,yn
	logical wrt,found
	akcess='sequential'
c
c  there are only four possiblities ...
c    1. file exists and its to be written
c    2. file does not exist and its to be written
c    3. file exists and its to be read
c    4. file does not exist and its to be read
c
10	inquire(file=file1,exist=found)
c
c  Case 4:
c  if a file is to be read and doesn't exist, ask user to correct the name
c
	if((.not. wrt).and.(.not. found)) then
	   call bell
	   print *,'  File not found, try again --- '
	   read *,file1
           go to 10
	endif
c
c  Case 3:
c  if file exists and is to be read ...
c
	if(found .and. .not. wrt) ftype='old'
c
c  Case 1:
c  if file is to written to and it already exists, see if user wants to
c   overwrite or append
c
	if(found .and. wrt) then
20	  print *,'  This file already exists, do you wish to '
          print *,'  Overwrite, Append, or Quit (O,A,Q)? --- '
	  read *,doit
	  call ucase(doit,1)
	  if((doit.ne.'O').and.(doit.ne.'A').and.(doit.ne.'Q')) then
	  call bell
	  print *,'  Invalid response'
	  go to 20
	  endif
	  if(doit .eq. 'O') then 
202		print *,'  Are you sure (y/n)? --- '
		read *,yn
		call ucase(yn,1)
			if((yn.ne.'Y').and.(yn.ne.'N')) then
			call bell
			print *,'  Invalid response'
			go to 202
			endif
		if(yn .eq. 'N') go to 20
		ftype='unknown'
	  call backup(file1)
	  endif
c
	  if(doit .eq. 'A') then
		call getfilt(iunit,file1,type,filt,numfil,errflg)
c getfilt leaves file open and rewound
		if(errflg .eq. 1) stop
			    akcess='append'
101		read(iunit,*,end=102)
		go to 101
102		return
	   endif
	   if(doit .eq. 'Q') stop
	endif
c
c  Case 2:
c  if file does not exist and it is to be written to...
c
	if(.not. found .and. wrt) ftype='new'
c
c  open the file
c
	open(unit=iunit,file=file1,status=ftype,err=100) 
	return
100	print *,'  Error opening file.'
	call bell
	stop
	end
c
c
	subroutine backup(file1)
c
c  this routine creates a backup copy of a file before overwriting it
c
	character*30 file1,fx,input
	character*80 str
	do j=1,30
	fx(j:j)= ' '
	enddo
c
c  locate start of file extension
c
	i=index(file1,'.')
	if(i .ne. 0) then
		do k=1,i
		fx(k:k)=file1(k:k)
		enddo
		fx(i+1:i+3)='BAK'
	else
		do k=1,30
		if(file1(k:k) .eq. ' ') go to 100
		fx(k:k)=file1(k:k)
		enddo
100		fx(k:k+3)='.BAK'
	endif
c
c
	print *,'  Enter the backup file name [',fx(1:nblank(fx)),']: '
	read(5,'(a30)') input
	if(input .ne. ' ') fx=input
c
	str='copy '//file1//' '//fx
	call system(str)
	return
	end
