	subroutine getfilt(iunit,file1,type,filt,numfil,errflg)
c
c  this routine reads the four different file types and returns the filter
c  designations used
c
c  input: iunit = unit number for file
c         file1 = file name
c         type  = expected file id code
c
c  output: filt = filter codes
c          numfil = number of filters (or mag + colors)
c          errflg = error flag (1 if error found)
c
	character*30 id,type,file1
	character*3 filt(5)
	errflg=0
c
c  open file and read its ID code
c
	open(unit=iunit,file=file1,status='old',err=999)
	read (iunit,90) id
90	format(a30)
	if(type .ne. id) then
	call bell
	print *,'  This is the wrong type of file.'
	errflg=1
	return
	endif
c
c  process by id type
c
	if(type .eq. 'RAW DATA FILE') go to 101
	if(type .eq. 'INSTRUMENTAL MAGNITUDES') go to 102
	if(type .eq. 'STAR LIST FILE') go to 103
	if(type .eq. 'STANDARD MAGNITUDES AND COLORS') go to 104
	if(type .eq. 'DIFFERENTIAL MAGNITUDES') go to 105
c
	call bell
	print *,'  This is the wrong type of file.'
	errflg=1
	return
c
101	do 70 i=1,4
70	read(iunit,*)
	call filters(iunit,numfil,filt)
	go to 999
c
102	read(iunit,*)
	call filters(iunit,numfil,filt)
	go to 999
c
103	call labels(iunit,numfil,filt)
	go to 999
c
104	read(iunit,*)
	call labels(iunit,numfil,filt)
	go to 999
c
105	read(iunit,*)
	read(iunit,*)
	call labels(iunit,numfil,filt)
c
999	rewind (iunit)
	return
	end	

	subroutine filters(iunit,numfil,filt)
c
c  See how many and which filters are used 
c
	character*3 filt(5),fil1
	character*10 star
c
	ii=1
	read(iunit,*) star,filt(1)
	fil1=filt(1)	
500	ii=ii+1
	if(ii .eq. 6) go to 499
	read(iunit,*,end=499) star,filt(ii)
	if(fil1 .ne. filt(ii)) go to 500
499	numfil=ii-1
	return
	end

	subroutine labels(iunit,numfil,filt)
c
c  get the filter labels
c
	character*3 filt(5)
	character text(40),c1,c2,c3
c
	read (iunit,1000) text
1000	format(40x,40a1)
	kl=1
	j=1
2500	if(text(j) .eq. ' ') go to 2000
		c1=text(j)
		c2=text(j+1)
		c3=text(j+2)
		filt(kl)=c1//c2//c3
		j=j+2
		kl=kl+1
2000	j=j+1
	if(j .le. 40) go to 2500
	numfil=kl-1
	return
	end
