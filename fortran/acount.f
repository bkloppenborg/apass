      program acount
c
c look for missing files in epp
c kludge, hardcoded string sizes
c written 2016-11-03 aah
c
      integer icnt(1000),ifl
      integer i,lastfile,ifile
      character file1*6,file2*80,file3*80,infile*80,outfile*80
      character basefile*14,extfile*4
      logical ok
c
      print *,'Enter filelist: '
      read (5,'(a)') file3
      open (unit=3,file=file3,status='old')
      print *,'Enter output file: '
      read (5,'(a)') outfile
      open (unit=2,file=outfile,status='new')
      write (2,901)
901   format ('  night',2x,'afile',2x,'nfile',2x,'lfile')
10    continue
      read (3,910,end=800) basefile,extfile
910   format (a14,a4)
      print *,'processing night ',basefile
      j=0
      ifl=0
      do i=101,1000
        write (infile,900) basefile,i
900     format (a14,'a.',i4.4)
        inquire (file=infile,exist=ok)
        if (ok) then
          lastfile = i
          ifl = ifl+1
        endif
      enddo
      ifile = lastfile-100
      if ((ifile-ifl).gt.30) then
      write (2,9000) basefile(1:7),ifl,ifile,lastfile
9000  format (a7,2x,i5,2x,i5,2x,i5)
      endif
      goto 10
800   continue
      close(2)
      close(3)
      stop
      end
