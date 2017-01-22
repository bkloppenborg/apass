      program missing_files
c
c look for missing files in epp
c kludge, hardcoded string sizes
c written 2016-11-03 aah
c
      integer icnt(1000),ifl(1000)
      integer i,lastfile,ifile
      character file1*6,file2*80,file3*80,infile*80,outfile*80
      character basefile*14,extfile*4
      logical ok
c
      print *,'Enter filelist '
      read (5,'(a)') file3
      open (unit=3,file=file3,status='old')
10    continue
      read (3,910,end=800) basefile,extfile
910   format (a14,a4)
      print *,'processing file ',basefile,extfile
      do i=1,1000
        icnt(i) = 0
        ifl(i) = 0
      enddo
      j=0
      do i=101,1000
        write (infile,900) basefile,i
900     format (a14,'a.',i4.4)
        inquire (file=infile,exist=ok)
        if (ok) then
          lastfile=i
          ifl(i) = 1
        endif
      enddo
      infile = basefile//'.epp'
      outfile = basefile//'_log.txt'
      open (unit=2,file=outfile,status='new')
      open (file=infile,unit=1,status='old')
      read (1,*) 
      read (1,*) 
100   continue
        read (1,901,end=400,err=300) ifile
901     format (265x,i5)
        if (ifile.lt.101.or.ifile.gt.1000) goto 100
        icnt(ifile) = icnt(ifile)+1
        goto 100
300   continue
      print *," error in file, ",ifile
400   continue
      close(1)
c
c now make output file
c
      itest = 0
      do i=101,lastfile
        if (ifl(i).eq.1.and.icnt(i).eq.0) itest = 1
      enddo
      write (2,903) itest, basefile,extfile
903   format ('xxx ',i5,3x,a14,a4)
      do i=101,lastfile
        write(2,904) i,ifl(i),icnt(i)
904     format (2i5,i10)
      enddo
      close(2)
      goto 10
800   continue
      stop
      end
