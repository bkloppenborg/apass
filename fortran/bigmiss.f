      program bigmiss
      integer itot,ifld,icnt
      character file2*22
      open (unit=1,file='problem_files.txt',status='old')
      open (unit=3,file='rerun.txt',status='new')
10    continue
      read (1,900,end=300) file2
900     format (a22)
        open (unit=2,file=file2,status='old')
        read (2,*)
        itot=0
20      continue
        read (2,901,end=100) ifld,icnt
901     format (5x,i5,i10)
        if (ifld.eq.1.and.icnt.eq.0) itot = itot+1
        goto 20
100   continue
      close (2)
      if (itot.gt.30) then
        write (3,902) itot,file2
902     format (i5,2x,a22)
      endif
      goto 10
300   continue
      close (1)
      close (2)
      stop
      end
