        subroutine jday(id,im,iyr,uthr,date)
c this sub calculates julian date from ut date
c inputs are id,im,iy, uthr--day,month,year,and ut hour
c output is date, julian date-2,400,000
c Henden, 1977
c
	double precision date,jdo,iiy
        dimension mo(12)
        integer id,im,iyr,iy
        data mo /0,31,59,90,120,151,181,212,243,273,304,334/
c
c leap is number of leap days since 1900
c
        iy = iyr - 1900
        leap=iy/4
        if((4*leap-iy) .eq. 0 .and. im .lt. 3) leap=leap-1
        iiy=365.0*float(iy)
        jdo=15020.0+iiy+float(leap)+float(mo(im))+float(id)
        date=jdo+uthr/24.0-0.5
        return
        end
