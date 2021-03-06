filelist=system("ls *.dat")

set term png size 14000,7000

# configure the pallette
set cbrange [0:4]
set palette maxcolors 4
set palette defined (0 "red", 1 "red", 2 "orange", 3 "green", 4 "purple")

# configure the plot
set nokey
set grid
set format x "%3.0f"
set format y "%3.0f"
set xtics  border mirror -2,2,362 font ",20"
set x2tics border mirror -2,2,362 font ",20"
set ytics  border mirror -91,2,91 font ",20"
set y2tics border mirror -91,2,91 font ",20"
set xrange [0:360]
set yrange [-90:90]

set title "APASS B Coverage Map"
set output 'apass_B.png'
plot for [filename in filelist] filename using 2:4:12 with points palette pointtype 0

set title "APASS V Coverage Map"
set output 'apass_V.png'
plot for [filename in filelist] filename using 2:4:13 with points palette pointtype 0

set title "APASS su Coverage Map"
set output 'apass_su.png'
plot for [filename in filelist] filename using 2:4:14 with points palette pointtype 0

set title "APASS sg Coverage Map"
set output 'apass_sg.png'
plot for [filename in filelist] filename using 2:4:15 with points palette pointtype 0

set title "APASS sr Coverage Map"
set output 'apass_sr.png'
plot for [filename in filelist] filename using 2:4:16 with points palette pointtype 0

set title "APASS si Coverage Map"
set output 'apass_si.png'
plot for [filename in filelist] filename using 2:4:17 with points palette pointtype 0

set title "APASS sz Coverage Map"
set output 'apass_sz.png'
plot for [filename in filelist] filename using 2:4:18 with points palette pointtype 0

set title "APASS Ha Coverage Map"
set output 'apass_Ha.png'
plot for [filename in filelist] filename using 2:4:19 with points palette pointtype 0

set title "APASS Zs Coverage Map"
set output 'apass_Zs.png'
plot for [filename in filelist] filename using 2:4:20 with points palette pointtype 0

set title "APASS Y Coverage Map"
set output 'apass_Y.png'
plot for [filename in filelist] filename using 2:4:20 with points palette pointtype 0