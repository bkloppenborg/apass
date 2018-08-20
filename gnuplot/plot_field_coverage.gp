set term png size 14000,7000

# configure the palette
#set cbrange [1:10]
#set palette maxcolors 10
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
set output 'apass_B_field.png'
plot 'fred_field_summary.txt' using 2:3:4 with points palette pointtype 5 pointsize 5.0

set title "APASS V Coverage Map"
set output 'apass_V_field.png'
plot 'fred_field_summary.txt' using 2:3:5 with points palette pointtype 5 pointsize 5.0

set title "APASS su Coverage Map"
set output 'apass_su_field.png'
plot 'fred_field_summary.txt' using 2:3:6 with points palette pointtype 5 pointsize 5.0

set title "APASS sg Coverage Map"
set output 'apass_sg_field.png'
plot 'fred_field_summary.txt' using 2:3:7 with points palette pointtype 5 pointsize 5.0

set title "APASS sr Coverage Map"
set output 'apass_sr_field.png'
plot 'fred_field_summary.txt' using 2:3:8 with points palette pointtype 5 pointsize 5.0

set title "APASS si Coverage Map"
set output 'apass_si_field.png'
plot 'fred_field_summary.txt' using 2:3:9 with points palette pointtype 5 pointsize 5.0

set title "APASS sz Coverage Map"
set output 'apass_sz_field.png'
plot 'fred_field_summary.txt' using 2:3:10 with points palette pointtype 5 pointsize 5.0

set title "APASS Ha Coverage Map"
set output 'apass_Ha_field.png'
plot 'fred_field_summary.txt' using 2:3:11 with points palette pointtype 5 pointsize 5.0

set title "APASS Zs Coverage Map"
set output 'apass_Zs_field.png'
plot 'fred_field_summary.txt' using 2:3:12 with points palette pointtype 5 pointsize 5.0

set title "APASS Y Coverage Map"
set output 'apass_Y_field.png'
plot 'fred_field_summary.txt' using 2:3:13 with points palette pointtype 5 pointsize 5.0