# gnuplot script to plot whole sky antenna patterns

set pm3d at s hidden3d 100 transparent
set style line 100 lt 5 lw 0.5
unset hidden3d
unset surf

set view 50, 45, 1, 1

set palette model HSV rgbformulae 3,2,2 negative

set terminal x11 1
set title "Cross"
splot 'whole_sky_cross.txt' matrix

set terminal x11 2
set title "Plus"
splot 'whole_sky_plus.txt' matrix

set terminal x11 3
set title "Sum"
splot 'whole_sky_sum.txt' matrix


set terminal x11 4
set title "Int. Cross"
splot 'int_whole_sky_cross.txt' matrix

set terminal x11 5
set title "Int. Plus"
splot 'int_whole_sky_plus.txt' matrix

set terminal x11 6
set title "Int. Sum"
splot 'int_whole_sky_sum.txt' matrix




#set palette model CMY functions gray, 1, 1
set palette model CMY rgbformulae 22,32,3 negative
set terminal pdf


set output "whole_sky_cross.pdf"
set title "Cross"
splot 'whole_sky_cross.txt' matrix

set output "whole_sky_plus.pdf"
set title "Plus"
splot 'whole_sky_plus.txt' matrix

set output "whole_sky_sum.pdf"
set title "Sum"
splot 'whole_sky_sum.txt' matrix


set output "int_whole_sky_cross.pdf"
set title "Int. Cross"
splot 'int_whole_sky_cross.txt' matrix

set output "int_whole_sky_plus.pdf"
set title "Int. Plus"
splot 'int_whole_sky_plus.txt' matrix

set output "int_whole_sky_sum.pdf"
set title "Int. Sum"
splot 'int_whole_sky_sum.txt' matrix




set palette model HSV rgbformulae 3,2,2 negative
set terminal png
set output "whole_sky_cross.png"
set title "Cross"
splot 'whole_sky_cross.txt' matrix

set output "whole_sky_plus.png"
set title "Plus"
splot 'whole_sky_plus.txt' matrix

set output "whole_sky_sum.png"
set title "Sum"
splot 'whole_sky_sum.txt' matrix

set output "int_whole_sky_cross.png"
set title "Int. Cross"
splot 'int_whole_sky_cross.txt' matrix

set output "int_whole_sky_plus.png"
set title "Int. Plus"
splot 'int_whole_sky_plus.txt' matrix

set output "int_whole_sky_sum.png"
set title "Int. Sum"
splot 'int_whole_sky_sum.txt' matrix


set terminal x11 6
