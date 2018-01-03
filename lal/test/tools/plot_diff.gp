set pm3d at s hidden3d 100 transparent
set style line 100 lt 5 lw 0.5
unset hidden3d
unset surf

set view 50, 45, 1, 1

set terminal x11 1
set palette model HSV rgbformulae 3,2,2 negative
set title "Abs. difference between LAL and GRASP for det. @ North Pole - Cross pol."
splot 'ff_diff_cros.txt' matrix

set terminal x11 2
set title "Abs. difference between LAL and GRASP for det. @ North Pole - Plus pol."
splot 'ff_diff_plus.txt' matrix

set output "ff_diff_cros.pdf"
#set palette model CMY functions gray, 1, 1
set palette model CMY rgbformulae 22,32,3 negative
set terminal pdf
set title "Abs. difference between LAL and GRASP for det. @ North Pole - Cross pol."
splot 'ff_diff_cros.txt' matrix

set output "ff_diff_plus.pdf"
set title "Abs. difference between LAL and GRASP for det. @ North Pole - Plus pol."
splot 'ff_diff_plus.txt' matrix


set palette model RGB
set output "ff_diff_cros.png"
set terminal png
set title "Abs. difference between LAL and GRASP for det. @ North Pole - Cross pol."
splot 'ff_diff_cros.txt' matrix

set output "ff_diff_plus.png"
set title "Abs. difference between LAL and GRASP for det. @ North Pole - Plus pol."
splot 'ff_diff_plus.txt' matrix

set terminal x11
set palette model HSV color rgbformulae 3,2,2 negative
