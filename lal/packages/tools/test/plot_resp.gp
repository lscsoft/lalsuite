set pm3d at s hidden3d 100 transparent
set style line 100 lt 5 lw 0.5
unset hidden3d
unset surf

set view 50, 45, 1, 1

set terminal x11 1
set palette model HSV rgbformulae 3,2,2 negative
set title "LAL response for det. @ North Pole - Cross pol."
splot 'ff_lal_cros.txt' matrix

set terminal x11 2
set title "LAL response for det. @ North Pole - Plus pol."
splot 'ff_lal_plus.txt' matrix

set terminal x11 3
set title "Local response for det. @ North Pole - Cross pol."
splot 'ff_local_cros.txt' matrix

set terminal x11 4
set title "Local response for det. @ North Pole - Plus pol."
splot 'ff_local_plus.txt' matrix

set terminal x11 1
set palette model HSV color rgbformulae 3,2,2 negative
