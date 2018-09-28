#if(ARG1 eq '')\
#  print "usage: gnuplot [-p] -c ".ARG0." datafile";\
#  exit
#
#datafile=ARG1
datafile='bench'

set style fill solid
set boxwidth 0.5

set grid

set key autotitle columnhead
set key right bottom

set xlabel 'MPIxOMP'
set ylabel 'MMUPs/s'

set offsets 0.5, 0.5
set yrange[0:1.5]

set multiplot layout 2, 1 title "Gromacs Hybrid Performance" font ",14"
plot \
    datafile \
    index 1 \
    using 4:xtic(1) \
    title '32-Cores' \
    with linespoints

plot \
    datafile \
    index 2 \
    using 4:xtic(1) \
    title '16-Cores' \
    with linespoints
