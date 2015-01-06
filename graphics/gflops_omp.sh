#!/usr/bin/gnuplot

set macro

# start define plot styles ####################################################
# Palette URL:
# http://colorschemedesigner.com/#3K40zsOsOK-K-

red_000    = "#F9B7B0"
red_025    = "#F97A6D"
red_050    = "#E62B17"
red_075    = "#8F463F"
red_100    = "#6D0D03"

blue_000   = "#A9BDE6"
blue_025   = "#7297E6"
blue_050   = "#1D4599"
blue_075   = "#2F3F60"
blue_100   = "#031A49"

green_000  = "#A6EBB5"
green_025  = "#67EB84"
green_050  = "#11AD34"
green_075  = "#2F6C3D"
green_100  = "#025214"

brown_000  = "#F9E0B0"
brown_025  = "#F9C96D"
brown_050  = "#E69F17"
brown_075  = "#8F743F"
brown_100  = "#6D4903"

grid_color = "#6a6a6a"
text_color = "#6a6a6a"

my_line_width = "4"
my_axis_width = "1.5"
my_ps = "2"

# set the style for the set 1, 2, 3...
set style line 1 linetype 2 linecolor rgbcolor blue_050 linewidth @my_line_width
set style line 2 linetype 1 linecolor rgbcolor red_050  linewidth 2

# this is to use the user-defined styles we just defined.
set style increment user

# set the color and width of the axis border
set border 31 lw @my_axis_width lc rgb text_color

# set key options
#set key top left inside box width -1 height 1 enhanced spacing 1 samplen 3   \
#    font 'Helvetica, 24'
set key top left width -2 height 1

# set grid color
set grid lc rgb grid_color

# set text font and size
set terminal postscript eps font 'Helvetica,24'

set encoding utf8
# end define plot styles ######################################################


set xlabel "Number of cores"
set ylabel "Gflop/s"

datafile = "./samples.dat"
set title "MUrB performance, 100 000 bodies"

set output "gflops_omp.eps"
set xrange [1:24]
set xtics ("2" 2, "4" 4, "6" 6, "8" 8, "10" 10, "12" 12, "14" 14, "16" 16, "18" 18, "20" 20, "22" 22, "24" 24)

plot datafile using 1:4 i 2 with linespoint title 'MUrB', \
     x*73.6                 with lines      title 'peak performance'
