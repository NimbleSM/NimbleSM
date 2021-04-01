
set terminal pdf enhanced font "Times-Roman,24" size 12in, 8in
set output "wave_in_bar_end_displacement.pdf"
set xlabel "Time (s)" font "Times-Roman,32"
set ylabel "Displacement (cm)" font "Times-Roman,32"
#set xrange [0.0:0.015]
#set key at 0.01, 1.18
#set key font ",16"
plot "wave_in_bar_end_displacement.txt" using 1:2 with lines lw 4 lc 2 title "Simulation", \
     "analytic_solution_free_end.txt" using 1:2 with lines lw 4 lc 3 title "Analytical Result"
