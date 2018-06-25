
set terminal pdf enhanced font "Times-Roman,24" size 12in, 8in
set output "stress_strain.pdf"
set xlabel "Strain (cm/cm)" font "Times-Roman,32"
set ylabel "Stress (dyne/cm^2)" font "Times-Roman,32"
#set xrange [0.0:0.015]
#set key at 0.01, 1.18
#set key font ",16"
plot "quasistatic_tension_stress_strain.txt" using 1:2 with lines lw 4 lc 0 notitle
