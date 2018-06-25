
# specimen dimensions
length = 25.4
cross_section = 645.16

set terminal pdf enhanced font "Times-Roman,24" size 12in, 8in
set title "Macroscale Engineering Stress-Strain Curves"
set output "stress_strain.pdf"
set xlabel "Engineering Strain (m/m)" font "Times-Roman,32"
set ylabel "Engineering Stress (MPa)" font "Times-Roman,32"
set key bottom right
plot "force_displacement_data.txt" using (($7 - $6)/length):(($8 - $9)/(2.0*cross_section)) with points pt 7 ps 2 lc 2 title "Parallel to fiber direction", \
     "force_displacement_data.txt" using (($3 - $2)/length):(($4 - $5)/(2.0*cross_section)) with points pt 7 ps 2 lc 3 title "Perpendicular to fiber direction"


