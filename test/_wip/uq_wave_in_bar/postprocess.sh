#!/bin/bash
scrape_1D.py uq_wave_in_bar.serial.e
rm solution.dat
if [ ! -e solution.dat ]; then
analytic_solution.py x.dat t.dat > solution.log
fi 
paste solution.dat uq_wave_in_bar.dat > u.dat
