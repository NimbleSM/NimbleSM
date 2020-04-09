#!/bin/bash
../../../build_NimbleSM_UQ/src/NimbleSM_Serial uq_wave_in_bar.in
#scrape_1D.py uq_wave_in_bar.serial.e
#analytic_solution.py x.dat t.dat > solution.log
#paste solution.dat uq_wave_in_bar.dat > u.dat
exodiff -f uq_wave_in_bar.exodiff uq_wave_in_bar.serial.e uq_wave_in_bar.gold.e > test.log ; tail -n2 test.log
