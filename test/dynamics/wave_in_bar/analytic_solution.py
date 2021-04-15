#!/usr/bin/env python

from math import pi, sqrt, sin

if __name__ == "__main__":

    # Description of bar
    density = 7.8     # g/cm^3
    modulus = 3.0e12  # dyne/cm^2
    length  = 2.0     # cm

    # Initial velocity
    initialVelocity = 1000.0 # cm/sec

    numEntries = 1000
    startTime = 0.0 # sec
    endTime = 1.0e-5 # sec
    times = [i*endTime/numEntries for i in range(numEntries+1)];
    displacements = []

    numTermsInSeries = 1000

    # Loop from 1 to N
    loopRange = range(numTermsInSeries+1)[1:]

    #
    # Analytic solution for midpoint
    #
    
    # Position along bar
    x = 1.0 # cm

    for t in times:

        displacement = 0.0
        for n in loopRange:
            wn = (2*n-1) * (pi/(2.0*length)) * sqrt(modulus/density)
            An = 8.0 * length * initialVelocity * sqrt(density/modulus)/ ((2*n-1)*(2*n-1)*pi*pi)
            displacement = displacement + ( An * sin(wn*t) * sin( (2*n-1)*pi*x/(length*2.0) ) )
        displacements.append(displacement)

    outFile = open("analytic_solution_midpoint.txt", 'w')
    for i in range(len(displacements)):
        outFile.write(str(times[i]) + " " + str(displacements[i]) + "\n")
    outFile.close()

    displacements = []
    
    #
    # Analytic solution for free end
    #
    
    # Position along bar
    x = 2.0 # cm

    for t in times:

        displacement = 0.0
        for n in loopRange:
            wn = (2*n-1) * (pi/(2.0*length)) * sqrt(modulus/density)
            An = 8.0 * length * initialVelocity * sqrt(density/modulus)/ ((2*n-1)*(2*n-1)*pi*pi)
            displacement = displacement + ( An * sin(wn*t) * sin( (2*n-1)*pi*x/(length*2.0) ) )
        displacements.append(displacement)

    outFile = open("analytic_solution_free_end.txt", 'w')
    for i in range(len(displacements)):
        outFile.write(str(times[i]) + " " + str(displacements[i]) + "\n")
    outFile.close()

    print "\nAnalytic solutions written to analytic_solution_midpoint.txt and analytic_solution_free_end.txt\n"
    
