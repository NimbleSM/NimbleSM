#!/usr/bin/env python

import sys
import subprocess
import string

RELATIVE_TOLERANCE = 1.0e-6

def modify_tolerance(line):

    modified_line = line[:]

    # Modify the tolerance to be specified as an absolute tolerance
    # This is helpful when tweaking tolerances for specific tests, platforms, etc.

    if "min" in line and "max" in line and "@" in line and "TIME STEPS" not in line:

        # handle some specific cases where noise tends
        # to appear in values that should be zero

        min_tolerances = {}
        min_tolerances["stress"] = 10.0
        min_tolerances["def_grad"] = 1.0e-12
        min_tolerances["displacement"] = 1.0e-14
        min_tolerances["velocity"] = 1.0e-14
        min_tolerances["acceleration"] = 1.0e-7
        min_tolerances["force"] = 1.0e-9

        modify_tolerance = False
        for key in min_tolerances.keys():
            if key in line:
                modify_tolerance = True
                min_tol = min_tolerances[key]

        max_value = 0.0
        vals = string.splitfields(line)
        variable_name = vals[0]
        for i in range(len(vals)):
            if vals[i] == "max:":
                max_value = float(vals[i+1])

        absolute_tolerance = max_value*RELATIVE_TOLERANCE

        if modify_tolerance:
            if absolute_tolerance < min_tol:
                print "Modifying", variable_name, "tolerance from", absolute_tolerance, "to", min_tol
                absolute_tolerance = min_tol

        absolute_tolerance_string = "absolute {:12.12e}    ".format(absolute_tolerance)

        index = string.find(modified_line, '#')
        modified_line = line[:index] + absolute_tolerance_string + line[index:]

    return modified_line

if __name__ == "__main__":

    print "\nCreating exodiff file..."

    if len(sys.argv) != 2:
        print "Usage:  create_exodiff_file.py <exodus_file.e>\n"
        sys.exit(1)

    exodus_file_name = sys.argv[1]
    print "  Running \"exodiff -summary " + exodus_file_name + "\""
    exodiff_summary = subprocess.check_output(["exodiff", "-summary", exodus_file_name])

    exodiff_file_name = exodus_file_name[:string.rfind(exodus_file_name, '.')] + ".exodiff"
    exodiff_file = open(exodiff_file_name, 'w')

    lines = string.splitfields(exodiff_summary, '\n')
    for line in lines:
        modified_line = modify_tolerance(line)
        exodiff_file.write(modified_line + "\n")
    exodiff_file.close()

    print "  Exodiff tolerances written to", exodiff_file_name, "\n"
