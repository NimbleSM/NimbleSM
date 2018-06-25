#!/usr/bin/env python

# This script requires exodus.py
# It can be tricky to get exodus.py to work because it requires libnetcdf.so and libexodus.so
import sys
sys.path.append('/Users/djlittl/Software/seacas/GCC_4.9.4_THREAD_SAFE/lib')
import exodus
import string

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "\nUsage:  PostProcess.py <exodus_file_name>\n"
        sys.exit(1)

    in_file_name = sys.argv[1]
    in_file = exodus.exodus(in_file_name, mode='r')

    out_file_label = string.splitfields(in_file_name, '.')[0] + "_"

    # Print database parameters from in_file
    print " "
    print "Database version:         " + str(round(in_file.version.value,2))
    print "Database title:           " + in_file.title()
    print "Database dimensions:      " + str(in_file.num_dimensions())
    print "Number of nodes:          " + str(in_file.num_nodes())
    print "Number of elements:       " + str(in_file.num_elems())
    print "Number of element blocks: " + str(in_file.num_blks())
    print "Number of node sets:      " + str(in_file.num_node_sets())
    print "Number of side sets:      " + str(in_file.num_side_sets())
    print " "

    # Extract the time steps
    time_vals = in_file.get_times()

    # Extract the model coordinates
    x_coords, y_coords, z_coords = in_file.get_coords()

    # Find the bar dimensions
    x_min = 1.0e50
    x_max = -1.0e50
    y_min = 1.0e50
    y_max = -1.0e50
    z_min = 1.0e50
    z_max = -1.0e50
    for i in range(len(x_coords)):
        if x_coords[i] < x_min:
            x_min = x_coords[i]
        if x_coords[i] > x_max:
            x_max = x_coords[i]
        if y_coords[i] < y_min:
            y_min = y_coords[i]
        if y_coords[i] > y_max:
            y_max = y_coords[i]
        if z_coords[i] < z_min:
            z_min = z_coords[i]
        if z_coords[i] > z_max:
            z_max = z_coords[i]

    tol = 1.0e-8
    x_min_index = []
    x_max_index = []
    for i in range(len(x_coords)):
        if abs(x_coords[i] - x_min) < tol:
            x_min_index.append(i)
        if abs(x_coords[i] - x_max) < tol:
            x_max_index.append(i)

    initial_length = x_max - x_min
    cross_section = (y_max - y_min) * (z_max - z_min)

    print "Min x value =", x_min, "at nodes", x_min_index
    print "Max x value =", x_max, "at nodes", x_max_index
    print "Initial length =", initial_length
    print "Cross sectional area =", cross_section

    # Extract nodal displacements and forces
    num_time_steps = in_file.num_times()
    node_variable_names = in_file.get_node_variable_names()
    if 'displacement_x' not in node_variable_names:
        print "\nERROR:  Failed to extract displacement_x data\n"
        sys.exit(1)
    if 'internal_force_x' not in node_variable_names:
        print "\nERROR:  Failed to extract internal_force_x data\n"
        sys.exit(1)

    force = []
    displacement = []

    for time_step in range(num_time_steps):

        displacement_x = in_file.get_node_variable_values('displacement_x', time_step+1)
        force_x = in_file.get_node_variable_values('internal_force_x', time_step+1)

        net_force_min_face = 0.0
        net_force_max_face = 0.0
        ave_displacement_min_face = 0.0
        ave_displacement_max_face = 0.0

        for index in x_min_index:
            net_force_min_face += force_x[index]
            ave_displacement_min_face += displacement_x[index]
        ave_displacement_min_face /= len(x_min_index)

        for index in x_max_index:
            net_force_max_face += force_x[index]
            ave_displacement_max_face += displacement_x[index]
        ave_displacement_max_face /= len(x_min_index)


        force.append(0.5*(net_force_min_face - net_force_max_face))
        displacement.append(ave_displacement_max_face - ave_displacement_min_face)

    in_file.close()

    stress = []
    strain = []
    for i in range(len(force)):
        stress.append(force[i]/cross_section)
        strain.append(displacement[i]/initial_length)
    print

    out_file_name = out_file_label + 'stress_strain.txt'
    data_file = open(out_file_name, 'w')
    for i in range(num_time_steps):
        data_file.write(str(strain[i]) + "  " + str(stress[i]) + "\n")
    data_file.close()
    print "Stress-strain data written to", out_file_name, "\n"
