#!/usr/bin/env python

# This script requires exodus.py
# It can be tricky to get exodus.py to work because it requires libnetcdf.so and libexodus.so
import sys
sys.path.append('/Users/djlittl/Software/seacas/lib')
import exodus
import string
import numpy
import matplotlib.pyplot as plt

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

    # Find the node id for the end node
    x_min = 1.0e50
    x_min_node_index = -1
    for i in range(len(x_coords)):
        if x_coords[i] < x_min:
            x_min = x_coords[i]
            x_min_node_index = i

    print "Min x value =", x_min, "at node index", x_min_node_index, "\n"

    # Extract nodal displacements
    num_time_steps = in_file.num_times()
    node_variable_names = in_file.get_node_variable_names()
    if 'displacement_x' not in node_variable_names:
        print "\nERROR:  Failed to extract displacement_x data\n"
        sys.exit(1)

    bar_end_displacement = []
    for timeStep in range(num_time_steps):

        displacement_x = in_file.get_node_variable_values('displacement_x', timeStep+1)
        bar_end_displacement.append(displacement_x[x_min_node_index])

    in_file.close()

    print

    out_file_name = out_file_label + 'end_displacement.txt'
    data_file = open(out_file_name, 'w')
    for i in range(num_time_steps):
        data_file.write(str(time_vals[i]) + "  " + str(bar_end_displacement[i]) + "\n")
    data_file.close()
    print "Bar end displacement data for written to", out_file_name, "\n"

#    out_file_name = out_file_label + 'load_displacement.pdf' 
#    fig, ax = plt.subplots()
#    ax.plot(nodeset_displacement[:],nodeset_force[:],color='blue',marker='o',label='1 elem/block')
#    plt.xlabel('displacement (mm)')
#    plt.ylabel('force (N)')
#    lg = plt.legend(loc = 4)
#    lg.draw_frame(False)
#    plt.tight_layout()
#    plt.show()
#    fig.savefig('load_displacment.pdf')
