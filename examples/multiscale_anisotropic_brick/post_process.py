#!/usr/bin/env python

# This script requires exodus.py (distributed with SEACAS)
import sys
sys.path.append('/Users/djlittl/Software/seacas/GCC_5.4.0/lib')
import exodus

import numpy
import string

if __name__ == "__main__":

    # Specimen geometry
    initial_length = 25.4
    initial_cross_section = 25.4 * 25.4

    if len(sys.argv) != 2:
        print "\nUsage:  post_process.py <exodus_file_name>\n"
        sys.exit(1)

    in_file_name = sys.argv[1]
    in_file = exodus.exodus(in_file_name, mode='r')

    # Extract nodal displacements and forces
    num_time_steps = in_file.num_times()
    times = in_file.get_times()
    node_variable_names = in_file.get_node_variable_names()
    if 'displacement_x' not in node_variable_names:
        print "\nERROR:  Failed to extract displacement_x data\n"
        sys.exit(1)
    if 'displacement_y' not in node_variable_names:
        print "\nERROR:  Failed to extract displacement_y data\n"
        sys.exit(1)
    if 'displacement_z' not in node_variable_names:
        print "\nERROR:  Failed to extract displacement_z data\n"
        sys.exit(1)
    if 'internal_force_x' not in node_variable_names:
        print "\nERROR:  Failed to extract internal_force_x data\n"
        sys.exit(1)

    # Read node sets
    node_set_ids = in_file.get_node_set_ids()
    node_set_nodes = {}
    for node_set_id in node_set_ids:
        node_ids = in_file.get_node_set_nodes(node_set_id)
        node_set_nodes[node_set_id] = node_ids[:]

    # In this particular case, we want to plot forces and displacements
    # where the forces and displacements are on specific node sets

    left_bottom_node_set = node_set_nodes[1]
    left_top_node_set = node_set_nodes[2]
    right_rear_node_set = node_set_nodes[11]
    right_front_node_set = node_set_nodes[12]

    left_bottom_displacement = []
    left_top_displacement = []
    left_bottom_force = []
    left_top_force = []
    right_rear_displacement = []
    right_front_displacement = []
    right_rear_force = []
    right_front_force = []

    stress_parallel_to_fibers = []
    strain_parallel_to_fibers = []
    stress_perpendicular_to_fibers = []
    strain_perpendicular_to_fibers = []

    for time_step in range(num_time_steps):

        displacement_x = in_file.get_node_variable_values('displacement_x', time_step+1)
        displacement_z = in_file.get_node_variable_values('displacement_z', time_step+1)
        force_x = in_file.get_node_variable_values('internal_force_x', time_step+1)
        force_z = in_file.get_node_variable_values('internal_force_z', time_step+1)

        displacement = 0.0
        force = 0.0
        for node_id in left_bottom_node_set:
            displacement += displacement_x[node_id-1]
            force += force_x[node_id-1]
        displacement /= len(left_bottom_node_set)
        left_bottom_displacement.append(displacement)
        left_bottom_force.append(force)

        displacement = 0.0
        force = 0.0
        for node_id in left_top_node_set:
            displacement += displacement_x[node_id-1]
            force += force_x[node_id-1]
        displacement /= len(left_top_node_set)
        left_top_displacement.append(displacement)
        left_top_force.append(force)

        displacement = 0.0
        force = 0.0
        for node_id in right_rear_node_set:
            displacement += displacement_z[node_id-1]
            force += force_z[node_id-1]
        displacement /= len(right_rear_node_set)
        right_rear_displacement.append(displacement)
        right_rear_force.append(force)

        displacement = 0.0
        force = 0.0
        for node_id in right_front_node_set:
            displacement += displacement_z[node_id-1]
            force += force_z[node_id-1]
        displacement /= len(right_front_node_set)
        right_front_displacement.append(displacement)
        right_front_force.append(force)

        strain_parallel_to_fibers.append( (right_front_displacement[-1] - right_rear_displacement[-1])/initial_length )
        stress_parallel_to_fibers.append( (right_rear_force[-1] - right_front_force[-1])/(2.0*initial_cross_section) )
        strain_perpendicular_to_fibers.append( (left_top_displacement[-1] - left_bottom_displacement[-1])/initial_length )
        stress_perpendicular_to_fibers.append( (left_bottom_force[-1] - left_top_force[-1])/(2.0*initial_cross_section) )

    in_file.close()

    # Print the specified material properties, which should bound the RVE response

    bulk_modulus_fiber = 1.563e5
    shear_modulus_fiber = 8.929e4
    youngs_modulus_fiber = 9.0 * bulk_modulus_fiber * shear_modulus_fiber / (3.0 * bulk_modulus_fiber + shear_modulus_fiber)
    bulk_modulus_matrix = 8.800e3
    shear_modulus_matrix = 9.103e2
    youngs_modulus_matrix = 9.0 * bulk_modulus_matrix * shear_modulus_matrix / (3.0 * bulk_modulus_matrix + shear_modulus_matrix)
    print "\nMaterial properties of RVE materials (should bound the macroscale response):"
    print "  Young's modulus for fiber material:                {0:.2e} MPa".format(youngs_modulus_fiber)
    print "  Young's modulus for matrix material:               {0:.2e} MPa".format(youngs_modulus_matrix)

    # Compute Young's modulus parallel to and perpendular to the fiber direction

    x = numpy.array(strain_parallel_to_fibers)
    y = numpy.array(stress_parallel_to_fibers)
    coeff = numpy.polyfit(x, y, 1)
    youngs_modulus_parallel = coeff[0]

    x = numpy.array(strain_perpendicular_to_fibers)
    y = numpy.array(stress_perpendicular_to_fibers)
    coeff = numpy.polyfit(x, y, 1)
    youngs_modulus_perpendicular = coeff[0]

    print "\nAnisotropic macroscale material response:"
    print "  Young's modulus parallel to fiber direction:       {0:.2e} MPa".format(youngs_modulus_parallel)
    print "  Young's modulus perpendicular to fiber direction:  {0:.2e} MPa".format(youngs_modulus_perpendicular)

    print
    out_file_name = 'force_displacement_data.txt'
    out_file = open(out_file_name, 'w')
    for time_step in range(num_time_steps):
        out_file.write(str(times[time_step]) + " " + str(left_bottom_displacement[time_step]) + " " + str(left_top_displacement[time_step]) + " " + str(left_bottom_force[time_step]) + " " + str(left_top_force[time_step]) + " " + str(right_rear_displacement[time_step]) + " " + str(right_front_displacement[time_step]) + " " + str(right_rear_force[time_step]) + " " + str(right_front_force[time_step]) + "\n")
    out_file.close()
    print "Data for written to", out_file_name
    print
