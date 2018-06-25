#!/usr/bin/env python

# This script requires exodus.py

import sys
import os
import string
import sets

sys.path.append('/path/to/seacas/install')
import exodus

if __name__ == "__main__":

    if len(sys.argv) < 4:
        print "\nUsage:  parallel_preprocessing.py <input_mesh_file_base_name> <output_mesh_file_base_name> <num_ranks>\n"
        sys.exit(1)

    input_base_name = sys.argv[1]
    output_base_name = sys.argv[2]
    num_ranks = int(sys.argv[3])

    print "\n-- parallel_preprocessing.py --"
    print "\nGenesis input file base name =", input_base_name
    print "Genesis output file base name =", output_base_name
    print "Number of partitions", num_ranks

    rank_string_size = len(str(num_ranks))

    # create a list of all the file names
    input_file_names = []
    output_file_names = []
    for i_rank in range(num_ranks):
        rank_string = str(i_rank)
        while len(rank_string) < rank_string_size:
            rank_string = "0" + rank_string
        suffix = ".g." + str(num_ranks) + "." + rank_string
        input_file_names.append(input_base_name + suffix)
        output_file_names.append(output_base_name + suffix)

    # store the global node ids on each rank
    # may want to store this data in a Set structure (from the sets module) or something similar for fast searching
    # for now just put it in a list of lists
    global_node_ids = []

    # Read the global node ids for each input file
    for name in input_file_names:

        genesis_file = exodus.exodus(name, 'r')

        # default map from local to global ids
        # there may be additional auxiliary node maps present (created by decomp) but they can be ignored here
        node_id_map = genesis_file.get_node_id_map()

        global_node_ids.append(node_id_map)

        genesis_file.close()

    # PLACEHOLDER FOR CLIQUE CALCULATION
    # find nodes on partition boundaries, so that we have something to plot to test this script
    node_on_boundary = {}
    for node_ids in global_node_ids:
      for node_id in node_ids:
        if node_id in node_on_boundary.keys():
          node_on_boundary[node_id] = 1
        else:
          node_on_boundary[node_id] = 0

    node_var_name = 'clique'
    node_var_values = []
    for i_rank in range(len(global_node_ids)):
      values = []
      for node_id in global_node_ids[i_rank]:
        values.append(node_on_boundary[node_id])
      node_var_values.append(values)

    # write the modified genesis files
    for i_rank in range(len(input_file_names)):

      name = input_file_names[i_rank]
      pos = string.find(name, '.')
      new_name = output_base_name + name[pos:]

      if os.path.exists(new_name):
        os.remove(new_name)

      global_vars = []
      nodal_vars = [node_var_name]
      element_vars = []

      new_database = exodus.exodus(name, 'r').copy(new_name)
      exodus.add_variables(new_database, global_vars, nodal_vars, element_vars)
      new_database.put_node_variable_values(node_var_name, 1, node_var_values[i_rank])
      new_database.close()

