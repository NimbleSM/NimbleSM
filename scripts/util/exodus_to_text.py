#!/usr/bin/env python

import os
import sys

# This script requires exodus.py (distributed with SEACAS)
import sys
sys.path.append('/path/to/seacas/install')
import exodus

from ctypes import *
import os

NETCDF_SO = "/path/to/seacas/lib/libnetcdf.dylib"
if os.uname()[0] == 'Darwin':
  EXODUS_SO = "/path/to/seacas/lib/libexodus.dylib"
else:
  EXODUS_SO = "/path/to/seacas/lib/libexodus.so"
NETCDF_LIB = cdll.LoadLibrary(NETCDF_SO)
EXODUS_LIB = cdll.LoadLibrary(EXODUS_SO)

EX_MAPS_INT64_API    = 0x2000  # all maps (id, order, ...) store int64_t values

def ConvertGenesisFileToText(exodus_file_name):

  exodus_file = exodus.exodus(exodus_file_name)

  print " "
  print "Database version:         " + str(round(exodus_file.version.value,2))
  print "Database title:           " + exodus_file.title()
  print "Database dimensions:      " + str(exodus_file.num_dimensions())
  print "Number of nodes:          " + str(exodus_file.num_nodes())
  print "Number of elements:       " + str(exodus_file.num_elems())
  print "Number of element blocks: " + str(exodus_file.num_blks())
  print "Number of node sets:      " + str(exodus_file.num_node_sets())
  print "Number of side sets:      " + str(exodus_file.num_side_sets())
  print " "

  if exodus_file.num_side_sets() > 0:
    print "Warning:  Side sets are present in the genesis file, but will be ignored (not currently implemented).\n"

  # nodal coordinates
  x_coord, y_coord, z_coord = exodus_file.get_coords()

  # default maps from local to global ids
  node_id_map = exodus_file.get_node_id_map()
  elem_id_map = exodus_file.get_elem_id_map()

  # check for auxiliary node maps, which are used by decomp
  # currently no direct interface in exodus.py, hence some ugliness
  if exodus_file.EXODUS_LIB.ex_int64_status(exodus_file.fileId) & EX_MAPS_INT64_API:
    aux_node_map = (c_longlong * exodus_file.numNodes.value)()
    num_aux_node_maps = c_longlong(0)
    num_aux_elem_maps = c_longlong(0)
  else:
    aux_node_map = (c_int * exodus_file.numNodes.value)()
    num_aux_node_maps = c_longlong(0)
    num_aux_elem_maps = c_longlong(0)

  exodus_file.EXODUS_LIB.ex_get_map_param(exodus_file.fileId, byref(num_aux_node_maps), byref(num_aux_elem_maps))

  if num_aux_node_maps.value == 1:
    print "Applying auxiliary node map (this is probably a decomposed mesh file).\n"
    obj_type = exodus.ex_entity_type("EX_NODE_MAP")
    exodus_file.EXODUS_LIB.ex_get_num_map(exodus_file.fileId, obj_type, 1, byref(aux_node_map))
    node_id_map = aux_node_map

  # element block ids
  block_ids = exodus_file.get_elem_blk_ids()

  # node sets ids
  node_set_ids = exodus_file.get_node_set_ids()

  node_set_parameters = {}
  for node_set_id in node_set_ids:
    node_set_parameters[node_set_id] = exodus_file.get_node_set_params(node_set_id)

  # node set node list contains the node ids for each node set
  node_set_node_ids = {}
  for node_set_id in node_set_ids:
    node_ids = exodus_file.get_node_set_nodes(node_set_id)
    node_set_node_ids[node_set_id] = []
    for i in range(len(node_ids)):
      node_set_node_ids[node_set_id].append(node_ids[i])

  # loop over element blocks and record the connectivity
  num_elements = 0
  elem_conn = {}
  for block_id in block_ids:
    elem_conn[block_id] = []
    print "Working on element block with ID:  " + str(block_id)
    print "Number of elements in this block:  " + str(exodus_file.num_elems_in_blk(block_id))
    print "Type of element in this block:     " + exodus_file.elem_type(block_id)
    connect, num_elem_this_block, num_nodes_per_elem = exodus_file.get_elem_connectivity(block_id)
    num_elements += num_elem_this_block
    for element_index in range(num_elem_this_block):
      elem_conn[block_id].append([])
      for node_index in connect[element_index*num_nodes_per_elem:(element_index+1)*num_nodes_per_elem]:
        elem_conn[block_id][element_index].append(node_index)

  output_file_name = exodus_file_name + ".txt"
  if os.path.exists(output_file_name):
    os.remove(output_file_name)
  output_file = open(output_file_name, mode='w')
  output_file.write("number_of_nodes " + str(len(x_coord)) + "\n")
  output_file.write("number_of_elements " + str(num_elements) + "\n")
  output_file.write("number_of_blocks " + str(len(block_ids)) + "\n")
  output_file.write("number_of_node_sets " + str(len(node_set_ids)) + "\n")
  output_file.write("node_coordinates\n")
  for index in range(len(x_coord)):
    output_file.write(str(node_id_map[index]) + " " \
                      + "{0:.16g}".format(x_coord[index]) + " " \
                      + "{0:.16g}".format(y_coord[index]) + " " \
                      + "{0:.16g}".format(z_coord[index]) + "\n")
  global_elem_index = 0
  for block_id in block_ids:
    output_file.write("element_block\n" + str(block_id) + \
                      " block_"+str(block_id) + " " \
                      + str(len(elem_conn[block_id])) + " " \
                      + str(len(elem_conn[block_id][0])) + "\n")
    for elem_connectivity in elem_conn[block_id]:
      conn_str = str(elem_id_map[global_elem_index]) + " "
      global_elem_index += 1
      for i_node in range(len(elem_connectivity)):
        node_id = elem_connectivity[i_node]
        conn_str += str(node_id)
        if i_node < len(elem_connectivity) - 1:
          conn_str += " "
      output_file.write(conn_str + "\n")
  for node_set_id in node_set_ids:
    output_file.write("nodeset\n" + str(node_set_id) + " " \
                      + "nodelist_"+str(node_set_id) + " " \
                      + str(len(node_set_node_ids[node_set_id])) + "\n")
    for node_id in node_set_node_ids[node_set_id]:
      output_file.write(str(node_id) + "\n")
  output_file.close()

  print "\nMesh written to", output_file_name, "\n"

  return

if __name__ == "__main__":

    print "\n-- exodus_to_text --\n"

    usage_message = "Usage:  exodus_to_text.py <exodus_file.e> [-p] [num_partitions] [-convert_test_directory]\n"

    if len(sys.argv) != 2 and len(sys.argv) != 4:
        print usage_message
        sys.exit(1)

    convert_test_directory = False
    if len(sys.argv) == 2:
      if sys.argv[1] == "-convert_test_directory":
        convert_test_directory = True
      else:
        exodus_file_name = sys.argv[1]

    partitioned_mesh = False
    if len(sys.argv) == 4:
      num_partitions = int(sys.argv[3])
      if num_partitions > 1:
        partitioned_mesh = True

    if partitioned_mesh == False and convert_test_directory == False:
      ConvertGenesisFileToText(exodus_file_name)
    elif partitioned_mesh == True:
      print "Processing", num_partitions, "partitions"
      for i in range(num_partitions):
        partition_file_name = exodus_file_name + "." + str(num_partitions) + "." + str(i)
        print " ", partition_file_name
        ConvertGenesisFileToText(partition_file_name)
    elif convert_test_directory == True:
      # Convert all the .g files in the test repository
      files_to_convert = []
      root_dir = os.getcwd()
      directory_tree = os.walk("./test/")
      for item in directory_tree:
        directory = item[0]
        file_names = item[2]
        for file_name in file_names:
          if file_name[-2:] == ".g" or (file_name[-6:-3] == ".g." and file_name[-4:] != ".txt"):
            files_to_convert.append([directory, file_name])
      os.chdir(root_dir)
      for item in files_to_convert:
        directory = item[0]
        file_name = item[1]
        os.chdir(directory)
        ConvertGenesisFileToText(file_name)
        os.chdir(root_dir)
    else:
      print usage_message
      sys.exit(1)

