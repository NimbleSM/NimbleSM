#!/usr/bin/env python

import os
import sys
import string

# This script requires exodus.py (distributed with SEACAS)
import sys
sys.path.append('/path/to/seacas/install')
import exodus

def ConvertTextFileToExodus(partition_file_name):

  text_file = open(partition_file_name)
  lines = text_file.readlines()
  text_file.close()

  number_of_nodes = 0
  number_of_elements = 0
  number_of_blocks = 0
  number_of_node_sets = 0
  number_of_side_sets = 0
  number_of_global_data = 0
  number_of_node_data = 0
  global_node_ids = []
  global_element_ids = []
  node_coordinate_x = []
  node_coordinate_y = []
  node_coordinate_z = []
  global_data_names = []
  block_ids = []
  block_names = []
  block_connectivity = []
  number_of_elements_per_block = []
  number_of_nodes_per_element = []
  node_set_ids = []
  node_set_names = []
  node_sets = []
  times = []
  global_data = []
  node_data_names = []
  node_data = []
  element_data_names = []
  element_data = []

  number_of_lines = len(lines)
  i_line = 0

  time_step = -1

  while i_line < number_of_lines:

    vals = string.splitfields(lines[i_line])
    i_line += 1

    if len(vals) > 0:
      key = vals[0]

      if key == "number_of_nodes":
        number_of_nodes = int(vals[1])

      elif key == "number_of_elements":
        number_of_elements = int(vals[1])

      elif key == "number_of_blocks":
        number_of_blocks = int(vals[1])

      elif key == "number_of_node_sets":
        number_of_node_sets = int(vals[1])

      elif key == "node_coordinates":
        for i_node in range(number_of_nodes):
          vals = string.splitfields(lines[i_line])
          i_line += 1
          global_id = int(vals[0])
          x = float(vals[1])
          y = float(vals[2])
          z = float(vals[3])
          global_node_ids.append(global_id)
          node_coordinate_x.append(x)
          node_coordinate_y.append(y)
          node_coordinate_z.append(z)

      elif key == "element_block":
        vals = string.splitfields(lines[i_line])
        i_line += 1
        block_id = int(vals[0])
        block_ids.append(block_id)
        block_name = vals[1]
        block_names.append(block_name)
        num_elem_in_block = int(vals[2])
        number_of_elements_per_block.append(num_elem_in_block)
        num_node_per_elem = int(vals[3])
        number_of_nodes_per_element.append(num_node_per_elem)
        connectivity = []
        for i_elem in range(num_elem_in_block):
          vals = string.splitfields(lines[i_line])
          i_line += 1
          global_id = int(vals[0])
          global_element_ids.append(global_id)
          for i_node in range(num_node_per_elem):
            connectivity.append(int(vals[1 + i_node]))
        block_connectivity.append(connectivity)

      elif key == "nodeset":
        vals = string.splitfields(lines[i_line])
        i_line += 1
        nodeset_id = int(vals[0])
        node_set_ids.append(nodeset_id)
        nodeset_name = vals[1]
        node_set_names.append(nodeset_name)
        num_nodes_in_nodeset = int(vals[2])
        nodeset = []
        for i_node in range(num_nodes_in_nodeset):
          vals = string.splitfields(lines[i_line])
          i_line += 1
          nodeset.append(int(vals[0]))
        node_sets.append(nodeset)

      elif key == "global_data":
        vals = string.splitfields(lines[i_line])
        i_line += 1
        number_of_global_data = int(vals[0])
        for i_data in range(number_of_global_data):
          vals = string.splitfields(lines[i_line])
          i_line += 1
          global_data_names.append(vals[0])

      elif key == "node_data":
        vals = string.splitfields(lines[i_line])
        i_line += 1
        number_of_node_data = int(vals[0])
        for i_data in range(number_of_node_data):
          vals = string.splitfields(lines[i_line])
          i_line += 1
          node_data_names.append(vals[0])

      elif key == "element_data":
        vals = string.splitfields(lines[i_line])
        i_line += 1
        number_of_element_data = int(vals[0])
        for i_data in range(number_of_element_data):
          vals = string.splitfields(lines[i_line])
          i_line += 1
          element_data_names.append(vals[0])

      elif key == "time":
        vals = string.splitfields(lines[i_line])
        i_line += 1
        time = float(vals[0])
        times.append(time)
        time_step += 1
        global_data.append([])
        node_data.append([])
        element_data.append([])
        for i_block in range(number_of_blocks):
          element_data[-1].append([])
          for i_elem_data in range(number_of_element_data):
            element_data[-1][-1].append([])

      elif key == "global_data_values":
        for i_data in range(number_of_global_data):
          vals = string.splitfields(lines[i_line])
          i_line += 1
          global_data[time_step].append(float(vals[0]))

      elif key == "node_data_values":
        for i_data in range(number_of_node_data):
          node_data[time_step].append([])
          for i_node in range(number_of_nodes):
            vals = string.splitfields(lines[i_line])
            i_line += 1
            node_data[time_step][i_data].append(float(vals[0]))

      elif key == "element_data_values":
        number_of_data_lists = number_of_blocks * number_of_element_data
        for i_list in range(number_of_data_lists):
          vals = string.splitfields(lines[i_line])
          i_line += 1
          data_name = vals[0]
          element_variable_index = int(vals[1]) - 1
          block_id = int(vals[2])
          num_elem_in_block = int(vals[3])
          block_index = block_ids.index(block_id)
          for i in range(num_elem_in_block):
            vals = string.splitfields(lines[i_line])
            i_line += 1
            datum = float(vals[0])
            element_data[time_step][block_index][element_variable_index].append(datum)

    key = "none"

  exodus_file_name = partition_file_name[:-4]
  exodus_title = "text_to_exodus.py translation of " + partition_file_name

  print "Converting", partition_file_name, "to", exodus_file_name

  if os.path.exists(exodus_file_name):
    os.remove(exodus_file_name)

  number_of_dimensions = 3
  exodus_file = exodus.exodus( exodus_file_name,
                               'w',
                               'ctype',
                               exodus_title,
                               number_of_dimensions,
                               number_of_nodes,
                               number_of_elements,
                               number_of_blocks,
                               number_of_node_sets,
                               number_of_side_sets )

  coord_names = ["X", "Y", "Z"]
  exodus_file.put_coord_names(coord_names)
  exodus_file.put_coords(node_coordinate_x, node_coordinate_y, node_coordinate_z)

  for i_block in range(len(block_ids)):
    block_id = block_ids[i_block]
    element_type = "HEX"
    number_of_attributes = 0
    exodus_file.put_elem_blk_info(block_id, element_type, number_of_elements_per_block[i_block], number_of_nodes_per_element[i_block], number_of_attributes)
    exodus_file.put_elem_blk_name(block_id, block_names[i_block])
    exodus_file.put_elem_connectivity(block_id, block_connectivity[i_block])

  for i_nodeset in range(number_of_node_sets):
    node_set_id = node_set_ids[i_nodeset]
    exodus_file.put_node_set_params(node_set_id, len(node_sets[i_nodeset]))
    exodus_file.put_node_set_name(node_set_id, node_set_names[i_nodeset])
    exodus_file.put_node_set(node_set_id, node_sets[i_nodeset])

  exodus_file.put_node_id_map(global_node_ids)
  exodus_file.put_elem_id_map(global_element_ids)

  exodus_file.set_global_variable_number(number_of_global_data)
  for i_name in range(len(global_data_names)):
    exodus_file.put_global_variable_name(global_data_names[i_name], i_name + 1)

  exodus_file.set_node_variable_number(number_of_node_data)
  for i_name in range(len(node_data_names)):
    exodus_file.put_node_variable_name(node_data_names[i_name], i_name + 1)

  exodus_file.set_element_variable_number(number_of_element_data)
  for i_name in range(len(element_data_names)):
    exodus_file.put_element_variable_name(element_data_names[i_name], i_name + 1)

  for i_time in range(len(times)):
    step = i_time + 1
    exodus_file.put_time(step, times[i_time])
    for i_data in range(len(global_data[i_time])):
      exodus_file.put_global_variable_values(global_data_names[i_data], step, global_variable_values[i_time][i_data])
    for i_data in range(len(node_data[i_time])):
      exodus_file.put_node_variable_values(node_data_names[i_data], step, node_data[i_time][i_data])
    for i_block in range(len(block_ids)):
      block_id = block_ids[i_block]
      for i_data in range(number_of_element_data):
        data_name = element_data_names[i_data]
        exodus_file.put_element_variable_values(block_id, data_name, step, element_data[i_time][i_block][i_data])

  exodus_file.close()

  return

if __name__ == "__main__":

    print "\n-- text_to_exodus --\n"

    usage_message = "Usage:  text_to_exodus.py <exodus_file.e> [-p] [num_partitions] [-convert_test_directory]\n"

    if len(sys.argv) != 2 and len(sys.argv) != 4:
        print usage_message
        sys.exit(1)

    convert_test_directory = False
    if len(sys.argv) == 2:
      if sys.argv[1] == "-convert_test_directory":
        convert_test_directory = True
      else:
        text_file_name = sys.argv[1]

    partitioned_mesh = False
    if len(sys.argv) == 4:
      num_partitions = int(sys.argv[3])
      if num_partitions > 1:
        partitioned_mesh = True

    if partitioned_mesh == False and convert_test_directory == False:
      ConvertTextFileToExodus(text_file_name)
    elif partitioned_mesh == True:
      print "Processing", num_partitions, "partitions"
      for i in range(num_partitions):
        partition_file_name = text_file_name + "." + str(num_partitions) + "." + str(i)
        print " ", partition_file_name
        ConvertTextFileToExodus(partition_file_name)
    elif convert_test_directory == True:
      # Convert all the .g.txt files in the test repository
      files_to_convert = []
      root_dir = os.getcwd()
      directory_tree = os.walk("./test/")
      for item in directory_tree:
        directory = item[0]
        file_names = item[2]
        for file_name in file_names:
          if file_name[-6:] == ".e.txt" or (file_name[-10:-7] == ".e." and file_name[-4:] == ".txt"):
            files_to_convert.append([directory, file_name])
      os.chdir(root_dir)
      for item in files_to_convert:
        directory = item[0]
        file_name = item[1]
        os.chdir(directory)
        ConvertTextFileToExodus(file_name)
        os.chdir(root_dir)
    else:
      print usage_message
      sys.exit(1)
