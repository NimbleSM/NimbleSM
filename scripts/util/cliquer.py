#!/usr/bin/env python

# This script requires exodus.py

import sys
import os
import string
import sets
import copy
import time
import random
from contextlib import contextmanager

sys.path.append('/path/to/seacas/install')
import exodus

class dev_null_output:
    def write(self, _): pass

@contextmanager
def run_silently():
    hold_stdout = sys.stdout
    sys.stdout = dev_null_output()
    yield
    sys.stdout = hold_stdout

@contextmanager
def time_operation(name):
    start = time.time()
    yield
    total = time.time() - start
    print name, "(", total, "seconds)"

def print_usage():
    print "\nUsage:  parallel_preprocessing.py <input_mesh_file_base_name> <output_mesh_file_base_name> <num_ranks>\n"

def read_io_info(args):
    return args[1], args[2], int(args[3])

def print_io_info(input_base_name, output_base_name, num_ranks):
    print "\n-- parallel_preprocessing.py --"
    print "\nGenesis input file base name =", input_base_name
    print "Genesis output file base name =", output_base_name
    print "Number of partitions = ", num_ranks

def genesis_files_from_base(basename, num_ranks):
    rank_string_size = len(str(num_ranks))
    filenames = []
    for i_rank in range(num_ranks):
        rank_string = str(i_rank)
        while len(rank_string) < rank_string_size:
            rank_string = "0" + rank_string
        suffix = ".g." + str(num_ranks) + "." + rank_string
        filenames.append(basename + suffix)
    return filenames

def get_node_id_map_from_file(filename):
    with run_silently():
        genesis_file = exodus.exodus(filename, 'r')
        node_id_map = genesis_file.get_node_id_map()
        genesis_file.close()
    return node_id_map

def copy_exodus_file(oldname, newname):
    with run_silently():
        old_database = exodus.exodus(oldname, 'r')
        new_database = old_database.copy(newname)
        old_database.close()
    return new_database

def map_clique_ids_to_nodes(all_global_node_ids, global_node_ids):
    # This is the intuitive thing to do:
    # ranks_sharing_each_node = dict.fromkeys(all_global_node_ids, [])
    # it breaks because rather than each dictionary element getting it's own array,
    # all dictionary elements end up with references to the same array.
    clique_id_lookup = {nodeid: 0 for nodeid in all_global_node_ids}
    clique_size_lookup = [0]
    clique_member_count = [len(all_global_node_ids)]
    clique_count = 1
    for node_ids in global_node_ids:
        clique_reassignment_dict = {}
        for node_id in node_ids:
            old_clique_id = clique_id_lookup[node_id]
            clique_member_count[old_clique_id] -= 1
            if old_clique_id in clique_reassignment_dict:
                new_clique_id = clique_reassignment_dict[old_clique_id]
                clique_id_lookup[node_id] = new_clique_id
                clique_member_count[new_clique_id] += 1
            else:
                new_clique_id = clique_count
                clique_size_lookup.append(clique_size_lookup[old_clique_id] + 1)
                clique_member_count.append(1)
                clique_reassignment_dict[old_clique_id] = new_clique_id
                clique_id_lookup[node_id] = new_clique_id
                clique_count = clique_count + 1

    def CliqueDowngradeGenerator():
        new_clique_num = 1
        for clique in xrange(clique_count):
            if clique_size_lookup[clique] > 1 and clique_member_count[clique] > 0:
                yield new_clique_num
                new_clique_num = new_clique_num + 1
            else:
                yield 0

    clique_downgrade_dict = list(CliqueDowngradeGenerator())
    clique_id_lookup = { nodeid: clique_downgrade_dict[cliqueid] for (nodeid, cliqueid) in clique_id_lookup.iteritems() }


    print "Total number of nodes = ", len(all_global_node_ids)
    print "Boundary node len = ", (len(clique_id_lookup) - clique_id_lookup.values().count(0))
    print "Number of cliques = ", (clique_count - clique_downgrade_dict.count(0))
    return clique_id_lookup;

def assign_clique_colors(node_clique_ids):
    # for each rank gets a set of the cliques on that rank
    rank_clique_ids = map(set, node_clique_ids)
    # removes the null clique
    for clique_ids in rank_clique_ids:
        clique_ids.remove(0)

    rank_clique_ids = map(list, rank_clique_ids)
    num_clique_ids = 1 + max(map(max, rank_clique_ids))

    # for each clique id gets the ranks with that clique
    clique_rank_ids = map(copy.copy,  num_clique_ids * [[]])
    for rank in range(len(rank_clique_ids)):
        for clique_id in rank_clique_ids[rank]:
            clique_rank_ids[clique_id].append(rank);

    # calculates the color number for each clique
    # the color number of a clique is the order it needs to be created
    # MPI_Comm_split initializes all cliques sharing the same color number
    # simultaneously, but it can only do that if no two cliques with the same
    # color number share a rank
    clique_color_number = num_clique_ids * [0]
    colors = []
    for clique in range(1, num_clique_ids):
        rank_list = clique_rank_ids[clique]
        found_color = False
        color_num = 1
        for color in colors:
            if color.isdisjoint(rank_list):
                color.update(rank_list)
                found_color = True
                break
            color_num += 1
        clique_color_number[clique] = color_num
        if not found_color:
            colors.append(set(rank_list))

    print 'color_number =', len(colors)
    return [[clique_color_number[c] for c in clique_ids_on_rank] for clique_ids_on_rank in node_clique_ids]


def main():
    if len(sys.argv) < 4:
        print_usage()
        sys.exit(1)

    with time_operation("Attaining lists of input files..."):
        (input_base_name, output_base_name, num_ranks) = read_io_info(sys.argv)

        print_io_info(input_base_name, output_base_name, num_ranks)

        input_files = genesis_files_from_base(input_base_name, num_ranks);
        output_files = genesis_files_from_base(output_base_name, num_ranks);

    with time_operation("Reading global node ids..."):
        global_node_ids = map(get_node_id_map_from_file, input_files)

    with time_operation("Finding the set of all global node ids..."):
        all_global_node_ids = set.union(*map(set, global_node_ids))

    with time_operation("Mapping node ids to cliques..."):
        clique_id_lookup = map_clique_ids_to_nodes(all_global_node_ids, global_node_ids)

    with time_operation("Transfering node ids to cliques..."):
        node_clique_ids = [[clique_id_lookup[node_id] for node_id in rank_id_list] for rank_id_list in global_node_ids]

    with time_operation("Calculating clique colors..."):
        node_clique_colors = assign_clique_colors(node_clique_ids)

    with time_operation("Cleaning existing genesis files..."):
        for file in output_files:
            if os.path.exists(file):
                os.remove(file)

    with time_operation("Writing modified genesis files..."):
        # write the modified genesis files
        for i_rank in range(len(input_files)):

            name = input_files[i_rank]
            new_name = output_files[i_rank]

            global_vars = []
            nodal_vars = ['clique', 'clique_color']
            element_vars = []

            #with run_silently():
            new_database = copy_exodus_file(name, new_name)
            print "DEBUGGING", name, i_rank, new_database.num_nodes(), new_database.num_times()
            #with run_silently():
            exodus.add_variables(new_database, global_vars, nodal_vars, element_vars)
            new_database.put_node_variable_values('clique', 1, node_clique_ids[i_rank])
            new_database.put_node_variable_values('clique_color', 1, node_clique_colors[i_rank])
            new_database.close()

    print "Finished"

#I separated it out because otherwise other functions could access valuables declared in the if block representing main.
if __name__ == "__main__":
    main()
