#!/usr/bin/env python

import os
import sys

# This script requires exodus.py (distributed with SEACAS)
import sys
sys.path.append('/path/to/seacas/install')
import exodus

import string

if __name__ == "__main__":

    print "\n-- compute_volume_averaged_rve_properties --\n"

    usage_message = "Usage:  compute_volume_averaged_rve_properties.py <input_deck.in> <rve_exodus_file.e>\n"

    if len(sys.argv) != 3:
        print usage_message
        sys.exit(1)

    input_deck_name = sys.argv[1]
    print "\nReading input deck:", input_deck_name
    input_deck = open(input_deck_name)
    input_deck_lines = input_deck.readlines()
    input_deck.close()

    material_density = {}
    material_bulk_modulus = {}
    material_shear_modulus = {}
    rve_block_material = {}
    for line in input_deck_lines:
        if "material parameters" in line:
            vals = string.splitfields(line, ":")
            mat_props = string.splitfields(vals[-1])
            mat_name = mat_props[0]
            for i in range(len(mat_props)):
                entry = mat_props[i]
                if entry == "density":
                    material_density[mat_name] = float(mat_props[i+1])
                if entry == "bulk_modulus":
                    material_bulk_modulus[mat_name] = float(mat_props[i+1])
                if entry == "shear_modulus":
                    material_shear_modulus[mat_name] = float(mat_props[i+1])
        if "microscale block" in line:
            vals = string.splitfields(line, ":")
            props = string.splitfields(vals[-1])
            rve_block_name = props[0]
            rve_material_name = props[1]
            rve_block_material[rve_block_name] = rve_material_name

    rve_exodus_file_name = sys.argv[2]
    print "Reading RVE exodus file:", rve_exodus_file_name
    exodus_file = exodus.exodus(rve_exodus_file_name, mode='r')

    elem_var_names = exodus_file.get_element_variable_names()
    if "volume" not in elem_var_names:
        print "\n**** Error:  Element variable \"volume\" not found!\n"
        sys.exit(1)

    step = 1
    block_ids = exodus_file.get_elem_blk_ids()
    block_names = exodus_file.get_elem_blk_names()
    block_volume = {}
    total_volume = 0.0
    for block_id in block_ids:
        block_volume[block_id] = 0.0
        volumes = exodus_file.get_element_variable_values(block_id, "volume", step)
        for volume in volumes:
            block_volume[block_id] += volume
            total_volume += volume

    print "\nTotal RVE volume", total_volume, "\n"
    for key in block_volume.keys():
        print "Block", key, "volume", block_volume[key], "fraction", block_volume[key]/total_volume
    print

    vol_ave_density = 0.0
    vol_ave_bulk_modulus = 0.0
    vol_ave_shear_modulus = 0.0
    for i in range(len(block_ids)):
        block_id = block_ids[i]
        block_name = block_names[i]
        fraction = block_volume[block_id]/total_volume
        material_name = rve_block_material[block_name]
        vol_ave_density += fraction * material_density[material_name]
        vol_ave_bulk_modulus += fraction * material_bulk_modulus[material_name]
        vol_ave_shear_modulus += fraction * material_shear_modulus[material_name]

    print "Volume averaged material properties:"
    print "  Density", vol_ave_density
    print "  Bulk modulus", vol_ave_bulk_modulus
    print "  Shear modulus", vol_ave_shear_modulus
    print
