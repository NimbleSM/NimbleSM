#!/usr/bin/env python

# This script requires exodus.py (distributed with SEACAS)
import sys
sys.path.append('/path/to/seacas/install')
import exodus

def coincident(node_1, node_2, tol):

    for dim in range(3):
        if abs(node_1[dim] - node_2[dim]) > tol:
            return False

    return True

if __name__ == "__main__":

    tol = 1.0e-12

    if len(sys.argv) != 2:
        print "\nUsage:  check_periodic_mesh.py <genesis_file_name>\n"
        sys.exit(1)

    in_file_name = sys.argv[1]
    in_file = exodus.exodus(in_file_name, mode='r')
    num_nodes = in_file.num_nodes()
    x_coord, y_coord, z_coord = in_file.get_coords()
    in_file.close()

    # determine the bounding box
    big = 1.0e50
    x_min = big
    x_max = -big
    y_min = big
    y_max = -big
    z_min = big
    z_max = -big
    for i in range(num_nodes):
        x = x_coord[i]
        y = y_coord[i]
        z = z_coord[i]
        if x < x_min:
            x_min = x
        if x > x_max:
            x_max = x
        if y < y_min:
            y_min = y
        if y > y_max:
            y_max = y
        if z < z_min:
            z_min = z
        if z > z_max:
            z_max = z
    print "\nBounding box:"
    print "  x min:", x_min
    print "  x max:", x_max
    print "  y min:", y_min
    print "  y max:", y_max
    print "  z min:", z_min
    print "  z max:", z_max, "\n"

    min_x_surface = []
    max_x_surface = []
    min_y_surface = []
    max_y_surface = []
    min_z_surface = []
    max_z_surface = []

    for i in range(num_nodes):
        x = x_coord[i]
        y = y_coord[i]
        z = z_coord[i]
        if abs(x_coord[i] - x_min) < tol:
            min_x_surface.append((x, y, z))
        if abs(x_coord[i] - x_max) < tol:
            max_x_surface.append((x, y, z))
        if abs(y_coord[i] - y_min) < tol:
            min_y_surface.append((x, y, z))
        if abs(y_coord[i] - y_max) < tol:
            max_y_surface.append((x, y, z))
        if abs(z_coord[i] - z_min) < tol:
            min_z_surface.append((x, y, z))
        if abs(z_coord[i] - z_max) < tol:
            max_z_surface.append((x, y, z))

    print "Number of surface nodes:"
    print "  x min:", len(min_x_surface)
    print "  x max:", len(min_x_surface)
    print "  y min:", len(min_y_surface)
    print "  y max:", len(min_y_surface)
    print "  z min:", len(min_z_surface)
    print "  z max:", len(min_z_surface)

    print "\nChecking nodes on x face."
    for node_1 in min_x_surface:
        found_match = False
        for node_2 in max_x_surface:
            if coincident(node_1, node_2, tol):
                found_match = True
                continue
        if not found_match:
            print "Failed to match node on x face:", node_1
            sys.exit(1)

    print "Checking nodes on y face."
    for node_1 in min_y_surface:
        found_match = False
        for node_2 in max_y_surface:
            if coincident(node_1, node_2, tol):
                found_match = True
                continue
        if not found_match:
            print "Failed to match node on y face:", node_1
            sys.exit(1)

    print "Checking nodes on z face."
    for node_1 in min_z_surface:
        found_match = False
        for node_2 in max_z_surface:
            if coincident(node_1, node_2, tol):
                found_match = True
                continue
        if not found_match:
            print "Failed to match node on z face:", node_1
            sys.exit(1)

    print "\nMesh is periodic.\n"
