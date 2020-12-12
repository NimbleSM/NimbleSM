#!/usr/bin/env python
import os,sys
sys.path.append('../../tools/')
import exodus

fname = sys.argv[1]
oname = fname.split(".")[0]+".dat"
db = exodus.exodus(fname)
times = list(db.get_times())
nsteps = len(times)
print("  {0} times {1}".format(nsteps, [ float("{0:7.4g}".format(t)) for t in times]))
o = open("t.dat","w")
for t in times: o.write("{0}\n".format(t))
o.close()
blks = list(db.get_elem_blk_ids())
nelms = db.num_elems()
nnodes = db.num_nodes()
print("  blocks {0} elements {1}".format(blks,nelms))
print("  nodes {0} coordinate names {1}".format(nnodes,db.get_coord_names()))
print("  nodal variables {0}".format(db.get_node_variable_names()))

xs,ys,zs = db.get_coords()
Xs = [ xs[i] for i in range(0,nnodes,4) ]
#nnodes = nnodes//4
o = open("x.dat","w")
for x in Xs: o.write("{0}\n".format(x))
o.close()

ntraj = 3
o = open(oname,"w")
for i in range(nsteps):
  step = i + 1
  us = db.get_node_variable_values("displacement_x",step)
  vs = db.get_node_variable_values("velocity_x",step)
  wws = [ db.get_node_variable_values("off_nom_displacement_"+str(j+1)+"_x",step) for j in range(ntraj) ]
  o.write("# time {0}\n".format(times[i]))
  for j in range(0,nnodes,4):
    o.write("{0:9.4f} {1:12.8e} {2:9.4f} ".format(xs[j],us[j],vs[j]))
    for ws in wws:
      o.write("{0:12.8e} ".format(ws[j]))
    o.write("\n")
  o.write("\n\n")
o.close()

db.close()
