#!/usr/bin/env python
import os,sys
sys.path.append('./')
import exodus

def amax(vs):
  return max(max(vs),-min(vs))

##################################################################
class exoreader:
  sscale = 1.0e-9
  def __init__(self,fname):
    self.oname = fname.split(".")[0]+".dat"
    self.db = exodus.exodus(fname)
    self.times = list(self.db.get_times())
    self.nsteps = len(self.times)
    print("  {0} times {1}".format(self.nsteps, [ float("{0:7.4g}".format(t)) for t in self.times]))
#   o = open("t.dat","w")
#   for t in self.times: o.write("{0}\n".format(t))
#   o.close()
    self.nnodes = self.db.num_nodes()
    xs,ys,zs = self.db.get_coords()
    self.xs = xs
    self.Xs = [ xs[i] for i in range(0,self.nnodes,4) ]
#   o = open("x.dat","w")
#   for x in self.Xs: o.write("{0}\n".format(x))
#   o.close()
    print("  nodes {0} coordinate names {1}".format(self.nnodes,self.db.get_coord_names()))
    print("  nodal variables {0}".format(self.db.get_node_variable_names()))
    print("  element variables {0}".format(self.db.get_element_variable_names()))

    self.nexact = 0
    for name in self.db.get_node_variable_names():
      if name[:19] == "exact_displacement_" and name[-2:] == "_x":
        print(name)
        self.nexact += 1
    print(" number of exact trajectories {0}".format(self.nexact))
    self.ntraj = 0
    for name in self.db.get_node_variable_names():
      if name[:20] == "sample_displacement_" and name[-2:] == "_x":
        print(name)
        self.ntraj += 1
    print(" number of sample trajectories {0}".format(self.ntraj))
    
    self.blk = 1 # NOTE
    elem_conn, num_blk_elems, num_elem_nodes = self.db.get_elem_connectivity(self.blk)
    self.conn = []
    k = 0
    for i in range(num_blk_elems):
      c = []
      for j in range(num_elem_nodes):
        c.append(elem_conn[k]-1)
        k += 1
      self.conn.append(c)
    counts = self.nnodes*[0]
    for c in self.conn:
      for i in c:  counts[i] += 1
    self.weights = self.nnodes*[0.0]
    for i,c in enumerate(counts): 
      if c > 0: self.weights[i] = 1.0/c
    
  def process(self):
    o = open(self.oname,"w")
    o.write("# {0} sample {1} exact trajectories\n".format(self.ntraj,self.nexact))
    o.write("# x u f u_1 u_2 ... f_1 f_2 ....\n")
    maxs = 0.0
    for i,time in enumerate(self.times):
      step = i + 1
      us = self.db.get_node_variable_values("displacement_x",step)
      vs = self.db.get_node_variable_values("velocity_x",step)
      fs = self.db.get_node_variable_values("internal_force_x",step)
      wxs = [ self.db.get_node_variable_values("exact_displacement_"+str(j+1)+"_x",step) for j in range(self.nexact) ]
      fxs = [ self.db.get_node_variable_values("exact_force_"+str(j+1)+"_x",step) for j in range(self.nexact) ]
      wws = [ self.db.get_node_variable_values("sample_displacement_"+str(j+1)+"_x",step) for j in range(self.ntraj) ]
      ffs = [ self.db.get_node_variable_values("sample_force_"+str(j+1)+"_x",step) for j in range(self.ntraj) ]
      print("time: {0:10.8f} u:{1:9.6f}".format(time,amax(us)),end=" ")
      for ws in wws:
        print("{0:9.6f}".format(amax(ws)),end=" ")
      print("")
      stress = self.db.get_element_variable_values(self.blk,'stress_xx',step)
      for s in stress: maxs = max(s,maxs)
      ss = self.nnodes*[0.0]
      for j,c in enumerate(self.conn):
        for k in c:
          ss[k-1] += stress[j]
      for j,w in enumerate(self.weights): ss[j] *= w
      o.write("# time {0:12.9f}\n".format(time))
      for j in range(0,self.nnodes,4):
        o.write("{0:7.4f} {1:9.6f} {2:9.6f} ".format(self.xs[j],us[j],fs[j]))
        for ws in wxs: o.write("{0:9.6f} ".format(ws[j]))
        for fs in fxs: o.write("{0:9.6f} ".format(fs[j]))
        for ws in wws: o.write("{0:9.6f} ".format(ws[j]))
        for fs in ffs: o.write("{0:9.6f} ".format(fs[j]))
        o.write("\n")
      o.write("\n\n")
    o.write("# {0} steps\n".format(len(self.times)))
    o.close()
    self.db.close()

##################################################################
if __name__ == "__main__":
##################################################################
  if len(sys.argv) < 2:
    print("!! usage: extract.py <db.e> !!")
    exit()
  fname = sys.argv[1]
  db = exoreader(fname)
  db.process()
