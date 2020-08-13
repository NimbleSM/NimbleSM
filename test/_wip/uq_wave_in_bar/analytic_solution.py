#!/usr/bin/env python
import os,sys
from math import pi, sqrt, sin, cos

# 1/2 (Cos[a - b] - Cos[a + b]) == Sin[a] Sin[b]

if __name__ == "__main__":
  density = 7.8     # g/cm^3
  modulus = 3.0e12  # dyne/cm^2 <<< NOTE what modulus
  initialVelocity = 1000.0 # cm/sec
 
  args = sys.argv[1:]

  if len(args) > 0:
    xfile = args[0]
    xs = [ float(line.split()[0]) for line in open(xfile).readlines() ]
    x0 = xs[-1]
    length = xs[-1]-xs[0]
  else:
    length  = 2.0     # cm
    nnodes = 101
    dx = 0.02
    x0 = 0.0
    xs = [x0+i*dx for i in range(nnodes)];

  if len(args) > 1:
    tfile = args[1]
    times = [ float(line.split()[0]) for line in open(tfile).readlines() ]
    nsteps = len(times)
  else:
    nsteps = 101
    tf = 1.0e-5 # sec
    dt = tf/(nsteps-1)
    times = [i*dt for i in range(nsteps)];

  nterms = 1000
  c = sqrt(modulus/density)
  k = pi/(2.0*length)
  w = k*c
  A = 8.0*length*initialVelocity/(c*pi*pi)
  
  o = open("solution.dat",'w')
  for j,t in enumerate(times):
    print(">> {0:4d}/{1:4d} time {2:8.4g}".format(j+1,nsteps,t))
    o.write("# time {0:9.6g}\n".format(t))
    for x in xs:
      z = x - x0
      u = 0.0
      v = 0.0
      for i in range(nterms):
        n = i + 1
        m = 2*n-1
        wn = m*w
        kn = m*k
        An = A/(m*m)
        u = u - An*sin(wn*t)*sin(kn*z) 
        v = v - An*cos(wn*t)*sin(kn*z)*wn 
      o.write("{0:7.4f} {1:9.6g} {2:9.6g} \n".format(x,u,v))
    o.write("\n\n")
  o.close()
