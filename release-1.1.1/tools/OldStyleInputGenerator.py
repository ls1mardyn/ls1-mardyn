#!/usr/bin/env python
#
# Input file generator for the Molecular Dynamics code LS1
#
# Author: Christoph Niethammer, 2009
# Mail:   christoph.niethammer@web.de
#

from math import *
from optparse import OptionParser
import random

parser = OptionParser()

# Default values
N=1000
T=0.7
Lx=97.0
Ly=97.0
Lz=97.0
t0=0.0
m=10000.
eps=1.0
sigm=1.0
rc=5.0
example="LJ"

parser.add_option("-N", action="store", type="int", dest="N", help="Number of particles")
parser.add_option("-m", action="store", type="float", dest="m", help="Particle mass for LJ only")
parser.add_option("-T", action="store", type="float", dest="T", help="Temperature")
parser.add_option("-n", action="store", type="float", dest="n", help="Particle density")
parser.add_option("-L", action="store", type="float", dest="L", help="Box length")
parser.add_option("-x", "--Lx", action="store", type="float", dest="Lx", help="Box length in x direction")
parser.add_option("-y", "--Ly", action="store", type="float", dest="Ly", help="Box length in y direction")
parser.add_option("-z", "--Lz", action="store", type="float", dest="Lz", help="Box length in z direction")
parser.add_option("-r", "--rc", action="store", type="float", dest="rc", help="Cutoff radius. Should match value from '*.cfg' file!")
parser.add_option("-t", "--example", action="store", type="string", dest="example", help="available options: LJ(default), EOX, MIX (=2CLJD and EOX mixed")

(options, args) = parser.parse_args()

if (options.N > 0):
  N=options.N
if (options.T > 0):
  T=options.T
if (options.L > 0):
  Lx=Ly=Lz=options.L
if (options.Lx > 0):
  Lx=options.Lx
if (options.Ly > 0):
  Lx=options.Ly
if (options.Lz > 0):
  Lx=options.Lz
if (options.n):
  Lx=Ly=Lz=pow((N/options.n),(1./3.))
if (options.m > 0):
  m=options.m
if (options.rc):
  rc=options.rc
if (options.example):
  example=options.example
else :
  example=""
  
vmax=sqrt((2.*T)/(3.*m))

print "MDProject trunk 20090721"
print "currentTime     %(t0).2f" % vars()
print "#rho   %f" % ( N/(Lx*Ly*Lz) )
print "# Vmax %f" % vmax
print "Temperature     %f" % T
print 'Length          %(Lx).4f %(Ly).4f %(Lz).4f' % vars()
if (example == "EOX"):
  print """NumberOfComponents      1
3 0 1 0 0
0.0000000000000000        0.0000000000000000        1.4681756272422370       1.59989999999999993E-002  1.96741662167147779E-004   5.8447340512780981   38.36144        0
-1.4739864075776508        0.0000000000000000      -0.83729029230229368       1.40269999999999995E-002  2.68352891066251406E-004   6.6643082884145430  38.36144        0
1.4739864075776508        0.0000000000000000      -0.83729029230229368       1.40269999999999995E-002  2.68352891066251406E-004   6.6643082884145430  38.36144        0
0.0000000000000000        0.0000000000000000       7.79229859719781925E-002   0.0000000000000000        0.0000000000000000      -1.00000000000000000       0.96744548351221815
0.0000000000000000        0.0000000000000000        0.0000000000000000
10000000000.0000000
""" % vars()
elif (example == "LJ"):
  print """NumberOfComponents      1
1 0 0 0 0
  0. 0. 0.       %(m).4f %(eps).4f %(sigm).4f %(rc).4f 0.0
  0. 0. 0.
1.e10""" % vars()
elif (example == "MIX"):
  print """NumberOfComponents      2
3 0 1 0 0
0.0000000000000000        0.0000000000000000        1.4681756272422370       1.59989999999999993E-002  1.96741662167147779E-004   5.8447340512780981   38.36144        0
-1.4739864075776508        0.0000000000000000      -0.83729029230229368       1.40269999999999995E-002  2.68352891066251406E-004   6.6643082884145430  38.36144        0
1.4739864075776508        0.0000000000000000      -0.83729029230229368       1.40269999999999995E-002  2.68352891066251406E-004   6.6643082884145430  38.36144        0
0.0000000000000000        0.0000000000000000       7.79229859719781925E-002   0.0000000000000000        0.0000000000000000      -1.00000000000000000       0.96744548351221815
0.0000000000000000        0.0000000000000000        0.0000000000000000
2 0 1 0 0
0.0 0.0 -2.245185	0.015035 0.000433822 6.59439    5.0    0
0.0 0.0 2.245185	0.015035 0.000433822 6.59439    5.0    0
0.0 0.0 0.0	0.0 0.0 1.0	-0.61537
0.0 0.0 0.0
0.99 0.98
10000000000.0000000
""" % vars()


print "NumberOfMolecules       " + `N`
if (example == "EOX" or example == "MIX"):
  print "MoleculeFormat  ICRVQD"
elif (example == "LJ"):
  print "MoleculeFormat  ICRV"

continuous_arc=True
cid=1
v2 = 0.
L=(Lx,Ly,Lz)

nmax=int(ceil(pow(N, 1./3.))) 

for n in range (N):
  
  if (example == "MIX"):
    cid = random.randint(1,2);

  myn=n
  nz=myn%nmax
  myn=n-nz
  ny=(myn/nmax)%nmax
  myn=myn-ny*nmax
  nx=myn/(nmax*nmax)
  
  rx=Lx/nmax*(0+nx)
  ry=Ly/nmax*(0+ny)
  rz=Lz/nmax*(0+nz)

  vx=vy=vz=0
  if (continuous_arc):
    a1=random.uniform(0.,2*pi)
    a2=random.uniform(0.,pi)
    vx=vmax*cos(a1)*sin(a2)
    vy=vmax*cos(a1)*cos(a2)
    vz=vmax*sin(a1)
  else:
    vx=random.choice((-vmax,vmax))
    vy=random.choice((-vmax,vmax))
    vz=random.choice((-vmax,vmax))

  q1=q2=q3=q4=0
  D1=D2=D3=0
  if(example == "EOX" or example == "MIX"):
    q1=random.choice((-0.5,0.5))
    q2=random.choice((-0.5,0.5))
    q3=random.choice((-0.5,0.5))
    q4=random.choice((-0.5,0.5))

  if(example == "EOX" or example == "MIX"):
    print "%i %i %f %f %f %f %f %f %f %f %f %f %f %f %f" % (n+1, cid, rx,ry,rz, vx,vy,vz, q1,q2,q3,q4, D1,D2,D3)
  elif(example == "LJ"):
    print "%i %i %f %f %f %f %f %f" % (n+1, cid, rx,ry,rz, vx,vy,vz)
