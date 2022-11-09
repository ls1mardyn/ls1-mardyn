#! /usr/bin/env python
# -*- coding: utf-8 -*-

## pm2xml.py
## convert the ITT pm files to the new XML format
## Martin Bernreuther <bernreuther@hlrs.de>

##**************************************************************************
##   Copyright (C) 2005 - 2011 by Martin Bernreuther                       *
##   author: Martin Bernreuther <bernreuther@hlrs.de>                      *
##                                                                         *
##   This program is free software; you can redistribute it and/or modify  *
##   it under the terms of the GNU General Public License as published by  *
##   the Free Software Foundation; either version 2 of the License, or     *
##   (at your option) any later version.                                   *
##                                                                         *
##   This program is distributed in the hope that it will be useful,       *
##   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
##   GNU General Public License for more details.                          *
##                                                                         *
##   You should have received a copy of the GNU General Public License     *
##   along with this program; if not, write to the                         *
##   Free Software Foundation, Inc.,                                       *
##   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
##**************************************************************************/

import sys
import numpy as np

def readKeyword(inpfhdl, keywd):
	for line in inpfhdl:
		tokens=str.split(line, '=')
		if len(tokens)>0 and str.strip(tokens[0])==str.strip(keywd):
			if len(tokens)>1:
				return str.strip(tokens[1])
			else:
				return ""
	return None

def readKeywords(inpfhdl, keywds):
	result={}
	for line in inpfhdl:
		tokens=str.split(line, '=')
		if len(tokens)==0: return result
		tokens[0]=str.strip(tokens[0])
		for k in keywds:
			if tokens[0]==str.strip(k):
				if len(tokens)>1:
					result[tokens[0]]=str.strip(tokens[1])
				else:
					result[tokens[0]]=None
				keywds.remove(k)
				if len(keywds)==0: return result
				break
	return result

def processInput():

	sites=[]
	numlj126=0
	numdipole=0
	numquadrupole=0
	totalmass=0.
	smx=0.
	smy=0.
	smz=0.

	nsitetypes=int(readKeyword(inpfhdl, "NSiteTypes"))
	for st in range(nsitetypes):
		props=readKeywords(inpfhdl, ["SiteType", "NSites"])
		sitetype=props["SiteType"]
		nsites=int(props["NSites"])
		if sitetype=="LJ126":
			for s in range(nsites):
				props=readKeywords(inpfhdl, ["x", "y", "z", "sigma", "epsilon", "mass"])
				x=float(props["x"])
				y=float(props["y"])
				z=float(props["z"])
				mass=float(props["mass"])
				sites.append(SiteLJ126(x, y, z, mass, float(props["sigma"]), float(props["epsilon"])))
				numlj126+=1
				totalmass+=mass
				smx+=mass*x
				smy+=mass*y
				smz+=mass*z
		elif sitetype=="Dipole":
			for s in range(nsites):
				props=readKeywords(inpfhdl, ["x", "y", "z", "theta", "phi", "dipole", "mass", "shielding"])
				x=float(props["x"])
				y=float(props["y"])
				z=float(props["z"])
				mass=float(props["mass"])
				sites.append(SiteDipole(x, y, z, mass, float(props["theta"]), float(props["phi"]), float(props["dipole"]), float(props["shielding"])))
				numdipole+=1
				totalmass+=mass
				smx+=mass*x
				smy+=mass*y
				smz+=mass*z
		elif sitetype=="Quadrupole":
			for s in range(nsites):
				props=readKeywords(inpfhdl, ["x", "y", "z", "theta", "phi", "quadrupole", "mass", "shielding"])
				x=float(props["x"])
				y=float(props["y"])
				z=float(props["z"])
				mass=float(props["mass"])
				sites.append(SiteQuadrupole(x, y, z, mass, float(props["theta"]), float(props["phi"]), float(props["quadrupole"]), float(props["shielding"])))
				numquadrupole+=1
				totalmass+=mass
				smx+=mass*x
				smy+=mass*y
				smz+=mass*z


	## adjusting coordinates to the center of gravity
	cx=smx/totalmass
	cy=smy/totalmass
	cz=smz/totalmass
	
	## moving center of gravity [cx, cy, cz] to origin
	
	A=np.zeros((3,3),float)
	for s in sites:
		s.x-=cx
		s.y-=cy
		s.z-=cz
		A[0, 0]+=(s.y*s.y+s.z*s.z)*s.mass
		A[0, 1]-=s.x*s.y*s.mass
		A[0, 2]-=s.x*s.z*s.mass
		A[1, 1]+=(s.x*s.x+s.z*s.z)*s.mass
		A[1, 2]-=s.y*s.z*s.mass
		A[2, 2]+=(s.x*s.x+s.y*s.y)*s.mass
	A[1,0] = A[0,1]
	A[2,0] = A[0,2]
	A[2,1] = A[1,2]
	(eigval,eigvec)=np.linalg.eig(A)	## Eigenvalues and row-wise Eigenvectors


	## rotating to principal axes

	A=np.zeros((3,3),float)
	for s in sites:
		x=s.x
		y=s.y
		z=s.z
		s.x=eigvec[0, 0]*x+eigvec[0, 1]*y+eigvec[0, 2]*z
		s.y=eigvec[1, 0]*x+eigvec[1, 1]*y+eigvec[1, 2]*z
		s.z=eigvec[2, 0]*x+eigvec[2, 1]*y+eigvec[2, 2]*z
		A[0, 0]+=(s.y*s.y+s.z*s.z)*s.mass
		A[0, 1]-=s.x*s.y*s.mass
		A[0, 2]-=s.x*s.z*s.mass
		A[1, 1]+=(s.x*s.x+s.z*s.z)*s.mass
		A[1, 2]-=s.y*s.z*s.mass
		A[2, 2]+=(s.x*s.x+s.y*s.y)*s.mass
	A[1,0] = A[0,1]
	A[2,0] = A[0,2]
	A[2,1] = A[1,2]
	I=A.diagonal()

	print("<components>")
	print("  <moleculetype id=\"1\">")
	for s in sites:
		s.prntXML()
	print("    <momentsofinertia rotaxes=\"xyz\" >")
	print("      <Ixx>%.16g</Ixx>"%I[0])
	print("      <Iyy>%.16g</Iyy>"%I[1])
	print("      <Izz>%.16g</Izz>"%I[2])
	print("    </momentsofinertia>")
	print("  </moleculetype>")
	print("</components>")


class Site:
	"molecule site"
	
	x=0.
	y=0.
	z=0.
	mass=0.
	
	def __init__(self, x, y, z, mass):
		self.x=x
		self.y=y
		self.z=z
		self.mass=mass
	
	def prnt(self):
		print("Site [%g,%g,%g]; mass: %g"%(self.x, self.y, self.z, self.mass))

class SiteLJ126(Site):
	"Lennard Jones 12-6 site"
	
	sigma=0.
	epsilon=0.
	
	def __init__(self, x, y, z, mass, sigma, epsilon):
		Site.__init__(self, x, y, z, mass)
		self.sigma=sigma
		self.epsilon=epsilon
	
	def prnt(self):
		print("LJ126 [%.16g,%.16g,%.16g]; mass: %.16g; sigma: %.16g; epsilon: %.16g"%(self.x, self.y, self.z, self.mass, self.sigma, self.epsilon))
	
	def prntXML(self):
		print(" <site type=\"LJ126\" >")
		print("  <coord> <x>%.16g</x> <y>%.16g</y> <z>%.16g</z> </coord>"%(self.x,self.y,self.z))
		print("  <mass>%.16g</mass>"%self.mass)
		print("  <sigma>%.16g</sigma>"%self.sigma)
		print("  <epsilon>%.16g</epsilon>"%self.epsilon)
		print(" </site>")

class SiteDipole(Site):
	"Dipole site"
	
	theta=0.
	phi=0.
	dipole=0.
	shielding=0.
	
	def __init__(self, x, y, z, mass, theta, phi, dipole, shielding):
		Site.__init__(self, x, y, z, mass)
		self.theta=theta
		self.phi=phi
		self.dipole=dipole
		self.shielding=shielding
	
	def prnt(self):
		print("Dipole [%.16g,%.16g,%.16g]; mass: %.16g; theta,phi: %.16g,%.16g; dipole,shielding: %.16g,%.16g"%(self.x, self.y, self.z, self.mass, self.theta, self.phi, self.dipole, self.shielding))

	def prntXML(self):
		print(" <site type=\"Dipole\" >")
		print("  <coord> <x>%.16g</x> <y>%.16g</y> <z>%.16g</z> </coord>"%(self.x,self.y,self.z))
		print("  <orientation theta=\"%.16g\" phi=\"%.16g\" />"%(self.theta,self.phi))
		print("  <mass>%.16g</mass>"%self.mass)
		print("  <shielding>%.16g</shielding>"%self.shielding)
		print("  <dipole>%.16g</dipole>"%self.dipole)
		print(" </site>")

class SiteQuadrupole(Site):
	"Quadrupole site"
	
	theta=0.
	phi=0.
	quadrupole=0.
	shielding=0.
	
	def __init__(self, x, y, z, mass, theta, phi, quadrupole, shielding):
		Site.__init__(self, x, y, z, mass)
		self.theta=theta
		self.phi=phi
		self.quadrupole=quadrupole
		self.shielding=shielding
	
	def prnt(self):
		print("Quadrupole [%.16g,%.16g,%.16g]; mass: %.16g; theta,phi: %.16g,%.16g; quadrupole,shielding: %.16g,%.16g"%(self.x, self.y, self.z, self.mass, self.theta, self.phi, self.quadrupole, self.shielding))

	def prntXML(self):
		print(" <site type=\"Quadrupole\" >")
		print("  <coord> <x>%.16g</x> <y>%.16g</y> <z>%.16g</z> </coord>"%(self.x,self.y,self.z))
		print("  <orientation theta=\"%.16g\" phi=\"%.16g\" />"%(self.theta,self.phi))
		print("  <mass>%.16g</mass>"%self.mass)
		print("  <shielding>%.16g</shielding>"%self.shielding)
		print("  <quadrupole>%.16g</quadrupole>"%self.quadrupole)
		print(" </site>")


if len(sys.argv)>1:
	try:
		with open(sys.argv[1], 'r') as inpfhdl:
			processInput()
	except IOError as e:
		(errno, strerror) = e.args
		print("Error opening %s (%s): %s" % (sys.argv[1],errno, strerror))
else:
	inpfhdl=sys.stdin
	processInput()
