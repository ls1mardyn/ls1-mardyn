#!/usr/bin/env python
# -*- coding: utf-8 -*-
## printMPIrestart.py
## print binary MPIrestart file content in a ASCII/UTF8 format (e.g. for debugging purposes)
## Martin Bernreuther <bernreuther@hlrs.de>, 2015

from __future__ import print_function

import sys
if sys.hexversion < 0x02070000:	sys.stderr.write("WARNING: running an old Python version <2.7: {0}\n".format(sys.version.replace("\n","\t")))

import argparse
import struct

import os
import time

import signal
# kill program silently, if SIGINT is received (e.g. due to Ctrl-C keyboard input)
signal.signal(signal.SIGINT, signal.SIG_DFL)
# kill program silently, if SIGPIPE is received (e.g. if stdout is piped to head)
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


argparser=argparse.ArgumentParser(description="print MPIrestart file content")
argparser.add_argument('inpfile', nargs='?', type=argparse.FileType('rb'), default=sys.stdin, help="ICS input file (default: stdin)")
argparsergroup_print = argparser.add_mutually_exclusive_group()
argparsergroup_print.add_argument("--onlyheader", action="store_true", help="print only header information")
argparsergroup_print.add_argument("--onlydata", action="store_true", help="print only data")
args=argparser.parse_args()

printheader=True
printdata=True

if args.onlyheader: printdata=False
if args.onlydata: printheader=False

inpfile_position=0
inpfile_buffer=""
inpfile_posrange=[0,0]

def read_inpfile_chunk(length):
	global args,inpfile_position,inpfile_buffer,inpfile_posrange
	inpfile_posrange[0]=inpfile_position
	inpfile_buffer=args.inpfile.read(length)
	if inpfile_buffer: inpfile_position+=length
	inpfile_posrange[1]=inpfile_position
	return inpfile_buffer

def read_inpfile_struct(format):
	length=struct.calcsize(format)
	datastring=read_inpfile_chunk(length)
	if datastring:
		data=struct.unpack(format,datastring)
	else:
		data=None
	return data

def read_inpfile_struct0(format):
	data=read_inpfile_struct(format)
	if data is not None and len(data)>=1: data=data[0]
	return data

def read_inpfile_endofstring(eos="\00"):
	global args,inpfile_position,inpfile_buffer,inpfile_posrange
	inpfile_buffer=""
	inpfile_posrange[0]=inpfile_position
	while args.inpfile:
		readchar=args.inpfile.read(1)
		if readchar:
			inpfile_position+=1
			if readchar==eos: break
		else:
			break
		inpfile_buffer+=readchar
		inpfile_posrange[1]=inpfile_position
	return inpfile_buffer


print("# {0}".format(args.inpfile.name))
#if args.inpfile.name!="<stdin>":
try:
	fileinfo = os.stat(args.inpfile.name)
	print(time.strftime("# modification time:\t%d.%m.%Y %H:%M:%S", time.localtime(fileinfo.st_mtime)))
	print("# file size        :\t{0} Bytes".format(fileinfo.st_size))
except OSError as e:
	#print ("Error {0} stat {1}: {2}".format(e.errno, args.inpfile.name, e.strerror))
	pass

magicVersion = read_inpfile_chunk(64-8-4).rstrip("\00")
if printheader: print("[{0},{1}]\tmagic versionstring:\t{2}".format(inpfile_posrange[0],inpfile_posrange[1]-1,magicVersion))

endiannesschar='<'
if magicVersion=="MarDyn20140817" or magicVersion=="MarDyn20140819" or magicVersion=="MarDyn20150211trunk":
	
	if magicVersion=="MarDyn20150211trunk":
		endiannesstest=read_inpfile_struct0(">i")
		if endiannesstest==0x0a0b0c0d:
			# big endian
			endiannesschar='>'
			if printheader: print("[{0},{1}]\tendiannesstest: big endian".format(inpfile_posrange[0],inpfile_posrange[1]-1))
		elif endiannesstest==0x0d0c0b0a:
			# little endian
			endiannesschar='<'
			if printheader: print("[{0},{1}]\tendiannesstest: little endian".format(inpfile_posrange[0],inpfile_posrange[1]-1))
		else:
			if printheader: print("WARNING: unknown endianness - assuming little endian")
	else:
		read_inpfile_chunk(4)
	gap=read_inpfile_struct0(endiannesschar+"Q")
	if printheader:
		print("[{0},{1}]\tgap:\t{2}".format(inpfile_posrange[0],inpfile_posrange[1]-1,gap))
		
		datalayout=read_inpfile_endofstring()
		print("[{0},{1}]\tdata layout:\t{2}".format(inpfile_posrange[0],inpfile_posrange[1]-1,datalayout))
		
		token=read_inpfile_endofstring()
		if token=="BB":
			print("[{0},{1}]\ttoken {2} found".format(inpfile_posrange[0],inpfile_posrange[1]-1,token))
			numbb=read_inpfile_struct0(endiannesschar+"Q")
			print("[{0},{1}]\tnumber of bounding boxes:\t{2}".format(inpfile_posrange[0],inpfile_posrange[1]-1,numbb))
			nummoleculessum=0
			#minstartidx=None
			#maxlastidx=None
			for i in range(numbb):
				pos=inpfile_position
				minx=read_inpfile_struct0(endiannesschar+"d")
				miny=read_inpfile_struct0(endiannesschar+"d")
				minz=read_inpfile_struct0(endiannesschar+"d")
				maxx=read_inpfile_struct0(endiannesschar+"d")
				maxy=read_inpfile_struct0(endiannesschar+"d")
				maxz=read_inpfile_struct0(endiannesschar+"d")
				startidx=read_inpfile_struct0(endiannesschar+"Q")
				#if minstartidx is None or startidx<minstartidx: minstartidx=startidx
				nummolecules=read_inpfile_struct0(endiannesschar+"Q")
				#lastidx=startidx+nummolecules-1
				#if maxlastidx is None or lastidx>maxlastidx: maxlastidx=lastidx
				nummoleculessum+=nummolecules
				print("[{0},{1}]\tbb{2}:\t{3},{4},{5}\t{6},{7},{8}\t{9}\t{10}".format(pos,inpfile_posrange[1]-1,i,minx,miny,minz,maxx,maxy,maxz,startidx,nummolecules))
		else:
			print("ERROR: read token {0} instead of \"BB\"".format(token))
		print("# {0} molecules altogether".format(nummoleculessum))
	else:
		read_inpfile_chunk(gap)
	if printdata:
		i=0
		while args.inpfile:
			pos=inpfile_position
			id=read_inpfile_struct0(endiannesschar+"Q")
			if id is None: break
			#componentid=read_inpfile_struct0("<i")
			#read_inpfile_chunk(4)
			componentid=read_inpfile_struct0(endiannesschar+"ixxxx")
			rx=read_inpfile_struct0(endiannesschar+"d")
			ry=read_inpfile_struct0(endiannesschar+"d")
			rz=read_inpfile_struct0(endiannesschar+"d")
			vx=read_inpfile_struct0(endiannesschar+"d")
			vy=read_inpfile_struct0(endiannesschar+"d")
			vz=read_inpfile_struct0(endiannesschar+"d")
			qw=read_inpfile_struct0(endiannesschar+"d")
			qx=read_inpfile_struct0(endiannesschar+"d")
			qy=read_inpfile_struct0(endiannesschar+"d")
			qz=read_inpfile_struct0(endiannesschar+"d")
			Dx=read_inpfile_struct0(endiannesschar+"d")
			Dy=read_inpfile_struct0(endiannesschar+"d")
			Dz=read_inpfile_struct0(endiannesschar+"d")
			i+=1
			print("[{0},{1}]\tm{2}:\t{3}\t{4}\t{5},{6},{7}\t{8},{9},{10}\t{11};{12},{13},{14}\t{15},{16},{17}".format(pos,inpfile_posrange[1]-1,i,id,componentid,rx,ry,rz,vx,vy,vz,qw,qx,qy,qz,Dx,Dy,Dz))
else:
	print("Unknown Version {0}".format(magicVersion))

args.inpfile.close()

