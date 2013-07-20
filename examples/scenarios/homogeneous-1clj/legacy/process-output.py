#!/usr/bin/python
# -*- coding: utf-8 -*-

# TODO: gnuplot?

import os
import re
import subprocess
import sys

# INFO: Simulating
# Simulating 10 steps
numTimeStepsLine = re.compile(r".*?Simulating (\d*) steps")

# INFO: Computation in main loop took: 
computationTimeLine = re.compile(r".*?INFO: Computation in main loop took: (-?\d+\.\d+)")

totalTimeLine = re.compile(r".*?main: used (-?\d+\.\d+) s")

def parse(ls1filename, datFileName):
	print "parsing ls1filename", ls1filename, " datFile", datFileName
	numTimeSteps = 0
	computationTime = 0
	totalTime = 0
	datfile = open(datFileName, 'a')
	with open(ls1filename) as f:
		for l in f:
			if numTimeSteps == 0:
				m = numTimeStepsLine.match(l)
				if m:
					numTimeSteps = float(m.group(1))
					print "numTimeSteps=", numTimeSteps
			m = computationTimeLine.match(l)
			if m:
				computationTime = float(m.group(1))
				print "computationTime=", computationTime
			m = totalTimeLine.match(l)
			if m:
				totalTime = float(m.group(1))
				print "totalTime=", totalTime
	if numTimeSteps == 0:
		print "Error: count == 0! Apparently there were no usable lines in the output log file.\n"
		sys.exit(1)
	datfile.write(str(totalTime) + "\n")
	datfile.close()
	
if __name__ == "__main__":
	#print "Hi there :)"
	if len(sys.argv) < 3:
		print "usage: process-output.py <ls1 output file> <data file>"
		sys.exit(1)
	parse(sys.argv[1], sys.argv[2])
