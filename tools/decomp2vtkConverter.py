#!/usr/bin/python

import sys
import StringIO
print ''
if len(sys.argv) < 3:
    print 'ERROR: Not enough arguments given, need at least:'
    print '[1] input file name'
    print '[2] output file name'
    
inputFileName = str(sys.argv[1])
outputFileName = str(sys.argv[2])

print 'Input file name:', inputFileName
print 'Output file name:', outputFileName

if inputFileName.find('.decomp') == -1:
    print 'ERROR: wrong input file format, should be *.decomp' 
    sys.exit()
if outputFileName.find('.vtk') == -1:
    print 'ERROR: wrong output file format, should be *.vtk' 
    sys.exit()

f_in = open(inputFileName, 'r')
f_out = open(outputFileName, 'w')

l = f_in.readline()
while l.find('decompData Regions') != -1:
    if l != '':
        print 'ERROR: reached end of file, but require "decompData Regions"'
        sys.exit()
    print l
    l = f_in.readline()
f_in.readline()  # skip decompData Regions

class Cell:
    def __init__(self, minmaxline, cellid):
        parts = minmaxline.split(' ')
        self.xmin = parts[0]
        self.ymin = parts[1]
        self.zmin = parts[2]
        self.xmax = parts[3]
        self.ymax = parts[4]
        self.zmax = parts[5]
        self.cellid = cellid
    def getPointCoordinates(self):
        output = StringIO.StringIO()
        print >> output, self.xmin, self.ymin, self.zmin
        print >> output, self.xmin, self.ymin, self.zmax
        print >> output, self.xmin, self.ymax, self.zmin
        print >> output, self.xmin, self.ymax, self.zmax
        print >> output, self.xmax, self.ymin, self.zmin
        print >> output, self.xmax, self.ymin, self.zmax
        print >> output, self.xmax, self.ymax, self.zmin
        print >> output, self.xmax, self.ymax, self.zmax
        returnstring = output.getvalue()
        output.close()
        return returnstring
    def getCellPointIndices(self):
        output = StringIO.StringIO()
        for i in range(8):
            print >> output, str(i + 8 * self.cellid),
        print >> output, ''
        
        returnstring = output.getvalue()
        output.close()
        return returnstring
        
numcells = 0
l = f_in.readline()
celllist = []
while l.find('particleData') == -1 and l != "":
    celllist.append(Cell(l.rstrip(), numcells))
    numcells += 1
    l = f_in.readline()

print 'Found', numcells, 'cells.'
f_out.write('# vtk DataFile Version 2.0\n\
MarDyn decomposition output\n\
ASCII\n\
DATASET UNSTRUCTURED_GRID\n')

f_out.write('POINTS ' + str(8 * numcells) + ' float\n')  # 8 points per cell
for cell in celllist:
    f_out.write(cell.getPointCoordinates())

f_out.write('CELLS ' + str(numcells) + ' ' + str(numcells * 9) + '\n')  # CELLS requires numcells and number of total list entries (=numcells + 8*numcells) 
for cell in celllist:
    f_out.write('8 ' + cell.getCellPointIndices())

f_out.write('CELL_TYPES ' + str(numcells) + '\n')
for cell in celllist:
    f_out.write('11\n')    

f_in.close()
f_out.close()
