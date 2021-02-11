#!/usr/bin/python3

import sys
import io
import re

################################# Documentation #################################

# This script turns a list of .decomp files into .vtk files.
# The file format is expected to be like this:
# xmin ymin zmin xmax ymax zmax {configKey1: configVal1 , configKeyN: configValN}
#
# Files like this can be generated via the DecompWriter when using the
# general domain decomposiotion

############################# classes and functions #############################
class Cell:
    def __init__(self, minMaxString, cellid, configString):
        parts = minMaxString.split(' ')
        self.extents = {
                'xmin' : parts[0],
                'ymin' : parts[1],
                'zmin' : parts[2],
                'xmax' : parts[3],
                'ymax' : parts[4],
                'zmax' : parts[5],
                }
        self.cellid = cellid
        self.config = {} # init dict
        # extract configuration data which is quasi-json
        configString = re.sub('[{}]', '', configString)
        for elem in configString.split(' , ') :
            [key, value] = elem.split(': ')
            self.config[key] = value
        
    def getPointCoordinates(self):
        output = io.StringIO()
        print(self.extents['xmin'], self.extents['ymin'], self.extents['zmin'], file=output)
        print(self.extents['xmin'], self.extents['ymin'], self.extents['zmax'], file=output)
        print(self.extents['xmin'], self.extents['ymax'], self.extents['zmin'], file=output)
        print(self.extents['xmin'], self.extents['ymax'], self.extents['zmax'], file=output)
        print(self.extents['xmax'], self.extents['ymin'], self.extents['zmin'], file=output)
        print(self.extents['xmax'], self.extents['ymin'], self.extents['zmax'], file=output)
        print(self.extents['xmax'], self.extents['ymax'], self.extents['zmin'], file=output)
        print(self.extents['xmax'], self.extents['ymax'], self.extents['zmax'], file=output)
        returnstring = output.getvalue()
        output.close()
        return returnstring

    def getCellPointIndices(self):
        output = io.StringIO()
        for i in range(8):
            print(str(i + 8 * self.cellid), file=output, end=' ')
        print('', file=output)

        returnstring = output.getvalue()
        output.close()
        return returnstring

################################ start of script ################################

#### input parsing
print('')
if len(sys.argv) < 2:
    print('ERROR: Not enough arguments given, need at least:')
    print('[1] input file name')

for inputFileName in sys.argv[1:]:
    outputFileName = inputFileName.replace('.decomp', '.vtk')

    if inputFileName.find('.decomp') == -1:
        print('ERROR: wrong input file format, should be *.decomp')
        sys.exit()

    f_in = open(inputFileName, 'r')
    f_out = open(outputFileName, 'w')
    ### data validation
    currentLine = f_in.readline()
    while currentLine.find('decompData Regions') != -1:
        if currentLine != '':
            print('ERROR: reached end of file, but require "decompData Regions"')
            sys.exit()
        currentLine = f_in.readline()
    f_in.readline()  # skip decompData Regions

    ### data parsing
    numcells = 0
    line = f_in.readline()
    celllist = []
    while line.find('particleData') == -1 and line != "":
        match=re.search(r'([0-9. ]+)(.*)', line)
        minMaxString=match.group(1).rstrip()
        configString=match.group(2)

        celllist.append(Cell(minMaxString, numcells, configString))
        numcells += 1
        line = f_in.readline()

    print(inputFileName + ' -> ' + outputFileName + ' (' + str(numcells) + ' cells)')
    # Write VTK header and dataset information
    f_out.write('# vtk DataFile Version 2.0\n'
            + 'MarDyn decomposition output\n'
            + 'ASCII\n'
            + 'DATASET UNSTRUCTURED_GRID\n'
            )

    # write cell corner points as point list
    f_out.write('POINTS ' + str(8 * numcells) + ' float\n')  # 8 points per cell
    for cell in celllist:
        f_out.write(cell.getPointCoordinates())

    # write cells. Every line is a cell
    f_out.write('CELLS ' + str(numcells) + ' ' + str(numcells * 9) + '\n')
    # CELLS requires numcells and number of total list entries (=numcells + 8*numcells)
    for cell in celllist:
        f_out.write('8 ' + cell.getCellPointIndices())

    # write vtk cell types
    f_out.write('CELL_TYPES ' + str(numcells) + '\n')
    for cell in celllist:
        f_out.write('11\n')  # type 11 = VTK_VOXEL = cuboid

    # write ids for every cell. They are not necessarily the same as the MPI Rank ID
    f_out.write('CELL_DATA ' + str(numcells) + '\n'
            + 'SCALARS cell_scalars int 1\n'
            + 'LOOKUP_TABLE default\n'
            )
    i=0
    for cell in celllist:
        f_out.write(str(i) + '\n')
        i = i + 1

    # write all configuration information as field data
    f_out.write('FIELD ConfigurationData '+ str(len(celllist[0].config.keys()) + 1) +'\n')
    # write one combined string
    f_out.write('FullConfiguration 1 ' + str(numcells) + ' string\n')
    for cell in celllist:
        f_out.write('+'.join(cell.config.values()) + '\n')
    # write each config property individually
    for configKey in celllist[0].config.keys():
        f_out.write('\n' + configKey.replace(' ', '') + ' 1 ' + str(numcells) + ' string\n')
        for cell in celllist:
            f_out.write(cell.config[configKey] + '\n')

    # finish
    f_in.close()
    f_out.close()

