import vtk

import sys
import os
import multiprocessing
import functools

"""Combines all vtu files referenced by one pvtu to one single vtu.
This script processes up to N pvtu files in parallel, where N is the number of logical cores on the system.
"""

convertedPrefix = 'combined_'

################### VTU ####################

def convertVtu(vtufilename, appender):

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtufilename)
    reader.Update()
    readerOutput = reader.GetOutput()

    appender.AddInputData(readerOutput)

################### PVTU ###################

def convertPVtu(pvtuFilename, outputDir):
    pvtuDir = os.path.dirname(pvtuFilename)
    outputFile = os.path.join(outputDir, convertedPrefix + os.path.splitext(os.path.basename(pvtuFilename))[0] + '.vtu')
    print(pvtuFilename + ' -> ' + outputFile)
    # This thing can collect data from multiple datasets
    appender = vtk.vtkAppendFilter()

    with open(pvtuFilename, 'r') as pvtuFileIn:
        for line in pvtuFileIn:
            if '<Piece Source' in line:
                # extract sub vtu filenames
                vtuFilename = line.split('"')[1]
                convertVtu(os.path.join(pvtuDir, vtuFilename), appender)

    appender.Update()
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(outputFile)
    appenderOutput = appender.GetOutput()
    writer.SetInputData(appenderOutput)
    # Create the output dir if it doesn't exist yet
    os.makedirs(outputDir, exist_ok=True)
    writer.Write()

############################################ Script ############################################
if __name__ == '__main__':
    # Require at least two arguments.
    if len(sys.argv) < 3:
         raise RuntimeError(f"""\
 Not enough input arguments.

 Usage: {os.path.basename(sys.argv[0])} outputDir files.pvtu...
 """)

    outputDir = os.path.normpath(sys.argv[1])
    files = sys.argv[2:]
    # bind the output dir argument to the function
    convertPVtuWithOutputDir = functools.partial(convertPVtu, outputDir=outputDir)
    # parallelize over the file list
    multiprocessing.Pool().map(convertPVtuWithOutputDir, files)
