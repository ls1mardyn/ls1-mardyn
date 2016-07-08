#!/usr/bin/python

# This script tests that no changes to simulation output were introduced between two MarDyn binaries
# E.g. we change the structure of the cache-s to SoA-s and we want that the program gives the same
# physical output as before.
#
# I'm starting a first basic version, many options will be added as this script is turned into a
# more general tool which Jenkins automatically tests.
#
# For the moment, this will test the output of a ResultWriter MarDyn output-plugin.
#
# Author: Nikola Tchipev
# 8. Mai 2016


from getopt import getopt
from sys import argv
from subprocess import call
from glob import glob
from ntpath import basename
import os
import ntpath
import cmd
from subprocess import Popen, PIPE
from shlex import split
import compareHelpers

mpi = '-1'
newMarDyn = ''
oldMarDyn = '-1'
cfgFilename = ''
inpFilename = ''
comparePlugin = 'ResultWriter'
numIterations = '25'


options, remainder = getopt(argv[1:], 'm:n:o:c:i:p:I:h', 
                            ['mpi=',
                             'newMarDyn=',
                             'oldMarDyn=',
                             'cfgFilename=',
                             'inpFilename=',
                             'plugin=',
                             'numIterations=',
                             'help'
                             ])

for opt, arg in options:
    if opt in ('-n', '--newMarDyn'):
        newMarDyn = arg
    elif opt in ('-o', '--oldMarDyn'):
        oldMarDyn = arg
    elif opt in ('-c', '--cfgFilename'):
        cfgFilename = arg
    elif opt in ('-i', '--inpFilename'):
        inpFilename = arg
    elif opt in ('-m', '--mpi'):
        mpi = arg
    elif opt in ('-p', '--plugin'):
        comparePlugin = arg
    elif opt in ('-I', '--numIterations'):
        numIterations = arg
    elif opt in ('-h', '--help'):
        print "Make sure two versions of mardyn produce identical simulation results. Sample usage:"
        print """ ./vr -m 4 -n MarDyn.PAR_RELEASE_AVX2 -o MarDyn.PAR_RELEASE_AOS -c ../examples/surface-tension_LRC/C6H12_500/C6H12_500_1R.cfg -i ../examples/surface-tension_LRC/C6H12_500/C6H12_500.inp -p GammaWriter -I 10 """
        exit(1)
    else:
        print "unknown option: " + opt
        exit(1)

SEQ = (mpi == '-1')
PAR = not SEQ   
        
noReferenceRun = (oldMarDyn == '-1')
doReferenceRun = not noReferenceRun

MPI_START='mpirun' # e.g. I need to set it to mpirun.mpich locally

if comparePlugin == 'ResultWriter':
    comparePostfix = '.res'
elif comparePlugin == 'GammaWriter':
    comparePostfix = '.gamma'
else:
    print "Plugin " + comparePlugin + " not supported yet."
    print "Have a look whether you can add it yourself."
    exit(1)
    


if noReferenceRun:
    print "no old version given. Will try to reuse existing output, by not erasing it at start."
        
# JUMP to validationRuns - extract path to validationRuns from argv[0]!
pathToValidationRuns = ntpath.dirname(os.path.realpath(__file__))
pathToValidationRuns = os.path.realpath(pathToValidationRuns)

print pathToValidationRuns
        
# first clean all the folders
cleanUpCommand = ['rm']
cleanUpCommand.extend(glob(pathToValidationRuns + '/*.cfg'))
cleanUpCommand.extend(glob(pathToValidationRuns + '/*.inp'))
cleanUpCommand.extend(glob(pathToValidationRuns + '/new/*'))
if doReferenceRun:
    cleanUpCommand.extend(glob(pathToValidationRuns + '/reference/*'))
cleanUpCommand.extend(glob(pathToValidationRuns + '/MarDyn*'))
call(cleanUpCommand)

# copy all there
call(['cp', newMarDyn, cfgFilename, inpFilename, pathToValidationRuns])
if doReferenceRun:
    call(['cp', oldMarDyn, pathToValidationRuns])

# go there
os.chdir(pathToValidationRuns)

# get the basenames
cfgBase = ntpath.basename(cfgFilename)
inpBase = ntpath.basename(inpFilename)
oldMarDynBase = ntpath.basename(oldMarDyn)
newMarDynBase = ntpath.basename(newMarDyn)

print "append ComparisonWriter here"
with open(cfgBase, "a") as myfile:
    myfile.write("output " + comparePlugin + " 1 val.comparison")
    
comparisonFilename = 'val.comparison' + comparePostfix

if doReferenceRun:
    call(['cp', cfgBase, 'reference/'])
    call(['cp', inpBase, 'reference/'])
    call(['cp', oldMarDynBase, 'reference/'])
call(['cp', cfgBase, 'new/'])
call(['cp', inpBase, 'new/'])
call(['cp', newMarDynBase, 'new/'])

# first run
os.chdir('new')
cmd = []
if PAR:
    cmd.extend([MPI_START, '-n', str(mpi)])
cmd.extend(['./' + newMarDynBase, cfgBase, numIterations]); 
print cmd
p = Popen(cmd, stdout=PIPE, stderr=PIPE)
out, err = p.communicate()

p = Popen(split("sed -i.bak /[Mm]ar[Dd]yn/d " + comparisonFilename))  # deletes lines containing MarDyn, mardyn, Mardyn or marDyn. 
# These are the lines containing timestamps, and have to be removed for proper comparison.
p.wait()
os.chdir('..')

## second run
if doReferenceRun:
    os.chdir('reference')
    cmd = []
    if PAR:
        cmd.extend([MPI_START, '-n', str(mpi)])
    cmd.extend(['./' + oldMarDynBase, cfgBase, numIterations])
    print cmd
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    
    p = Popen(split("sed -i.bak /[Mm]ar[Dd]yn/d " + comparisonFilename))
    p.wait()
    os.chdir('..')

#call(['diff' 'reference/val.comparison.res' 'new/val.comparison.res'])
returnValue=compareHelpers.compareFiles("reference/" + comparisonFilename, "new/" + comparisonFilename)

if returnValue == 0:
    print ""
    print "Identical values!"
    print ""
    exit(0)
else:
    print ""
    print "mismatches"
    print ""
    exit(1)
