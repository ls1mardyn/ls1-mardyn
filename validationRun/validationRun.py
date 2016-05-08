#!/usr/bin/python

# This script test that no changes to simulation output were introduced between two MarDyn binaries
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
from subprocess import Popen
from shlex import split

newMarDyn = ''
oldMarDyn = ''
cfgFilename = ''
inpFilename = ''
comparePlugin = 'ResultWriter'


options, remainder = getopt(argv[1:], '', ['newMarDyn=',
                                                      'oldMarDyn=',
                                                      'cfgFilename=',
                                                      'inpFilename='
                                                      ])
# Options will come eventually
#print 'OPTIONS   :', options
#
# for opt, arg in options:
#     if opt in ('-o', '--output'):
#         output_filename = arg
#     elif opt in ('-v', '--verbose'):
#         verbose = True
#     elif opt == '--version':
#         version = arg

if len(options) == 0:
    print "No options given:"
    print "    trying format <newMarDyn> <oldMarDyn> <cfgFilename> <inpFilename>"
    print "    with default values ResultWriter, frequency 1, 50 iterations\n\n\n"
    if len(argv) != 5:
        print "did not specify 4 arguments"
        exit(1)
    else:
        newMarDyn = argv[1]
        oldMarDyn = argv[2]
        cfgFilename = argv[3]
        inpFilename = argv[4]
        
        
# JUMP to validationRuns - extract path to validationRuns from argv[0]!
pathToValidationRuns = ntpath.dirname(argv[0])
        
# first clean all the folders
cleanUpCommand = ['rm']
cleanUpCommand.extend(glob(pathToValidationRuns + '/*.cfg'))
cleanUpCommand.extend(glob(pathToValidationRuns + '/*.inp'))
cleanUpCommand.extend(glob(pathToValidationRuns + '/new/*'))
cleanUpCommand.extend(glob(pathToValidationRuns + '/reference/*'))
cleanUpCommand.extend(glob(pathToValidationRuns + '/MarDyn*'))
call(cleanUpCommand)

# copy all there
call(['cp', oldMarDyn, newMarDyn, cfgFilename, inpFilename, pathToValidationRuns])

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

call(['cp', cfgBase, 'reference/'])
call(['cp', inpBase, 'reference/'])
call(['cp', oldMarDynBase, 'reference/'])
call(['cp', cfgBase, 'new/'])
call(['cp', inpBase, 'new/'])
call(['cp', newMarDynBase, 'new/'])

# first run
os.chdir('new')
call(['./' + newMarDynBase, cfgBase, '50'])
Popen(split("sed -i.bak /MarDyn/d val.comparison.res"))
os.chdir('..')

## second run
os.chdir('reference')
call(['./' + oldMarDynBase, cfgBase, '50'])
Popen(split("sed -i.bak /MarDyn/d val.comparison.res"))
os.chdir('..')

#call(['diff' 'reference/val.comparison.res' 'new/val.comparison.res'])
child = Popen(split("diff reference/val.comparison.res new/val.comparison.res"))
child.communicate()[0]
returnValue = child.returncode


if returnValue == 0:
    print "Identical values!"
    exit(0)
else:
    print "mismatches"
    print returnValue
    exit(1)
