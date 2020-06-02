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
import os
import ntpath
import cmd
from subprocess import Popen, PIPE
from shlex import split
import compareHelpers
import time

# from twisted.internet.defer import returnValue

mpi = '-1'
newMarDyn = ''
oldMarDyn = '-1'
xmlFilename = ''
inpFilename = ''
additionalFilenames = []
comparePlugins = ['ResultWriter', 'GammaWriter', 'RDF']
disabledPlugins = []
numIterations = '25'
baseisnormal = 0
remote = ''
remoteprefix = '/scratch'
# shortopts: if they have an argument, then add : after shortcut
options, remainder = getopt(argv[1:], 'M:m:n:o:c:i:p:I:hbr:R:LB:a:AS',
                            ['mpicmd=',
                             'mpi=',
                             'newMarDyn=',
                             'oldMarDyn=',
                             'xmlFilename=',
                             'inpFilename=',
                             'plugin=',
                             'numIterations=',
                             'help',
                             'baseisnormal',
                             'remote=',
                             'remoteprefix=',
                             'baseIsLocal',
                             'baseRemote=',
                             'additionalFile=',
                             'allMPI',
                             'legacy-cell-processor',
                             'srunFix',
                             'disablePlugin='
                             ])
nonDefaultPlugins = False
baseIsLocal = False
legacyCellProcessor = False

allMPI = False
MPI_START = 'mpirun'  # e.g. I need to set it to mpirun.mpich locally
print(options)
baseRemote = ""
for opt, arg in options:
    if opt in ('-n', '--newMarDyn'):
        newMarDyn = arg
    elif opt in ('-o', '--oldMarDyn'):
        oldMarDyn = arg
    elif opt in ('-c', '--xmlFilename'):
        xmlFilename = arg
    elif opt in ('-i', '--inpFilename'):
        inpFilename = arg
    elif opt in ('-m', '--mpi'):
        mpi = arg
    elif opt in ('-M', '--mpicmd'):
        MPI_START = arg
    elif opt in ('-p', '--plugin'):
        if (not nonDefaultPlugins):  # first encounter of "-p" -> clear plugin list
            nonDefaultPlugins = True
            comparePlugins = []
        comparePlugins.append(arg)
    elif opt in ('--disablePlugin'):
        disabledPlugins.append(arg)
    elif opt in ('-I', '--numIterations'):
        numIterations = arg
    elif opt in ('-h', '--help'):
        print("Make sure two versions of mardyn produce identical simulation results. Sample usage:")
        print(" multiple -p are possible. Currently ResultWriter, GammaWriter and RDF are supported.")
        print(
            """ ./vr -m 4 -n MarDyn.PAR_RELEASE_AVX2 -o MarDyn.PAR_RELEASE_AOS -c ../examples/surface-tension_LRC/C6H12_500/C6H12_500_1R.xml -i ../examples/surface-tension_LRC/C6H12_500/C6H12_500.inp -p GammaWriter -I 10 """)
        print(" -b specifies, that the base (old file) is assumed to be sequential")
        print(" -r specifies the remote host ")
        print(
            " --remoteprefix changes the prefix of the directory, that is used on the remote host. as relative path to $HOME or absolute path if path starts with /")
        exit(1)
    elif opt in ('-b', '--baseisnormal'):
        baseisnormal = 1
    elif opt in ('-r', '--remote'):
        remote = arg
        print("remote", remote)
    elif opt in ('-R', '--remoteprefix'):
        remoteprefix = arg
    elif opt in ('--baseIsLocal'):
        baseIsLocal = True
    elif opt in ('-B', '--baseRemote'):
        baseRemote = arg
    elif opt in ('-a', '--additionalFile'):
        additionalFilenames.append(arg)
    elif opt in ('-A', '--allMPI'):
        allMPI = True
    elif opt in ('', '--legacy-cell-processor'):
        legacyCellProcessor = True
    elif opt in ('-S', '--srunFix'):
        os.environ["I_MPI_PMI_LIBRARY"] = "/usr/lib64/libpmi.so"
    else:
        print("unknown option: " + opt)
        exit(1)

if baseIsLocal and baseRemote:
    print("defined baseIsLocal and defined a base remote host. this contradicts itself. exiting...")
    exit(1)

# disable disabled plugins:
for dPlugin in disabledPlugins:
    comparePlugins.remove(dPlugin)

SEQ = (mpi == '-1')
PAR = not SEQ

noReferenceRun = (oldMarDyn == '-1')
doReferenceRun = not noReferenceRun

comparePostfixes = []
for i in range(len(comparePlugins)):
    comparePlugin = comparePlugins[i]
    if comparePlugin == 'Resultwriter' or comparePlugin == 'ResultWriter':
        comparePlugins[i] = comparePlugin = 'ResultWriter'
        comparePostfixes.append('.res')
    elif comparePlugin == 'GammaWriter':
        comparePostfixes.append('.gamma')
    elif comparePlugin == 'RDF':
        comparePostfixes.append('.rdf')
    else:
        print("Plugin " + comparePlugin + " not supported yet.")
        print("Have a look whether you can add it yourself.")
        exit(1)

if noReferenceRun:
    print("no old version given. Will try to reuse existing output, by not erasing it at start.")

# JUMP to validationRuns - extract path to validationRuns from argv[0]!
pathToValidationRuns = ntpath.dirname(os.path.realpath(__file__))
pathToValidationRuns = os.path.realpath(pathToValidationRuns)

print(pathToValidationRuns)

# first clean all the folders
cleanUpCommand = ['rm']
cleanUpCommand.extend(glob(pathToValidationRuns + '/*.xml'))
cleanUpCommand.extend(glob(pathToValidationRuns + '/*.cfg'))  # remains for cleanup
cleanUpCommand.extend(glob(pathToValidationRuns + '/*.inp'))
cleanUpCommand.extend(glob(pathToValidationRuns + '/new/*'))
if doReferenceRun:
    cleanUpCommand.extend(glob(pathToValidationRuns + '/reference/*'))
cleanUpCommand.extend(glob(pathToValidationRuns + '/MarDyn*'))
p = Popen(cleanUpCommand, stdout=PIPE, stderr=PIPE)
p.communicate()  # suppresses possible errors if nothing there yet, as we don't want them for rm

# copy all there
call(['cp', newMarDyn, xmlFilename, inpFilename] + additionalFilenames + [pathToValidationRuns])
if doReferenceRun:
    call(['cp', oldMarDyn, pathToValidationRuns])

# go there
os.chdir(pathToValidationRuns)

# get the basenames
xmlBase = ntpath.basename(xmlFilename)
inpBase = ntpath.basename(inpFilename)
additionalFileBases = []
additionalFileBases[:] = [ntpath.basename(fn) for fn in additionalFilenames]
oldMarDynBase = ntpath.basename(oldMarDyn)
newMarDynBase = ntpath.basename(newMarDyn)

# print "append ComparisonWriter here"
# print "append ComparisonWriter here"
with open(xmlBase, "r") as prev_file:
    with open("tmp.xml", "w") as new_file:
        contents = prev_file.readlines()
        # Now contents is a list of strings and you may add the new line to this list at any position
        # contents.insert(4, "\n This is a new line \n ")

        for i in range(len(contents)):
            line = contents[i]
            if 'RDF' in comparePlugins and line.find("<run>") != -1:
                contents.insert(i + 1, """<equilibration><steps>0</steps></equilibration>\n""")
                i += 1
                continue

            if line.find("<output>") != -1:
                for comparePlugin in comparePlugins:
                    if comparePlugin == 'RDF':  # configuring RDF within the xml is different...
                        contents.insert(i + 1, """
        <outputplugin name="RDF">
            <writefrequency>10</writefrequency>
            <outputprefix>val.comparison</outputprefix>
            <intervallength>0.003</intervallength>
            <bins>1000</bins>
        </outputplugin>\n""")
                        i += 1
                        # myfile.write("initStatistics 0\nRDF 0.003 1000\nRDFOutputTimesteps 10\nRDFOutputPrefix val.comparison\n")
                    else:
                        contents.insert(i + 1, """
        <outputplugin name=\"""" + comparePlugin + """\">
            <writefrequency>1</writefrequency>
            <outputprefix>val.comparison</outputprefix>
        </outputplugin>""")
                        i += 1
                    # myfile.write("output " + comparePlugin + " 1 val.comparison\n")
        new_file.write("".join(contents))
call(['mv', 'tmp.xml', xmlBase])

comparisonFilenames = []
for comparePostfix in comparePostfixes:
    comparisonFilenames.append('val.comparison' + comparePostfix)

if doReferenceRun:
    call(['mkdir', '-p', 'reference/'])
    call(['cp', xmlBase, 'reference/'])
    call(['cp', inpBase, 'reference/'])
    if len(additionalFileBases):
        call(['cp'] + additionalFileBases + ['reference/'])
    call(['cp', oldMarDynBase, 'reference/'])
call(['mkdir', '-p', 'new/'])
call(['cp', xmlBase, 'new/'])
call(['cp', inpBase, 'new/'])
if len(additionalFileBases):
    call(['cp'] + additionalFileBases + ['new/'])
call(['cp', newMarDynBase, 'new/'])


def doRun(directory, MardynExe):
    # first run
    if baseRemote and directory == "reference":
        localRemote = baseRemote
    else:
        localRemote = remote
    os.chdir(directory)
    call(['chmod', '+x', MardynExe])
    cmd = []

    doRemote = localRemote and (directory == 'new' or not baseIsLocal)

    if doRemote:
        rsyncremote = localRemote
        if localRemote.endswith('-mic0') or localRemote.endswith('-mic1'):
            rsyncremote = localRemote[:-5]
        command = "mkdir -p " + remoteprefix
        mkdircmd = []
        mkdircmd.extend(['ssh', rsyncremote, command])
        p = Popen(mkdircmd, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if p.returncode:
            print("error on mkdir -p:")
            print(out, err)
            exit(1)
        remotedirectory = remoteprefix + "/" + directory
        command = "rsync --delete-before -r ../" + directory + " " + rsyncremote + ":" + remoteprefix
        print(command)
        p = Popen(split(command))
        p.wait()
        if p.returncode:
            print("error on rsync")
            exit(1)
        command = "cd " + remotedirectory + " && pwd && "
        cmd.extend(['ssh', localRemote, command])

    if allMPI:
        cmd.extend(split(MPI_START))
        if directory == 'new' or not baseisnormal:
            cmd.extend(['-n', str(mpi)])
        else:
            cmd.extend(['-n', '1'])
    else:
        if PAR and (directory == 'new' or not baseisnormal):
            cmd.extend(split(MPI_START))
            cmd.extend(['-n', str(mpi)])

    if legacyCellProcessor and directory == "new":
        cmd.extend(
            ['./' + MardynExe, "--legacy-cell-processor", "--final-checkpoint=0", xmlBase, "--steps", numIterations])
    else:
        cmd.extend(['./' + MardynExe, "--final-checkpoint=0", xmlBase, "--steps", numIterations])
    # cmd.extend(['/work_fast/tchipevn/SDE/sde-external-7.41.0-2016-03-03-lin/sde64', '-knl', '--', './' + MardynExe, "--final-checkpoint=0", xmlBase, numIterations]);
    print(cmd)
    print("================")
    t = time.time()
    while True:
        # repeatedly try this if srun was not working
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if p.returncode == 1 and (
                "Job violates accounting/QOS policy" in err or "Socket timed out on send/recv" in err):
            print("srun submit limit reached or socket timed out error, trying again in 60s")
            time.sleep(60)
            continue
        break
    t = time.time() - t
    print("elapsed time:", t)
    if p.returncode:
        print("error while executing program:")
        print(out, err)
        exit(1)
    print(out, err)
    if doRemote:  # sync back
        command = "rsync " + rsyncremote + ":" + remotedirectory + "/* ./"
        print(command)
        p = Popen(split(command))
        p.wait()

    if "RDF" in comparePlugins:
        p = Popen(['ls', '-r'] + glob("val.comparison*.rdf"), stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        p = Popen(split("cp " + split(out)[0] + " val.comparison.rdf"))
        # Copy newest rdf file to val.comparison.rdf
        p.wait()
    for comparisonFilename in comparisonFilenames:
        # possible switch/if statements if other comparison plugins require different output.
        p = Popen(split(
            "sed -i.bak '/^#/d; s/[[:blank:]]*$//; /^$/d' " + comparisonFilename))  # deletes lines starting with #.
        # These are the lines containing timestamps, and have to be removed for proper comparison.
        p.wait()
    os.chdir('..')


print("new run:")
# first run
doRun('new', newMarDynBase)

# second run
if doReferenceRun:
    print("reference run:")
    doRun('reference', oldMarDynBase)
returnValue = 0
# call(['diff' 'reference/val.comparison.res' 'new/val.comparison.res'])
print("")
for i in range(len(comparePlugins)):
    localReturn = compareHelpers.compareFiles("reference/" + comparisonFilenames[i], "new/" + comparisonFilenames[i])
    returnValue += localReturn
    if localReturn == 0:
        print("Identical values! for ", comparePlugins[i])
    else:
        print("mismatches for ", comparePlugins[i])

if returnValue == 0:
    print("")
    print("Identical values!")
    print("")
    exit(0)
else:
    print("")
    print("mismatches")
    print("")
    exit(1)
