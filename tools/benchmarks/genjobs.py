#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# genjobs.py
# generate jobs from a template job, subsituting parameters
#
# Martin Bernreuther <bernreuther@hlrs.de>, November 2011
# license: GPL (see http://www.gnu.org/licenses/gpl.html)

import sys
import os
import optparse
import ConfigParser
#import xml.dom.minidom
import shutil
import shlex
import subprocess
import time
#import re

desc="generate&start benchmark runs"
version="2011.11.25"
author="Martin Bernreuther <bernreuther@hlrs.de>"

print time.strftime("genjobs.py starting at %a, %d.%m.%Y %H:%M:%S %Z")

label=time.strftime("%Y%m%dT%H%M%S")

configfile="config"
rootdir="."	# relative path definitions will be relative to the configfiledir
dstroot=label
delimiter='\f'	# '\f','\t','\000',' ',':'
gentemplate="gentemplate"
genparareplfiles=[]
gencondition=None
gencommand=None
pptemplate=None
ppscript=None
ppcommand=None


if sys.version_info < (2, 6):
	#raise "ERROR: python version too old: "+str(sys.version_info)+"<2.6"
	sys.stderr.write("\nWARNING: using old python version "+str(sys.version_info)+" < 2.6\n\n")

progname=os.path.basename(sys.argv[0])
optparser = optparse.OptionParser(usage="usage: {0} [options] [configfile]\n{1}\n(version: {2}, author: {3}, license: GPL)".format(progname,desc,version,author))	## %prog
optparser.add_option("-V", "--version", action="store_true", dest="print_version", default=False, help="print version and exit [default: {0}]".format("False"))	## %default
#optparser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="be more verbose [default: {0}]".format("False"))
#optparser.add_option("-c", "--config", type="string", dest="configfile", help="configuration file [default: {0}]".format(configfilename))
optparser.add_option("-r", "--root", type="string", dest="rootdir", help="root directory for templates [default: <directory of configuration file>]")
optparser.add_option("-d", "--dstroot", type="string", dest="dstroot", help="destination root directory [default: {0}]".format(dstroot))
optparser.add_option("-t", "--delimiter", type="string", dest="delimiter", help="delimiter for lists [default (hex): {0}]".format(''.join([hex(ord(c))[2:].rjust(2,'0') for c in delimiter])))
optparser.add_option("-g", "--gentemplate", type="string", dest="gentemplate", help="generator job template directory [default: {0}]".format(gentemplate))
optparser.add_option("-p", "--pptemplate", type="string", dest="pptemplate", help="postprocessing template directory [default: {0}]".format(pptemplate))
(options, args) = optparser.parse_args()

if options.print_version:
	## print version and exit
	print version
	sys.exit(0)

parasubs={}

print "options ----------------------------------------------------------------"
#if options.configfile is not None:
#	configfile=options.configfile
if len(args)>1: configfile=args[1]
if not os.path.isfile(configfile):
	sys.stderr.write("ERROR: configuration file \"{0}\" not found!\n".format(configfile))
	sys.exit(1)
configfiledir=os.path.dirname(os.path.realpath(configfile))
print "configuration file:\t{0}".format(configfile)
#print "         directory:\t{0}".format(configfiledir)

cfgparser = ConfigParser.SafeConfigParser()
cfgparser.optionxform = str
cfgparser.read([configfile])

if not os.path.isabs(rootdir):
	rootdir=os.path.join(configfiledir,rootdir)
if cfgparser.has_option("general","root"):
	rootdir=cfgparser.get("general","root")
	if not os.path.isabs(rootdir):
		rootdir=os.path.join(configfiledir,rootdir)
if options.rootdir is not None:
	rootdir=options.rootdir
if not os.path.isdir(rootdir):
	sys.stderr.write("ERROR: root directory {0} does not exist\n".format(rootdir))
	sys.exit(2)
rootdir=os.path.realpath(rootdir)
print "root directory:\t{0}".format(rootdir)

if cfgparser.has_option("general","dstroot"):
	dstroot=cfgparser.get("general","dstroot")
if options.dstroot is not None:
	dstroot=os.path.join(rootdir,options.dstroot)
parasubs["$DSTROOT"]=dstroot
parasubs["$DSTROOTPATH"]=os.path.realpath(dstroot)
print "destination directory:\t{0}".format(dstroot)

if cfgparser.has_option("general","delimiter"):
	delimiter=cfgparser.get("general","delimiter")
if options.delimiter is not None:
	delimiter=options.delimiter
parasubs["$DELIMITER"]=delimiter
print "list seperator (hex):\t{0}".format(''.join([hex(ord(c))[2:].rjust(2,'0') for c in delimiter]))

if cfgparser.has_option("generator","template"):
	gentemplate=cfgparser.get("generator","template")
if options.gentemplate is not None:
	gentemplate=options.gentemplate
gentemplate=os.path.join(rootdir,gentemplate)
if not os.path.isdir(gentemplate):
	sys.stderr.write("ERROR: template directory {0} does not exist\n".format(gentemplate))
	sys.exit(3)
parasubs["$GENTEMPLATE"]=gentemplate
parasubs["$GENTEMPLATEPATH"]=os.path.realpath(gentemplate)
print "generator template directory:\t{0}".format(gentemplate)

if cfgparser.has_option("generator","parafiles"):
	genparareplfiles=cfgparser.get("generator","parafiles").split()
print "generator files to replace parameters:\t{0}".format(genparareplfiles)

if cfgparser.has_option("generator","condition"):
	gencondition=cfgparser.get("generator","condition")
	print "parameters condition:\t",gencondition

if cfgparser.has_option("generator","command"):
	gencommand=cfgparser.get("generator","command")
	print "generator command:\t",gencommand

if cfgparser.has_option("postproc","template"):
	pptemplate=cfgparser.get("postproc","template")
if options.pptemplate is not None:
	pptemplate=options.pptemplate
if pptemplate is not None:
	pptemplate=os.path.join(rootdir,pptemplate)
	if not os.path.isdir(pptemplate):
		sys.stderr.write("ERROR: postprocessing template directory {0} does not exist\n".format(pptemplate))
		sys.exit(4)
	parasubs["$PPTEMPLATE"]=pptemplate
	parasubs["$PPTEMPLATEPATH"]=os.path.realpath(pptemplate)
	print "postprocessing template:\t{0}".format(pptemplate)
	
	if cfgparser.has_option("postproc","parafiles"):
		ppparareplfiles=cfgparser.get("postproc","parafiles").split()
	print "postprocessing files to replace parameters:\t{0}".format(genparareplfiles)
	
	if cfgparser.has_option("postproc","command"):
		ppcommand=cfgparser.get("postproc","command")
		print "postprocessing command:\t",ppcommand


print "parameters:"
print "name\tvalues\t(width,int)"
parameters=[]
numvar=0
for i in cfgparser.items("parameters"):
	paraname,paravalue=i
	paravalues=paravalue.split()
	maxlen=0
	onlydigits=True
	for paravalue in paravalues:
		if len(paravalue)>maxlen: maxlen=len(paravalue)
		if not paravalue.isdigit(): onlydigits=False
	parameters.append([paraname,paravalues,[maxlen,onlydigits]])
	print "{0}\t{1}\t({2},{3})".format(paraname,paravalues,maxlen,onlydigits)
	if not numvar:
		numvar=len(paravalues)
	else:
		numvar*=len(paravalues)



def execmd(cmd,wd="."):
	#print "execute command \"{0}\" (working directory: {1})".format(cmd,wd)
	rc=-1
	stdoutdata=""
	stderrdata=""
	try:
		try:
			#shell scripts have to be started using sh
			stdoutdata=subprocess.check_output(shlex.split(cmd),cwd=wd)
			# with shell=True it should also be possible to start shell scripts directly, but this option might introduce a security risk
			#stdoutdata=subprocess.check_output(shlex.split(cmd),cwd=wd,shell=True)
		except subprocess.CalledProcessError,(rc,stdoutdata):
			print "ERROR ({0}) executing {1}: {2}".format(rc,cmd,stdoutdata)
		else:
			rc=0
	except AttributeError:
		try:
			pid=subprocess.Popen(shlex.split(cmd),cwd=wd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			stdoutdata,stderrdata=pid.communicate()
			rc=int(pid.returncode)
		except ValueError:
			print "Popen ERROR executing "+cmd
	return rc,stdoutdata,stderrdata

def adjrlinks(dstdir,srcdir):
	#print "adapt relative links contained in directory {0} copied from {1}".format(dstdir,srcdir)
	dstdirqueue=["."]
	while len(dstdirqueue)>0:
		for dirname, dirnames, filenames in os.walk(os.path.join(dstdir,dstdirqueue.pop())):
			for subdirname in dirnames:
				dstdirqueue.append(subdirname)
			for filename in filenames:
				dstfile=os.path.join(dstdir,filename)
				if os.path.islink(dstfile):
					dstlinkfile=os.path.realpath(dstfile)
					templatelinkfile=os.path.realpath(os.path.join(srcdir,filename))
					if dstfile!=templatelinkfile:
						os.remove(dstfile)
						os.symlink(os.path.relpath(templatelinkfile,dstdir),dstfile)



if os.path.exists(dstroot):
	sys.stderr.write("ERROR: destination directory {0} already exists\n".format(dstroot))
	sys.exit(5)

if pptemplate is not None:
	print "creating directory {0} copying {1}".format(dstroot,pptemplate)
	shutil.copytree(pptemplate,dstroot,symlinks=True)
	adjrlinks(dstroot,pptemplate)
else:
	print "creating directory {0}".format(dstroot)
	#os.mkdir(dstroot)
	os.makedirs(dstroot)
if not (os.path.isdir(dstroot) and os.access(dstroot, os.W_OK)):
	sys.stderr.write("ERROR: destination directory {0} not accessible\n".format(dstroot))
	sys.exit(6)

print

print "generate jobs ##########################################################"
createdjobs=[]
cmd_output=[]
v=[0]*len(parameters)
for i in range(numvar):
	cond=gencondition
	cmd=gencommand
	jobname=""
	for p in range(len(v)):
		paraname=parameters[p][0]
		jobname+=paraname
		paravalues=parameters[p][1]
		if v[p]>=len(parameters[p][1]):
			v[p+1]+=1
			v[p]=0
		paravalue=paravalues[v[p]]
		parasubs['$'+paraname]=paravalue
		maxlen,onlydigits=parameters[p][2]
		if onlydigits:
			jobname+=paravalue.zfill(maxlen)
		else:
			jobname+=str(paravalue)
	parasubs["$JOBNAME"]=jobname
	# substitute all parameters within condition & cmd
	for paraname,paravalue in parasubs.items():
		if cond is not None:
			cond=cond.replace(paraname,str(paravalue))
		if cmd is not None:
			cmd=cmd.replace(paraname,str(paravalue))
	if cond is not None:
		try:
			evalcond=eval(cond)
		except Exception, errno:
			evalcond=False
	else:
		evalcond=True
	if evalcond:
		jobdir=os.path.join(dstroot,jobname)
		parasubs["$JOBDIR"]=jobdir
		parasubs["$JOBDIRPATH"]=os.path.realpath(jobdir)
		print '=',i+1,'/',numvar," benchmark job directory:",jobdir,'='*(39-len(jobdir)-len(str(i+1))-len(str(numvar)))
		if os.path.exists(jobdir):
			print "destination directory {0} already exists ... skipping".format(jobdir)
			continue
		shutil.copytree(gentemplate,jobdir,symlinks=True)
		adjrlinks(jobdir,gentemplate)
		# substitute parameters within all parafiles
		for parareplacefile in genparareplfiles:
			substfile=os.path.join(jobdir,parareplacefile)
			try:
				substfhdl=open(substfile, 'r')
			except IOError, (errno, strerror):
				print "Error ({0}) opening {1} for reading: {2}".format(errno,substfile,strerror)
				sys.exit(7)
			substfilecontent=substfhdl.read()	# substfile must not be too large to be read at once!
			# (optional size argument could limit the portion to be read, but then substfilecontent might be truncated...)
			substfhdl.close()
			# alternatively loop over lines (with e.g. readline(), xreadlines or fileinput) and read-replace-write line-by-line
			# replace parameters within file content
			for paraname,paravalue in parasubs.items():
				substfilecontent=substfilecontent.replace(paraname,str(paravalue))
			substfile=os.path.join(jobdir,parareplacefile)
			try:
				substfhdl=open(substfile, 'w')
			except IOError, (errno, strerror):
				print "Error ({0}) opening {1} for writing: {2}".format(errno,substfile,strerror)
				sys.exit(8)
			substfhdl.write(substfilecontent)
			substfhdl.close()
		createdjobs.append(jobname)
		#
		if cmd is not None:
			print "- command:",cmd,'-'*(60-len(cmd))
			rc,stdoutdata,stderrdata=execmd(cmd,jobdir)
			if rc!=0:
				print stdoutdata
				print stderrdata
				print '-'*64-len(str(rc)),rc,"failed"
			else:
				print stdoutdata
				print '-'*67,"done"
				cmd_output.append(stdoutdata.rstrip("\n"))
	else:
		print "condition",cond,"does not hold"
	v[0]+=1

for p in range(len(parameters)):
	paraname=parameters[p][0]
	del parasubs['$'+paraname]
parasubs["$CREATEDJOBS"]=delimiter.join(createdjobs)
parasubs["$GENCMDOUTPUT"]=delimiter.join(cmd_output)

print "################################################################### done"


if pptemplate is not None:
	print "postprocessing job #####################################################"
	# substitute parameters within all parafiles
	for parareplacefile in ppparareplfiles:
		substfile=os.path.join(dstroot,parareplacefile)
		try:
			substfhdl=open(substfile, 'r')
		except IOError, (errno, strerror):
			print "Error ({0}) opening {1} for reading: {2}".format(errno,substfile,strerror)
			sys.exit(9)
		substfilecontent=substfhdl.read()	# substfile must not be too large to be read at once!
		substfhdl.close()
		# replace parameters within file content
		for paraname,paravalue in parasubs.items():
			substfilecontent=substfilecontent.replace(paraname,str(paravalue))
		substfile=os.path.join(dstroot,parareplacefile)
		try:
			substfhdl=open(substfile, 'w')
		except IOError, (errno, strerror):
			print "Error ({0}) opening {1} for writing: {2}".format(errno,substfile,strerror)
			sys.exit(10)
		substfhdl.write(substfilecontent)
		substfhdl.close()
	if ppcommand is not None:
		ppcmd=ppcommand
		for paraname,paravalue in parasubs.items():
				ppcmd=ppcmd.replace(paraname,str(paravalue))
		print "- command:",ppcmd,'-'*(60-len(ppcmd))
		rc,stdoutdata,stderrdata=execmd(ppcmd,dstroot)
		if rc!=0:
			print stdoutdata
			print stderrdata
			print '-'*64-len(str(rc)),rc,"failed"
		else:
			print stdoutdata
			print '-'*67,"done"
			cmd_output.append(stdoutdata)
	print "################################################################### done"

print time.strftime("finished at %a, %d.%m.%Y %H:%M:%S %Z")
