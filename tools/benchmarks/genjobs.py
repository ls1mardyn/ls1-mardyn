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
version="2011.12.28"
author="Martin Bernreuther <bernreuther@hlrs.de>"

print time.strftime("genjobs.py starting at %a, %d.%m.%Y %H:%M:%S %Z")

label=time.strftime("%Y%m%dT%H%M%S")

configfile="config"
rootdir="."	# relative path definitions will be relative to the configfiledir
dstroot=label
delimiter='\f'	# '\f','\t','\000',' ',':'
gentemplate="gentemplate"
genparareplfiles=[]
gencmdcondition=None
gencommand=None
breakcondition=None
pptemplate=None
ppscript=None
ppcommand=None


if sys.version_info < (2, 6):
	#raise "ERROR: python version too old: "+str(sys.version_info)+"<2.6"
	sys.stderr.write("\n!!! WARNING: using old python version "+str(sys.version_info)+" < 2.6 !!!\n\n")
#else
#	print "using Python version {0}".format(sys.version_info)

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
if len(args)>0: 
	configfile=args[0]
if not os.path.isfile(configfile):
	sys.stderr.write("ERROR: configuration file \"{0}\" not found!\n".format(configfile))
	sys.exit(1)
print "configuration file:\t{0}".format(configfile)
configfiledir=os.path.dirname(os.path.realpath(configfile))
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
	gencmdcondition=cfgparser.get("generator","condition")
	print "generator command condition:\t",gencmdcondition

if cfgparser.has_option("generator","command"):
	gencommand=cfgparser.get("generator","command")
	print "generator command:\t",gencommand

if cfgparser.has_option("generator","break"):
	breakcondition=cfgparser.get("generator","break")
	print "break condition:\t",breakcondition


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



class Parameter:
	"""Parameter set"""
	
	def __init__(self,parameter,value=None):
		if isinstance(parameter,Parameter):
			self.__name=parameter.name()
			self.__values=parameter.values()
		elif isinstance(parameter,str):
			self.__name=parameter
		if isinstance(value,str):
			self.__values=value.split()
		self.__initialize()
	
	def __initialize(self):
		self.__maxlen=0
		self.__onlydigits=True
		for v in self.__values:
			if len(v)>self.__maxlen: self.__maxlen=len(v)
			if not v.isdigit(): self.__onlydigits=False
		
	
	def name(self):
		return self.__name
	
	def values(self):
		return self.__values
	
	def maxlen(self):
		return self.__maxlen
		
	def onlydigits(self):
		return self.__onlydigits
	
	def numvalues(self):
		return len(self.__values)
	
	def strvalues(self,delim=" "):
		return delim.join(self.__values)
	
	def __str__(self):
		return self.__name
	
	def __int__(self):
		return self.numvalues()
	


print "parameters:"
print "name\tvalues\t(width,int)"
parameters=[]
numvar=0
jobname=""
for i in cfgparser.items("parameters"):
	if numvar==0: numvar=1
	p=Parameter(i[0],i[1])
	parasubs['$'+p.name()]=p.strvalues()
	if p.numvalues()>1:
		parameters.append(Parameter(p))
		numvar*=p.numvalues()
		print "{0}\t{1}\t({2},{3})\t*{4}".format(p.name(),p.values(),p.maxlen(),p.onlydigits(),p.numvalues())
	else:
		#parasubs['$'+p.name()]=p.values()[0]
		jobname+=p.name()+p.values()[0]
		print "{0}\t{1}\t(substitution)".format(p.name(),p.values())
print "number of jobs:\t{0}".format(numvar)


def execmd(cmd,wd="."):
	#print "execute command \"{0}\" (working directory: {1})".format(cmd,wd)
	status=-1
	stdoutdata=""
	stderrdata=""
	try:
		try:
			#stdoutdata=subprocess.check_output(cmd,cwd=wd,shell=True)	# security risk?
			stdoutdata=subprocess.check_output(shlex.split(cmd),cwd=wd)	# shell scripts have to be started using sh
			#stdoutdata=subprocess.check_output(shlex.split(cmd),cwd=wd,shell=True)
		except subprocess.CalledProcessError,(status,stdoutdata):
			print "ERROR ({0}) executing {1}: {2}".format(status,cmd,stdoutdata)
		else:
			status=0
	except AttributeError:
		try:
			pid=subprocess.Popen(cmd,cwd=wd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
			#                         executable=shlex.split(cmd)[0]
			#pid=subprocess.Popen(shlex.split(cmd),cwd=wd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
			stdoutdata,stderrdata=pid.communicate()
			status=int(pid.returncode)
			#if status: print "ERROR {0} executing \"{1}\"".format(status,cmd)	# return status and let caller deal with it
		except ValueError:
			print "Popen ERROR executing "+cmd
	return status,stdoutdata,stderrdata

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

def replaceparameters(src,paradefs):
	s=src
	for paraname,paravalue in paradefs.items():
		s=s.replace(paraname,str(paravalue))
	return s


if os.path.exists(dstroot):
	sys.stderr.write("ERROR: destination directory {0} already exists\n".format(dstroot))
	sys.exit(5)

print "creating directory {0}".format(dstroot)
if pptemplate is not None:
	print "\t(copying {0})".format(pptemplate)
	shutil.copytree(pptemplate,dstroot,symlinks=True)
	adjrlinks(dstroot,pptemplate)
else:
	#os.mkdir(dstroot)
	os.makedirs(dstroot)
if not (os.path.isdir(dstroot) and os.access(dstroot, os.W_OK)):
	sys.stderr.write("ERROR: destination directory {0} not accessible\n".format(dstroot))
	sys.exit(6)

print

print "generate jobs ##########################################################"
createdjobs=[]
cmd_output=[]
parasubs["$COMMANDSTATUS"]=0
v=[0]*max(1,len(parameters))
for i in range(numvar):
	if numvar>1:
		jobname=""
		for p in range(len(v)):
			paraname=parameters[p].name()
			jobname+=paraname
			paravalues=parameters[p].values()
			if v[p]>=parameters[p].numvalues():
				v[p+1]+=1
				v[p]=0
			paravalue=paravalues[v[p]]
			parasubs['$'+paraname]=paravalue
			if parameters[p].onlydigits():
				jobname+=paravalue.zfill(parameters[p].maxlen())
			else:
				jobname+=str(paravalue)
	else:
		if jobname=="": jobname="NO_parameters"
	parasubs["$JOBNAME"]=jobname
	# substitute all parameters within cmdcond & cmd
	#cmdcond=replaceparameters(gencmdcondition,parasubs)
	#cmd=replaceparameters(gencommand,parasubs)
	cmdcond=gencmdcondition
	cmd=gencommand
	for paraname,paravalue in parasubs.items():
		if cmdcond is not None:
			cmdcond=cmdcond.replace(paraname,str(paravalue))
		if cmd is not None:
			cmd=cmd.replace(paraname,str(paravalue))
	if cmdcond is not None:
		try:
			evalcond=eval(cmdcond)
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
			#substfilecontent=replaceparameters(substfilecontent,parasubs)
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
			print "-",time.strftime("%H:%M:%S"),"command:",cmd,'-'*(50-len(cmd))
			status,stdoutdata,stderrdata=execmd(cmd,jobdir)
			parasubs["$COMMANDSTATUS"]=status
			if status!=0:
				print stdoutdata
				print stderrdata
				print '-'*(56-len(str(status))),status,"failed",time.strftime("%H:%M:%S")
			else:
				print stdoutdata
				print '-'*58,"done",time.strftime("%H:%M:%S")
				cmd_output.append(stdoutdata.rstrip("\n"))
	else:
		print "condition",gencmdcondition,"does not hold after substitution:",cmdcond
	# check break condition
	if breakcondition is not None:
		breakcond=replaceparameters(breakcondition,parasubs)
		try:
			if eval(breakcond):
				print "break condition",breakcondition," fulfilled after substituion:",breakcond
				break
		except Exception, errno:
			print "ERROR evaluating break condition",breakcondition,"after substitution:",breakcond
			#break
			pass
	v[0]+=1

print "",time.strftime("%H:%M:%S"),"created",len(createdjobs),"jobs"

for p in parameters:
	#del parasubs['$'+p.name()]
	parasubs['$'+p.name()]=p.strvalues()
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
		print "-",time.strftime("%H:%M:%S"),"command:",ppcmd,'-'*(50-len(ppcmd))
		status,stdoutdata,stderrdata=execmd(ppcmd,dstroot)
		#parasubs["$COMMANDSTATUS"]=status
		if status!=0:
			print stdoutdata
			print stderrdata
			print '-'*(56-len(str(status))),status,"failed",time.strftime("%H:%M:%S")
		else:
			print stdoutdata
			print '-'*58,"done",time.strftime("%H:%M:%S")
			cmd_output.append(stdoutdata)
	print "################################################################### done"

print time.strftime("finished at %a, %d.%m.%Y %H:%M:%S %Z")
