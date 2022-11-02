#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# genjobs.py
# generate jobs from a template job, subsituting parameters
#
# Martin Bernreuther <bernreuther@hlrs.de>, November 2011

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
#
#import logging

desc="generate&start benchmark runs"
version="2013.02.26"
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
logfile=None

logfhdl=None

if sys.version_info < (2, 6):
	#raise "ERROR: python version too old: "+str(sys.version_info)+"<2.6"
	sys.stderr.write("\n!!! WARNING: using old python version "+str(sys.version_info)+" < 2.6 !!!\n\n")
#else
#	print "using Python version {0}".format(sys.version_info)

progname=os.path.basename(sys.argv[0])
optparser = optparse.OptionParser(usage="usage: {0} [options] [configfile]\n{1}\n(version: {2}, author: {3})".format(progname,desc,version,author))	## %prog
optparser.add_option("-V", "--version", action="store_true", dest="print_version", default=False, help="print version and exit [default: {0}]".format("False"))	## %default
#optparser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="be more verbose [default: {0}]".format("False"))
#optparser.add_option("-c", "--config", type="string", dest="configfile", help="configuration file [default: {0}]".format(configfilename))
optparser.add_option("-r", "--root", type="string", dest="rootdir", help="root directory for templates [default: <directory of configuration file>]")
optparser.add_option("-d", "--dstroot", type="string", dest="dstroot", help="destination root directory [default: {0}]".format(dstroot))
optparser.add_option("-t", "--delimiter", type="string", dest="delimiter", help="delimiter for lists [default (hex): {0}]".format(''.join([hex(ord(c))[2:].rjust(2,'0') for c in delimiter])))
optparser.add_option("-g", "--gentemplate", type="string", dest="gentemplate", help="generator job template directory [default: {0}]".format(gentemplate))
optparser.add_option("-p", "--pptemplate", type="string", dest="pptemplate", help="postprocessing template directory [default: {0}]".format(pptemplate))
optparser.add_option("-l", "--logfile", type="string", dest="logfile", help="log file [default: {0}]".format(logfile))
(options, args) = optparser.parse_args()

if options.print_version:
	## print version and exit
	print version
	sys.exit(0)

parasubs={}

def replaceparameters(src,paradefs):
	s=src
	for paraname,paravalue in paradefs.items():
		s=s.replace(paraname,str(paravalue))
	return s

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

parasubs["$DSTROOTNAME"]=dstroot
if cfgparser.has_option("general","dstroot"):
	dstroot=cfgparser.get("general","dstroot")
	parasubs["$DSTROOTNAME"]=dstroot
if options.dstroot is not None:
	dstroot=options.dstroot
	parasubs["$DSTROOTNAME"]=dstroot
	dstroot=os.path.join(rootdir,dstroot)
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
parasubs["$GENTEMPLATENAME"]=gentemplate
gentemplate=os.path.join(rootdir,gentemplate)
if not os.path.isdir(gentemplate):
	sys.stderr.write("ERROR: template directory {0} does not exist\n".format(gentemplate))
	sys.exit(3)
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

if cfgparser.has_option("generator","logfile"):
	logfile=cfgparser.get("generator","logfile")
if options.logfile:
	logfile=options.logfile
if logfile:
	logfile=replaceparameters(logfile,parasubs)
	parasubs["$LOGFILENAME"]=logfile
	logfile=os.path.join(rootdir,logfile)
	try:
		logfhdl=open(logfile, 'w')
		parasubs["$LOGFILEPATH"]=os.path.realpath(logfile)
		print "log file:\t{0}".format(logfile)
	except IOError, (errno, strerror):
		print "Error ({0}) opening {1} for writing: {2}".format(errno,logfile,strerror)
		#logfhdl=None
		print "\tlogging deactivated"
		del parasubs["$LOGFILENAME"]
		#parasubs["$LOGFILENAME"]=""
		#parasubs["$LOGFILEPATH"]=""
#else:
#	parasubs["$LOGFILENAME"]=""
#	parasubs["$LOGFILEPATH"]=""
if logfhdl is not None:
	logfhdl.write("# log file:\t{0}\n".format(parasubs["$LOGFILEPATH"]))
	logfhdl.write("# label:\t{0}\n".format(label))
	logfhdl.write("# generator template:\t{0}\n".format(parasubs["$GENTEMPLATEPATH"]))
	logfhdl.write("# destination root:\t{0}\n".format(parasubs["$DSTROOTPATH"]))

if cfgparser.has_option("postproc","template"):
	pptemplate=cfgparser.get("postproc","template")
if options.pptemplate is not None:
	pptemplate=options.pptemplate
if pptemplate is not None:
	parasubs["$PPTEMPLATENAME"]=pptemplate
	pptemplate=os.path.join(rootdir,pptemplate)
	if not os.path.isdir(pptemplate):
		sys.stderr.write("ERROR: postprocessing template directory {0} does not exist\n".format(pptemplate))
		sys.exit(4)
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
	
	def __init__(self,parameter,values=None):
		if isinstance(parameter,Parameter):
			self.__name=parameter.__name
			#self.__name=parameter.name()
			self.__values=parameter.__values
			#self.__values=parameter.values()
			self.__maxlen=parameter.__maxlen
			#self.__maxlen=parameter.maxlen()
			self.__onlydigits=parameter.__onlydigits
			#self.__onlydigits=parameter.onlydigits()
			self.__index=parameter.__index
			self.__value=parameter.__value
		elif isinstance(parameter,str):
			self.__name=parameter
			self.__initialize(values)
		else:
			self.__name=str(parameter)
			self.__initialize(values)
	
	def __initialize(self,values=None):
		self.__values=[""]
		if isinstance(values,str) and values:
			self.__values=values.split()
		self.__value=self.__values[0]
		self.__maxlen=0
		self.__index=0
		self.__onlydigits=True
		for v in self.__values:
			if len(v)>self.__maxlen: self.__maxlen=len(v)
			if not v.isdigit(): self.__onlydigits=False
		
	
	def name(self):
		return self.__name
	
	def maxlen(self):
		return self.__maxlen
		
	def onlydigits(self):
		return self.__onlydigits
	
	def numvalues(self):
		return len(self.__values)
	
	def index(self):
		return self.__index
	
	def setindex(self,idx):
		self.__index=(idx)%self.numvalues()
	
	def values(self):
		return self.__values
	
	def value(self):
		#self.__value=self.__values[max(self.__index,0)]
		return self.__value
	
	def nextvalue(self):
		self.setindex(self.__index+1)
		self.__value=self.__values[self.__index]
		return self.value()
	
	def strvalues(self,delim=" "):
		return delim.join(self.__values)
	
	def strvalue(self,padding=False):
		if padding:
			if self.onlydigits:
				#return self.value().rjust(self.maxlen(),'0')
				return self.value().zfill(self.maxlen())
			else:
				return self.value().rjust(self.maxlen(),'_')
		else:
			return str(self.value())
	
	#def eval_value(self,parasubs):
	#	paravalue=self.value()
	#	if paravalue[0]=='=': paravalue=paravalue[1:]
	#	if parasubs is not None:
	#		paravalue=replaceparameters(paravalue,parasubs)
	#	return str(eval(paravalue))
	
	def __str__(self):
		return self.__name
	
	def __int__(self):
		return self.numvalues()
	


print "parameters:"
print "name\tvalues\t(width,int)"
parameters=[]
paraformulae=[]
numvar=0
jobname=""
if cfgparser.has_section("parameters"):
	for i in cfgparser.items("parameters"):
		if numvar==0: numvar=1
		p=Parameter(i[0],i[1])
		parasubs['$'+p.name()]=p.strvalues()
		if p.numvalues()>1:
			parameters.append(Parameter(p))
			numvar*=p.numvalues()
			print "{0}\t{1}\t({2},{3})\t*{4}".format(p.name(),p.values(),p.maxlen(),p.onlydigits(),p.numvalues())
		else:
			if p.value()[0]=='=':
				paraformulae.append(Parameter(p))
				print "{0}\t{1}\t(formula substitution)".format(p.name(),p.value())
				# initial assignment should not be the formula itself, so do something useful(?) here...
				#del parasubs['$'+p.name()]
				#parasubs['$'+p.name()]=""
				parasubs['$'+p.name()]="("+p.value()[1:]+")"	# works for formulas without variables
			else:
				jobname+=p.name()+p.value()
				print "{0}\t{1}\t(substitution only)".format(p.name(),p.value())
	print "number of jobs:\t{0}".format(numvar)
else:
	print "WARNING: ",configfile," does not have a \"parameters\" section"
	if logfhdl is not None:
		logfhdl.write("WARNING: configuration file {0} does not contain a \"parameters\" section\n".format(configfile))
	


def execmd(cmd,wd="."):
	global logfhdl
	#print "execute command \"{0}\" (working directory: {1})".format(cmd,wd)
	status=-1
	stdoutdata=""
	stderrdata=""
	try:
		try:
			#stdoutdata=subprocess.check_output(cmd,cwd=wd,shell=True)	# security risk?
			stdoutdata=subprocess.check_output(shlex.split(cmd),cwd=wd)	# shell scripts have to be started using sh
			#stdoutdata=subprocess.check_output(shlex.split(cmd),cwd=wd,shell=True)
		except subprocess.CalledProcessError,e:
			print "ERROR ({0}) executing {1}: {2}".format(e.returncode,cmd,e.output)
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
	if logfhdl is not None:
			#logfhdl.write("cd {0}; {1}; cd $OLDPWD # cmd rc={2}\n".format(wd,cmd,status))
			#logfhdl.write("cd {0}; {1}; cd - # cmd rc={2}\n".format(wd,cmd,status))
			logfhdl.write("(cd {0}; {1}) # cmd rc={2}\n".format(wd,cmd,status))
	return status,stdoutdata,stderrdata

def adjrlinks(dstdir,srcdir):
	""" adjust relative links
	"""
	global logfhdl
	#print "adapt relative links contained in directory {0} copied from {1}".format(dstdir,srcdir)
	dstdirqueue=["."]
	while len(dstdirqueue)>0:
		qdir=dstdirqueue.pop()
		dstqdir=os.path.normpath(os.path.join(dstdir,qdir))
		for dirname, dirnames, filenames in os.walk(dstqdir):
			for subdirname in dirnames:
				#dstsubdir=os.path.normpath(os.path.join(dstqdir,subdirname))
				dstdirqueue.append(subdirname)
			for filename in filenames:	# filenames also contains links to directories
				dstfile=os.path.normpath(os.path.join(dstqdir,filename))
				if os.path.islink(dstfile):	# correct (relative) links
					#dstlinkfile=os.path.realpath(dstfile)
					templatelinkfile=os.path.realpath(os.path.join(srcdir,qdir,filename))
					if dstfile!=templatelinkfile:
						os.remove(dstfile)
						if logfhdl is not None:
							logfhdl.write("rm {0}\n".format(dstfile))
						symlinkpath=os.path.relpath(templatelinkfile,dstqdir)
						os.symlink(symlinkpath,dstfile)
						if logfhdl is not None:
							logfhdl.write("ln -s {0} {1}\n".format(symlinkpath,dstfile))


if os.path.exists(dstroot):
	sys.stderr.write("ERROR: destination directory {0} already exists\n".format(dstroot))
	sys.exit(5)

print "creating directory {0}".format(dstroot)
if pptemplate is not None:
	print "\t(copying {0})".format(pptemplate)
	shutil.copytree(pptemplate,dstroot,symlinks=True)
	if logfhdl is not None:
		logfhdl.write("cp -a {0} {1}\n".format(pptemplate,dstroot))
	adjrlinks(dstroot,pptemplate)
else:
	#os.mkdir(dstroot)
	os.makedirs(dstroot)
	if logfhdl is not None:
		logfhdl.write("mkdir {0}\n".format(dstroot))
if not (os.path.isdir(dstroot) and os.access(dstroot, os.W_OK)):
	sys.stderr.write("ERROR: destination directory {0} not accessible\n".format(dstroot))
	sys.exit(6)

print

print "generate jobs ##########################################################"
createdjobs=[]
cmd=None
cmd_output=[]
#cmd_status=[]
cmd_failed=0
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
	# substitution&evaluation of formula parameters
	for p in paraformulae:
		paravalue=replaceparameters(p.value()[1:],parasubs)
		parasubs['$'+p.name()]=str(eval(paravalue))
		#print "${0} = {1}".format(p.name(),parasubs['$'+p.name()])
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
		print '=',i+1,'/',numvar," benchmark job directory:",jobdir,'='*(39-len(str(i+1))-len(str(numvar))-len(jobdir))
		if os.path.exists(jobdir):
			print "destination directory {0} already exists ... skipping".format(jobdir)
			continue
		shutil.copytree(gentemplate,jobdir,symlinks=True)
		if logfhdl is not None:
			logfhdl.write("cp -a {0} {1}\n".format(gentemplate,jobdir))
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
			print "-",time.strftime("%H:%M:%S"),"command:",cmd,'-'*(51-len(cmd))
			status,stdoutdata,stderrdata=execmd(cmd,jobdir)
			parasubs["$COMMANDSTATUS"]=status
			#cmd_status.append(status)
			if status!=0:
				print stdoutdata
				print stderrdata
				print '-'*(51-len(str(status))),"failed (",status,")",time.strftime("%H:%M:%S")
			else:
				print stdoutdata
				print '-'*58,"done",time.strftime("%H:%M:%S")
				cmd_output.append(stdoutdata.rstrip("\n"))
	else:
		print '-',i+1,'/',numvar," benchmark job:",jobname,'-'*(49-len(str(i+1))-len(str(numvar))-len(jobname))
		print "condition",gencmdcondition,"does not hold after substitution:",cmdcond,"is false"
	# check break condition
	if breakcondition is not None:
		breakcond=replaceparameters(breakcondition,parasubs)
		try:
			if eval(breakcond):
				print "break condition",breakcondition," fulfilled after substituion:",breakcond,"is true"
				break
		except Exception, errno:
			print "ERROR evaluating break condition",breakcondition,"after substitution:",breakcond,"is erroneous"
			#break
	v[0]+=1

print "-"*72
print "",time.strftime("%H:%M:%S"),"created",len(createdjobs),"jobs"
print "",len(cmd_output),"commands executed without errors"
if cmd is not None:
	#  cmd_failed=len(filter(lambda s: s>0, cmd_status))
	if cmd_failed>0:
		print " Errors detected for executed commands:",cmd_failed,"failed"
	else:
		print " No errors detected for executed commands"

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
		print "-",time.strftime("%H:%M:%S"),"command:",ppcmd,'-'*(51-len(ppcmd))
		status,stdoutdata,stderrdata=execmd(ppcmd,dstroot)
		#parasubs["$COMMANDSTATUS"]=status
		if status!=0:
			print stdoutdata
			print stderrdata
			print '-'*(51-len(str(status))),"failed (",status,")",time.strftime("%H:%M:%S")
		else:
			print stdoutdata
			print '-'*58,"done",time.strftime("%H:%M:%S")
			cmd_output.append(stdoutdata)
	print "################################################################### done"

	if logfhdl is not None: logfhdl.close()

print time.strftime("finished at %a, %d.%m.%Y %H:%M:%S %Z")
