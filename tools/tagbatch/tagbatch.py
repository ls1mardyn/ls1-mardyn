#!/usr/bin/env python3
# coding=utf-8

"""
   tagbatch.py
  
	Created on: 09.11.2015
		Author: Konstantin Kroschewski (kokr)
		Version: 0.4
"""

#TODO: Dateinamen bei -m ueberpruefen, ggf. postfix rein; 

import sys
#import re
import odf.opendocument
from odf.table import Table, TableRow, TableCell
from odf.text import P
#import copy

#for marking files as executable:
import os
import stat

import argparse

#ODSReader by Marco Conti

class GrowingList(list):
	def __setitem__(self, index, value):
		if index >= len(self):
			self.extend([None]*(index + 1 - len(self)))
		list.__setitem__(self, index, value)

class ODSReader:

	# loads the file
	def __init__(self, file, clonespannedcolumns=None):
		self.clonespannedcolumns = clonespannedcolumns
		self.doc = odf.opendocument.load(file)
		self.SHEETS = {}
		for sheet in self.doc.spreadsheet.getElementsByType(Table):
			self.readSheet(sheet)

	# reads a sheet in the sheet dictionary, storing each sheet as an
	# array (rows) of arrays (columns)
	def readSheet(self, sheet):
		name = sheet.getAttribute("name")
		rows = sheet.getElementsByType(TableRow)
		arrRows = []

		# for each row
		for row in rows:
			row_comment = ""
			arrCells = GrowingList()
			cells = row.getElementsByType(TableCell)

			# for each cell
			count = 0
			for cell in cells:
				# repeated value?
				repeat = cell.getAttribute("numbercolumnsrepeated")
				if(not repeat):
					repeat = 1
					spanned = int(cell.getAttribute('numbercolumnsspanned') or 0)
					# clone spanned cells
					if self.clonespannedcolumns is not None and spanned > 1:
						repeat = spanned

				ps = cell.getElementsByType(P)
				textContent = ""

				# for each text/text:span node
				for p in ps:
					for n in p.childNodes:
						if (n.nodeType == 1 and n.tagName == "text:span"):
							for c in n.childNodes:
								if (c.nodeType == 3):
									textContent = u'{}{}'.format(textContent, n.data)

						if (n.nodeType == 3):
							textContent = u'{}{}'.format(textContent, n.data)

				if(textContent):
					if(textContent[0] != "#"):  # ignore comments cells
						for _rr in range(int(repeat)):  # repeated?
							arrCells[count]=textContent
							count+=1
					else:
						row_comment = row_comment + textContent + " "
				else:
					for _rr in range(int(repeat)):
						count+=1

			# if row contained something
			if(len(arrCells)):
				arrRows.append(arrCells)

			else:
				arrRows.append([None])
			#	print("Empty or commented row (", row_comment, ")")

		self.SHEETS[name] = arrRows

	# returns a sheet as an array (rows) of arrays (columns)
	def getSheet(self, name):
		return self.SHEETS[name]


#main program

#read input parameters
#tagOpen = u"@{"
#tagClose = u"}"
#tagFile = ""
#outputName = ""
#calcFile = ""
#useSheet = ""
#firstTag = ""
#verbose = False
#showHelp = False
#oneOutputFile = False
#makeBashScript = False


def fileExist(file):
	fExist = False
	try:
		fh = open(file, 'r')
	except IOError:
		fExist = False
	else:
		fh.close()
		fExist = True
	return fExist

def optionParse():
	description = "tagbatch v0.4"
	usage = "%(prog)s -t <input-tagged-config-file> -c <ods-table-file> [-s <ods-table-sheet>] [-f <ods-table-first-tag>] [-o <output-config-file-prefix>] [-v] [--help] [-h]"
	epilog= "Output file numeration will always have 4 digits. It starts at 0001. A tag looks like this: @{tagname}"
	parser = argparse.ArgumentParser(usage=usage, description=description, epilog=epilog)
	
	parser.add_argument('-t', '--tagfile', required=True, default="", help='specify tagged input config file')
	parser.add_argument('-c', '--tablefile', required=True, default="", help='specify the ods table file holding the batch values')
	parser.add_argument('-s', '--tablesheet', required=False, default="", help='specify the sheet of the ods table to get the data from [optional: only needed if file has multiple sheets]')
	parser.add_argument('-f', '--firsttag', required=False, default="", help='specify the first tag-name in the ods-table file [optional: otherwise first non-empty text cell will be considered as first-tag]')
	parser.add_argument('-o', '--outputname', required=False, default="", help='specify output name prefix [optional: otherwise input config file name will be used]')
	parser.add_argument('-a', '--onefile', required=False, action="store_true", help='Do not create multiple files, write output to one file. Output filename will be like specified at -o.')
	parser.add_argument('-x', '--bashscript', required=False, action="store_true", help='Same behaviour like -a, but add bash-script-header to output file and mark file as executable.')
	parser.add_argument('-ot', '--customopentag', required=False, default=u"@{", help='Set custom open tag. Standard is @{.')
	parser.add_argument('-ct', '--customclosetag', required=False, default=u"}", help='Set custom close tag. Standard is }.')
	parser.add_argument('-v', '--verbose', required=False, action="store_true", help='verbose output')
	#parser.add_argument('-h', '--help', required=False, action="store_true", help='show this help-text')
	
	return parser.parse_args()

args = optionParse()

if args.bashscript:
	args.onefile = True

#assert args.tagfile != "", "input-tagged-config-file not specivied"
#assert calcFile != "", "ods-table-file not specivied"

if args.outputname == "":
	if args.onefile:
		args.outputname = args.tagfile[0:args.tagfile.rfind(".")] + "_out"
	else:
		args.outputname = args.tagfile[0:args.tagfile.rfind(".")]
	

if args.verbose:
	print("tagbatch v0.4 by kokr\n")
	print("input parameters:")
	print("-> tag-file:", args.tagfile)
	print("-> output-name:", args.outputname)
	print("-> table-file:", args.tablefile)
	print("-> Tagformat:", args.customopentag+"tagname"+args.customclosetag, "\n")

#check file existance
assert fileExist(args.tagfile), "input-tagged-config-file not existant"
assert fileExist(args.tablefile), "ods-table-file not existant"

#read ods file to find valid tags and extract the tag values
doc = ODSReader(args.tablefile, clonespannedcolumns=True)

#find right sheet
if args.tablesheet != "":
	try:
		table = doc.getSheet(args.tablesheet)
	except KeyError:
		raise AssertionError("sheet not found in ods-file:", args.tablesheet)
else:
	for key in doc.SHEETS.iterkeys():
		if len(doc.SHEETS[key]) > 1:
			if args.tablesheet != "":
				raise AssertionError("multiple non-empty sheets found, please see --help")
			args.tablesheet = key
	if args.tablesheet == "":
		raise AssertionError("ods table file is empty")
	table = doc.getSheet(args.tablesheet)

#find start cell
startCell = [None, None]
for i in range(len(table)):
	for j in range(len(table[i])):
		if args.firsttag == "":
			if table[i][j] != None:
				alpha = False
				for x in range(len(table[i][j])):
					c = table[i][j][x]
					if c != u'0' and c != u'1' and c != u'2' and c != u'3' and c != u'4' and c != u'5' and c != u'6' and c != u'7' and c != u'8' and c != u'9' and c != u',' and c != u'.':
						alpha = True
						break
				if alpha == True:
					startCell = [i, j]
					break
		else:
			if table[i][j] == args.firsttag:
				startCell = [i, j]
				break
	else:
		continue
	break
	
if args.verbose:
	print("processing ods-table...")
	print("-> using ods-table start-cell: {0}".format(startCell))
assert startCell != [None, None], "table sheet has no tags or start-cell not found"

#get values of inner data-table
tableTags = []
for tidx in range(startCell[1], len(table[startCell[0]])):
	if (table[startCell[0]][tidx] == None):
		break
	else:
		tableTags.append([table[startCell[0]][tidx]])
		for tidy in range(startCell[0]+1, len(table)):
			if tidx >= len(table[tidy]):
				break
			elif table[tidy][tidx] == None:
				break
			else:
				tableTags[-1].append(table[tidy][tidx])

#look for duplicate tags
for dupi in range(len(tableTags)):
	for dupi2 in range(dupi,len(tableTags)):
		if(tableTags[dupi][0] == tableTags[dupi2][0] and id(tableTags[dupi]) != id(tableTags[dupi2])):
			raise AssertionError("Duplicate Tag found in ods-file:", tableTags[dupi][0])
			
if args.verbose:
	print("-> found {0} tags".format(len(tableTags)))


#determine batch-length and ignore tags having not enought entries
blength = len(tableTags[0])
deleteTableTags = []
if len(tableTags) > 1:
	for bli in range(1,len(tableTags)):
		blength2 = len(tableTags[bli])
		if blength < blength2:
			blength = blength2
	for bli in range(0,len(tableTags)):
		if len(tableTags[bli]) < blength:
			print("-> Warning: Tag \'{0}\' has less entries ({1}) than the other tag(s) ({2}). It will be ignored.".format(tableTags[bli][0], len(tableTags[bli])-1, blength-1))
			deleteTableTags.append(tableTags[bli])
	for di in deleteTableTags:
		tableTags.remove(di)

if args.verbose:
	print("-> there are {0} rows which will be processed".format(blength-1))

#find and read tags in tag-file 
try:
	ftag = open(args.tagfile, "r")
	stag = ftag.readlines()
	
except IOError as err:
	errno, strerror = err.args
	print("I/O error({0}): {1}".format(errno, strerror))
	raise
except:
	print("Unexpected error:", sys.exc_info()[0])
	raise

ftag.close()
fileTags = []
for li in range(len(stag)):
	topen = stag[li].find(args.customopentag)
	while topen != -1:
		tclose = stag[li].find(args.customclosetag, topen+len(args.customopentag))
		if tclose == -1:
			raise AssertionError("Tag not closed. Line: {0}".format(li+1))
		elif stag[li].find(args.customopentag, topen+len(args.customopentag), tclose) != -1:
			raise AssertionError("Tag openend before last tag closed. Line: {0}".format(li+1))
		else:
			fileTags.append([stag[li][topen+len(args.customopentag):tclose],li,topen,tclose])
		topen = stag[li].find(args.customopentag, tclose+len(args.customclosetag))
		
assert len(fileTags) != 0, "The tag-file does not contain any tags. A tag looks like this: "+args.customopentag+"tagname"+args.customclosetag

if args.verbose:
	print("\nprocessing tag-file:")
	print("-> found {0} tags".format(len(fileTags)))
	print("\ncompare table and tag-file:")

#give warning if tags are missing
mTableTags = set()
for mi in tableTags:
	mTableTags.add(mi[0])
	
mFileTags = set()
for mi in fileTags:
	mFileTags.add(mi[0])

mInter = mTableTags.intersection(mFileTags)
assert len(mInter) != 0, "No tag in table belongs to any tag in tag-file."

subTable = mTableTags.difference(mInter)
subFile = mFileTags.difference(mInter)

deleteTableTags = []
for sT in subTable:
	print("-> Warning: Tag \'{0}\' only exists in table-file. It will be ignored.".format(sT))
	for ti in tableTags:
		if ti[0] == sT:
			deleteTableTags.append(ti) 
	
deleteFileTags = []
for sF in subFile:
	print("-> Warning: Tag \'{0}\' only exists in tag-file. It will be ignored.".format(sF))
	for fi in fileTags:
		if fi[0] == sF:
			deleteFileTags.append(fi)

for dT in deleteTableTags:
	tableTags.remove(dT)
	
for dF in deleteFileTags:
	fileTags.remove(dF)
	
if args.verbose:
	print("-> {0} tags are in both, the tag-file and the table-file: ".format(len(fileTags)), end='')
	for ctags in fileTags:
		if ctags != fileTags[-1]:
			print(ctags[0] + ", ", end='')
		else:
			print(ctags[0])
	print("\nwrite output files:")



#create batches
oSuffix = args.tagfile[args.tagfile.rfind("."):len(args.tagfile)]
affectedLines = set()
for aI in fileTags:
	affectedLines.add(aI[1])
	
if args.onefile:
	if args.bashscript:
		bName = args.outputname + ".sh"
		
	else:
		bName = args.outputname + oSuffix
		
	try:
		ftag = open(bName, "w")
		
		if args.bashscript:
			ftag.write("#!/bin/bash\n\n")
			
		for bI in range(1,blength):
			for lI in range(len(stag)):
				if lI in affectedLines:
					rstag = stag[lI][:]
					for tI in tableTags:
						rstag = rstag.replace(args.customopentag + tI[0] + args.customclosetag, tI[bI])
					ftag.write(rstag)
				else:
					ftag.write(stag[lI])
			ftag.write("\n")
		ftag.close()
		
		#mark file as executable if needed
		if args.bashscript:
			st = os.stat(bName)
			os.chmod(bName, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH )
			
		if args.verbose:
			print("-> {0} written".format(bName))
			
	except IOError as err:
		errno, strerror = err.args
		print("I/O error({0}): {1}".format(errno, strerror))
		raise
	except:
		print("Unexpected error:", sys.exc_info()[0])
		raise
	
else:	
	for bI in range(1,blength):
		bName = args.outputname + str(bI).zfill(4) + oSuffix
		try:
			ftag = open(bName, "w")
			for lI in range(len(stag)):
				if lI in affectedLines:
					rstag = stag[lI][:]
					for tI in tableTags:
						rstag = rstag.replace(args.customopentag + tI[0] + args.customclosetag, tI[bI])
					ftag.write(rstag)
				else:
					ftag.write(stag[lI])
			ftag.close()
			if args.verbose:
				print("-> {0} written".format(bName))
				
		except IOError as err:
			errno, strerror = err.args
			print("I/O error({0}): {1}".format(errno, strerror))
			raise
		except:
			print("Unexpected error:", sys.exc_info()[0])
			raise
		



