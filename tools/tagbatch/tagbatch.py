#!/usr/bin/env python
# coding=utf-8

"""
   tagbatch.py
  
    Created on: 09.11.2015
        Author: Konstantin Kroschewski (kokr)
        Version: 0.3
"""

#TODO: Dateinamen bei -m ueberpruefen, ggf. postfix rein; 

import sys
import re
import odf.opendocument
from odf.table import Table, TableRow, TableCell
from odf.text import P
import copy

#for marking files as executable:
import os
import stat

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
                        for rr in range(int(repeat)):  # repeated?
                            arrCells[count]=textContent
                            count+=1
                    else:
                        row_comment = row_comment + textContent + " "
                else:
                    for rr in range(int(repeat)):
                        count+=1

            # if row contained something
            if(len(arrCells)):
				arrRows.append(arrCells)

            else:
				arrRows.append([None])
            #    print ("Empty or commented row (", row_comment, ")")

        self.SHEETS[name] = arrRows

    # returns a sheet as an array (rows) of arrays (columns)
    def getSheet(self, name):
        return self.SHEETS[name]


#main program

#read input parameters
tagOpen = u"@{"
tagClose = u"}"
tagFile = ""
outputName = ""
calcFile = ""
useSheet = ""
firstTag = ""
verbose = False
showHelp = False
oneOutputFile = False
makeBashScript = False


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


argc = len(sys.argv)

if argc > 1:
	
	ignoreNext = False
	for i in range(1, argc):
		
		if ignoreNext:
			ignoreNext = False
			continue

		if sys.argv[i] == "-h" or sys.argv[i] == "--help":
			showHelp = True
		
		elif sys.argv[i] == "-v":
			verbose = True
			
		elif sys.argv[i] == "-a":
			oneOutputFile = True
			
		elif sys.argv[i] == "-x":
			oneOutputFile = True
			makeBashScript = True
		
		elif sys.argv[i] == "-t":
			if i < (argc - 1):
				tagFile = sys.argv[i+1]
				ignoreNext = True
			else:
				raise AssertionError("no Parameter given after \"-t\"")
			
		elif sys.argv[i] == "-o":
			if i < (argc - 1):
				outputName = sys.argv[i+1]
				ignoreNext = True
			else:
				raise AssertionError("no Parameter given after \"-o\"")
			
		elif sys.argv[i] == "-c":
			if i < (argc - 1):
				calcFile = sys.argv[i+1]
				ignoreNext = True
			else:
				raise AssertionError("no Parameter given after \"-c\"")
				
		elif sys.argv[i] == "-s":
			if i < (argc - 1):
				useSheet = sys.argv[i+1]
				ignoreNext = True
			else:
				raise AssertionError("no Parameter given after \"-s\"")
				
		elif sys.argv[i] == "-f":
			if i < (argc - 1):
				firstTag = sys.argv[i+1]
				ignoreNext = True
			else:
				raise AssertionError("no Parameter given after \"-f\"")
				
		elif sys.argv[i] == "-ot":
			if i < (argc - 1):
				tagOpen = sys.argv[i+1]
				ignoreNext = True
			else:
				raise AssertionError("no Parameter given after \"-ot\"")
				
		elif sys.argv[i] == "-ct":
			if i < (argc - 1):
				tagClose = sys.argv[i+1]
				ignoreNext = True
			else:
				raise AssertionError("no Parameter given after \"-ct\"")
				
		else:
			raise AssertionError("Unknown Parameter:", sys.argv[i])
			
			
else:
	showHelp = True
	
if showHelp:
	print "tagbatch v0.3\n"
	print "usage:"
	print sys.argv[0], "-t <input-tagged-config-file> -c <ods-table-file> [-s <ods-table-sheet>]\n[-f <ods-table-first-tag>] [-o <output-config-file-prefix>] [-v] [--help] [-h]\n"
	print "Parameters:"
	print "-t <input-tagged-config-file>\tspecify tagged input config file"
	print "-c <ods-table-file>\t\tspecify the ods table file holding the batch values"
	print "-s <ods-table-sheet>\t\tspecify the sheet of the ods table to get the data from\n\t\t\t\t[optional: only needed if file has multiple sheets]"
	print "-f <ods-table-first-tag>\tspecify the first tag-name in the ods-table file\n\t\t\t\t[optional: otherwise first non-empty text cell will be\n\t\t\t\tconsidered as first-tag]"
	print "-o <output-file-prefix>\t\tspecify output name prefix\n\t\t\t\t[optional: otherwise input config file name will be used]"
	print "-a \t\t\t\tDo not create multiple files, write output to one file.\n\t\t\t\tOutput filename will be like specified at -o."
	print "-x \t\t\t\tSame behaviour like -a, but add bash-script-header to output file\n\t\t\t\tand mark file as executable."
	print "-ot <custom-open-tag> \t\tSet custom open tag. Standard is @{."
	print "-ct <custom-close-tag> \t\tSet custom close tag. Standard is }."
	print "-v\t\t\t\tverbose output"
	print "-h or --help\t\t\tshow this help-text\n"
	print "\nOutput file numeration will always have 4 digits. It starts at 0001.\n"
	print "A tag looks like this: @{tagname}\n"
	
else:
	
	assert tagFile != "", "input-tagged-config-file not specivied"
	assert calcFile != "", "ods-table-file not specivied"
	
	if outputName == "":
		if oneOutputFile:
			outputName = tagFile[0:tagFile.rfind(".")] + "_out"
		else:
			outputName = tagFile[0:tagFile.rfind(".")]
		
	
	if verbose:
		print "tagbatch v0.2 by kokr\n"
		print "input parameters:"
		print "-> tag-file:", tagFile
		print "-> output-name:", outputName
		print "-> table-file:", calcFile
		print "-> Tagformat:", tagOpen+"tagname"+tagClose, "\n"
	
	#check file existance
	assert fileExist(tagFile), "input-tagged-config-file not existant"
	assert fileExist(calcFile), "ods-table-file not existant"
	
	#read ods file to find valid tags and extract the tag values
	calcFile = unicode(calcFile)
	
	doc = ODSReader(calcFile, clonespannedcolumns=True)
	
	#find right sheet
	if useSheet != "":
		try:
			table = doc.getSheet(unicode(useSheet))
		except KeyError:
			raise AssertionError("sheet not found in ods-file:", useSheet)
	else:
		for key in doc.SHEETS.iterkeys():
			if len(doc.SHEETS[key]) > 1:
				if useSheet != "":
					raise AssertionError("multiple non-empty sheets found, please see --help")
				useSheet = key
		if useSheet == "":
			raise AssertionError("ods table file is empty")
		table = doc.getSheet(useSheet)
	
	#find start cell
	startCell = [None, None]
	for i in range(len(table)):
		for j in range(len(table[i])):
			if firstTag == "":
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
				if table[i][j] == unicode(firstTag):
					startCell = [i, j]
					break
		else:
			continue
		break
		
	if verbose:
		print "processing ods-table..."
		print "-> using ods-table start-cell: {0}".format(startCell) 
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
				
	if verbose:
		print "-> found {0} tags".format(len(tableTags))
	
	
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
				print "-> Warning: Tag \'{0}\' has less entries ({1}) than the other tag(s) ({2}). It will be ignored.".format(tableTags[bli][0], len(tableTags[bli])-1, blength-1)
				deleteTableTags.append(tableTags[bli])
		for di in deleteTableTags:
			tableTags.remove(di)
	
	if verbose:
		print "-> there are {0} rows which will be processed".format(blength-1)
	
	#find and read tags in tag-file 
	try:
		ftag = open(tagFile, "r")
		stag = ftag.readlines()
		
	except IOError as (errno, strerror):
		print "I/O error({0}): {1}".format(errno, strerror)
		raise
	except:
		print "Unexpected error:", sys.exc_info()[0]
		raise
	
	ftag.close()
	fileTags = []
	for li in range(len(stag)):
		topen = stag[li].find(tagOpen)
		while topen != -1:
			tclose = stag[li].find(tagClose, topen+len(tagOpen))
			if tclose == -1:
				raise AssertionError("Tag not closed. Line: {0}".format(li+1))
			elif stag[li].find(tagOpen, topen+len(tagOpen), tclose) != -1:
				raise AssertionError("Tag openend before last tag closed. Line: {0}".format(li+1))
			else:
				fileTags.append([stag[li][topen+len(tagOpen):tclose],li,topen,tclose])
			topen = stag[li].find(tagOpen, tclose+len(tagClose))
			
	assert len(fileTags) != 0, "The tag-file does not contain any tags. A tag looks like this: {tagname}"
	
	if verbose:
		print "\nprocessing tag-file:"
		print "-> found {0} tags".format(len(fileTags))
		print "\ncompare table and tag-file:"
	
	#give warning if tags are missing
	mTableTags = set()
	for mi in tableTags:
		mTableTags.add(mi[0])
		
	mFileTags = set()
	for mi in fileTags:
		mFileTags.add(unicode(mi[0]))

	mInter = mTableTags.intersection(mFileTags)
	assert len(mInter) != 0, "No tag in table belongs to any tag in tag-file."
	
	subTable = mTableTags.difference(mInter)
	subFile = mFileTags.difference(mInter)
	
	deleteTableTags = []
	for sT in subTable:
		print "-> Warning: Tag \'{0}\' only exists in table-file. It will be ignored.".format(sT)
		for ti in tableTags:
			if ti[0] == sT:
				deleteTableTags.append(ti) 
		
	deleteFileTags = []
	for sF in subFile:
		print "-> Warning: Tag \'{0}\' only exists in tag-file. It will be ignored.".format(sF)
		for fi in fileTags:
			if unicode(fi[0]) == sF:
				deleteFileTags.append(fi)
	
	for dT in deleteTableTags:
		tableTags.remove(dT)
		
	for dF in deleteFileTags:
		fileTags.remove(dF)
		
	if verbose:
		print "-> {0} tags are in both, the tag-file and the table-file:".format(len(fileTags)),
		for mIs in range(len(fileTags)):
			if mIs == 0:
				print " " + fileTags[mIs][0] + ", ",
			elif mIs == (len(fileTags)-1):
				print fileTags[mIs][0]
			else:
				print fileTags[mIs][0] + ", ",	
		print "\nwrite output files:"
	
	
	
	#create batches
	oSuffix = tagFile[tagFile.rfind("."):len(tagFile)]
	affectedLines = set()
	for aI in fileTags:
		affectedLines.add(aI[1])
		
	if oneOutputFile:
		if makeBashScript:
			bName = outputName + ".sh"
			
		else:
			bName = outputName + oSuffix
			
		try:
			ftag = open(bName, "w")
			
			if makeBashScript:
				ftag.write("#!/bin/bash\n\n")
				
			for bI in range(1,blength):
				for lI in range(len(stag)):
					if lI in affectedLines:
						rstag = stag[lI][:]
						rstag = unicode(rstag)
						for tI in tableTags:
							rstag = rstag.replace(tagOpen + tI[0] + tagClose, tI[bI])
						ftag.write(rstag.encode('utf-8'))
					else:
						ftag.write(stag[lI])
				ftag.write("\n")
			ftag.close()
			
			#mark file as executable if needed
			if makeBashScript:
				st = os.stat(bName)
				os.chmod(bName, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH )
				
			if verbose:
				print "-> {0} written".format(bName)
				
		except IOError as (errno, strerror):
			print "I/O error({0}): {1}".format(errno, strerror)
			raise
		except:
			print "Unexpected error:", sys.exc_info()[0]
			raise
		
	else:	
		for bI in range(1,blength):
			bName = outputName + str(bI).zfill(4) + oSuffix
			try:
				ftag = open(bName, "w")
				for lI in range(len(stag)):
					if lI in affectedLines:
						rstag = stag[lI][:]
						rstag = unicode(rstag)
						for tI in tableTags:
							rstag = rstag.replace(tagOpen + tI[0] + tagClose, tI[bI])
						ftag.write(rstag.encode('utf-8'))
					else:
						ftag.write(stag[lI])
				ftag.close()
				if verbose:
					print "-> {0} written".format(bName)
					
			except IOError as (errno, strerror):
				print "I/O error({0}): {1}".format(errno, strerror)
				raise
			except:
				print "Unexpected error:", sys.exc_info()[0]
				raise
			

	
	
