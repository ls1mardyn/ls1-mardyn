'''
Created on 19.05.2017

@author: kokr
'''

'''
This file contains classes to manage the config file I/O.
To add other configfile formats beside xml, add a FileHandlerYourFileformat and a
FileInformationYourFileformat class and link to the first in the FileHandler class.
'''

from . import SettingStruct
from . import Editors
from xml.dom.minidom import Document
import xml.dom.minidom
from PyQt5.QtCore import (QFile, QIODevice)


'''
The FileHandler class stores information about the opened configuration file and instances specific classes to read and write them.
'''
class FileHandler(object):
	def __init__(self, fileformat, parent=None):
		self.specificHandler = None
		self.fileformat = fileformat
		self.fileInformationObjects = []
		self.fileName = ""
		
		#this links the given fileformat to a specific file handler class to handle it
		if self.fileformat == "xml":
			self.specificHandler = FileHandlerXML(self)
	
	def isValid(self):
		if self.specificHandler != None:
			return self.specificHandler.isValid()
		else:
			return False
		
	def readFile(self, targetModel, libraryModel, plugin):
		return self.specificHandler.readFile(targetModel, libraryModel, self.fileName, plugin)
	
	def saveFile(self, sourceModel, filename=""):
		if filename == "":
			filename = self.fileName
		return self.specificHandler.saveFile(sourceModel, filename)
	
	def newFileInformationObject(self):
		return self.specificHandler.newFileInformationObject()
	
	def getFileExtension(self):
		if self.specificHandler != None:
			return self.specificHandler.fileExtension + ";;All files (*.*)"
		else:
			return "All files (*.*)"
	
	fileExtension = property(getFileExtension)

'''
The FileHandlerXML class reads and stores config file in the xml format
'''
class FileHandlerXML(object):
	def __init__(self, parent=None):
		self.fileExtension = "XML (*.xml)"
	
	def isValid(self):
		return True
	
	def newFileInformationObject(self):
		return FileInformationXML()
	
	def saveFile(self, sourceModel, filename):
		saveResult = None
		search = SettingStruct.DeepSearchSettingTree(sourceModel)
		doc = Document()
		for childTreeEntryTupel in search.getElementGenerator():
			settings = childTreeEntryTupel[1].data(1).getDependSettings()
			for setting in settings:
				fileInfo = setting.fileInformation
				element = doc
				for path in fileInfo.xmlpath:
					found = False
					for node in reversed(element.childNodes):
						if node.nodeType == node.ELEMENT_NODE:
							if node.tagName == path:
								if not path is fileInfo.xmlpath[-1]:
									found = True
									element = node
								else:
									if fileInfo.attribute != "":
										if not node.hasAttribute(fileInfo.attribute):
											found = True
											element = node
									if len(node.childNodes) == 0:
										found = True
										element = node
								break
					if not found:
						newElement = doc.createElement(path)
						element.appendChild(newElement)
						element = newElement
				
				if fileInfo.attribute != "":
					element.setAttribute(fileInfo.attribute, setting.value)
				else:
					if setting.value != "":
						textNode = doc.createTextNode(setting.value)
						element.appendChild(textNode)
		prettyXML = doc.toprettyxml(indent="  ", newl='\n', encoding="UTF-8")
		try:
			file = QFile(filename)
			if file.open(QIODevice.WriteOnly):
				file.write(prettyXML)
				file.close()
			else:
				saveResult = ""
		except IOError as e:
			saveResult = str(e)
		finally:
			doc.unlink()
		return saveResult
	
	def readFile(self, targetModel, libraryModel, fileName, plugin):
		errorMessageLog = []
		error = False
		libModelLookupList = []
		libModelSearch = SettingStruct.DeepSearchSettingTree(libraryModel)
		for nextLibModelElement in libModelSearch.getElementGenerator():
			libModelLeaf = nextLibModelElement[1]
			for settingItemIndex, settingItem in enumerate(libModelLeaf.data(1).settings):
				if settingItem.fileInformation != None:
					if settingItem.fileInformation.isValid():
						libModelLookupList.append([settingItem.fileInformation, settingItemIndex, libModelLeaf, settingItem])
		parsedXML = xml.dom.minidom.parse(fileName)
		rootElement = parsedXML.documentElement
		nodeList = [[0, rootElement]]
		currentPath = []
		lastAddedTargetModelChildIndex = None
		editorProbe = Editors.EditorProbe()
			
		while len(nodeList) != 0:
			level, XMLLeaf = nodeList.pop()
			maxPathIdx = len(currentPath)-1
			if maxPathIdx < level:
				currentPath.append(XMLLeaf.tagName)
			elif maxPathIdx == level:
				currentPath[level] = XMLLeaf.tagName
			elif maxPathIdx > level:
				for i in range(maxPathIdx-level):
					del currentPath[maxPathIdx-i]
				currentPath[level] = XMLLeaf.tagName
					
			for node in reversed(XMLLeaf.childNodes):
				if node.nodeType == node.ELEMENT_NODE:
					nodeList.append([level+1, node])
			
			#XMLLeafCheckList Entry: [currentPath, "<if any: attribute name>", "value"]
			XMLLeafCheckList = []
			
			#check if the leaf has child elements
			hasChildElements = False
			for item in XMLLeaf.childNodes:
				if item.nodeType == XMLLeaf.ELEMENT_NODE:
					hasChildElements = True
					break
				
			#consider this element as value element if it has no child elements
			if not hasChildElements:
				if XMLLeaf.hasChildNodes():
					XMLLeafCheckList.append([currentPath, "", XMLLeaf.childNodes[0].data])
			
			
			if XMLLeaf.hasAttributes():
				attr = XMLLeaf.attributes
				for i in range(attr.length):
					attribute = attr.item(i)
					XMLLeafCheckList.append([currentPath, attribute.name, attribute.value])
					
			targetModelSearch = SettingStruct.DeepSearchSettingTree(targetModel)
			for XMLLeafCheckListEntry in XMLLeafCheckList:
				
				found = False
				for nextTargetModelElement in targetModelSearch.getElementGenerator():
					targetModelLeafData = nextTargetModelElement[1].data(1)
					if targetModelLeafData.settingMarks != None:
						for settingItemIndex, settingItem in enumerate(targetModelLeafData.settings):
							if settingItem.fileInformation != None:
								if settingItem.fileInformation.isValid():
									if self.isPathSimilar(XMLLeafCheckListEntry[0], settingItem.fileInformation.xmlpath):
										if XMLLeafCheckListEntry[1] == settingItem.fileInformation.attribute:
											settingFits = False
											if editorProbe.isReadOnly(settingItem):
												if XMLLeafCheckListEntry[2] == settingItem.defaultvalue:
													settingFits = True
											else:
												settingFits = True
											
											if settingFits:
												if targetModelLeafData.settingMarks[settingItemIndex]:
													#this setting already has been read
													errorMessageLog.append(self.generateErrorMessage(XMLLeafCheckListEntry[0], XMLLeafCheckListEntry[1], XMLLeafCheckListEntry[2], targetModelLeafData, "child setting not complete and double read of setting, other values will be reset to default."))
													targetModelLeafData.settingMarks = None	
													found = True
													break
												else:
													#this setting has not been read
													targetModelLeafData.settingMarks[settingItemIndex] = True
													
													#do checks here if value is valid
													settingItem.value = XMLLeafCheckListEntry[2]
													probeResult = editorProbe.isValidValue(settingItem)
													if probeResult != None:
														defaultvalue = settingItem.defaultvalue
														errorMessageLog.append(self.generateErrorMessage(XMLLeafCheckListEntry[0], XMLLeafCheckListEntry[1], XMLLeafCheckListEntry[2], targetModelLeafData, probeResult[0] + ' Value will be set to default value: "' + defaultvalue + '".'))
														settingItem.value = defaultvalue
													
													#mark this child as complete if complete
													found = True
													if self.allSettingsRead(targetModelLeafData):
														targetModelLeafData.settingMarks = None
														break
													
						if found:
							break
				
				if not found:
					#copy child with setting from the completemodel to the activemodel
					for libModelLookupListEntry in libModelLookupList:						
						if self.isPathSimilar(XMLLeafCheckListEntry[0], libModelLookupListEntry[0].xmlpath):
							if XMLLeafCheckListEntry[1] == libModelLookupListEntry[0].attribute:
								
								settingFits = False
								if editorProbe.isReadOnly(libModelLookupListEntry[3]):
									if XMLLeafCheckListEntry[2] == libModelLookupListEntry[3].defaultvalue:
										settingFits = True
								else:
									settingFits = True
								
								if settingFits:
								
									targetModelToAddTreePath = libModelLookupListEntry[2].makeTreePath()
									
									insertResult = targetModel.insertPath(targetModelToAddTreePath, lastAddedTargetModelChildIndex, False)
									if type(insertResult) is str:
										errorMessageLog.append(self.generateErrorMessage(XMLLeafCheckListEntry[0], XMLLeafCheckListEntry[1], XMLLeafCheckListEntry[2], libModelLookupListEntry[2].data(1), insertResult, True))
										error = True
										break
									else:
										lastAddedTargetModelChildIndex = insertResult
										found = True
										newChildDataEntry = insertResult.internalPointer().data(1)
										newChildDataEntry.settingMarks = [False] * len(newChildDataEntry.settings)
										
										#workaround for having xml silbings which cant be seperated by a special attribute value
										for settingIndex, settingElement in enumerate(newChildDataEntry.settings):
											if settingElement.editor == "" and settingElement.defaultvalue == "":
												newChildDataEntry.settingMarks[settingIndex] = True
										
										newChildDataEntry.settings[libModelLookupListEntry[1]].value = XMLLeafCheckListEntry[2]
										newChildDataEntry.settingMarks[libModelLookupListEntry[1]] = True
										probeResult = editorProbe.isValidValue(newChildDataEntry.settings[libModelLookupListEntry[1]])
										if probeResult != None:
											defaultvalue = newChildDataEntry.settings[libModelLookupListEntry[1]].defaultvalue
											errorMessageLog.append(self.generateErrorMessage(XMLLeafCheckListEntry[0], XMLLeafCheckListEntry[1], XMLLeafCheckListEntry[2], libModelLookupListEntry[2], probeResult[0] + ' Value will be set to default value: "' + defaultvalue + '".'))
											newChildDataEntry.settings[libModelLookupListEntry[1]].value = defaultvalue
										break
				
				if error:
					break
				
				if not found:
					errorMessageLog.append(self.generateErrorMessage(XMLLeafCheckListEntry[0], XMLLeafCheckListEntry[1], XMLLeafCheckListEntry[2], None, "No corresponding SettingTree entry has been found. Please add it to "+plugin[1].absoluteFilePath("content.xml")+". This setting will be discarded."))
				
			if error:
				break
		
		# check if all settings were written to the active tree and no setting is missing
		targetModelSearch = SettingStruct.DeepSearchSettingTree(targetModel)
		for targetModelItem in targetModelSearch.getElementGenerator():
			childData = targetModelItem[1].data(1)
			
			if childData.settingMarks != None:
				# distribute value data to all unwritten setting entries with the same xmlpath
				for idx, mark in enumerate(childData.settingMarks):
					if not mark:
						if childData.settings[idx].fileInformation.isValid():
							xmlpath = childData.settings[idx].fileInformation.xmlpath
							for mIdx, mSetting in enumerate(childData.settings):
								if childData.settingMarks[mIdx]:
									if self.isPathSimilar(xmlpath, mSetting.fileInformation.xmlpath):
										childData.settings[idx].value = mSetting.value
										childData.settingMarks[idx] = True
				
			# get settings whose dependencies are met
			try:
				sparedDependSettings = childData.getDependSettings(True)
				
				if childData.settingMarks == None:
					childData.settingMarks = [True] * len(childData.settings)
				
				for idx, mark in enumerate(childData.settingMarks):
					
					# reset values of settings whose dependencies are not met
					if sparedDependSettings[idx] == None:
						childData.settings[idx].value = childData.settings[idx].defaultvalue
					
					# thow error message if important setting is missing
					if not mark:
						settingItem = sparedDependSettings[idx]
						if settingItem != None:
							errorMessageLog.append(self.generateErrorMessage(None, None, None, childData, 'The element at '+settingItem.fileInformation.__str__()+'Default value: "'+settingItem.defaultvalue+'" does not exist in the XML-File. Default value will be used.'))
				childData.settingMarks = None
				
			except RecursionError:
				errorMessageLog.append(self.generateErrorMessage(None, None, None, childData, "This SettingTreeElement has settings which depend on each other in a loop. This is not possible. Please edit "+plugin[1].absoluteFilePath("content.xml"), True))
				error = True
				break

		
		return errorMessageLog, error
	
	def generateErrorMessage(self, xmlpath, attribute, value, settingTreeItem, errorString, isError=False):
		if isError:
			message = 'Error: '
		else:
			message = 'Warning: '
		if xmlpath != None:
			if len(xmlpath) != 0:
				message += 'At XML-Path: "'
				for element in xmlpath:
					message += element
					if element is not xmlpath[-1]:
						message += '/'
				message += '" '
		if attribute != None:
			if attribute != '':
				message += 'Attribute: "' + attribute + '" '
		if value != None:
			if value != '':
				message += 'Value: "' + value + '" '
		if isinstance(settingTreeItem, SettingStruct.SettingData):
			message += 'SettingTreeElement: "' + settingTreeItem.name + '" '
		elif isinstance(settingTreeItem, SettingStruct.SettingEntry):
			message += 'SettingTreeElement: "' + settingTreeItem.__str__() + '" '
		message += errorString
		return message
	
	def isPathSimilar(self, path1, path2):
		similar = True
		if len(path1) == len(path2):
			for i in range(len(path1)):
				if path1[i] != path2[i]:
					similar = False
					break
		else:
			similar = False
		return similar
	
	def allSettingsRead(self, childData):
		read = True
		for mark in childData.settingMarks:
			if mark == False:
				read = False
				break
		return read

'''
The FileInformationXML class is used to store file format specific content to the SettingData to map
the setting to an xml element of the config file.
'''
class FileInformationXML(object):
	def __init__(self, parent=None):
		self.xmlpath = []
		self.attribute = ""
		
	def isValid(self):
		return len(self.xmlpath) > 0
	
	def parseSettingElement(self, child):
		message = None
		fittingEntry = child.nodeName == "filexmlpath"
		if fittingEntry and len(self.xmlpath) == 0 and child.hasChildNodes():
			validEntry = True
			filexmlpath = child.childNodes[0].data
			pathList = filexmlpath.split("/")
			for pathElement in pathList:
				spaceCount = pathElement.count("@")
				if spaceCount != 0:
					if pathElement is not pathList[-1] or spaceCount != 1:
						validEntry = False
						message = 'Invalid filexmlpath content: "' + filexmlpath + '". Found "@" character in path more than once. This should only occur once at the last path element to declare a attribute.'
			if validEntry:
				lastElements = pathList[-1].split("@")
				if len(lastElements) == 2:
					self.attribute = lastElements[1]
					pathList[-1] = lastElements[0]
					self.xmlpath = pathList
				else:
					self.xmlpath = pathList
		return fittingEntry, message
	
	def __str__(self):
		if self.isValid():
			std::string = 'XML-Path: "'
			for element in self.xmlpath:
				std::string += element
				if element is not self.xmlpath[-1]:
					std::string += '/'
			std::string += '" '
			if self.attribute != "":
				std::string += 'Attribute: "' + self.attribute + '" '
			return string
		else:
			return ""
	