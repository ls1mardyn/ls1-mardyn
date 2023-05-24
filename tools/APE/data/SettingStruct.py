'''
Created on 28.04.2017

@author: kokr
'''

from PyQt5.QtCore import (QAbstractItemModel, QModelIndex, Qt, QSortFilterProxyModel)
import copy

'''
The SettingEntry and the SettingTree classes are used to meet the needs for the model based treeView widgets.
It is used to organize the stored data of the plugin content.
The childData, settingData and headerData classes are used to store the actual internal plugin content data.
'''

class SettingEntry(object):

	def __init__(self, data, parent=None):
		self.parentItem = parent
		self.itemData = data
		self.childItems = []
		
	def appendChild(self, item, row=-1):
		if row == -1:
			self.childItems.append(item)
		else:
			self.childItems.insert(row, item)
		
	def child(self, row):
		return self.childItems[row]

	def childCount(self):
		return len(self.childItems)

	def columnCount(self):
		return len(self.itemData)

	def data(self, column):
		try:
			return self.itemData[column]
		except IndexError:
			return None

	def parent(self):
		return self.parentItem

	def row(self):
		if self.parentItem:
			return self.parentItem.childItems.index(self)

		return 0
	
	def setData(self, data, column):
		try:
			self.itemdata[column] = data
			return True
		except IndexError:
			return False
		
	def removeChild(self, i):
		del self.childItems[i]
		
	def addChildData(self, childData, row=-1):
		self.appendChild(SettingEntry([childData.__str__(), childData], self), row)
		
	def getSortedInsertRow(self, newItem):
		insertRow = 0
		sourceChildren = newItem.parent().childItems
		newItemRow = newItem.row()
		for sRow, sItem in enumerate(sourceChildren):
			if sRow > newItemRow:
				break
			else:
				try:
					while self.childItems[insertRow].__str__() == sItem.__str__():
						insertRow = insertRow + 1
						
				except IndexError:
					insertRow = self.childCount()
					break
					
			
		return insertRow
		
	def checkPlurality(self, newChild):
		result = None
		nPlural = newChild.data(1).plurality
		if nPlural == "alone" and self.childCount() != 0:
			result = "You try to add a setting entry that can only exist alone under the parent setting entry. Remove all silblings and try again."
		elif nPlural == "single":
			found = False
			for leaf in self.childItems:
				if leaf.__str__() == newChild.__str__():
					found = True
					break
			if found:
				result = "You try to add a second instance of an already existing setting entry. This setting entry does not support multiple instances. Remove the other instances and try again."
		return result
		
	def makeTreePath(self):
		treePath = [self]
		if treePath[-1].parent() != None:
			while treePath[-1].parent().parent() != None:
				treePath.append(treePath[-1].parent())
		return treePath
		
		
	def __str__(self):
		return self.itemData[0]
	
class SettingTree(QAbstractItemModel):
	def __init__(self, header=None, parent=None):
		super(SettingTree, self).__init__(parent)
		self.rootItem = SettingEntry(["Settings", header])
		#self.setupModelData(self.rootItem)

	def columnCount(self, parent):
		if parent.isValid():
			return parent.internalPointer().columnCount()
		else:
			return self.rootItem.columnCount()

	def data(self, index, role):
		if not index.isValid():
			return None

		if role != Qt.DisplayRole:
			
			return None

		item = index.internalPointer()

		return item.data(index.column())
	
	def childData(self, index):
		if not index.isValid():
			return None
		return index.internalPointer().data(1)

	def flags(self, index):
		if not index.isValid():
			return Qt.NoItemFlags

		return Qt.ItemIsEnabled | Qt.ItemIsSelectable

	def headerData(self, section, orientation, role):
		if orientation == Qt.Horizontal and role == Qt.DisplayRole:
			return self.rootItem.data(section)

		return None

	def index(self, row, column, parent):
		if not self.hasIndex(row, column, parent):
			return QModelIndex()

		if not parent.isValid():
			parentItem = self.rootItem
		else:
			parentItem = parent.internalPointer()

		childItem = parentItem.child(row)
		if childItem:
			return self.createIndex(row, column, childItem)
		else:
			return QModelIndex()

	def parent(self, index):
		if not index.isValid():
			return QModelIndex()

		childItem = index.internalPointer()
		parentItem = childItem.parent()

		if parentItem == self.rootItem:
			return QModelIndex()

		return self.createIndex(parentItem.row(), 0, parentItem)

	def rowCount(self, parent):
		if parent.column() > 0:
			return 0

		if not parent.isValid():
			parentItem = self.rootItem
		else:
			parentItem = parent.internalPointer()

		return parentItem.childCount()
	
	def removeRows(self, row, count, parent = QModelIndex()):
		returnVal = False
		rows = self.rowCount(parent)
		if (row + count) <= rows:
			if parent.isValid():
				parentItem = parent.internalPointer()
			else:
				parentItem = self.rootItem
			self.beginRemoveRows(parent, row, row+count-1)
			for i in reversed(range(row, row + count)):
				parentItem.removeChild(i)
			self.endRemoveRows()
			returnVal = True
		return returnVal
	
	def indexOfEqual(self, index):
		treePath = self.makeTreePath(index)
		localItem = self.rootItem
		while len(treePath) != 0:
			leaf = treePath.pop()
			found = False
			for child in localItem.childItems:
				if leaf.__str__() == child.__str__():
					localItem = child
					found = True
					break
			if found == False:
				break
		return self.createIndex(localItem.row(), 0, localItem)
	
	def insertPath(self, treePathSource, index, addAllFollowingSingleChilds=True):
		treePath = []
		for treePathElement in treePathSource:
			treePath.append(treePathElement)
		if addAllFollowingSingleChilds:
			while treePath[0].childCount() == 1:
				treePath.insert(0, treePath[0].child(0))
		dItem = self.rootItem
		while len(treePath) != 0:
			sItem = treePath.pop()
			found = False
			if len(treePath) != 0:
				possibleChilds = []
				for child in dItem.childItems:
					if child.__str__() == sItem.__str__():
						possibleChilds.append(child)
				if index != None:
					if len(possibleChilds) != 0:
						for child in possibleChilds:
							foundIndex = False
							treePathInsert = self.makeTreePath(index)
							for leaf in treePathInsert:
								if child is leaf:
									foundIndex = True
									break
							
							if foundIndex:
								dItem = child
								found = True
						if found == False:
							dItem = possibleChilds[-1]
							found = True
						
				else:
					if len(possibleChilds) != 0:
						dItem = possibleChilds[-1]
						found = True
			if found == False:
				pluralReturn = dItem.checkPlurality(sItem)
				if pluralReturn == None:
					if dItem is not self.rootItem:
						index = self.createIndex(dItem.row(), 0, dItem)
					else:
						index = QModelIndex()
					insertRow = dItem.getSortedInsertRow(sItem)
					self.beginInsertRows(index, insertRow, insertRow)
					dItem.addChildData(copy.deepcopy(sItem.data(1)), insertRow)
					self.endInsertRows()
					dItem = dItem.childItems[insertRow]
				else:
					return pluralReturn
		return self.createIndex(dItem.row(), 0, dItem)
		
	def makeTreePath(self, index):
		if index.isValid():
			return index.internalPointer().makeTreePath()
		else: 
			return self.rootItem.makeTreePath()
		
	def makeIndexTreePath(self, index):
		treePath = self.makeTreePath(index)
		indextreePath = []
		for item in treePath:
			newIndex = self.createIndex(item.row(), 0, item)
			indextreePath.append(newIndex)
		return indextreePath
	
	def checkImportance(self, compareModel):
		importanceLog = []
		error = False
		search = DeepSearchSettingTree(self)
		for elementTuple in search.getElementGenerator(True):
			element = elementTuple[1]
			if element is not self.rootItem:
				elementIndex = self.createIndex(element.row(), 0, element)
				compareElementIndex = compareModel.indexOfEqual(elementIndex)
				treePath = self.makeTreePath(elementIndex)
				treePathCompare = compareModel.makeTreePath(compareElementIndex)
				
				assert len(treePath) == len(treePathCompare), "There are elements in the active model, which can't be there, because they do not exist in the complete model. Spooky. (Or indexofEqual fails.)"
				compareElement = compareElementIndex.internalPointer()
			else:
				elementIndex = None
				treePathCompare = []
				compareElement = compareModel.rootItem
				
			for compareElementChild in compareElement.childItems:
				if compareElementChild.data(1).importance == "mandatory":
					found = False
					for elementChild in element.childItems:
						if elementChild.data(1).name == compareElementChild.data(1).name:
							found = True
							break
					if not found:
						treePathCompareLocal = copy.copy(treePathCompare)
						treePathCompareLocal.insert(0, compareElementChild)
						
						insertResult = self.insertPath(treePathCompareLocal, elementIndex)
						if type(insertResult) == str:
							message = 'Error: The Element "'
							for treePathCompareLocalElement in reversed(treePathCompareLocal):
								message += treePathCompareLocalElement.data(1).name
								if treePathCompareLocalElement is not treePathCompareLocal[0]:
									message += " -> "
							message += '" could not be added. It is marked as importance=mandatory, but plurality marks prevent integration: ' + insertResult
							error = True
						else:
							message = 'Warning: The Element "'
							for treePathCompareLocalElement in reversed(treePathCompareLocal):
								message += treePathCompareLocalElement.data(1).name
								if treePathCompareLocalElement is not treePathCompareLocal[0]:
									message += " -> "
							message += '" has been added. It is marked as importance=mandatory, but was not added by the input file. Values of it are set to default values.'
						importanceLog.append(message)
						
		return importanceLog, error
		
		
class ChildData(object):
	def __init__(self, name, plurality, importance, visibility, parent=None):
		self.plurality = plurality
		self.importance = importance
		self.name = name
		self.visibility = (visibility != "invisible")
		self.settings = []
		self.settingMarks = None
		
	def addSetting(self, settingObject):
		self.settings.append(settingObject)
		
	# only return settings whose dependencies are met
	def getDependSettings(self, withSpare=False):
		dSettings = []
		
		# go through every setting entry and check its dependencies
		for sEntry in self.settings:
			if sEntry.settingDependOk(self.settings):
				# copy over positive tested entries
				dSettings.append(sEntry)
			elif withSpare:
				dSettings.append(None)
				
		return dSettings
	
	def __str__(self):
		return self.name

class SettingData(object):
	def __init__(self, name, defaultvalue, varType, editor, editorOptions, description, dependList=None, fileInformation=None, parent=None):
		self.name = name
		self.defaultvalue = defaultvalue
		self.value = defaultvalue
		self.varType = varType
		self.editor = editor
		self.editorOptions = editorOptions
		self.description = description
		self.fileInformation = fileInformation
		if dependList != None:
			self.dependList = dependList
		else:
			self.dependList = []
		
	def settingDependOk(self, settings):
		for dependEntry in self.dependList:
			if not dependEntry.dependOk(settings):
				return False
		return True
		
	def __str__(self):
		return self.name
		
class HeaderData(object):
	def __init__(self, isDictionary = False, pluginname = ""):
		self.isDictionary = isDictionary
		self.pluginname = pluginname
		
class SettingElementDependEntry(object):
	def __init__(self):
		self.dependSettingName = ""
		self.dependSettingValues = []
		
	def dependOk(self, settings):
		for sEntry in settings:
			if sEntry.name == self.dependSettingName:
				if sEntry.settingDependOk(settings):
					for value in self.dependSettingValues:
						if sEntry.value == value:
							return True
		return False
		
	def isValid(self):
		return self.dependSettingName != "" and len(self.dependSettingValues) > 0
	
	def __str__(self):
		return self.dependSettingName + "\n" + str(self.dependSettingValues) + "\n"
'''
The XMLSettingElementHelper class reads setting elements of the plugin content xml file.
'''
class XMLSettingElementHelper(object):
	def __init__(self, element, fileHandler, parent=None):
		self.isError = False
		self.parent = parent
		self.element = element
		self.name = ""
		self.defaultvalue = ""
		self.type = ""
		self.editor = ""
		self.editoroptions = ""
		self.description = ""
		self.fileHandler = fileHandler
		self.fileInformation = fileHandler.newFileInformationObject()
		self.dependList = []
		
		
		for child in self.element.childNodes:
			if child.nodeType == child.ELEMENT_NODE:
				if child.nodeName == "name" and self.name == "" and child.hasChildNodes():
					self.name = child.childNodes[0].data
				elif child.nodeName == "defaultvalue" and self.defaultvalue == "" and child.hasChildNodes():
					self.defaultvalue = child.childNodes[0].data
				elif child.nodeName == "type" and self.type == "" and child.hasChildNodes():
					self.type = child.childNodes[0].data
				elif child.nodeName == "editor" and self.editor == "" and child.hasChildNodes():
					self.editor = child.childNodes[0].data
				elif child.nodeName == "editor-options" and self.editoroptions == "" and child.hasChildNodes():
					self.editoroptions = child.childNodes[0].data
				elif child.nodeName == "description" and self.description == "" and child.hasChildNodes():
					self.description = child.childNodes[0].data
				elif child.nodeName == "depends" and child.hasChildNodes():
					dependEntry = SettingElementDependEntry()
					for dependChild in child.childNodes:
						if dependChild.nodeType == dependChild.ELEMENT_NODE:
							if dependEntry.dependSettingName == "" and dependChild.nodeName == "settingname" and dependChild.hasChildNodes():
								dependEntry.dependSettingName = dependChild.childNodes[0].data
							elif dependChild.nodeName == "value" and dependChild.hasChildNodes():
								dependEntry.dependSettingValues.append(dependChild.childNodes[0].data)
					if dependEntry.isValid():
						self.dependList.append(dependEntry)
				else:
					fitting, message = self.fileInformation.parseSettingElement(child)
					self.isError = message != None
					if self.isError:
						if self.parent != None:
							parent.makeMessage(element, message, True)
							break
					if not fitting:
						if parent != None:
							parent.makeMessage(element, 'The Element "' + child.nodeName + '" is unknown.')
			
		if self.name == "":
			self.name = self.parent.settingNameGenerator()
			
	def makeSettingData(self):
		return SettingData(self.name, self.defaultvalue, self.type, self.editor, self.editoroptions, self.description, self.dependList, self.fileInformation)

'''
The XMLSettingElementHelper class reads child elements of the plugin content xml file.
'''
class XMLChildElementHelper(object):
	def __init__(self, element, fileHandler, parent=None):
		optionChecker = XMLOptionChecker()
		self.isError = False
		self.parent = parent
		self.element = element
		self.childs = []
		self.settings = []
		self.name = ""
		self.plurality = ""
		self.importance = ""
		self.description = ""
		self.visibility = ""
		self.fileHandler = fileHandler
		
		if element.hasAttribute("plurality"):
			self.plurality = element.getAttribute("plurality")
			checkMessage = optionChecker.checkOption("plurality", self.plurality)
			self.isError = checkMessage != None
			if self.isError:
				if self.parent != None:
					self.parent.makeMessage(element, checkMessage, True)
			
		if element.hasAttribute("importance"):
			self.importance = element.getAttribute("importance")
			checkMessage = optionChecker.checkOption("importance", self.importance)
			self.isError = checkMessage != None
			if self.isError:
				if self.parent != None:
					self.parent.makeMessage(element, checkMessage, True)
					
		if element.hasAttribute("visibility"):
			self.visibility = element.getAttribute("visibility")
			checkMessage = optionChecker.checkOption("visibility", self.visibility)
			self.isError = checkMessage != None
			if self.isError:
				if self.parent != None:
					self.parent.makeMessage(element, checkMessage, True)
		
		if not self.isError:
		
			for child in self.element.childNodes:
				if child.nodeType == child.ELEMENT_NODE:
					if child.nodeName == "name" and self.name == "" and child.hasChildNodes():
						self.name = child.childNodes[0].data
					elif child.nodeName == "child":
						self.childs.append(child)
					elif child.nodeName == "setting":
						self.settings.append(child)
					elif child.nodeName == "description" and child.hasChildNodes():
						self.description = child.childNodes[0].data
					else:
						if parent != None:
							parent.makeMessage(element, 'The Element "' + child.nodeName + '" is unknown.')
					
			if self.name == "":
				self.name = self.parent.childNameGenerator()
		
	def hasChilds(self):
		return len(self.childs) != 0
	
	def hasSettings(self):
		return len(self.settings) != 0
	
	def makeChildData(self):
		childData = ChildData(self.name, self.plurality, self.importance, self.visibility)
		if self.description != "":
			childData.addSetting(SettingData("", "", "", "textshow", "", self.description, None, None))
		if self.hasSettings():
			for settingXML in self.settings:
				setting = XMLSettingElementHelper(settingXML, self.fileHandler, self.parent)
				self.isError = setting.isError
				if not self.isError:
					childData.addSetting(setting.makeSettingData())
				else:
					childData = None
		return childData

class XMLOptionChecker(object):
	def checkOption(self, optionType, xmlString):
		if optionType == "plurality":
			options = "single", "alone", "multi"
		elif optionType == "importance":
			options = "mandatory", "optional"
		elif optionType == "visibility":
			options = "visible", "invisible"
		else:
			raise Exception('Unknown optionType!')
		
		found = False
		for option in options:
			if xmlString == option:
				found = True
		
		if found:
			return None
		else:
			message = 'The attribute "' + optionType + '" has a invalid value: "' + xmlString + '". Valid values are: "'
			for option in options:
				message += option + '"'
				if option != options[-1]:
					message += ', '
				else:
					message += '. '
			return message

'''
The FilterProxyModel class is used to hide as  invisible marked childs in the treeView widgets.
'''
class FilterProxyModel(QSortFilterProxyModel):
	def __init__(self, parent=None):
		super(FilterProxyModel, self).__init__(parent)
		
	def filterAcceptsRow(self, sourceRow, sourceParent):
		return self.sourceModel().childData(self.sourceModel().index(sourceRow, 0, sourceParent)).visibility
	
'''
The ContentXMLImporter class is used to transform the content xml file to an internal datastructure to work with
'''
class ContentXMLImporter(object):
	def __init__(self, settingtree, fileHandler, parent=None):
		self.settingtree = settingtree
		self.fileHandler = fileHandler
		self.tree = []
		self.childNameGeneratorIndex = 0
		self.settingNameGeneratorIndex = 0
		self.messages = []
		self.isError = False
		
		#prepare the initial tree stack content that is used to perform the deep search over the content xml file
		rootHelper = XMLChildElementHelper(self.settingtree, self.fileHandler, self)
		self.isError = rootHelper.isError
		
		if not self.isError:
			for child in reversed(rootHelper.childs):
				helper = XMLChildElementHelper(child, self.fileHandler, self)
				self.isError = helper.isError
				if not self.isError:
					self.tree.append([0, helper])
				else:
					break
			
	'''
	The childNameGenerator method will be used by the XMLChildElementHelper and XMLSettingElementHelper classes to enumerate unnamed childs entries.
	'''
	def childNameGenerator(self):
		name = "Unnamed_Child_"+str(self.childNameGeneratorIndex)+"_APE_Internal_Name"
		self.childNameGeneratorIndex += 1
		return name
	
	'''
	The childNameGenerator method will be used by the XMLChildElementHelper and XMLSettingElementHelper classes to enumerate unnamed setting entires.
	'''
	def settingNameGenerator(self):
		name = "Unnamed_Setting_"+str(self.settingNameGeneratorIndex)+"_APE_Internal_Name"
		self.settingNameGeneratorIndex += 1
		return name
			
	def isValid(self):
		return len(self.tree) != 0 and not self.isError
	

	'''
	The importXML class performs a deep search over the content xml entries and populates the internal SettingTree datastructure
	'''
	def importXML(self, model):
		item = model.rootItem
		oldLevel = 0	
		while len(self.tree) != 0:
			level, leaf = self.tree.pop()
			
			#item is root item
			if item.parent() == None:
				childData = leaf.makeChildData()
				self.isError = leaf.isError
				if not self.isError:
					item.addChildData(childData)
					item = item.child(0)
				else:
					break
			else:
				#this item is a sibling of the previous item
				if level == oldLevel:
					childData = leaf.makeChildData()
					self.isError = leaf.isError
					if not self.isError:
						item.parent().addChildData(childData)
						item = item.parent().child(-1)
					else:
						break
					
				#this item is a child of the previous item
				elif level > oldLevel:
					oldLevel = level
					childData = leaf.makeChildData()
					self.isError = leaf.isError
					if not self.isError:
						item.addChildData(childData)
						item = item.child(0)
					else:
						break
				
				#this item is a parent or grandparent or grandgrandparent or ... of the previous item
				elif level < oldLevel:
					while level < oldLevel:
						oldLevel = oldLevel - 1
						item = item.parent()
					childData = leaf.makeChildData()
					self.isError = leaf.isError
					if not self.isError:
						item.parent().addChildData(childData)
						item = item.parent().child(-1)
					else:
						break
					
			#append all childs of the current item to the tree stack
			for child in reversed(leaf.childs):
				helper = XMLChildElementHelper(child, self.fileHandler, self)
				self.isError = helper.isError
				if not self.isError:
					self.tree.append([level+1, helper])
				else:
					self.tree = []
					break
	
	'''
	The makeMessage method will be used by the XMLChildElementHelper and XMLSettingElementHelper classes to collect error and warning messages
	'''
	def makeMessage(self, node, text, error=False):
		#append text with path of node to self.messages
		currentNode = node
		nodeNameTree = []
		while currentNode.parentNode != None:
			nodeName2 = ""
			for element in currentNode.childNodes:
				if element.nodeName == "name" and element.hasChildNodes():
					nodeName2 = element.childNodes[0].data
					break
					
			if nodeName2 != "":
				nodeNameTree.append(currentNode.nodeName+"("+nodeName2+")")
			else:
				nodeNameTree.append(currentNode.nodeName)
			currentNode = currentNode.parentNode
			
		if error:
			string = "Error: "
		else:
			std::string = "Warning: "
			
		std::string += 'At XML-Path: "'
		for nodeName in reversed(nodeNameTree):
			std::string += nodeName
			if nodeName is not nodeNameTree[0]:
				std::string += '/'
		std::string += '" '
		std::string += text
		self.messages.append(std::string)
'''
The DeepSearchSettingTreeIndex class uses the DeepSearchSettingTree class to do the deep search with indexes as return values
'''
class DeepSearchSettingTreeIndex(object):
	def __init__(self, model, rootIndex=None):
		self.model = model
		rootItem = None
		if isinstance(rootIndex, QModelIndex):
			if rootIndex.isValid():
				rootItem = rootIndex.internalPointer()
		self.deepSearch = DeepSearchSettingTree(self.model, rootItem)
		
	def getElementGenerator(self, includeRoot=False):
		for element in self.deepSearch.getElementGenerator(includeRoot):
			yield self.model.createIndex(element[1].row(), 0, element[1])
		
'''
The DeepSearchSettingTree class performs a iterated deep search of the setting tree entries
'''
class DeepSearchSettingTree(object):
	def __init__(self, model, rootItem=None):
		self.model = model
		self.rootItem = rootItem
		assert isinstance(self.model, SettingTree), "DeepSearchSettingTree constructor needs an SettingTree object."
		if self.rootItem == None:
			self.rootItem = self.model.rootItem
	
	def getElementGenerator(self, includeRoot=False):
		if self.model != None:
			batch = []
			if includeRoot:
				batch.append([-1, self.rootItem])
			else:
				for child in reversed(self.rootItem.childItems):
					batch.append([0, child])
			while len(batch) != 0:
				level, leaf = batch.pop()
				yield level, leaf
				for child in reversed(leaf.childItems):
					batch.append([level+1, child])
				