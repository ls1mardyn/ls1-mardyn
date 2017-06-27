#!/usr/bin/env python3
# coding=utf-8

'''
Assisted Parameter Editor (APE)

Created on 29.03.2017

@author: kokr
'''

from PyQt5.QtWidgets import (QDialog, QMainWindow, QApplication, QMessageBox, QFileDialog, QLabel, QFrame, QTreeView, QDialogButtonBox, QFormLayout, QMenu)
from PyQt5.QtCore import (QDir, QCoreApplication, QFile, QIODevice, pyqtSlot, QObject, QItemSelection, QItemSelectionModel, Qt, QFileInfo, QModelIndex)
from PyQt5.QtGui import (QIcon, QPixmap)
import sys
from data import apeui, selectpluginform, addsettingform, importlog, about
from data.Editors import SettingEditorContainer
from data import FileStruct
import os
from types import MethodType
import importlib
from data.SettingStruct import (SettingEntry, SettingTree, HeaderData, ChildData, SettingData, XMLSettingElementHelper, XMLChildElementHelper, FilterProxyModel, ContentXMLImporter, DeepSearchSettingTreeIndex)
from xml.dom.minidom import parse
import xml.dom.minidom
from platform import node

'''
This function replaces the method from the main treeView to deselect an item on click on whitespace.
'''
def mousePressEvent(self, event):
	if not self.indexAt(event.pos()).isValid():
		self.clearSelection()
	QTreeView.mousePressEvent(self, event)

'''
The AboutWidget class shows a widget with general program information.
'''
class AboutWidget(QDialog, about.Ui_AboutWidget):
	def __init__(self, icon, parent=None):
		super(AboutWidget, self).__init__(parent)
		self.setupUi(self)
		self.pushButton.clicked.connect(self.accept)
		self.label_6.setPixmap(QPixmap(icon))

'''
The LogShow class is the parent class that shows a log viewer for the import process.
'''
class LogShow(QDialog, importlog.Ui_ImportLog):
	def __init__(self, parent=None):
		super(LogShow, self).__init__(parent)
		self.setupUi(self)
		self.buttonBox.accepted.connect(self.accept)
		self.label.setWordWrap(True)

'''
The ImportLog class inherits the LogShow class for displaying config file import warnings and errors.
'''
class ImportLog(LogShow):
	def __init__(self, errorLog, error, fileName, parent=None):
		super(ImportLog, self).__init__(parent)
		
		if error:
			labelText = 'The import of the file "' + fileName + '" failed because at least one error occured. '
		else:
			labelText = 'The import of the file "' + fileName + '" succeeded but at least one warning occured. '
		
		labelText += "Check and correct the first warning entry, it may produce following warning and error entries."
		self.label.setText(labelText)
			
		self.plainTextEdit.clear()
		text = ""
		for entry in errorLog:
			text += entry + os.linesep + os.linesep
		self.plainTextEdit.insertPlainText(text)
		
'''
The ImportLog class inherits the LogShow class for displaying plugin import warnings and errors.
'''
class LoadPluginLog(LogShow):
	def __init__(self, errorLog, error, pluginName, contentXMLFile, parent=None):
		super(LoadPluginLog, self).__init__(parent)
		
		self.setWindowTitle("Plugin Import Log")
		
		if error:
			self.label.setText('The import of the plugin "' + pluginName + '" from File "' + contentXMLFile + '" failed because at least one error occured:')
		else:
			self.label.setText('The import of the plugin "' + pluginName + '" from File "' + contentXMLFile + '"  succeeded but at least one warning occured:')
			
		self.plainTextEdit.clear()
		text = ""
		for entry in errorLog:
			text += entry + os.linesep + os.linesep
		self.plainTextEdit.insertPlainText(text)

'''
The SelectPluginForm class shows a form to select a plugin when creating a new config file
'''
class SelectPluginForm(QDialog, selectpluginform.Ui_PluginSelect):
	def __init__(self, pluginContainer, parent=None):
		super(SelectPluginForm, self).__init__(parent)
		self.setupUi(self)
		
		self.buttonBox.rejected.connect(self.reject)
		self.buttonBox.accepted.connect(self.accept)
		
		self.buttonBox.button(QDialogButtonBox.Ok).setEnabled(False)
		self.comboBox.currentIndexChanged.connect(self.setIndex)
		
		self.selectedIndex = None
		
		
		#check the existing plugins and show them in the drop down menu
		if len(pluginContainer) != 0:
			for pluginEntry in pluginContainer:
				plugin = pluginEntry[0]
				contentXMLFile = pluginEntry[1].absoluteFilePath(plugin.getContentXML())
				contentXML = xml.dom.minidom.parse(contentXMLFile)
				ape = contentXML.documentElement
				for node in ape.childNodes:
					if node.nodeType == node.ELEMENT_NODE:
						if node.tagName == "general":
							for generalnode in node.childNodes:
								if generalnode.nodeType == generalnode.ELEMENT_NODE:
									if generalnode.tagName == "pluginname":
										pluginname = generalnode
										self.comboBox.addItem(pluginname.childNodes[0].data, pluginEntry)
										self.buttonBox.button(QDialogButtonBox.Ok).setEnabled(True)
		else:
			self.comboBox.addItem("<no plugins available>")
			
	def setIndex(self, idx):
		self.selectedIndex = idx


'''
The AddSettingsForm class provides a widget with a treeview to add new setting elements.
'''
class AddSettingsForm(QDialog, addsettingform.Ui_AddSettingForm):
	def __init__(self, parent=None):
		super(AddSettingsForm, self).__init__(parent)
		self.setupUi(self)
		
		self.buttonBox.rejected.connect(self.reject)
		self.buttonBox.accepted.connect(self.accept)
		
		self.buttonBox.button(QDialogButtonBox.Ok).setEnabled(False)
		
		self.model = None
		self.elementIndex = None
		
	def initModel(self, model):
		self.model = model
		self.proxyModel = FilterProxyModel(self)
		self.proxyModel.setSourceModel(model)
		self.treeView.setModel(self.proxyModel)
		self.treeView.setColumnHidden(1, True)
		self.treeView.selectionModel().selectionChanged.connect(self.setSelectedElement)
		
	def setSelectedElement(self, itemTo):
		indexList = itemTo.indexes()
		if len(indexList) != 0:
			self.elementIndex = self.proxyModel.mapToSource(indexList[0])
			self.buttonBox.button(QDialogButtonBox.Ok).setEnabled(True)
		else:
			self.buttonBox.button(QDialogButtonBox.Ok).setEnabled(False)
			self.elementIndex = None
	
	def setTreeState(self, scrollIndex, selectIndex):
		scrollIndex = self.proxyModel.mapFromSource(scrollIndex)
		selectIndex = self.proxyModel.mapFromSource(selectIndex)
		self.treeView.scrollTo(scrollIndex)
		self.treeView.selectionModel().select(QItemSelection(selectIndex, selectIndex), QItemSelectionModel.ClearAndSelect)


'''
The MainApp class shows the main window and uses all other classes.
'''
class MainApp(QMainWindow, apeui.Ui_MainWindow):
	def __init__(self, parent=None):
		super(MainApp, self).__init__(parent)
		self.setupUi(self)
		
		#init statusbar
		self.setStyleSheet("QStatusBar::item { border: 0px solid black }; ")
		self.statusBar = QLabel()
		self.statusbar.addWidget(self.statusBar)
		self.statusBar.setText("No configuration file loaded.")
		
		#object variables
		self.currentPlugin = False
		self.pluginContainer = None
		self.selectedElement = None
		self.unsavedContent = False
		self.editorContainer = None
		self.ignoreWidgetsToSettingsData = False
		self.fileHandler = None
		self.lastFile = ""
		
		#buttons
		self.button_add.setEnabled(False)
		self.button_del.setEnabled(False)
		
		#signals
		self.actionLoad_Configuration.triggered.connect(self.loadFileDialog)
		self.actionSave_Configuration.triggered.connect(self.saveFile)
		self.actionSave_Configuration_File_As.triggered.connect(self.saveFileAs)
		self.actionNew_Configuration_File.triggered.connect(self.newFile)
		self.actionExit.triggered.connect(self.closeEventButton)
		self.actionAbout_This_Program.triggered.connect(self.showAboutWidget)
		self.button_add.clicked.connect(self.showAddSettingWidget)
		self.button_del.clicked.connect(self.deleteSetting)
		
		#treemodel
		self.treeView.setContextMenuPolicy(Qt.CustomContextMenu)
		self.treeView.customContextMenuRequested.connect(self.treeViewOpenMenu)
		self.treeView.mousePressEvent = MethodType(mousePressEvent, self.treeView)
		self.activeModel = SettingTree(self)
		self.proxyModel = FilterProxyModel(self)
		self.proxyModel.setSourceModel(self.activeModel)
		self.treeView.setModel(self.proxyModel)
		self.treeView.setColumnHidden(1, True)
		self.completeModel = SettingTree(HeaderData(True), self)
		
		#plugin content
		self.contentXML = None
		
		#icon
		self.icon = ""
		dataDir = QDir(os.path.dirname(os.path.realpath(__file__)))
		dataDir.cd("data")
		if dataDir.exists("icon.png"):
			self.icon = dataDir.absoluteFilePath("icon.png")
			self.setWindowIcon(QIcon(self.icon))
			
		#programm configuration for window size and last filename
		self.programConfig = ProgramConfiguration(self)
		self.programConfig.loadConfiguration()
		
	def showAboutWidget(self):
		aboutWidget = AboutWidget(self.icon, self)
		aboutWidget.exec()
		
	def checkCollapsed(self):
		if not self.isSelectionVisible():
			self.treeView.selectionModel().clearSelection()
		
	def isSelectionVisible(self):
		isVisible = True
		if self.selectedElement != None:
			parentIndex = self.activeModel.parent(self.selectedElement)
			if parentIndex.isValid():
				for index in self.activeModel.makeIndexTreePath(parentIndex):
					if not self.treeView.isExpanded(self.proxyModel.mapFromSource(index)):
						isVisible = False
						break
		else: 
			isVisible = False
		return isVisible
		
	def treeViewOpenMenu(self, position):
		menu = QMenu()
		hasChilds = False
		isVisible = True
		if self.selectedElement != None:
			isVisible = self.isSelectionVisible()
			if isVisible:
				hasChilds = self.activeModel.rowCount(self.selectedElement) > 0
		else:
			isVisible = False
		if not isVisible:
			self.treeView.selectionModel().clearSelection()
			rowCount = self.activeModel.rowCount(QModelIndex())
			if rowCount > 0:
				for i in range(rowCount):
					if self.activeModel.rowCount(self.activeModel.index(i, 0, QModelIndex())) > 0:
						hasChilds = True
						break
		
		actionExpand = menu.addAction("Expand Tree", self.expandTree)
		actionExpand.setEnabled(hasChilds)
		actionCollapse = menu.addAction("Collapse Tree", self.collapseTree)
		actionCollapse.setEnabled(hasChilds)
		menu.addSeparator()
		actionAdd = menu.addAction("Add Setting", self.showAddSettingWidget)
		actionAdd.setEnabled(self.button_add.isEnabled())
		actionDel = menu.addAction("Remove Setting", self.deleteSetting)
		actionDel.setEnabled(self.button_del.isEnabled())
		menu.exec_(self.treeView.viewport().mapToGlobal(position))
		
	def selectTreeViewIndex(self, index):
		self.treeView.selectionModel().select(QItemSelection(index, index), QItemSelectionModel.ClearAndSelect)
		
	def expandCollapseTree(self, expand=True):
		if self.fileHandler != None:
			deepSearch = DeepSearchSettingTreeIndex(self.activeModel, self.selectedElement)
			for index in deepSearch.getElementGenerator(self.selectedElement != None):
				indexView = self.proxyModel.mapFromSource(index)
				if indexView.isValid():
					if self.activeModel.rowCount(index) > 0:
						if expand:
							self.treeView.expand(indexView)
						else:
							self.treeView.collapse(indexView)
		
	def expandTree(self):
		self.expandCollapseTree(True)
		
	def collapseTree(self):
		self.expandCollapseTree(False)
		
	def closeEventButton(self):
		self.closeEvent()
		
	def closeEvent(self, event=None):
		doExit = True
		if self.unsavedContent:
			result = QMessageBox.warning(self, "Unsaved Content", "There are unsaved changes. Do you want to save them?", QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel, QMessageBox.Save)
			if result == QMessageBox.Save:
				if not self.saveFile():
					doExit = False
			elif result == QMessageBox.Cancel:
				doExit = False
				
		if doExit:
			self.programConfig.writeConfiguration()
			if event != None:
				event.accept()
			else:
				QCoreApplication.instance().quit()
		else:
			if event != None:
				event.ignore()
	def newFile(self):
		ps = SelectPluginForm(self.pluginContainer)
		if ps.exec() == QDialog.Accepted:
			if ps.selectedIndex != None:
				self.closeConfiguration()
				self.currentPlugin = ps.comboBox.itemData(ps.selectedIndex)
				if self.loadPluginModel():
					self.activeModel.beginResetModel()
					self.activeModel.checkImportance(self.completeModel)
					self.activeModel.endResetModel()
					self.button_add.setEnabled(True)
					self.label.setText("Select a setting to edit it.")
					self.windowIndicator()
					
					
		
	def widgetToSettingData(self):
		retVal = True
		if self.editorContainer != None:
			for item in self.editorContainer.editorList:
				returnMessage = item.setValueToSettingData()
				if returnMessage != None:
					if type(returnMessage[0]) is str:
						QMessageBox.critical(self, "Error", returnMessage[0], QMessageBox.Ok, QMessageBox.NoButton)
						returnMessage[1].setFocus(Qt.PopupFocusReason)
						retVal = False
		return retVal
		
	def cleanEditorWidgets(self):
		for i in range(self.formLayout_2.count()):
			lItem = self.formLayout_2.takeAt(0)
			lItem.widget().setVisible(False)
			self.formLayout_2.removeItem(lItem)
		self.editorContainer = None
		
	def showEditorLabel(self):
		if self.formLayout_2.count() != 0:
			self.cleanEditorWidgets()
		self.formLayout_2.setWidget(0, QFormLayout.FieldRole, self.label)
		self.label.setVisible(True)
		
	def setSelectedElement(self, itemTo):
		indexList = itemTo.indexes()
		newIndex = None
		if len(indexList) != 0:
			newIndex = self.proxyModel.mapToSource(indexList[0])
		block = False
		if newIndex != None and self.selectedElement != None:
			if newIndex.internalPointer() is self.selectedElement.internalPointer():
				block = True
		
		if not block:
			transferResult = True
			if self.ignoreWidgetsToSettingsData:
				self.ignoreWidgetsToSettingsData = False
			else:
				transferResult = self.widgetToSettingData()
			if transferResult:
				if newIndex != None:
					self.selectedElement = newIndex
					self.button_del.setEnabled(True)
					self.cleanEditorWidgets()
					self.editorContainer = SettingEditorContainer(self.activeModel.childData(self.selectedElement), self.scrollAreaWidgetContents_4, self)
					for editor in self.editorContainer.editorList:
						for item in editor.widgetList:
							if item.label == None:
								self.formLayout_2.addRow(item.widget)
							else:
								self.formLayout_2.addRow(item.label, item.widget)
				else:
					self.selectedElement = None
					self.button_del.setEnabled(False)
					self.showEditorLabel()
			else:
				if self.selectedElement != None:
					proxyIndex = self.proxyModel.mapFromSource(self.selectedElement)
					self.treeView.scrollTo(proxyIndex)
					self.treeView.selectionModel().select(QItemSelection(proxyIndex, proxyIndex), QItemSelectionModel.ClearAndSelect)
		
	def deleteSetting(self):
		if self.selectedElement != None:
			self.ignoreWidgetsToSettingsData = True
			self.setUnsavedContent(True)
			parent = self.activeModel.parent(self.selectedElement)
			row = self.selectedElement.internalPointer().row()
			self.activeModel.removeRow(row, parent)
		
	def showAddSettingWidget(self):
		asw = AddSettingsForm()
		asw.initModel(self.completeModel)
		indexList = self.treeView.selectedIndexes()
		if len(indexList) != 0:
			indexActive = self.proxyModel.mapToSource(indexList[0])
			indexComplete = self.completeModel.indexOfEqual(indexActive)
			if self.completeModel.rowCount(indexComplete) != 0:
				scrollIndex = self.completeModel.index(0, 0, indexComplete)
			else:
				scrollIndex = indexComplete
			selectIndex = indexComplete
			asw.setTreeState(scrollIndex, selectIndex)
		else:
			indexActive = None
		
		if asw.exec() == QDialog.Accepted:
			selectedIndex = asw.elementIndex
			if selectedIndex != None:
				treePath = self.completeModel.makeTreePath(selectedIndex)
				newIndex = self.activeModel.insertPath(treePath, indexActive)
				if type(newIndex) is str:
					QMessageBox.critical(self, "Error", newIndex, QMessageBox.Ok, QMessageBox.NoButton)
				else:
					newIndex = self.proxyModel.mapFromSource(newIndex)
					self.setUnsavedContent(True)
					self.treeView.scrollTo(newIndex)
					self.treeView.selectionModel().select(QItemSelection(newIndex, newIndex), QItemSelectionModel.ClearAndSelect)
	
	def saveFileAs(self):
		retVal = False
		filename = QFileDialog.getSaveFileName(self, 'Save configuration file to', self.FilePath(self.lastFile))
		if filename[0] != "":
			retVal = self.saveFile(False, filename[0])
			self.fileHandler.fileName = filename[0]
			self.lastFile = self.fileHandler.fileName
			self.unsavedContent = True
			self.setUnsavedContent(False)
			self.windowIndicator()
		return retVal
	
	def saveFile(self, checked=False, fileName=""):
		saveOK = False
		if self.fileHandler != None:
			if fileName == "":
				fileName = self.fileHandler.fileName
			
			if fileName == "":
				return self.saveFileAs()
				
			if self.widgetToSettingData():
				if fileName != None:
					saveResult = self.fileHandler.saveFile(self.activeModel, fileName)
					if saveResult == None:
						self.setUnsavedContent(False)
						saveOK = True
					else:
						QMessageBox.critical(self, "Error", 'Unable to save to file "' + fileName + '". ' + saveResult, QMessageBox.Ok, QMessageBox.NoButton)
		return saveOK
	
	def closeConfiguration(self):
		self.showEditorLabel()
		self.completeModel.beginResetModel()
		self.completeModel = SettingTree(HeaderData(True), self)
		self.completeModel.endResetModel()
		self.activeModel.beginResetModel()
		self.activeModel = SettingTree(self)
		self.proxyModel = FilterProxyModel(self)
		self.proxyModel.setSourceModel(self.activeModel)
		self.treeView.setModel(self.proxyModel)
		self.treeView.setColumnHidden(1, True)
		self.activeModel.endResetModel()
		self.button_del.setEnabled(False)
		self.button_add.setEnabled(False)
		self.label.setText("To start open or create a configuration file.")
		self.unsavedContent = False
		self.fileHandler = None
		self.setWindowTitle("APE")
		self.statusBar.setText("No configuration file loaded.")
	
	def setUnsavedContent(self, status):
		self.unsavedContent = status
		self.windowIndicator()
				
	def windowIndicator(self):
		filename = ""
		if self.fileHandler != None:
			filename = self.fileHandler.fileName
		statusbarText = ""
		pluginname = self.completeModel.headerData(1, Qt.Horizontal, Qt.DisplayRole).pluginname
		if pluginname != "":
			statusbarText += "Plugin: " + pluginname + " "
			
		if self.fileHandler != None:
			if filename != "":
				statusbarText += "File: " + filename
				
		self.statusBar.setText(statusbarText)
		if filename != "":
			info = QFileInfo(filename)
			fn = info.fileName()
		else:
			fn = ""
		if fn != "":
			if self.unsavedContent:
				self.setWindowTitle("* " + fn + " - APE")
			else:
				self.setWindowTitle(fn + " - APE")
		else:
			if self.unsavedContent:
				self.setWindowTitle("* APE")
			else:
				self.setWindowTitle("APE")
		
	def FilePath(self, filename=""):
		path = QDir.currentPath()
		if filename == "":
			if self.fileHandler != None:
				filename = self.fileHandler.fileName
		if filename != "":
			fInfo = QFileInfo(filename)
			qPath = fInfo.absoluteDir()
			if qPath.exists():
				path = qPath.absolutePath()
		return path
	
	def loadFileDialog(self):
		doLoad = False
		if self.unsavedContent:
			reply = QMessageBox.question(self, "Warning", "There are unsaved changes. Do you want to save them?", QMessageBox.Yes|QMessageBox.No|QMessageBox.Abort, QMessageBox.NoButton)
			if reply == QMessageBox.Yes:
				self.saveFile()
				doLoad = not self.unsavedContent
			elif reply == QMessageBox.No:
				self.closeConfiguration()
				doLoad = True
		else:
			doLoad = True
		if doLoad:
			filename = QFileDialog.getOpenFileName(self, 'Open configuration file', self.FilePath(self.lastFile))
			if filename[0]:
				self.closeConfiguration()
				file = QFile(filename[0])
				file.open(QIODevice.ReadOnly | QIODevice.Text)
				currentPluginTmp = False
				
				for item in self.pluginContainer:
					file.seek(0)
					if item[0].checkMagic(file):
						currentPluginTmp = item
						break
				
				file.close()
				
				if currentPluginTmp == False:
					QMessageBox.critical(self, "Error", 'Error: Selected file does not seem to fit to any installed plugin. Check the "checkMagic" method inside the "__init__.py" file of the desired plugins folder.', QMessageBox.Ok, QMessageBox.NoButton)
				else:
					self.currentPlugin = currentPluginTmp
					if self.loadPluginModel():
						self.fileHandler.fileName = filename[0]
						self.lastFile = self.fileHandler.fileName
						self.activeModel.beginResetModel()
						importLog, importError = self.fileHandler.readFile(self.activeModel, self.completeModel, self.currentPlugin)
						importanceError = False
						if not importError:
							importanceLog, importanceError = self.activeModel.checkImportance(self.completeModel)
							for logItem in importanceLog:
								importLog.append(logItem)
						
						self.activeModel.endResetModel()
						if importanceError or importError:
							self.closeConfiguration()
						else:
							self.button_add.setEnabled(True)
							self.label.setText("Select a setting to edit it.")
							self.windowIndicator()
							
						if len(importLog) != 0:
							importLogWindow = ImportLog(importLog, importError, filename[0], self)
							importLogWindow.exec()
					else:
						QMessageBox.critical(self, "Error", "Error: Selected file can't be read.", QMessageBox.Ok, QMessageBox.NoButton)
						self.closeConfiguration()
				
	def loadPluginModel(self):
		settingTreeOk = False
		returnVal = False
		contentXMLFile = self.currentPlugin[1].absoluteFilePath(self.currentPlugin[0].getContentXML())
		try:
			self.contentXML = xml.dom.minidom.parse(contentXMLFile)
		except BaseException as e:
			QMessageBox.critical(self, 'Error', 'Error: Plugin XML-File at "' + contentXMLFile + '" can not be read: '+ str(e), QMessageBox.Ok, QMessageBox.NoButton)
			return False
		ape = self.contentXML.documentElement
		
		#check general information of the content xml file
		if ape.nodeName == "ape" and ape.getAttribute("version") == "1.0":
			general = None
			settingtree = None
			for node in ape.childNodes:
				if node.nodeName == "general":
					general = node
				if node.nodeName == "settingtree":
					settingtree = node
			if general != None and settingtree != None:
				pluginname = None
				fileformat = None
				for node in general.childNodes:
					if node.nodeName == "pluginname":
						pluginname = node
					if node.nodeName == "fileformat":
						fileformat = node
				if pluginname != None and fileformat != None:
					pluginname = pluginname.childNodes[0].data
					fileformat = fileformat.childNodes[0].data
					
					self.fileHandler = FileStruct.FileHandler(fileformat, self)
					if self.fileHandler.isValid():
						importer = ContentXMLImporter(settingtree, self.fileHandler, self)
						if importer.isValid():
							#setup treeview models
							self.completeModel.beginResetModel()
							self.completeModel = SettingTree(HeaderData(True, pluginname), self)
							self.activeModel.beginResetModel()
							self.activeModel = SettingTree(self)
							self.proxyModel = FilterProxyModel(self)
							self.proxyModel.setSourceModel(self.activeModel)
							self.treeView.setModel(self.proxyModel)
							self.treeView.setColumnHidden(1, True)
							self.treeView.selectionModel().selectionChanged.connect(self.setSelectedElement)
							self.treeView.collapsed.connect(self.checkCollapsed)
							self.activeModel.endResetModel()
							importer.importXML(self.completeModel)
							self.completeModel.endResetModel()
						returnVal = not importer.isError
						if len(importer.messages) != 0:
							loadPluginLogWindow = LoadPluginLog(importer.messages, importer.isError, pluginname, contentXMLFile, self)
							loadPluginLogWindow.exec()
		return returnVal

class PluginManager(object):
	def initPlugins(self):
		#scan for plugins and import them
		pgmPath = os.path.dirname(os.path.realpath(__file__))
		self.pluginContainer = []
		pluginMainPath = QDir(pgmPath)
		if pluginMainPath.cd("plugins"):
			pluginEntries = pluginMainPath.entryList(QDir.Dirs, QDir.NoSort)
			if "." in pluginEntries: pluginEntries.remove(".")
			if ".." in pluginEntries: pluginEntries.remove("..")
			for item in pluginEntries:
				itemImport = "plugins." + item
				module = importlib.import_module(itemImport)
				pluginPath = QDir(pluginMainPath)
				pluginPath.cd(item)
				newPlugin = module.Plugin()
				self.pluginContainer.append([newPlugin,pluginPath])
			return True
		else:
			return False

class ProgramConfiguration(object):
	def __init__(self, parent=None):
		self.parent = parent
		self.confFile = "ape.cfg"
		
	def loadConfiguration(self):
		if self.parent != None:
			try:
				lastfile = ""
				windowpositionx = ""
				windowpositiony = ""
				windowsizex = ""
				windowsizey = ""
				
				configDOM = xml.dom.minidom.parse(self.confFile)
				configElement = configDOM.documentElement
				if configElement.tagName == "apeconfig":
					for node in configElement.childNodes:
						if node.nodeType == node.ELEMENT_NODE:
							if node.tagName == "lastfile" and node.hasChildNodes():
									lastfile = node.childNodes[0].data
							elif node.tagName == "windowpositionx" and node.hasChildNodes():
								windowpositionx = node.childNodes[0].data
							elif node.tagName == "windowpositiony" and node.hasChildNodes():
								windowpositiony = node.childNodes[0].data
							elif node.tagName == "windowsizex" and node.hasChildNodes():
								windowsizex = node.childNodes[0].data
							elif node.tagName == "windowsizey" and node.hasChildNodes():
								windowsizey = node.childNodes[0].data
					
					if lastfile != "":
						self.parent.lastFile = lastfile
					if windowpositionx != "" and windowpositiony != "" and windowsizex != "" and windowsizey != "":
						self.parent.setGeometry(int(windowpositionx), int(windowpositiony), int(windowsizex), int(windowsizey))
			except Exception:
				pass
	
	def writeConfiguration(self):
		if self.parent != None:
			try:
				lastfile = self.parent.lastFile
				windowRect = self.parent.geometry()
				windowpositionx = str(windowRect.x())
				windowpositiony = str(windowRect.y())
				windowsizex = str(windowRect.width())
				windowsizey = str(windowRect.height())
				
				configDOM = xml.dom.minidom.Document()
				configElement = configDOM.createElement("apeconfig")
				configDOM.appendChild(configElement)
				
				element = configDOM.createElement("lastfile")
				configElement.appendChild(element)
				configText = configDOM.createTextNode(lastfile)
				element.appendChild(configText)
				
				element = configDOM.createElement("windowpositionx")
				configElement.appendChild(element)
				configText = configDOM.createTextNode(windowpositionx)
				element.appendChild(configText)
				
				element = configDOM.createElement("windowpositiony")
				configElement.appendChild(element)
				configText = configDOM.createTextNode(windowpositiony)
				element.appendChild(configText)
				
				element = configDOM.createElement("windowsizex")
				configElement.appendChild(element)
				configText = configDOM.createTextNode(windowsizex)
				element.appendChild(configText)
				
				element = configDOM.createElement("windowsizey")
				configElement.appendChild(element)
				configText = configDOM.createTextNode(windowsizey)
				element.appendChild(configText)
				
				prettyXML = configDOM.toprettyxml(indent="  ", newl='\n', encoding="UTF-8")
				file = QFile(self.confFile)
				if file.open(QIODevice.WriteOnly):
					file.write(prettyXML)
					file.close()
				configDOM.unlink()
				
			except Exception:
				pass

def main():
	app = QApplication(sys.argv)
	mw = MainApp()
	mw.show()
	pi = PluginManager()
	if pi.initPlugins():
		mw.pluginContainer = pi.pluginContainer
		
		sys.exit(app.exec_())
	else:
		QMessageBox.critical(mw, "Error", "Error: Plugin path not found. Directory 'plugins' is missing in script path.", QMessageBox.Ok, QMessageBox.NoButton)

if __name__ == '__main__':
	main()
