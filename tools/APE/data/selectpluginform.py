# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'SelectPluginForm.ui'
#
# Created by: PyQt5 UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_PluginSelect(object):
    def setupUi(self, PluginSelect):
        PluginSelect.setObjectName("PluginSelect")
        PluginSelect.resize(443, 115)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(PluginSelect)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.label = QtWidgets.QLabel(PluginSelect)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.comboBox = QtWidgets.QComboBox(PluginSelect)
        self.comboBox.setObjectName("comboBox")
        self.verticalLayout.addWidget(self.comboBox)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self.buttonBox = QtWidgets.QDialogButtonBox(PluginSelect)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)
        self.verticalLayout_2.addLayout(self.verticalLayout)

        self.retranslateUi(PluginSelect)
        QtCore.QMetaObject.connectSlotsByName(PluginSelect)

    def retranslateUi(self, PluginSelect):
        _translate = QtCore.QCoreApplication.translate
        PluginSelect.setWindowTitle(_translate("PluginSelect", "Select A Plugin"))
        self.label.setText(_translate("PluginSelect", "Select a plugin for creating an new configuration file:"))

