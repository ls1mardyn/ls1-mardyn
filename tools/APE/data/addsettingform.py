# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'AddSettingForm.ui'
#
# Created by: PyQt5 UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_AddSettingForm(object):
    def setupUi(self, AddSettingForm):
        AddSettingForm.setObjectName("AddSettingForm")
        AddSettingForm.resize(468, 591)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(AddSettingForm.sizePolicy().hasHeightForWidth())
        AddSettingForm.setSizePolicy(sizePolicy)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(AddSettingForm)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(AddSettingForm)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.treeView = QtWidgets.QTreeView(AddSettingForm)
        self.treeView.setObjectName("treeView")
        self.verticalLayout.addWidget(self.treeView)
        self.buttonBox = QtWidgets.QDialogButtonBox(AddSettingForm)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)
        self.verticalLayout_2.addLayout(self.verticalLayout)

        self.retranslateUi(AddSettingForm)
        QtCore.QMetaObject.connectSlotsByName(AddSettingForm)

    def retranslateUi(self, AddSettingForm):
        _translate = QtCore.QCoreApplication.translate
        AddSettingForm.setWindowTitle(_translate("AddSettingForm", "Add Settings"))
        self.label.setText(_translate("AddSettingForm", "Please select an entry or a tree to add to the configuration:"))

