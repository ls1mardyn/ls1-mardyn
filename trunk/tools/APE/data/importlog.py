# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ImportLog.ui'
#
# Created by: PyQt5 UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_ImportLog(object):
    def setupUi(self, ImportLog):
        ImportLog.setObjectName("ImportLog")
        ImportLog.resize(806, 489)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(ImportLog)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(ImportLog)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.plainTextEdit = QtWidgets.QPlainTextEdit(ImportLog)
        self.plainTextEdit.setReadOnly(True)
        self.plainTextEdit.setObjectName("plainTextEdit")
        self.verticalLayout.addWidget(self.plainTextEdit)
        self.buttonBox = QtWidgets.QDialogButtonBox(ImportLog)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)
        self.verticalLayout_2.addLayout(self.verticalLayout)

        self.retranslateUi(ImportLog)
        QtCore.QMetaObject.connectSlotsByName(ImportLog)

    def retranslateUi(self, ImportLog):
        _translate = QtCore.QCoreApplication.translate
        ImportLog.setWindowTitle(_translate("ImportLog", "File Import Log"))
        self.label.setText(_translate("ImportLog", "Explain"))

