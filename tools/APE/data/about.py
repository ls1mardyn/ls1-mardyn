# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'About.ui'
#
# Created by: PyQt5 UI code generator 5.5.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_AboutWidget(object):
    def setupUi(self, AboutWidget):
        AboutWidget.setObjectName("AboutWidget")
        AboutWidget.resize(530, 188)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(AboutWidget.sizePolicy().hasHeightForWidth())
        AboutWidget.setSizePolicy(sizePolicy)
        AboutWidget.setMinimumSize(QtCore.QSize(530, 188))
        AboutWidget.setMaximumSize(QtCore.QSize(530, 188))
        self.label = QtWidgets.QLabel(AboutWidget)
        self.label.setGeometry(QtCore.QRect(100, 20, 391, 21))
        self.label.setStyleSheet("font-weight: bold")
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(AboutWidget)
        self.label_2.setGeometry(QtCore.QRect(100, 40, 391, 17))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(AboutWidget)
        self.label_3.setGeometry(QtCore.QRect(100, 110, 391, 17))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(AboutWidget)
        self.label_4.setGeometry(QtCore.QRect(100, 130, 391, 17))
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(AboutWidget)
        self.label_5.setGeometry(QtCore.QRect(100, 70, 391, 17))
        self.label_5.setObjectName("label_5")
        self.pushButton = QtWidgets.QPushButton(AboutWidget)
        self.pushButton.setGeometry(QtCore.QRect(430, 150, 85, 27))
        self.pushButton.setObjectName("pushButton")
        self.label_6 = QtWidgets.QLabel(AboutWidget)
        self.label_6.setGeometry(QtCore.QRect(20, 30, 64, 64))
        self.label_6.setText("")
        self.label_6.setObjectName("label_6")

        self.retranslateUi(AboutWidget)
        QtCore.QMetaObject.connectSlotsByName(AboutWidget)

    def retranslateUi(self, AboutWidget):
        _translate = QtCore.QCoreApplication.translate
        AboutWidget.setWindowTitle(_translate("AboutWidget", "Advanced Parameter Editor"))
        self.label.setText(_translate("AboutWidget", "Advanced Parameter Editor"))
        self.label_2.setText(_translate("AboutWidget", "Version 0.9"))
        self.label_3.setText(_translate("AboutWidget", "Lehrstuhl für Thermodynamik und Energietechnik"))
        self.label_4.setText(_translate("AboutWidget", "Universität Paderborn"))
        self.label_5.setText(_translate("AboutWidget", "APE is a modular and customizable visual config file editor."))
        self.pushButton.setText(_translate("AboutWidget", "OK"))

