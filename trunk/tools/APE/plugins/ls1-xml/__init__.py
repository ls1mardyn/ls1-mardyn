'''
Created on 27.04.2017

@author: kokr
'''

class Plugin(object):
	
	#check if the selected file matches this plugin
	def checkMagic(self, file):
		returnVal = False
		line = file.readLine()
		if line[:14] == "<?xml version=":
			line = file.readLine()
			returnVal = line[:16] == "<mardyn version="
				
		return returnVal
	
	#get content xml filename
	def getContentXML(self):
		return "content.xml"