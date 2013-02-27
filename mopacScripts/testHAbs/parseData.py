import os

import xlwt
import glob
from datetime import datetime

os.chdir("2013.02.26 - training/")
wbk = xlwt.Workbook()
sheet = wbk.add_sheet('ts auto-generator')
sheet.write(0,0,'Rxn Number')
sheet.write(0,1,'Reactant 1')
sheet.write(0,2,'Reactant 2')
sheet.write(0,3,'Product 1')
sheet.write(0,4,'Product 2')
sheet.write(0,5,'Activation Energy')
sheet.write(0,6,'TS vib')           
sheet.write(0,7,'*1')               
sheet.write(0,8,'*2')               
sheet.write(0,9,'*3')               
sheet.write(0,10,'*1 to *2')         
sheet.write(0,11,'*2 to *3')         
sheet.write(0,12,'*1 to *3')
sheet.write(0,13,'Notes')

i = 1
for files in glob.glob("activationER*"):
	# get the number of the reaction in this list
	listNum = files.split('.')[0].split('R')[-1]
	sheet.write(i,0,listNum)
	loadFile = file(files)
	loadFile = loadFile.readlines()
	loadFile.pop(0)
	infoLines = []
	for line in loadFile:
		infoLines.append(line.split()[-1])
	for j in range(len(infoLines)):
		sheet.write(i, j+1,infoLines[j])
	i+=1

wbk.save('reactionData.xls')