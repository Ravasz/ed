'''
Created on 31 Jan 2019

@author: mate

So HGNC IDs are only for human genes, not for mouse. So this script that aims to map mouse protein identifiers to HGNC IDs promptly returns empty lists. 
This script is discontinued, starting a new one instead that will map mouse IDs to human ones first.
'''

from tools import prot_id_converter

protName = "/home/mate/code/ed/src/data/bob/24H_T4_original_proteinGoups_avg.csv"

outName = protName.split(".")[0] + "_hgnc.csv"
fileList = []
convL = []

outF = open(outName,"w")

with open(protName,"r") as inpF:
  headerLine = next(inpF)
  
  lineCount = 1

  for inpLine in inpF:
    if lineCount % 100 == 0:
      resL = prot_id_converter(protList = convL, outDB = "hgncid")
      print(resL)
      convL = []
      outCount = 0
      for fileItem in fileList:
        outF.write(fileItem[0] + ",")
        multiFlag = False
        for nameItem in fileItem[0].split(";"):
          if multiFlag: outF.write(";")
          for itemD in resL:
            if itemD["InputValue"] == nameItem:
              outF.write(itemD["HGNC ID"])
              multiFlag = True
              break

        outF.write("," + fileItem[1] + "," + fileItem[2] + "\n")
        outCount += 1
        
      fileList = []

    inpList = inpLine.rstrip().split(",")
    fileList.append(inpList)
    if ";" in inpList[0]:
      for geneI in inpList[0].split(";"):
        convL.append(geneI)
    else: convL.append(inpList[0])

      
    lineCount += 1
    
outF.close()
    
