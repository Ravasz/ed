'''
Created on 26 May 2016

@author: mate
There is something wrong with hte way I process the ProteinGroups file. This script will test it.


'''

inpFile = open("/home/mate/workspace/katamari/src/root/ed/bob/processed/abs_quant_CTL-11-15-26-05-20162.txt","rU")

headerFlag = True

for inpLine in inpFile:
  if headerFlag:
    headerFlag = False
    continue
  

  inpL = inpLine.split(",")
  
  if inpL[1] == "A2BH40":
    print inpLine
  
  if int(inpL[5].split(";")[0]) < int(inpL[6].split(";")[0]):
    print inpLine

