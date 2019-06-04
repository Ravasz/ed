'''
Created on 10 May 2019

@author: mate


Take in a list of uniprot identifiers and generate a file with the gene names and the GO terms.

requires my tools.py to work 
'''

def master():
  """call all the functions in order, write out results"""
  from tools import prot_id_converter
  
  fileDict = file_processor()
  
  fileList = list(fileDict.keys())
  
  fileL = []
  for fileI in fileList:
    if "-" in fileI: fileL.append(fileI.split("-")[0])
    else: fileL.append(fileI)

  
  nameL = prot_id_converter(fileL, "10090", "uniprotaccession", "genesymbol")
  
  for indexNum in range(len(fileList)):
    fileDict[fileList[indexNum]].insert(0,nameL[indexNum])
    
  
  with open("/home/mate/code/ed/src/data/processed/Ago2_mass_spec_2019_GO.csv","w") as outF:
    outF.write("Uniprot Accession,Gene symbol,FASTA header,GO Molecular function,GO Biological process, GO Cellular component\n")
    
    for uniS in fileList:
      if "-" in uniS: uniUp = uniS.split("-")[0]
      else: uniUp = uniS
      uniD = go_finder(uniID = uniUp)
      outF.write(uniS + ",")

      for outI in fileDict[uniS]:
        outF.write(outI + ",")
      
      if len(uniD["Function"]) > 0:
        for funI in uniD["Function"][:-1]: outF.write(funI + ";")
        outF.write(uniD["Function"][-1] + ",")
      else: outF.write("-,")
      
      if len(uniD["Process"]) > 0:     
        for funI in uniD["Process"][:-1]: outF.write(funI + ";")
        outF.write(uniD["Process"][-1] + ",")
      else: outF.write("-,")
      
      if len(uniD["Component"]) > 0:    
        for funI in uniD["Component"][:-1]: outF.write(funI + ";")
        outF.write(uniD["Component"][-1] + "\n")
      else: outF.write("-\n")




def file_processor():
  """open .csv file and extract list of uniprot IDs from Kat's file. return them as dict with uniprot ID as key, and a list with the gene description as single element"""
  
  with open("/home/mate/code/ed/src/data/Ago2_mass_spec_2019.csv","r") as inpF:
    next(inpF) # skip header
    
    inpDict = {}
    
    for inpLine in inpF:
      inpList = inpLine.split(",")
      inpDict[inpList[0]] = [inpList[1].rstrip()]
      
  return inpDict
      


def go_finder(uniID, goFolder = "/home/mate/workspace/katamari/src/ed/datafiles/go_terms/"):
  """take in the folder where GO terms are located and a uniprot ID, and return any associated GO terms in a dict where keys are strings of the type of GO term, and values are a list of terms"""
  
  from tools import go_term_advanced_lookup
  
  go_term_advanced_lookup(uniID) # download go term file for protein if not already present

  with open(goFolder+uniID+".txt","r") as goF:
    goFun = []
    goProc = []
    goComp = []
    
    
    try:
      next(goF) # skip header
    except StopIteration:
      print(uniID)
      raise
    
    for goLine in goF:
  
      goL = goLine.split("\t")
      currGo = goL[7].replace(",","-")
      if currGo == "molecular_function" or currGo == "biological_process" or currGo == "cellular_component":
        continue
      
      if goL[11] == "Function" or goL[11] == "molecular_function":
        if currGo not in goFun:
          goFun.append(currGo)
      elif goL[11] == "Process" or goL[11] == "biological_process":
        if currGo not in goProc:
          goProc.append(currGo)
      elif goL[11] == "Component" or goL[11] == "cellular_component":
        if currGo not in goComp:
          goComp.append(currGo)        
      else:
        print("got something odd:")
        print(goLine)
        print(goL[11])
        raise ValueError
  
#   funS = ""
#   for funI in goFun:
#     funS += funI
#     if funI is not goFun[-1]:
#       funS += ";"
#   
#   procS = ""
#   for procI in goProc:
#     procS += procI
#     if procI is not goProc[-1]:
#       procS += ";"
#       
#   compS = ""
#   for compI in goComp:
#     compS += compI
#     if compI is not goComp[-1]:
#       compS += ";"                  
      
  return {"Function":goFun, "Process":goProc, "Component": goComp}
  
master()     