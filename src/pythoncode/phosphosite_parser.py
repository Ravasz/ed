'''
Created on 30 Jun 2016

@author: mate
explore phosphorylation sites outputted by maxquant in Phospho (STY)sites.txt
'''

def phospho_gene_counter():
  """count how many proteins there are in the phospho dataset"""
  with open("/home/mate/workspace/katamari/src/ed/datafiles/Phospho-STY-Sites.txt","r") as inpF:
    with open("/home/mate/workspace/katamari/src/ed/datafiles/phospho-histones.txt","w") as outF:
      geneL = []
      startFlag = True
      for inpLine in inpF:
        if startFlag:
          startFlag = False
          print(inpLine)
          continue
        inpL = inpLine.split("\t")
        curL = []
        if "Hist" in inpL[5]:
          print(inpL)
          for histI in inpL:
            outF.write(histI.split(";")[0])
            if histI != inpL[-1]:
              outF.write(",")
  
        for inpGene in inpL[5].split(";"):
          if inpGene == "":
            curL.append(inpL[2].split(";")[0])
            break
          else: curL.append(inpGene)
        curG = min(curL, key=len)
        if curG not in geneL and "REV__" not in curG:
          geneL.append(curG)
        
      print(geneL)
      print(len(geneL))
    
def phospho_to_regular():
  """compare the STY dataset proteins to the original 24H T4 calulcated dataset from Dundee. THey should mostly overlap"""
  with open("/home/mate/workspace/katamari/src/ed/datafiles/proteinGroups_phospho_24H_T4.txt","r") as inpPhos:
    with open("/home/mate/workspace/katamari/src/ed/datafiles/24h proteinGroups.txt","rU") as inpOld:
      oldL = []
      for inpLine in inpOld:
        inpL = inpLine.split("\t")
        inpIL = inpL[1].split(";")
        inpI = min(inpIL, key=len)
        if "REV__" in inpI or "CON__" in inpI: continue
        oldL.append(inpI)
      print("old:")
      print(len(oldL))
      
      phosL = []
      for inpLine in inpPhos:
        inpL = inpLine.split("\t")
        inpIL = inpL[1].split(";")
        inpI = min(inpIL, key=len)
        if "REV__" in inpI or "CON__" in inpI: continue
        phosL.append(inpI)
      print("phos:")
      print(len(phosL))      
      
      diffL = []  
      for oI in oldL:
        if oI not in phosL:
          diffL.append(oI)
      
      print(len(diffL))
      print(diffL)
      
      diff2L = []
      for oI in phosL:
        if oI not in oldL:
          diff2L.append(oI)
      
      print(len(diff2L))
      
    
    


phospho_gene_counter()