'''
Created on 29 Jun 2016

@author: mate

Take in a list of hits from my OST mass spectrometry prepared by Thierry, and return Gene names for each uniprot ID in them.
'''

from ed.tools import prot_id_converter

with open("/home/mate/workspace/katamari/src/ed/datafiles/ot1-16062016.csv","r") as inpF:
  idL = []
  headerFlag = True
  for inpLine in inpF:
    if headerFlag:
      headerFlag = False
      continue
    firstTab = 0
    secTab = 0
    counterN = 0
    for inpC in inpLine:
      counterN += 1
      if inpC == "|" and firstTab == 0:
        firstTab = counterN - 1

        continue
      elif inpC == "|" and secTab == 0:
        secTab = counterN - 1

        break
    inpI = inpLine[firstTab + 1:secTab]
    idL.append(inpI)
    
  geneL = prot_id_converter(idL, outDB="genesymbol")
  
  with open("/home/mate/workspace/katamari/src/ed/datafiles/ot1-16062016-gene-names.csv","w") as outF:
    for geneI in geneL:
      outF.write(geneI + "\n")
      
print "done"
        