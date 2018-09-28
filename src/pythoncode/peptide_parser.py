'''
Created on 7 Aug 2015

@author: mate

'''


def main():
  print("this is peptide parser")
  # prot_sequence_drawer("/home/mate/code/ed/src/data/eif5a/eif5-29-08-2018-peptides.txt", "/home/mate/code/ed/src/data/eif5a/eif5-28-08-2018-Phospho (STY)Sites.txt")
  # peptide_diff_plotter("Ptpn22", "/home/mate/code/ed/src/data/r619w/txt/peptides.txt","/home/mate/code/ed/src/data/r619w/txt/phosphosites.txt", ["-OST-1","-OST-2","-OST-3"], ["R-619W-1","-R619W-2","-R619W-3"])
  peptide_plotter("Eif5a", "/home/mate/code/ed/src/data/eif5a/eif5-29-08-2018-peptides.txt","/home/mate/code/ed/src/data/eif5a/eif5-28-08-2018-Phospho (STY)Sites.txt", ["E1_CT_LY","E1_CT_TR","E2_CT_LY","E2_CT_TR","E3_CT_LY","E3_CT_TR"])
  # peptide_plotter("Ptpn22", "/home/mate/code/ed/src/data/r619w/txt/peptides.txt","/home/mate/code/ed/src/data/r619w/txt/phosphosites.txt", ["R-619W-1","-R619W-2","-R619W-3"])
  # ["-OST-1","-OST-2","-OST-3"]
  

def peptide_diff_plotter(protName, peptidesFile, phosphoSitesFile, sample1List, sample2List):
  """take a protein name, download its sequence, and create a bar plot showing the difference in peptide intensities between two groups
  """
  
  from tools import prot_id_converter, prot_entrez_fetch
  # from copy import deepcopy
  import os.path
  import sys
  from collections import defaultdict
  import matplotlib.pyplot as plt
  from itertools import cycle
  
  print("working on: " + protName)
  
  print("\nprocessing files: ")
  if os.path.isfile(peptidesFile): 
    print(peptidesFile)
  else: 
    print("file %s not found" % (peptidesFile,))
    sys.exit(0)
    
  if os.path.isfile(phosphoSitesFile): 
    print(phosphoSitesFile)
  else: 
    print("file %s not found" % (phosphoSitesFile,))
    sys.exit(0)   

  
  with open(peptidesFile,"r") as inpF: 
    pep1D, uniId = peptide_collector(targetFile=inpF, targetS = protName, sampleList = sample1List, roundN = 1) # find peptides from peptides.txt for the protein name @UnusedVariable
  print(pep1D)
  print("peptides found for first hit")
  
  with open(peptidesFile,"r") as inpF:
    pep2D, uniId = peptide_collector(targetFile=inpF, targetS = protName, sampleList = sample2List, roundN = 1)  # @UnusedVariable

  print(pep2D)
  print("peptides found for second hit")  
  targetL = [protName]
  idList = prot_id_converter(targetL, "10090", inpDB = "genesymbol",outDB="refseqproteingi")
  seqL = prot_entrez_fetch(idList, retM="gb", retT="fasta")
  for seqItem in seqL:
    seqS = seqItem.split("\n")[1]
    print(seqS)
  # seqS = "MDQREILQQLLKEAQKKKLNSEEFASEFLKLKRQSTKYKADKIYPTTVAQRPKNIKKNRYKDILPYDHSLVELSLLTSDEDSSYINASFIKGVYGPKAYIATQGPLSTTLLDFWRMIWEYRILVIVMACMEFEMGKKKCERYWAEPGETQLQFGPFSISCEAEKKKSDYKIRTLKAKFNNETRIIYQFHYKNWPDHDVPSSIDPILQLIWDMRCYQEDDCVPICIHCSAGCGRTGVICAVDYTWMLLKDGIIPKNFSVFNLIQEMRTQRPSLVQTQEQYELVYSAVLELFKRHMDVISDNHLGREIQAQCSIPEQSLTVEADSCPLDLPKNAMRDVKTTNQHSKQGAEAESTGGSSLGLRTSTMNAEEELVLHSAKSSPSFNCLELNCGCNNKAVITRNGQARASPVVGEPLQKYQSLDFGSMLFGSCPSALPINTADRYHNSKGPVKRTKSTPFELIQQRKTNDLAVGDGFSCLESQLHEHYSLRELQVQRVAHVSSEELNYSLPGACDASCVPRHSPGALRVHLYTSLAEDPYFSSSPPNSADSKMSFDLPEKQDGATSPGALLPASSTTSFFYSNPHDSLVMNTLTSFSPPLNQETAVEAPSRRTDDEIPPPLPERTPESFIVVEEAGEPSPRVTESLPLVVTFGASPECSGTSEMKSHDSVGFTPSKNVKLRSPKSDRHQDGSPPPPLPERTLESFFLADEDCIQAQAVQTSSTSYPETTENSTSSKQTLRTPGKSFTRSKSLKIFRNMKKSVCNSSSPSKPTERVQPKNSSSFLNFGFGNRFSKPKGPRNPPSAWNM"
  print("protein sequence found")
  # annotS = deepcopy(seqS)  
  
  pepPos1D = defaultdict(list)
  pepPos2D = defaultdict(list)
  
  for pepI in pep1D: # make new dict with each peptide log2 intensity, start position in protein sequence, and peptide length
    pepPos1D[pepI].append(pep1D[pepI])
    pepPos1D[pepI].append(seqS.index(pepI))
    pepPos1D[pepI].append(len(pepI))
  
  for key, value in pepPos1D.items():
    print(key,value)
  print(len(pepPos1D))
  
  for pepI in pep2D:
    pepPos2D[pepI].append(pep2D[pepI])
    pepPos2D[pepI].append(seqS.index(pepI))
    pepPos2D[pepI].append(len(pepI))
  
  for key, value in pepPos2D.items():
    print(key,value)
  print(len(pepPos2D))  
  
  uniD = defaultdict(list)
  
  for pepN, pepO in pepPos1D.items():
    if pepN in pepPos2D:
      if pepO[1] == pepPos2D[pepN][1]:
        uniD[pepN] = [float(pepO[0]) - float(pepPos2D[pepN][0]),pepO[1],pepO[2]]
        if uniD[pepN][0] == 0:  uniD[pepN][0] = 0.1
      else:
        uniD[pepN] = [int(pepO[0]),pepO[1],pepO[2]]
    else:
      uniD[pepN] = [int(pepO[0]),pepO[1],pepO[2]]
     
  for pepN, pepO in pepPos2D.items():
    if pepN in pepPos1D:
      if pepO[1] == pepPos1D[pepN][1]:
        pass
      else:
        uniD[pepN] = [int(pepO[0]),pepO[1],pepO[2]]
    else:
      uniD[pepN] = [int(pepO[0]),pepO[1],pepO[2]]     
  
  cycolS = cycle('bgrcmk')
  
  fig = plt.figure()
  ax = fig.add_subplot(111) 
  
  # rects1 = ax.bar(ind, yvals, width, color='r')
  
  for pepN,pepO in uniD.items():
    currBar = ax.bar(pepO[1],pepO[0],pepO[2], alpha = 0.5, color = next(cycolS))
    plt.draw()
  
  
  plt.show()
  

  
def peptide_plotter(protName, peptidesFile, phosphoSitesFile, sampleList):
  """take a protein name, download its sequence, and create a bar plot with each peptide from it highlighted and quantified above the sequence. 
  Mark phosphorylation sites in the sequence too"""
  
  from tools import html_creator, prot_id_converter, prot_entrez_fetch
  # from copy import deepcopy
  import os.path
  import sys
  from collections import defaultdict
  import matplotlib.pyplot as plt
  from itertools import cycle
  
  print("working on: " + protName)
  
  print("\nprocessing files: ")
  if os.path.isfile(peptidesFile): 
    print(peptidesFile)
  else: 
    print("file %s not found" % (peptidesFile,))
    sys.exit(0)
    
  if os.path.isfile(phosphoSitesFile): 
    print(phosphoSitesFile)
  else: 
    print("file %s not found" % (phosphoSitesFile,))
    sys.exit(0)   

  
  with open(peptidesFile,"r") as inpF: 
    pepD, uniId = peptide_collector(targetFile=inpF, targetS = protName, sampleList = sampleList) # find peptides from peptides.txt for the protein name
  print(pepD)
  print("peptides found")
  targetL = [protName]
  idList = prot_id_converter(targetL, "10090", inpDB = "genesymbol",outDB="refseqproteingi")
  seqL = prot_entrez_fetch(idList, retM="gb", retT="fasta")
  for seqItem in seqL:
    seqS = seqItem.split("\n")[1]
    print(seqS)
  # seqS = "MDQREILQQLLKEAQKKKLNSEEFASEFLKLKRQSTKYKADKIYPTTVAQRPKNIKKNRYKDILPYDHSLVELSLLTSDEDSSYINASFIKGVYGPKAYIATQGPLSTTLLDFWRMIWEYRILVIVMACMEFEMGKKKCERYWAEPGETQLQFGPFSISCEAEKKKSDYKIRTLKAKFNNETRIIYQFHYKNWPDHDVPSSIDPILQLIWDMRCYQEDDCVPICIHCSAGCGRTGVICAVDYTWMLLKDGIIPKNFSVFNLIQEMRTQRPSLVQTQEQYELVYSAVLELFKRHMDVISDNHLGREIQAQCSIPEQSLTVEADSCPLDLPKNAMRDVKTTNQHSKQGAEAESTGGSSLGLRTSTMNAEEELVLHSAKSSPSFNCLELNCGCNNKAVITRNGQARASPVVGEPLQKYQSLDFGSMLFGSCPSALPINTADRYHNSKGPVKRTKSTPFELIQQRKTNDLAVGDGFSCLESQLHEHYSLRELQVQRVAHVSSEELNYSLPGACDASCVPRHSPGALRVHLYTSLAEDPYFSSSPPNSADSKMSFDLPEKQDGATSPGALLPASSTTSFFYSNPHDSLVMNTLTSFSPPLNQETAVEAPSRRTDDEIPPPLPERTPESFIVVEEAGEPSPRVTESLPLVVTFGASPECSGTSEMKSHDSVGFTPSKNVKLRSPKSDRHQDGSPPPPLPERTLESFFLADEDCIQAQAVQTSSTSYPETTENSTSSKQTLRTPGKSFTRSKSLKIFRNMKKSVCNSSSPSKPTERVQPKNSSSFLNFGFGNRFSKPKGPRNPPSAWNM"
  print("protein sequence found")
  # annotS = deepcopy(seqS)  
  
  pepPosD = defaultdict(list)
  
  for pepI in pepD: # make new dict with each peptide log2 intensity, start position in protein sequence, and peptide length
    
    if pepD[pepI] > 1: pepPosD[pepI].append(pepD[pepI])
    else: pepPosD[pepI].append(1.0) # if peptide intensity is 0, show as an intensity of 1 to make it visible
    pepPosD[pepI].append(seqS.index(pepI))
    pepPosD[pepI].append(len(pepI))
  
  for key, value in pepPosD.items():
    print(key,value)
  print(len(pepPosD))
  cycolS = cycle('bgrcmk')
  
  fig = plt.figure()
  ax = fig.add_subplot(111) 
  
  # rects1 = ax.bar(ind, yvals, width, color='r')
  
  for pepN,pepO in pepPosD.items():
    print(pepO[1],pepO[0],pepO[2])
    currBar = ax.bar(pepO[1],pepO[0],pepO[2], align = "edge", alpha = 0.5, color = next(cycolS))
    plt.draw()
  
  
  plt.show()
  


def prot_sequence_drawer(inpFileName = "/home/mate/code/ed/src/data/r619w/txt/peptides.txt", phosphoFileName = "/home/mate/code/ed/src/data/r619w/txt/phosphosites.txt"):
  """open peptides.txt and find peptides in there corresponding to the gene name specified in the queryS variable.
  Then find the protein refseq GI using db2db
  Then, download the full protein sequence from Entrez.
  Then, mark the found peptides in the protein sequence and output the resulting marked sequence in a html file as protein name.html
  Then, use the phosphorylation site dataset in the datafiles folder which is the complete phosphosite.org phosphosite map and mark phosphosites on the HTML"""
  
  from tools import html_creator, prot_id_converter, prot_entrez_fetch
  from copy import deepcopy
  import re, os
  print("this is peptide parser")
  
  queryS = "Eif5a" # this is the search term that will be worked on. It should be a protein name like "Ptpn22"
  trashList = []
  
  
  print("working on: " + queryS)
  with open(inpFileName,"r") as inpF: 
    pepL, uniId = peptide_finder(targetFile=inpF, targetS = queryS) # find peptides from peptides.txt for the protein name
  print(pepL)
  print("peptides found")
  targetL = [queryS]
  idList = prot_id_converter(targetL, "10090", inpDB = "genesymbol",outDB="refseqproteingi")
  seqL = prot_entrez_fetch(idList, retM="gb", retT="fasta")
  for seqItem in seqL:
    seqS = seqItem.split("\n")[1]
    print(seqS)
  print("protein sequence found")
  annotS = deepcopy(seqS)
  pStartL = []
  pEndL = []
  for pepItem in pepL: # locate peptides in full protein sequence and store positions for starts and ends. merge overlapping peptides.
    try:
      pepStart = seqS.index(pepItem)
    except ValueError:
      trashList.append(pepItem)
      continue
    pepEnd = pepStart + len(pepItem)
    
    startCount = 0 # handle starts
    for startItem in pStartL:
      if startItem <= pepStart:
        startCount += 1
    endCount = 0
    for endItem in pEndL:
      if endItem <= pepStart:
        endCount += 1
    if startCount == endCount and pepStart not in pEndL: # start new peptide
      pStartL.append(pepStart)
    elif startCount == endCount and pepStart in pEndL: # start a new peptide at the end of another peptide
      pEndL.remove(pepStart)
    
    overlapCount = 0
    for startItem in pStartL[:]: # handle ends
      if pepStart<startItem<=pepEnd:
        pStartL.remove(startItem) # remove extra starts
        overlapCount += 1
        
    for endItem in pEndL[:]:
      if pepStart<=endItem<=pepEnd: # remove extra ends
        pEndL.remove(endItem)
        overlapCount -= 1
    
    if pepEnd not in pEndL and overlapCount <= 0: # add end
      curStart = 500000
      for pSI in pStartL:
        if curStart > pSI > pepEnd:
          curStart = pSI
      curEnd = 500000
      for pEI in pEndL:
        if curEnd > pEI > pepEnd:
          curEnd = pEI
      if curStart <= curEnd:  # check if next tag is start or end. if start, add end. if end, do nothing
        pEndL.append(pepEnd)
     
  # print(uniId)
  phL = []
  with open(phosphoFileName, "r") as phInp: # now for the phosphosite data
    for phLine in phInp:
      phList = phLine.split("\t")
      # print(phList)
      
      phNameL = phList[0].split(";")
      
      if uniId in phNameL: 
        pepString = phList[67]
        while True:
          if pepString.find("(") == (-1): break
          
          modProb = float(pepString[pepString.index("(")+1:pepString.index(")")])
          # print(modProb)
          
          if modProb > 0.9:
            if seqS.find(re.sub("[^A-Z]+", "", pepString)) == -1: break # regex to select only text bits and leave out actual locations, check if the peptide is even in the target protein
            
            newSite = pepString.index("(") - 1 + seqS.index(re.sub("[^A-Z]+", "", pepString))
            # print(pepString[pepString.index("(") - 1])
            if newSite not in phL: phL.append(newSite)
            # print(seqS[phL[-1]])
            # print(pepString)
            
          pepString = pepString[pepString.index(")") + 1:]

        # phL.append(int(phList[1].split(";")[0]) - 1)
        
  # print(phL)
  # print(pStartL,pEndL)
  
  fullL = pStartL + pEndL
  for phItem in phL:
    if phItem not in fullL: fullL.append(phItem)
  fullL.sort()
  
  offsetN = 0
  for posI in fullL: # from the resulting intervals, create emphasis in html file, mark phosphosites in red
    if posI in pStartL:
      annotS = annotS[:posI+offsetN] + "<mark>" + annotS[posI+offsetN:]
      offsetN += 6
    elif posI in pEndL:
      annotS = annotS[:posI+offsetN] + "</mark>" + annotS[posI+offsetN:]
      offsetN += 7
    if posI in phL:
      annotS = annotS[:posI+offsetN] + r"""<strong style="color: red;">""" + annotS[posI+offsetN] + r"""</strong>""" + annotS[posI+offsetN + 1 :]
      offsetN += 37
        
  print(annotS)
  html_creator(queryS + " peptides", annotS, os.path.dirname(phosphoFileName) + "/" + queryS + ".html")
  print("found peptides marked in the file: ", end=" ")
  print(queryS + ".html")  
  print("\ntrash that did not match protein sequence:")
  for trashItem in trashList:
    print(trashItem)


  
def peptide_collector(targetFile, targetS, sampleList, roundN = 0):  
  """Find all peptides which are matched to the targetS gene name in peptides.txt
  return peptide sequences as a dict of strings and uniprot ID as string. The two together work as a tuple."""
  
  from collections import defaultdict
  import sys
  from math import log2

  headerFlag = True
  peptideD = defaultdict(list)
  uniprotId = ""
  namePos = -1
  seqPos = -1
  idPos = -1
  intensityPosL = []
  
  for lineS in targetFile:
    if headerFlag:
      headerFlag = False
      lineL = lineS.split("\t")
      posCount = 0
      for headerI in lineL:
        if headerI == "Gene names": namePos = posCount
        if headerI == "Sequence": seqPos = posCount
        if headerI == "Leading razor protein": idPos = posCount
        for sI in sampleList:
          if headerI == "Intensity " + sI: intensityPosL.append(posCount)
            
        posCount += 1
      if namePos == -1 or seqPos == -1 or idPos == -1 or len(intensityPosL) == 0:
        print("stuff is missing from the entry file")
        print(namePos,seqPos, idPos, intensityPosL)
        sys.exit(0)
      continue
    lineL = lineS.split("\t")
    # print(lineL)
    geneNames = lineL[namePos].split(";")
    if geneNames[0] == "": continue
    if geneNames[0] == targetS:
      if lineL[seqPos] not in peptideD:
        for intensityPos in intensityPosL:
          peptideD[lineL[seqPos]].append(lineL[intensityPos])
        uniprotId = lineL[idPos]
        if "-" in uniprotId:
          uniprotId = uniprotId[:uniprotId.index("-")]
    
  pepD = {} # calculate averages of intensity scores if multiple samples are given as input
  for pepItem in peptideD:
    pepSum = 0
    for pepNum in peptideD[pepItem]:
      pepSum += float(pepNum)
      if pepSum == 0: pepSum = len(peptideD[pepItem])
    pepAvg = round(log2(pepSum/len(peptideD[pepItem])),roundN)
    pepD[pepItem] = pepAvg
    
  
  return pepD, uniprotId  

def peptide_finder(targetFile, targetS = "Ptpn22"):
  """Find all peptides which are matched to the targetS gene name in peptides.txt
  return peptide sequences as a list of strings and uniprot ID as string. The two together work as a tuple.
  """
  headerFlag = True
  peptideL = []
  uniprotId = ""
  uniDict = {}
  for lineS in targetFile:
    if headerFlag:
      headerFlag = False
      continue
    lineL = lineS.split("\t")
    # print(lineL)
    geneNames = lineL[38].split(";")
    if geneNames[0] == "": continue
    # if geneNames[0] == targetS:
    if targetS in geneNames:
      if lineL[0] not in peptideL:
        peptideL.append(lineL[0])
        uniprotId = lineL[34]
        if ";" in uniprotId: 
          for uniI in uniprotId.split(";"):
            if "-" in uniI: 
              uniJ = uniI[:uniI.index("-")]
              if uniJ in uniDict: uniDict[uniJ] += 1
              else: uniDict[uniJ] = 0
            else:
              if uniI in uniDict: uniDict[uniI] += 1
              else: uniDict[uniI] = 0
            
        else: 
          if "-" in uniprotId: uniprotId = uniprotId[:uniprotId.index("-")]
          if uniI in uniDict: uniDict[uniI] += 1
          else: uniDict[uniI] = 0
        
        
  maxVal = 0
  maxStr = ""
  for keyS, valueN in uniDict.items():
    if valueN > maxVal:
      maxStr = keyS
      maxVal = valueN
      
  uniprotId = maxStr    
  print("\nbest match for protein ID:")
  print(uniprotId)      
        
  return peptideL, uniprotId
  
if __name__ == "__main__":
  main()