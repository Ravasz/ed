'''
Created on 7 Aug 2015

@author: mate

'''


def main():
  print("this is peptide parser")
  peptide_plotter("Ptpn22", "/home/mate/code/ed/src/data/r619w/txt/peptides.txt","/home/mate/code/ed/src/data/r619w/txt/phosphosites.txt", ["R-619W-1","-R619W-2","-R619W-3"])
  # ["-OST-1","-OST-2","-OST-3"]
  
  
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
    pepPosD[pepI].append(pepD[pepI])
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
    currBar = ax.bar(pepO[1],pepO[0],pepO[2], alpha = 0.5, color = next(cycolS))
    plt.draw()
  
  
  plt.show()
  


def prot_sequnce_drawer():
  """open peptides.txt and find peptides in there corresponding to the gene name specified in the queryS variable.
  Then find the protein refseq GI using db2db
  Then, download the full protein sequence from Entrez.
  Then, mark the found peptides in the protein sequence and output the resulting marked sequence in a html file as protein name.html
  Then, use the phosphorylation site dataset in the datafiles folder which is the complete phosphosite.org phosphosite map and mark phosphosites on the HTML"""
  
  from tools import html_creator, prot_id_converter, prot_entrez_fetch
  from copy import deepcopy
  print("this is peptide parser")
  
  queryS = "Ptpn22" # this is the search term that will be worked on. It should be a protein name like "Ptpn22"
  
  
  print("working on: " + queryS)
  with open("/home/mate/code/ed/src/data/r619w/txt/peptides.txt","r") as inpF: 
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
    pepStart = seqS.index(pepItem)
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
     
  print(uniId)
  phL = []
  with open("/home/mate/code/ed/src/data/r619w/txt/phosphosites.txt") as phInp: # now for the phosphosite data
    for phLine in phInp:
      phList = phLine.split("\t")
      try:
        if uniId == phList[0]: 
          phL.append(int(phList[1].split(";")[0]) - 1)
      except IndexError:
        continue
  print(phL)
  
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
  html_creator(queryS + " peptides", annotS, queryS + ".html")
  print("found peptides marked in the file: ", end=" ")
  print(queryS + ".html")  


  
def peptide_collector(targetFile, targetS, sampleList):  
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
    pepAvg = round(log2(pepSum/len(peptideD[pepItem])),0)
    pepD[pepItem] = pepAvg
    
  
  return pepD, uniprotId  

def peptide_finder(targetFile, targetS = "Ptpn22"):
  """Find all peptides which are matched to the targetS gene name in peptides.txt
  return peptide sequences as a list of strings and uniprot ID as string. The two together work as a tuple.
  """
  headerFlag = True
  peptideL = []
  uniprotId = ""
  for lineS in targetFile:
    if headerFlag:
      headerFlag = False
      continue
    lineL = lineS.split("\t")
    # print(lineL)
    geneNames = lineL[38].split(";")
    if geneNames[0] == "": continue
    if geneNames[0] == targetS:
      if lineL[0] not in peptideL:
        peptideL.append(lineL[0])
        uniprotId = lineL[34]
        if "-" in uniprotId:
          uniprotId = uniprotId[:uniprotId.index("-")]
  return peptideL, uniprotId
  
if __name__ == "__main__":
  main()