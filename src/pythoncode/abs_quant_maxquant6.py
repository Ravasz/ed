'''
Created on 5 Aug 2019

@author: mate


calculate absolute numbers of proteins using the proteinGroups.txt file. 
This module is written specifically for Bob's 24H T4 stimulation dataset, 
but should work with anything else in the same format.


based on abs_quant

'''

def master(fileName = "", outFolder = "", tagS="julia_t0"):
  """run all the other scripts in this module in the appropriate order"""
  
  # previous: /home/mate/workspace/katamari/src/root/ed/datafiles/24H_T4_recalc/MassSpec/txt_RS_24hr unique razor/proteinGroups.txt
  # previous: /home/mate/workspace/katamari/src/root/ed/datafiles/CTL_Bob_recalc/txt_CTL_unique razor/proteinGroups.txt
  
  import time 
  import os.path
  
  if fileName == "":
    fileName = os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "julia_time_course", "proteinGroups_Julia_t0.txt")
    
  if outFolder == "": 
    outFolder = os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "processed", "julia_time_course/")
  
  curDate = time.strftime("%H-%M-%d-%m-%Y")
  
  expTag = tagS + "-" + curDate
  
  file_parser(fileName, outFolder, expTag)
  entry_parser(outFolder + expTag + "_1.txt", outFolder, expTag)
  

def file_parser(fileName, outFolder, expTag):
  """From a proteinGroups.txt file, take out the following things:
  - numeric ID
  - protein uniprot ID
  - protein symbol
  - fasta header
  - molecular weight
  - peptides
  - unique + razor peptides
  - intensity 1
  - intensity 2
  - intensity 3
  - LFQ 1
  - LFQ 2
  - LFQ 3
         
  Do not select contaminants or reverse peptides""" 

  print("this is file parser in abs_quant")
  with open(fileName,"r") as inpF:
    cN = -1
    outF = open(outFolder + expTag + "_1.txt","w")
    outF.write("Numeric ID,")
    
    headerL = next(inpF).strip("\r\n").split("\t") # process header, and collect required columns
    colL = []
    conList = []
    colCount = 0
    foundFlag = False
    for headerI in headerL:
      # print(headerI)
      if headerI == "Potential contaminant":
        conList.append(colCount)
      elif headerI == "Reverse":
        conList.append(colCount)
      
      elif "Fasta header" in headerI:
        if foundFlag: outF.write(",")
        foundFlag = True
        colL.append(colCount)
        outF.write(headerI)      
      elif "Mol. weight" in headerI:
        if foundFlag: outF.write(",")
        foundFlag = True
        colL.append(colCount)
        outF.write(headerI)            
      elif "Intensity " in headerI:
        if foundFlag: outF.write(",")
        foundFlag = True
        colL.append(colCount)
        outF.write(headerI)
      elif "LFQ intensity" in headerI:
        if foundFlag: outF.write(",")
        foundFlag = True
        colL.append(colCount)
        outF.write(headerI)
      elif "Majority protein ID" in headerI:
        if foundFlag: outF.write(",")
        foundFlag = True
        colL.append(colCount)
        outF.write("Majority protein IDs, Protein names, Gene names")
      elif "Peptides" == headerI:
        if foundFlag: outF.write(",")
        foundFlag = True
        colL.append(colCount)
        outF.write(headerI)
      elif "Razor + unique peptides" == headerI:
        if foundFlag: outF.write(",")
        foundFlag = True
        colL.append(colCount)
        outF.write(headerI)
          
      colCount += 1
    if foundFlag: outF.write(",Protein Count 1, Protein count 2, Protein count 3, Mean protein count\n")
      
    # print(colL)
    for inpLine in inpF:
      cN += 1
      contFlag = False
      # if cN == 100: break
      inpLine = inpLine.rstrip("\r\n").replace(",","-")
      inpItems = inpLine.split("\t")
      # print(inpItems)
      for conI in conList:
        if inpItems[conI] == "+": 
          # print cN,
          contFlag = True # get rid of contaminants and reverse proteins
        
      if contFlag: continue
      
      
      outF.write(str(cN) + ",")
      for colNum in colL[:-1]:
        # print(inpItems[colNum])
        # if ";" in inpItems[colNum] and len(inpItems[colNum]) < 10: pass
          #print(inpItems[colNum])
        outF.write(inpItems[colNum] + ",")
      outF.write(inpItems[colL[-1]] + "\n")
   
    outF.close()
  print("\n", cN, "lines parsed successfully")      
  
def entry_parser(fileName, outFolder, expTag):
  """remove duplicate protein name and total peptide count cell entries"""
  from operator import add
  from tools import go_term_advanced_lookup
  
  goC = 0
  errC = 0
  
  print("this is entry parser in abs_quant")
  
  with open(fileName,"r") as inpF:
    outF = open(outFolder + expTag + "_2.txt","w")
    cN = 0
    outDict = {}
    
    for inpLine in inpF:
      cN += 1
      if cN == 1: 
        outF.write(inpLine) #.rstrip("\n"))
        print("header written")
        print(inpLine)
        continue # skip header
      inpLine = inpLine.rstrip()
      inpItem = inpLine.split(",")
      # print(inpItem)
      
      geneL = inpItem[1].split(";") # process uniprot ID here
      # print(geneL)
      geneNameL = []
      protNameL = []
      for geneI in geneL:
        geneNameL.append(geneI.split("|")[1])
        protNameL.append(geneI.split("|")[2].split("_")[0])
      # print(geneNameL)
      lenS = len(geneNameL[0])
      curGene = geneNameL[0]
      for geneI in geneNameL:
        if len(geneI) < lenS:
          lenS = len(geneI)
          curGene = geneI
      
      curProt = ""

      for protI in protNameL:
        if len(protI) == 6 and not protI[1].isalpha() or protI.startswith("A0A"): pass
        else: 
          curProt = protI[0] + protI[1:].lower()
          break
      
      if curProt == "": curProt = protNameL[0]
        

    
      if curGene[-2] == "-":
        curGene = curGene[:-2]
      if curGene[-3] == "-":
        curGene = curGene[:-3]
        
      protFullName = inpItem[2].split("OS=")[0].split("MOUSE ")[1]
      
      
      
      outL = [int(inpItem[0]), curGene, protFullName, curProt, inpItem[2]]
      outSum = 0
      for outI in inpItem[3:5]:
        outL.append(int(outI))
        outSum += int(outI)
      outL.append(float(inpItem[5]))
      for outI in inpItem[6:]:
        outL.append(int(outI))
        outSum += int(outI)      
      if outSum == 0: continue
      # print(outL)
          
        
      if curGene in outDict: # handle duplicate protein entries and merge them together
        print("%s is duplicate") % curGene
        # mergeL = outL[:4] # names and such will be from the newly found duplicate
        mergeL = outL[:2]
        mergeL.append(outL[2] + ";" + outDict[curGene][2])
        mergeL.append(outL[3])
        mergeL.extend([max(outL[4],outDict[curGene][4]), max(outL[5],outDict[curGene][5]), max(outL[6],outDict[curGene][6])]) # add the higher of the molecular weight and the total peptide count, and unique peptide count
        mergeL.extend(list(map(add, outL[7:], outDict[curGene][7:]))) # add intensities
        outL = mergeL
  
      outDict[curGene] = outL # assemble all the stuff in this dict
    
    finD = absolute_quant(outDict)
    
    outN = 0  
    print("looking up GO terms. This might take several minutes")
    for outDV in sorted(list(finD.items()), key=lambda x:x[1][0]): # sort the dict based on the ID they have
      outN += 1
      # print(outDV)
      
      # this bit looks up and downloads GO terms. Takes about 4s per proteins, so beware
      try:
        go_term_advanced_lookup(outDV[1][1])
        goC += 1
        if goC%100 == 0:
          print(".", end=" ")
        # print goC
      except:
        errC += 1
        print(outDV[1][1], " not found")
      
      # if outN == 100: break
      for outI in outDV[1][:-1]:
        outF.write(str(outI) + ",")
      outF.write(str(outDV[1][-1]) + "\n")
  
    print("\nunique proteins: ", outN)
    print("lines parsed: ", cN)
    print("GO terms found:", goC)
    print("GO_terms not found: ", errC)
    outF.close()

def absolute_quant(completeDict):
  """quantify all proteins using the proteomic ruler approach, detailed here: http://www.ncbi.nlm.nih.gov/pubmed/25225357
  to be called by entry_parser
  
  [Nprot] = ([protein intensity] * Navogadro * Mdna) / ([MWprot] * 1000 * sum[Histone intensity]) 
  
  """
  
  mDNA = 5.5225 * 10**-12 # weight of DNA in diploid murine cells, in grams
  nA = 6.022 * 10**23 # avogadro's constant
  histIntL = [0,0,0]
  totIntL = [0,0,0]
  
  histoneIDs = ["O35216", "P10922", "P43275", "P43277", "P15864", "P43274", "P43276", "Q8CJI4", "Q8VIK3", "Q8VIK3-2", 
                "Q07133", "Q8CGP7", "Q8CGP5", "Q8CGP6", "P22752", "Q6GSS7", "Q64522", "Q64523", "Q8BFU2", "Q9CQ70", 
                "Q8R1M2", "Q8CCK0", "Q3THW5", "P27661", "Q9QZQ8", "Q9QZQ8-2", "P0C0S6", "P70696", "Q64475", "Q6ZWY9", 
                "P10853", "Q64478", "Q8CGP1", "P10854", "Q8CGP2", "Q8CGP2-2", "Q64525", "Q64524", "Q9D2U9", "Q8CGP0", 
                "Q9D9Z7", "P68433", "P84228", "P84244", "P02301", "P62806", "Q9QYL0", "Q8CBB6", "Q497L1", "B2RVF0", 
                "Q3UA95", "Q149V4", "Q5M8Q2", "Q5SZA3", "B2RVP5", "B2RTM0", "A3KPD0", "Q9DAD9", "B9EI85", "B2RVD5", 
                "A9Z055", "Q5SWB1", "Q810S6", "B2RWH3", "Q149Z9", "A2AB79", "Q8R029", "S4R1E0", "J3QP08", "B2RTK3", 
                "Q8CGP4", "I7HFT9", "G3X9D5", "Q80ZM5", "A0A0G2JGI2", "A0A087WP11", "S4R1M3", "G3UX40", "A2BFR3", 
                "G3UWL7","S4R1G7", "A0A0G2JEV0", "F8WIX8", "L7MU04", "F8WI35", "E0CZ52", "E0CYL2", "E0CYN1", "E0CYR7", 
                "E0CZ27"]
  
  histCount = 0
  for keyS, itemL in completeDict.items(): # calculate the histone intensities in each lane first
    if keyS in histoneIDs:
      histCount += 1
      print(itemL)
      histIntL[0] += itemL[11]
      histIntL[1] += itemL[12]
      histIntL[2] += itemL[13]

    # print(itemL)
    totIntL[0] += itemL[11]
    totIntL[1] += itemL[12]
    totIntL[2] += itemL[13]
    

  # print histIntL
  # print totIntL
      
  # print len(histoneIDs), histCount
  
  
  resD = {} # this is the output will all data plus the absolute quantifications
  for keyS, itemL in completeDict.items(): # calculate protein abundances using the histone counts made above. write it all out to a new dict
    sampleCount = 0
    resD[keyS] = itemL
    for intI in itemL[11:]:
      Nprot = (float(intI) * nA * mDNA) / (float(itemL[7]) * 1000 * histIntL[sampleCount])
      sampleCount += 1
      resD[keyS].append(int(Nprot)) 
    resD[keyS].append(round(sum(resD[keyS][-3:])/3,0)) # add average

  
  print("protein numbers calculated")
  return resD
    
    

if __name__ == "__main__":
  master()