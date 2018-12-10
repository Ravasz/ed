'''
Created on 10 Dec 2018

@author: mate
'''


def proteome_list_maker():
  """take a complete proteome file used as mass spectrometry dataset and take out all the protein names. Store protein names in a file."""
  
  geneL = []
  
  outFile = open("/home/mate/code/ed/src/data/protein_lists/mouse_full_protein_list.txt","w")
  
  with open("/home/mate/code/ed/src/data/protein_lists/uniprot-proteome_3AUP000000589.fasta","r") as inpFile:
    for inpLine in inpFile:
      if inpLine.startswith(">"): # find fasta headers
        inpList = inpLine.split("|") # split header up into sections
        inpHeadL = inpList[2].split(" ") # split the protein description up
        for headI in inpHeadL:
          if headI.startswith("GN="): # find gene name in the protein description
            geneN = headI.split("=")[1].strip() 
            if geneN not in geneL: 
              geneL.append(geneN) # store just the gene name
            
              # print(geneN)
        
  outList = sorted(geneL) # sort gene names in alphabetical order
  
  for outI in outList:
    outFile.write(outI + "\n") # write them out to a file
  
  outFile.close()
  
def sublist_maker():
  """take out two columns from Thomas's .csv file and save them as separate lists"""

  outFileP = open("/home/mate/code/ed/src/data/protein_lists/mouse_PPP.txt","w")
  outFileG = open("/home/mate/code/ed/src/data/protein_lists/mouse_PPG.txt","w")
  
  with open("/home/mate/code/ed/src/data/protein_lists/Bootstrap_genelist.csv","r") as inpFile:  
    next(inpFile)
    for inpLine in inpFile:
      inpList = inpLine.split(",")
      if inpList[2] != "\n": outFileP.write(inpList[2].strip() + "\n")
      if inpList[3] != "\n": outFileG.write(inpList[3].strip() + "\n")
      
  outFileP.close()
  outFileG.close()
  
def list_combiner():
  """take in a full proteme gene list and two sublists, and do a bootstrap analysis on how samples of the proteome list overlap with the sublists"""
  
  from random import sample
  import statistics
  
  # read files in to memory and make lists:
  with open("/home/mate/code/ed/src/data/protein_lists/mouse_full_protein_list.txt","r") as proteomeF:
    proteomeL = []
    for protLine in proteomeF:
      protI = protLine.strip()
      proteomeL.append(protI)

  with open("/home/mate/code/ed/src/data/protein_lists/mouse_PPP.txt","r") as PFile:
    PL = []
    for protLine in PFile:
      protI = protLine.strip()
      PL.append(protI)  
    PFile.close()
    
  with open("/home/mate/code/ed/src/data/protein_lists/mouse_PPG.txt","r") as GFile:
    GL = []
    for protLine in GFile:
      protI = protLine.strip()
      GL.append(protI) 
    GFile.close() 
  
  rangeNum = 1000
  subSampleSize = 2000

  countD = {"P":[],"G":[], "Both":[],"No":[]} # store results
  
  print("iterating " + str(rangeNum) + " times...")
  
  for i in range(1, rangeNum  + 1): # set number of iterations 
    print(".", end="")
    if i % 100 == 0: print("\n")
    currL = sample(proteomeL,subSampleSize) # take random subsample with n members
    PCount = 0
    GCount = 0
    BothCount = 0
    NoCount = 0
    for protN in currL:
      foundFlag = False
      if protN in PL:
        if protN in GL:
          BothCount += 1
        else:
          PCount += 1
        foundFlag = True
        
      if protN in GL and not protN in PL: 
        GCount += 1
        foundFlag = True
      
      if not foundFlag: NoCount += 1
      
    countD["P"].append(PCount)
    countD["G"].append(GCount)
    countD["Both"].append(BothCount)
    countD["No"].append(NoCount)
  
  resultL = zip(countD["P"],countD["G"],countD["Both"],countD["No"])
  
  with open("/home/mate/code/ed/src/data/protein_lists/PPP_PPG_occurrences_it-" + str(rangeNum) + "_subsample-" + str(subSampleSize) + ".txt","w") as outF:
    outF.write("PPP,PPG,Both,Neither\n")
    
    for outS in resultL:
      for outI in outS:
        if outI is outS[-1]: outF.write(str(outI) + "\n") 
        else: outF.write(str(outI) + ",") 
  
  PAverage = round(sum(countD["P"])/len(countD["P"]),2)
  PSTD = statistics.stdev(countD["P"])
  
  GAverage = round(sum(countD["G"])/len(countD["G"]),2)
  GSTD = statistics.stdev(countD["G"])  

  BothAverage = round(sum(countD["Both"])/len(countD["Both"]),2)
  BothSTD = statistics.stdev(countD["Both"])  
    
  NoAverage = round(sum(countD["No"])/len(countD["No"]),2)
  NoSTD = statistics.stdev(countD["No"])  

  
  print(str(i) + " iterations, with a random sample size of " + str(subSampleSize) + "\n" )
  
  print("\nfound in only in PPP list: " + str(PAverage) + " +- " + str(round(PSTD,2)) + "\n")
  print("found in only in PPG list: " + str(GAverage) + " +- " + str(round(GSTD,2)) + "\n")  
  print("found in both lists: " + str(BothAverage) + " +- " + str(round(BothSTD,2)) + "\n")  
  print("not found in either list: " + str(NoAverage) + " +- " + str(round(NoSTD,2)) + "\n")
    
list_combiner()
    


