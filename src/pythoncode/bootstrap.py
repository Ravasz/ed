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
  
  proteomeF = open("/home/mate/code/ed/src/data/protein_lists/mouse_full_protein_list.txt","r")
  PFile = open("/home/mate/code/ed/src/data/protein_lists/mouse_PPP.txt","r")
  GFile = open("/home/mate/code/ed/src/data/protein_lists/mouse_PPG.txt","r")  
  
  # read files in to memory and make lists:
  
  proteomeL = []
  for protLine in proteomeF:
    protI = protLine.strip()
    proteomeL.append(protI)
    
  PL = []
  for protLine in PFile:
    protI = protLine.strip()
    PL.append(protI)  
    
  GL = []
  for protLine in GFile:
    protI = protLine.strip()
    GL.append(protI) # read files in to memory and make lists
  
  rangeNum = 1000
  subSampleSize = 1000

  countD = {"P":[],"G":[],"No":[]} # store results
  
  print("iterating " + str(rangeNum) + " times...")
  
  for i in range(1, rangeNum  + 1): # set number of iterations 
    print(".", end="")
    if i % 100 == 0: print("\n")
    currL = sample(proteomeL,subSampleSize) # take random subsample with n members
    PCount = 0
    GCount = 0
    NoCount = 0
    for protN in currL:
      foundFlag = False
      if protN in PL: 
        PCount += 1
        foundFlag = True
      if protN in GL: 
        GCount += 1
        foundFlag = True
      if not foundFlag: NoCount += 1
    countD["P"].append(PCount)
    countD["G"].append(GCount)
    countD["No"].append(NoCount)
  
  
  PAverage = round(sum(countD["P"])/len(countD["P"]),2)
  PSTD = statistics.stdev(countD["P"])
  
  GAverage = round(sum(countD["G"])/len(countD["G"]),2)
  GSTD = statistics.stdev(countD["G"])  
    
  NoAverage = round(sum(countD["No"])/len(countD["No"]),2)
  NoSTD = statistics.stdev(countD["No"])  
  
  print(str(i) + " iterations, with a random sample size of " + str(rangeNum) + "\n" )
  
  print("\nfound in PPP list: " + str(PAverage) + " +- " + str(round(PSTD,2)) + "\n")
  print("found in PPG list: " + str(GAverage) + " +- " + str(round(GSTD,2)) + "\n")  
  print("not present in either list: " + str(NoAverage) + " +- " + str(round(NoSTD,2)) + "\n")
    
list_combiner()
    


