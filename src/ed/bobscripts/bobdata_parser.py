'''
Created on 19 May 2015

@author: mate
find significantly enriched proteins in bob's data. 
- first parse groteinGroups.txt with file_parser to extract uniprot identifiers, gene names, unique+razor peptide counts, and ko1- ko2- ko3- wt1- wt2- wt3- LFQ intensitites.
- then, using entry_parser, remove duplicate uniprot IDs and only keep the shortest one, or the first one if tied, with the matching unique+razor count. copy that with everything else to a new file.
- then, using lfq_parser remove 0 (undetected) lfq values and replace them random whole numbers between 1 and 100
- then, use this dataset in R to calculate P values using the bob_ttest.R script there on the intensity values.
- finally, using the stat_parser function here, take the gene names that have a p value of at least 0.05 and pool them for a DAVID analysis.
- alternatively to stat_parser, volcano_plotter can be used to prepare data for the volcano plot R script
- yet another alternative is to use set_fdr on the R output to calculate false discovery rates
'''

def main():
  print "call a function here to start"
  lfq_parser()

def pipeline_runner():
  """do the analysis of proteomics datasets by calling the appropriate functions one after the other.
  Not written yet. """
  file_parser() # take raw data file and extract columns of interest. remove contaminants.
  entry_parser() # remove duplicates, faulty lines and format the whole thing normally.
  lfq_parser() # replace 0s in lfq reading with random small numbers for t testing purposes
  # open Rstudio and do T testing there
  from ed.tools import ROutputFormatter
  ROutputFormatter() # reformat R output to something more appealing, add FDR and fold change values
  
def kegg_converter():
  """process list of uniprot accessions for KEGG pathway analysis"""
  from ed.tools import prot_id_converter
  
  protList = []
  headerFlag = True
  with open("../bob/processed/24h_bobprots_up_full.csv","r") as inpF:
    for inpLine in inpF:
      if headerFlag:
        headerFlag = False
        continue
      inpList = inpLine.split(",")
      protList.append(inpList[1])
  print prot_id_converter(protList, outDB="kegggeneid")
  
def scrambler():
  """
  take in a csv with 6 columns of numbers, 
  and return a csv with 6 columns of the same numbers, 
  but all mixed up
  """
  import random
  with open("../bob/bob_lfq.csv", "r") as inpF:
    giantList = []
    headerFlag = True
    for inpLine in inpF:
      if headerFlag: 
        headerFlag = False
        continue
      inpL = inpLine.split(",")
      inpL[-1] = inpL[-1].strip()
      for inpI in inpL:
        giantList.append(float(inpI))

  shuffledList = random.sample(giantList, len(giantList)) # this does the randomization
  print "input scrambled"
  with open("../bob/bob_decoy_lfq.csv", "w") as outF:
    outF.write("LFQ scrambled 1,LFQ scrambled 2,LFQ scrambled 3,LFQ scrambled 4,LFQ scrambled 5,LFQ scrambled 6\n")
    lineCount = 0
    for listItem in shuffledList:
      if lineCount < 5:
        outF.write(str(listItem) + ",")
        lineCount += 1
      else:
        outF.write(str(listItem) + "\n")
        lineCount = 0
  print "output written to: " + str(outF.name)

def volcano_plotter():
  """take in a full list of genes and reformat them for the volcano plot R script. 
  output the reformatted data to a new file"""
  print "this is volcano plotter"
  from math import log
  with open("../bob/processed/24h_bobdata_ed2_volcano.csv", "w") as outF:
    outF.write("Gene log2FoldChange pvalue\n")
    with open("../bob/processed/24h_bobdata_ed2.csv", "r") as inpF:
      skipFlag = True
      missCount = 1
      for inpLine in inpF:
        if skipFlag:
          skipFlag = False
          continue
        inpLine = inpLine.split("\" \"")
        curLine = []
        for inpI in inpLine:
          try:
            curLine.append(float(inpI.strip("\"\n ")))
          except ValueError:
            curLine.append(inpI.strip("\"\n ")) # by this point, each line in the entry file is processed into a neat list
        if curLine[2] == "": # if no gene name is given, just add a placeholder
          curLine[2] = "Noname" + str(missCount)
          missCount += 1
        # calculate log2foldChange here:
        try:
          FAvg = (curLine[4] + curLine[5] + curLine[6])/3.0 # KO
          SAvg = (curLine[7] + curLine[8] + curLine[9])/3.0 # WT
        except TypeError:
          print curLine
          raise
        logFoldChange = log(SAvg/FAvg,2) # so positive numbers are more abundant in the wt cells, negatives number in the KO, at least for the 24H bobdata file
        outF.write(curLine[2] + " " + str(logFoldChange) + " " + str(curLine[10]) + "\n") # write out results to file

        
    


def set_fdr(fdrN = 0.05):
  """using the specified fdr number find all the p values (and corresponding genes in Bob's dataset)
  which meet the criteria outlined in the benjamini-hochberg FDR correction thing from 1995.
  Basically sort p values, and than start taking the smallest ones until one of them hits the corrected
  p value proposed by the equation in that paper:
  (i/m)Q, where i is the rank, m is the total number of tests, and Q is the false discovery rate you choose
  from here:  http://www.biostathandbook.com/multiplecomparisons.html
  This could be done in R too... I might do it there one day
  
  This reads the whole dataset into memory, and so it might be demanding.
  """
  print "this is set_fdr"
  def p_value_key(protItem):
    """mini function returning the last element of a list. just because I do not like unnamed functions"""
    return protItem[-1]
  
  protList = []
  curL = []
  headerFlag = True
  with open("../bob/processed/24h_bobdata_ed2.csv", "r") as inpF: # read and process the csv with protein names and p values
    for inpLine in inpF:
      if headerFlag: 
        headerFlag = False
        continue
      inpLine = (inpLine.rstrip().split("\" \""))
      for inpItem in inpLine:
        curL.append(inpItem.strip(" \""))
      if curL[-1] == "NaN":
        curL[-1] = 1
      curL[-1] = float(curL[-1])
      protList.append(curL)
      curL = []
      
  protList.sort(key = p_value_key) # sort the whole list on p value (lowest to highest) 
  i = 0.0 # see i and m in the function description
  m = float(len(protList))
  print "dataset length: ", m
  for protListI in protList:
    i += 1
    critVal = (i/m)*fdrN # this is the benjamini-hochberg defined critical value
    print "threshold: ", critVal # this is the adjusted p value the current measurement has to pass
    print "current p value: ", protListI[-1]
    if protListI[-1] < critVal:
      print protListI
    else:
      print "p value did not pass threshold. No other significant proteins in dataset."
      break


    

def interactor_finder():
  """take a list of protein names and check if they are in Bob's dataset"""
  from ed.tools import prot_id_converter

  proteinList = []
  with open("../datafiles/known_interactors.txt","r") as inpProt: # create list of gene names from hand-made text file with known ptp22 interactors
    for protLine in inpProt:
      if protLine != "\n":
        curName = protLine.strip().split("\t")[0]
        curName = curName[0] + curName[1:].lower()
        proteinList.append(curName)
  inpIdL = prot_id_converter(proteinList, "10090", "genesymbol", "uniprotaccession") # convert to uniprot accessions
  print inpIdL
  
  with open("../bob/processed/bobprots_all.csv","r") as targetF: # create list of all uniprot accessions in Bob's dataset (unique razor proteins only)
    targetD = {}
    for targetLine in targetF:
      targetD[targetLine.split(",")[0]] = targetLine.split(",")[1].strip()
  for inpIdItem in inpIdL:
    for queryI in inpIdItem:
      if queryI in targetD:
        print targetD[queryI]
        break
        
    
        
  


def stat_parser():
  """take protein names with a significant p value and out them to a result file"""
  from ed.tools import file_importer, file_outporter
  from math import log
  
  print "this is stat parser"
  
  relPath = "bob/processed/24h_bobdata_ed2.csv"
  outPathUp = "bob/processed/24h_bobprots_up_full.csv"
  outPathDown = "bob/processed/24h_bobprots_down_full.csv"
  inpF = file_importer(relPath)
  outFUp = file_outporter(outPathUp)
  outFDown = file_outporter(outPathDown)
  
  
  skipFlag = True
  
  for inpLine in inpF:
    if skipFlag:
      skipFlag = False
      outFDown.write("ID,Uniprot ID,Gene name,unique peptides (unique+razor),KO1,KO2,KO3,WT1,WT2,WT3,enrichment,P value\n")
      outFUp.write("ID,Uniprot ID,Gene name,unique peptides (unique+razor),KO1,KO2,KO3,WT1,WT2,WT3,enrichment,P value\n")
      continue
    inpLine = inpLine.split("\" \"")
    curLine = []
    for inpI in inpLine:
      curLine.append(inpI.strip("\"\n"))
    try: 
      curLine[-1] = float(curLine[-1])
    except ValueError:
      curLine[-1] = 1   
    if curLine[-1] < 0.05 and int(curLine[3]) > 1: # check if protein has at least 2 unique peptides and has a significant p value
      curLine[4:10] = [int(x) for x in curLine[4:10]]
      enrScore = log((sum(curLine[4:7]) / 3.0)/(sum(curLine[7:10]) / 3.0),2) # calculate log2 enrichment score
      # print int(sum(curLine[4:7]) / 3.0), int(sum(curLine[7:10]) / 3.0)
      if sum(curLine[4:7]) / 3.0 > sum(curLine[7:10]) / 3.0: # if the mean of the KO intensities is higher than the wt  
        for outI in curLine:
          outFDown.write(str(outI).strip(" "))
          if outI is not curLine[-1]:
            outFDown.write(",")
            if outI is curLine[-2]:
              outFDown.write(str(enrScore)+ ",")
          else:
            outFDown.write("\n")
        # outFDown.write(curLine[1] + "," + curLine[2] + "\n")
      else:
        # outFUp.write(curLine[1] + "," + curLine[2] + "\n")
        for outI in curLine:
          outFUp.write(str(outI).strip(" "))
          if outI is not curLine[-1]:
            outFUp.write(",")
            if outI is curLine[-2]:
              outFUp.write(str(enrScore)+ ",")
          else:
            outFUp.write("\n")
  
  inpF.close()
  outFUp.close()
  outFDown.close()
  print "stat_parser completed"

def protein_name_collector():
  """take all uniprot IDs (but only one per protein) from bob's dataset and return them as a single list"""
  resL = []
  with open("bob/processed/bobprots_down.csv", "r") as inpF:
    for inpLine in inpF:
      inpLine = inpLine.split(",")
      resL.append(inpLine[0].strip(" \n"))
  return resL

def lfq_parser():
  """remove 0 values from lfq measurements and replace them with a random number between 1 and 100
  This is needed for ttseting later in R, as each measurement there has to have some sort of noise in it"""
  from ed.tools import file_importer, file_outporter
  from random import randint
  
  print "this is lfq parser"
  
  relPath = "bob/processed/24h_bobdata_ed.csv"
  outPath = "bob/processed/24h_bobdata_no0_ed.csv"
  inpF = file_importer(relPath)
  outF = file_outporter(outPath)  
  headerFlag = True
  for inpLine in inpF:
    if headerFlag: 
      headerFlag = False
      outF.write(inpLine)
      continue
    inpLine = inpLine.strip()
    inpItems = inpLine.split(",")
    for inpI in inpItems[0:4]: # copy over gene name and such to new file
      outF.write(inpI)
      outF.write(",")
    
    commaCount = 0
    for inpJ in inpItems[4:]: # copy over lfq values while replacing 0-s with random values
      commaCount += 1
      if int(inpJ) == 0:
        randNum = randint(1,100)
        outF.write(str(randNum))
      else:
        outF.write(inpJ)
      if commaCount < 6:
          outF.write(",")
    outF.write("\n")
  inpF.close()
  outF.close()
  print "lfq parser finished successfully"

def entry_parser():
  """remove duplicate protein name and total peptide count cell entries from bob's dataset"""
  from ed.tools import file_importer, file_outporter
  from operator import add
  
  print "this is entry parser"
  
  relPath = "bob/processed/24h_bobdata.csv"
  outPath = "bob/processed/24h_bobdata_ed.csv"
  inpF = file_importer(relPath)
  outF = file_outporter(outPath)
  cN = 0
  outDict = {}
  for inpLine in inpF:
    cN += 1
    inpLine = inpLine.strip()
    inpItem = inpLine.split(",")
    geneL = inpItem[1].split(";")
    lenS = len(geneL[0])
    curGene = geneL[0]
    for geneI in geneL:
      if len(geneI) < lenS:
        lenS = len(geneI)
        curGene = geneI
    if "__" in curGene: continue
    try: # get rid of wonky lines introduced by excel
      int(curGene)
      continue
    except ValueError: 
      pass
    numL = inpItem[3].split(";")
    curNum = numL[geneL.index(curGene)]
    
    protL = inpItem[2].split(";") # process protein name here
    try:
      curProt = protL[geneL.index(curGene)]
    except IndexError:
      curProt = protL[0]
    if curProt == "":
      # print "no protein name found. adding the uniprot ID."
      curProt = curGene
  
    if curGene[-2] == "-":
      curGene = curGene[:-2]
    if curGene[-3] == "-":
      curGene = curGene[:-3]
    
    outL = [int(inpItem[0]),curGene,curProt,curNum] + inpItem[4:10] 
    
    try:
      for inpN in inpItem[4:10]:
        inpItem[inpItem.index(inpN)] = int(inpN)
      countFlag = True
    except ValueError:
      print inpItem[4:10]
      countFlag = False
    if countFlag:
      if sum(inpItem[4:10]) == 0: continue # there are some unexpressed proteins in there
      
    

    if curGene in outDict: # handle duplicate protein entries and merge them together
      # print "%s is duplicate" % curGene
      addL = []
      for i in outDict[curGene][3:]:
        addL.append(int(i))

      addL2 = []
      for j in outL[3:]:
        addL2.append(int(j))

      outL[3:] = map(add, addL, addL2) # admittedly this looks terrible

    outDict[curGene] = outL # assemble all the stuff in this dict
  
  outN = 0  
  for outDV in sorted(outDict.items(), key=lambda x:x[1][0]): # sort the dict based on the ID they have
    outN += 1
    # print outDV[1]
    # if outN == 100: break
    for outI in outDV[1][:-1]:
      outF.write(str(outI) + ",")
    outF.write(str(outDV[1][-1]) + "\n")
  

  print "unique proteins: ", outN
  print "lines parsed: ", cN
  inpF.close()
  outF.close()
  
def file_parser():
  """from bob"s proteinGroups.txt take: Majority protein IDs Peptide counts (razor+unique) ['LFQ intensity KO1', 'LFQ intensity KO2', 'LFQ intensity KO3', 'LFQ intensity WT1', 'LFQ intensity WT2', 'LFQ intensity WT3']
  and write them to a new file. do not select contaminants or reverse peptides"""

  from ed.tools import file_importer, file_outporter
  print "this is file parser"
  inpF = file_importer("bob/24h_proteingroups.csv")
  outF = file_outporter("bob/processed/24h_bobdata.csv")
  for inpLine in inpF:
    inpP = inpLine.split("\r")
    cN = 0
    print len(inpP)
    for inpI in inpP:
      inpItems = inpI.split("\t") 
      if inpItems[100] == "+" or inpItems[101] == "+": continue # get rid of contaminants and reverse proteins
      outF.write(str(cN) + "," + inpItems[1] + "," + inpItems[6] + "," + inpItems[3] + "," + inpItems[86] + "," + inpItems[87] + "," + inpItems[88] + "," + inpItems[89] + "," + inpItems[90] + "," + inpItems[91] + "\n")
      # print inpItems [1],"+++", inpItems [3],"+++", inpItems [6],"+++", inpItems[86:92]
      cN += 1
      # if cN == 40: break

    break

  inpF.close()
  outF.close()
  print cN

  
if __name__ == "__main__":
  main()