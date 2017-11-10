'''
Created on 12 Oct 2017

@author: mate

find significantly enriched proteins in maxquant output data. This script is a modification of bobdata_parser, with all the features there, but adapted for any proteingroups.txt file from maxquant.
 
- first parse groteinGroups.txt with file_parser to extract uniprot identifiers, gene names, unique+razor peptide counts, and ko1- ko2- ko3- wt1- wt2- wt3- LFQ intensitites.
- then, using entry_parser, remove duplicate uniprot IDs and only keep the shortest one, or the first one if tied, with the matching unique+razor count. copy that with everything else to a new file.
- then, using lfq_parser remove 0 (undetected) lfq values and replace them random whole numbers between 1 and 100
- then, use this dataset in R to calculate P values using the bob_ttest.R script there on the intensity values.
- finally, using the stat_parser function here, take the gene names that have a p value of at least 0.05 and pool them for a DAVID analysis.
- alternatively to stat_parser, volcano_plotter can be used to prepare data for the volcano plot R script
- yet another alternative is to use set_fdr on the R output to calculate false discovery rates


cfg file name: proteingroups_analyzer_params.cfg

cfg file layout:
<Files>
/home/mate/code/ed/src/data/cav1ko/txt_cav1ko-1-17082017/proteinGroups.txt
/home/mate/code/ed/src/data/cav1ko/txt_cav1ko_11_10_2017/proteinGroups.txt

</Files>

<Outputfolder>
/home/mate/code/ed/src/data/cav1ko/processed
</Outputfolder>

<Samples>
<Ungrouped>
</Ungrouped>
<Group1>
</Group1>
<Group2>
</Group2>
<Group3>
</Group3>
<Group4>
</Group4>
</Samples>

'''

def main():
  print("call a function here to start")
  #file_analyzer()
  file_combiner()
  
  

def file_analyzer():
  """
  extract column names from the input files in proteingroups_analyzer_params.cfg
  and extract samples names from them.
  Write sample names into the cfg file into the ungrouped section for the user to group.
  
  needs column_finder, empty_line_remover, and file_picker
  """
  import os.path, sys
  
  global fileCount
  fileCount = 1
  if os.path.isfile("proteingroups_analyzer_params.cfg"): # check if cfg file exists
    pass
  else:
    print("the config file named proteingroups_analyzer_params.cfg was not found.\n One needs to be written with the names and full paths of files to be analyzed and placed in the same directory as this script")
    sys.exit(0) # end the program if config file is not found
    
  with open("proteingroups_analyzer_params.cfg","r") as cfgF: # locate all the groups in cfgFile
    groupList = []
    samplesFlag = False
    for cfgLine in cfgF:
      if cfgLine == "<Samples>\n": 
        samplesFlag = True
        continue
      elif cfgLine == "</Samples>\n" : 
        samplesFlag = False
        continue
      
      if samplesFlag and "<Group" in cfgLine:
        groupList.append(cfgLine.rstrip("\n"))
  
  empty_line_remover("<Ungrouped>") # remove all empty lines from cfgFile groups section
  for groupI in groupList:
    empty_line_remover(groupI)

  for fileI in file_picker("proteingroups_analyzer_params.cfg"):
    column_finder(fileI)
    fileCount += 1

def file_combiner():
  """take all the files from the cfg file and combine them to a single file. 
  Read in the next file line by line.
  Remove contaminants and reverse proteins.
  Use uniprot ID as key and check if the protein is already present.
  add new entry to existing columns or a new line.
  """

  from collections import defaultdict
  from copy import copy
  import sys, os.path
  import pandas as pd
  
  pd.set_option("display.expand_frame_repr", False)
  
  print("this is file_combiner")
  

  inpFileList = []
  groupList = []
  
  with open("proteingroups_analyzer_params.cfg","r") as cfgFile:

    for cfgLine in cfgFile:
      
      if cfgLine == "<Files>\n": # collect file names here
        lineCount = 0
        while True:
          newLine = next(cfgFile)
          lineCount += 1
          if lineCount == 20: 
            print("something terrible has happened!")
            sys.exit(0)
          if newLine == "\n": continue
          elif newLine == "</Files>\n":
            break
          else: 
            inpFileList.append(newLine.rstrip("\n"))
        continue
      
      if cfgLine == "<Outputfolder>\n": # collect output folder here
        lineCount = 0
        while True:
          newLine = next(cfgFile)
          lineCount += 1
          if lineCount == 20: 
            print("something terrible has happened!")
            sys.exit(0)
          if newLine == "\n": continue
          else: 
            outputFolder = newLine.rstrip("\n")
            break
            
      
      if cfgLine == "<Samples>\n": # place all samples from all groups into a single list 
        lineCount = 0
        while True:
          newLine = next(cfgFile)

          if newLine.rstrip("\n") == "</Samples>": break  
          lineCount += 1
          if lineCount == 40: 
            print("something terrible has happened!")
            sys.exit(0)
          if newLine == "\n": continue
          if "<Group" in newLine:
            subLineCount = 0
            while True:
              newLine = next(cfgFile)
              if "</Group" in newLine:
                break
              subLineCount += 1
              if subLineCount == 10: 
                print("something terrible has happened!")
                sys.exit(0)
              if newLine == "\n": continue
              
              groupList.append(newLine.rstrip("\n"))

  if len(groupList) < 2:
    print("there should be at least 2 groups to compare. Please assign samples to groups in the cfg file.")
    sys.exit(0)
  
  print(inpFileList)
  print(groupList)
      
  if os.path.isdir(outputFolder):    # check if output folder is correctly extracted from cfg file
    print(outputFolder)
  else:
    print("output folder named %s not found. please check the outputfolder area in the config file and add a proper folder name in there" % (outputFolder,))
    sys.exit(0)

  # by this point everything is collected from the cfg file that we need      
  
  outF = open(os.path.join(outputFolder, "combined.csv"),"w")
  
  # write header of output file: 
  
  outF.write("Majority protein ID,Protein name,Gene name,") 
  for j in groupList:
    outF.write("Peptides " + j + ",")
  
  for j in groupList:
    outF.write("Razor + unique peptides " + j + ",")
    
  for j in groupList:
    outF.write("LFQ intensity " + j + ",")
    
  outF.write("Fold change (Group1/Group2),Pvalue(Group1/Group2)\n")
  
  # header written
    
  
  
  dataCounter = 0
  fileCount = 0
  colNums = []
  finDict = defaultdict(list)
 
  for dataFile in inpFileList:
    dataCounter += 1
    currGroupL = []   # select samples to analyze in each file
    currColCount = 0 # count how many columns there are to analyze in each file
    for groupI in groupList:
      groupN = int(groupI.split("-")[-1])
      if groupN == dataCounter:
        currGroupL.append(groupI[:groupI.rindex("-")])
        currColCount += 1
    colNums.append(currColCount)
    # print currGroupL
    
    colL = [] # list of columns that are to be collected
    colL += ["Majority protein IDs","Protein names","Gene names"]
    for sI in currGroupL:
      colL.append("Peptides "+ sI)
    for sI in currGroupL:
      colL.append("Razor + unique peptides "+ sI)
    for sI in currGroupL:
      colL.append("LFQ intensity "+ sI)

    
     
    with open(dataFile,"r") as inpF:
   
      cN = 0
      outDict = {}
      headerFlag = True
      neededColIndexes = []
      neededColNames = []
      fileCount += 1
          
      for inpLine in inpF:
        cN += 1
        if headerFlag:
          headerFlag = False
          inpItem = inpLine.rstrip("\r\n").split("\t")
          for inpI in inpItem: # collect the positions and names of columns that are to be copied into the results file
            if inpI in colL:
              neededColIndexes.append(inpItem.index(inpI))
              neededColNames.append(inpI)
          # print neededColIndexes
          continue
        
        inpLine = inpLine.rstrip("\r\n")
        inpItem = inpLine.split("\t")
        geneL = inpItem[neededColIndexes[0]].split(";")
        lenS = len(geneL[0])
        curGene = geneL[0]
        for geneI in geneL: # find gene name with the shortest length
          if len(geneI) < lenS:
            lenS = len(geneI)
            curGene = geneI
        if "__" in curGene: continue # get rid of contaminant lines
        corrGeneN = geneL.index(curGene)

 
        if curGene[-2] == "-":
          curGene = curGene[:-2]
        elif curGene[-3] == "-":
          curGene = curGene[:-3]
        
        workingL = []  
        for neededN in neededColIndexes:
          workingL.append(inpItem[neededN])
        
        # print workingL
        
        uniqueL = []
        for workingI in workingL[:3]: # remove ambiguities 
          workingIL = workingI.split(";")
          try:
            uniqueL.append(workingIL[corrGeneN])
          except IndexError:
            uniqueL.append(workingIL[0])
        uniqueL.extend(workingL[3:])

        if uniqueL[2] == "":
          # print "no protein name found. adding the uniprot ID."
          uniqueL[2] = curGene
        
        # print uniqueL

        if curGene in outDict: # handle duplicate protein entries and merge them together
          print("%s is duplicate" % curGene)
          mergedL = []
          for mI in outDict[curGene][:3]:
            mergedL.append(mI)
          
          posCount = 3  
          finN = 3+(len(currGroupL)*2)
          for mI in outDict[curGene][3:finN +1]:
            mergedL[posCount] = max(int(outDict[curGene][posCount]),int(uniqueL[posCount])) # take higher number of peptide and unique peptide counts
            posCount += 1
          for mI in outDict[curGene][finN +1:]:
            mergedL[posCount] = sum(int(outDict[curGene][posCount]),int(uniqueL[posCount])) # add up LFQ scores
            posCount += 1
          outDict[curGene] = mergedL # add in unified scores
          
        else: outDict[curGene] = uniqueL
      
      print(outDict["P20934"])
      print(fileCount)
      
      # so it collects the relevant data from each of the files. As a next step
      # these dicts have to be merged to a single dict and then written to a file
      # adding in pandas here instead of fussing with dicts more
      
    outDF = pd.DataFrame.from_dict(outDict, orient = "index")
    del outDF[0]
    # print(outDF)
    outDF.columns = neededColNames[1:]
    outDF.index.name = neededColNames[0]
    print(outDF)

      
#       for outKey,outValue in outDict.items(): 
#         if outKey not in finDict:
#           finDict[outKey] = outDict[outKey][:3]
#           for i in range(1,fileCount): # add zeroes to the line to represent previous iteration
#             for j in range(colNums[fileCount-1]):
#               finDict[outKey].append("0") # this bit isn't working yet

            
            
#         if outKey in finDict:
#           if len(finDict[outKey]) == fileCount -1: # if protein was present in previous file
            

#     for outKey,outValue in outDict.items(): 
#       if outKey in finDict: # add modified dicts together into single, unified dict
#         # print fileCount, finDict[outKey]
#         # print outValue
#         outIndex = 0
#         for outItem in outValue:
#           finDict[outKey][outIndex].append(outItem)
#           outIndex += 1
#         # print finDict[outKey]
#    
#       else:  # or just add new entries
#         if fileCount == 1:
#           for outItem in outValue:
#             finDict[outKey].append([outItem])
#          
#         else: # fill up entries that were not present in the previous cycle
#           loopCount = 0
#           while loopCount < fileCount - 1:
#             for i in range(len(outValue)):
#               if len(finDict[outKey]) == i:
#                 finDict[outKey].append([])
#               else:
#                 # print finDict[outKey]
#                 finDict[outKey][i].append("")
#             loopCount += 1
#           outIndex = 0
#           for outItem in outValue:
#             # print finDict[outKey]
#             finDict[outKey][outIndex].append(outItem)          
#             outIndex += 1
#             
# #     for fKey, fValue in finDict.items():
# #       pass
# #       print fValue
#    
#     for testKey in finDict: # fill up entries in result dict which were not present in previous file
#       if len(finDict[testKey][0]) < fileCount - 1:
#         for i in range(len(finDict[testKey])):
#           if i < 3:
#             finDict[testKey][i].append("")
#           else:
#             finDict[testKey][i].append("0")
#           
# #     for fKey, fValue in finDict.items():
# #       pass
# #       print fValue

     
         
       
#       outN = 0  
#       # prepare header for file:
#       headList = headerLine.strip("\n\r").split("\t")
#       if fileCount > 1:
#         for headerItem in headList[:-1]:
#           headerI = headerItem.replace(",",".")
#           headerCount = 1
#           while headerCount < fileCount:
#             outF.write(headerI + "-" + str(headerCount) + "|")
#             headerCount += 1  
#           outF.write(headerI + "-" + str(headerCount) + "\t")
#            
#         headerCount = 1
#         while headerCount < fileCount:
#           outF.write(headList[-1] + "-" + str(headerCount) + "|")
#           headerCount += 1
#          
#         outF.write(headList[-1] + "-" + str(headerCount) + "\n")
#      
#       elif fileCount == 1:
#         for headerItem in headList[:-1]:
#           headerI = headerItem.replace(",",".")    
#           outF.write(headerI + "\t")
#         outF.write(headList[-1].replace(",",".") + "\n")
#        
#       else:
#         print "number of input files should be at least one. Got less somehow"
#         raise ValueError
#     
#   
#   for outDK, outDV in finDict.items(): # write out assembled results to a file
#     outN += 1
#     if len(outDK) > 30: print "this line should not be displayed"
#     # print outDV[1]
#     # if outN == 100: break
#     nameCount = 0
#     for outI in outDV:
#       # if nameCount == 0: print outI
#       for outPiece in outI[:-1]:
#         outU = outPiece.replace(",",".")
#         if outU == "": outF.write("_|")
#         else: outF.write(str(outU) + "|")
#       if outI[-1] == "": # handle missing entries
#         if nameCount == 6: outF.write(outDV[0][0] + "\t") # replace missing gene names with their uniprot ID
#         else: outF.write("_\t")
#       else: outF.write(str(outI[-1]).replace(",",".") + "\t")
#       nameCount += 1
#     outF.write("\n")
#   
# 
#   print "unique proteins: ", outN
#   print "lines parsed: ", cN
#   # print headerLine
#   inpF.close()
#   outF.close()
  

    
def file_picker(cfgS):
  """open files for analysis from proteingroups_analyzer_params.cfg in the same directory"""
  
  
  print("this is file_picker")
  
  import os.path, sys
  
  fileL = []
  with open(cfgS,"r") as cfgFile:
    fileList = False
    fileListComplete = False
    for fileLine in cfgFile: # extract file names from cfg file here
      lineS = fileLine.rstrip("\r\n")
      if lineS == "<Files>":
        fileList = True
        continue
      if lineS == "": continue
      if lineS == "</Files>": 
        fileListComplete = True
        break
      if fileList:
        if not os.path.isfile(lineS):
          print("file %s not found. Please check proteingroups_analyzer_params.cfg and try again." % (lineS,))
          sys.exit(0)
        fileL.append(lineS)
    if not fileListComplete or not fileList:
      print("something is wrong in the cfg file. Could not find file block. Please check if <File> and </File> tags are present")
      sys.exit(0)
  
  
  
  return fileL

def empty_line_remover(markerS):
  """helper function to remove empty lines from cfg file. Needs markerS as argument, like "<Ungrouped>" """
  
  import sys
  
  extraLinesFlag = False
  with open("proteingroups_analyzer_params.cfg","r") as outF:
    contentsL = outF.readlines()
    ungroupedStart = contentsL.index(markerS + "\n")
    ungroupedEnd = contentsL.index("</" + markerS[1:] + "\n")
    if not ungroupedStart + 1 == ungroupedEnd:
      extraLinesFlag = True
      extraLines = contentsL[ungroupedStart + 1:ungroupedEnd]
      for extraLine in extraLines:
        if not extraLine.rstrip("\r\n") == "":
          print("there are some %s samples in the config file already. Please clear them before starting again" % (markerS[1:-1].lower()))
          sys.exit(0)
  
  if extraLinesFlag:
    with open("proteingroups_analyzer_params.cfg","w") as outF:
      firstHalf = contentsL[:ungroupedStart + 1]
      lastHalf = contentsL[ungroupedEnd:]
      outF.writelines(firstHalf)
      outF.writelines(lastHalf)
    
def column_finder(filePath):
  """identify columns for analysis in all other functions"""
  
  import os.path
  import sys

  print("\nchecking ",filePath, " ...\n")
    
  currF = open(filePath,"r")
  headerI = currF.readline().rstrip("\r\n")
  headerL = headerI.split("\t")
  # print headerL
  

  sampleL = []
  for headI in headerL:
    if headI[:9] == "Peptides ":
      # print headI[9:]
      sampleL.append(headI[9:])
  print("samples found:")
  sampleCount = 1
  for sampleS in sampleL:
    print("%d) %s-%d" % (sampleCount,sampleS,fileCount,))
    sampleCount += 1
  currF.close()
  
  with open("proteingroups_analyzer_params.cfg","r") as outF:
    contentsL = outF.readlines()
    global outputFolder
    outputFolder = contentsL[contentsL.index("<Outputfolder>\n") + 1].rstrip("\r\n")
    if not os.path.isdir(outputFolder):
      print(outputFolder, " was not found or is not a folder")
      sys.exit(0)
    for sampleI in sampleL:
      unGroupedN = contentsL.index("</Ungrouped>\n")
      contentsL.insert(unGroupedN, "%s-%d\n" % (sampleI,fileCount,))
  with open("proteingroups_analyzer_params.cfg","w") as outF:
    outF.writelines(contentsL)

def ROutputFormatter():
  """take a terrible output file from R and format it it in a more nice way, 
  like remove leftover spaces and commas in it then add fold change and FDR score"""
  from math import log
  from tools import file_importer, file_outporter
  
  fdrN = 0.05
  def p_value_key(protItem):
    """mini function returning the last element of a list. just because I do not like unnamed functions"""
    return protItem[-1]
  
  protList = []
  headerFlag = True
  missCount = 1
  inpF = file_importer("bob/processed/OST-24-05-2017_combined_ttest_2.csv")  # read and process the csv with protein names and p values
  for inpLine in inpF:
    if headerFlag: 
      headerFlag = False
      continue
    inpLine = inpLine.split("\" \"")
    curLine = []
    
    for inpI in inpLine:
      if inpI.strip("\"\n ") == "NaN":
        inpI = 1
      try:
        curLine.append(int(inpI.strip("\"\n ")))
      except ValueError:
        try:
          curLine.append(float(inpI.strip("\"\n ")))
        except ValueError:
          curLine.append(inpI.strip("\"\n ")) # by this point, each line in the entry file is processed into a neat list
    if curLine[2] == "" or curLine[2] == "_": # if no gene name is given, just add a placeholder
      curLine[2] = "Noname" + str(missCount)
      missCount += 1

    protList.append(curLine)
      
  protList.sort(key = p_value_key) # sort the whole list on p value (lowest to highest) 
  i = 0.0 
  m = float(len(protList))
  print("dataset length: ", int(m))
  outputF = file_outporter("bob/processed/OST-24-05-2017_combined_ttest_ed_2.csv")
  outputF.write("ID,UniprotID,Gene name,OST1,OST2,OST3,WT1,WT2,WT3,pValue-wilcox,FDR,log2FoldChange\n")
  for protListI in protList:
    i += 1
    critVal = (i/m)*fdrN # this is the benjamini-hochberg defined critical value
    protListI.append(critVal)
    try:
      FAvg = (protListI[3] + protListI[4] + protListI[5])/3.0 # OST
      SAvg = (protListI[6] + protListI[7] + protListI[8])/3.0 # OT1
    except TypeError:
      print(curLine)
      raise
    try:
      logFoldChange = log(FAvg/SAvg,2) # so positive numbers are more abundant in the OST cells, negatives number in the OT1 cells, at least for the OST IP mass spec file
    except ZeroDivisionError:
      logFoldChange = log(FAvg/0.5,2)
    protListI.append(logFoldChange)
    
    for outI in protListI:
      outputF.write(str(outI))
      if outI is protListI[-1]:
        outputF.write("\n")
      else:
        outputF.write(",")
      
  print("formatting complete")
  
def kegg_converter():
  """process list of uniprot accessions for KEGG pathway analysis"""
  from tools import prot_id_converter
  
  protList = []
  headerFlag = True
  with open("../bob/processed/24h_bobprots_up_full.csv","r") as inpF:
    for inpLine in inpF:
      if headerFlag:
        headerFlag = False
        continue
      inpList = inpLine.split(",")
      protList.append(inpList[1])
  print(prot_id_converter(protList, outDB="kegggeneid"))
  
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
  print("input scrambled")
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
  print("output written to: " + str(outF.name))

def volcano_plotter():
  """take in a full list of genes and reformat them for the volcano plot R script. 
  output the reformatted data to a new file"""
  print("this is volcano plotter")
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
          print(curLine)
          raise
        logFoldChange = log(SAvg/FAvg,2) # so positive numbers are more abundant in the wt cells, negatives number in the KO, at least for the 24H bobdata file
        outF.write(curLine[2] + " " + str(logFoldChange) + " " + str(curLine[10]) + "\n") # write out results 

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
  print("this is set_fdr")
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
  print("dataset length: ", m)
  for protListI in protList:
    i += 1
    critVal = (i/m)*fdrN # this is the benjamini-hochberg defined critical value
    print("threshold: ", critVal) # this is the adjusted p value the current measurement has to pass
    print("current p value: ", protListI[-1])
    if protListI[-1] < critVal:
      print(protListI)
    else:
      print("p value did not pass threshold. No other significant proteins in dataset.")
      break

def interactor_finder():
  """take a list of protein names and check if they are in Bob's dataset"""
  from tools import prot_id_converter

  proteinList = []
  with open("../datafiles/known_interactors.txt","r") as inpProt: # create list of gene names from hand-made text file with known ptp22 interactors
    for protLine in inpProt:
      if protLine != "\n":
        curName = protLine.strip().split("\t")[0]
        curName = curName[0] + curName[1:].lower()
        proteinList.append(curName)
  inpIdL = prot_id_converter(proteinList, "10090", "genesymbol", "uniprotaccession") # convert to uniprot accessions
  print(inpIdL)
  
  with open("../bob/processed/bobprots_all.csv","r") as targetF: # create list of all uniprot accessions in Bob's dataset (unique razor proteins only)
    targetD = {}
    for targetLine in targetF:
      targetD[targetLine.split(",")[0]] = targetLine.split(",")[1].strip()
  for inpIdItem in inpIdL:
    for queryI in inpIdItem:
      if queryI in targetD:
        print(targetD[queryI])
        break
        
def stat_parser():
  """take protein names with a significant p value and out them to a result file"""
  from tools import file_importer, file_outporter
  from math import log
  
  print("this is stat parser")
  
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
  print("stat_parser completed")

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
  # from tools import file_importer, file_outporter
  from random import randint
  import os.path
  # from math import log10
  
  print("this is lfq parser")
  
  """
  relPath = "bob/processed/OST-24-05-2017_combined.csv"
  outPath = "bob/processed/OST-24-05-2017_combined_no0_1000.csv"
  inpF = file_importer(relPath)
  outF = file_outporter(outPath)  
  """
  inpF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "cav1ko", "processed", "cav1ko-1.csv"),"r")
  outF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "cav1ko", "processed", "cav1ko-1_no0.csv"),"w")
  
  headerFlag = True
  rowCount = 0
  for inpLine in inpF:
    if headerFlag: 
      headerFlag = False
      lfqColCount = 0
      inpList = inpLine.split("\t")
      for headerI in inpList:
        if "LFQ intensity" in headerI:
          break
        else: lfqColCount += 1
      outF.write("ID,Protein ID, Gene name,") # write header into output file
      if len(inpList[lfqColCount].split("|")) > 1:
        for headI in inpList[lfqColCount].split("|"):
          outF.write(headI + ",")
        for headI in inpList[lfqColCount + 1].split("|")[:-1]:
          outF.write(headI + ",")
        outF.write(inpList[lfqColCount + 1].split("|")[-1] + "\n")
      else:
        outF.write(inpList[lfqColCount] + ",")
        outF.write(inpList[lfqColCount + 1] + ",")
        outF.write(inpList[lfqColCount + 2] + "\n")
      rowCount += 1
      continue
    
    outF.write(str(rowCount) + ",")
    
    inpLine = inpLine.strip()
    inpItems = inpLine.split("\t")
    inpName = max(inpItems[0].split("|"), key=len) # get unique protein ID
    inpGeneName = max(inpItems[6].split("|"), key=len) # and gene name
    outF.write(inpName + "," + inpGeneName + ",")
    
    inpLFQ = inpItems[lfqColCount].split("|") + inpItems[lfqColCount + 1].split("|") + inpItems[lfqColCount + 2].split("|")# get lfq intensity scores
    # print inpLFQ
    for lfqI in inpLFQ[:-1]: # write lfq values
      if lfqI == "_" or lfqI == "0":
        outF.write(str(randint(1,1000)) + ",") ################## try with log10 values this time
      else:
        try:
          outF.write(str(int(lfqI)) + ",")
        except ValueError:
          print(inpItems)
          raise
    
    if inpLFQ[-1] == "_" or inpLFQ[-1] == "0": outF.write(str(randint(1,1000)) + "\n")
    else: outF.write(str(inpLFQ[-1]) + "\n")
    
    
    
    rowCount += 1

    
    
    """
    
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
  
  """
  
  print("lfq parser finished successfully")
  
def lfq_parser_2x():
  """remove 0 values from lfq measurements and replace them with a random number between 1 and 100
  This is needed for ttseting later in R, as each measurement there has to have some sort of noise in it
  
  Only include hits which appear at least in two OST samples
  """
  from tools import file_importer, file_outporter
  # from random import random
  from math import log10
  
  print("this is lfq parser_2x")
  
  relPath = "bob/processed/OST-24-05-2017_combined.csv"
  outPath = "bob/processed/OST-24-05-2017_combined_no0_2.csv"
  inpF = file_importer(relPath)
  outF = file_outporter(outPath)  
  headerFlag = True
  rowCount = 0
  for inpLine in inpF:
    if headerFlag: 
      headerFlag = False
      lfqColCount = 0
      inpList = inpLine.split("\t")
      for headerI in inpList:
        if "LFQ intensity" in headerI:
          break
        else: lfqColCount += 1
      outF.write("ID,Protein ID, Gene name,") # write header into output file
      for headI in inpList[lfqColCount].split("|"):
        outF.write(headI + ",")
      for headI in inpList[lfqColCount + 1].split("|")[:-1]:
        outF.write(headI + ",")
      outF.write(inpList[lfqColCount + 1].split("|")[-1] + "\n")
      rowCount += 1
      continue
    
    inpLine = inpLine.strip()
    inpItems = inpLine.split("\t")    
    inpLFQ = inpItems[lfqColCount].split("|") + inpItems[lfqColCount + 1].split("|") # get lfq intensity scores
    
    procLFQ = []
    for lfqi in inpLFQ:
      if lfqi == "_": procLFQ.append(0)
      else: procLFQ.append(int(lfqi))
    if sum(procLFQ[:3])<=sum(procLFQ[3:]): continue # get rid of contaminants in control sample
    
    sumOne = inpLFQ[1] + inpLFQ[2]
    sumTwo = inpLFQ[1] + inpLFQ[3]
    sumThree = inpLFQ[2] + inpLFQ[3]
    
    if sumOne == "__" or sumTwo == "__" or sumThree == "__": continue # test if protein is being detected in at least 2 OST samples
        
    outF.write(str(rowCount) + ",")
    

    inpName = max(inpItems[0].split("|"), key=len) # get unique protein ID
    inpGeneName = max(inpItems[6].split("|"), key=len) # and gene name
    outF.write(inpName + "," + inpGeneName + ",")
    
    inpLFQ = inpItems[lfqColCount].split("|") + inpItems[lfqColCount + 1].split("|") 
    # print inpLFQ
    for lfqI in inpLFQ[:-1]: # write lfq values
      if lfqI == "_" or lfqI == "0":
        outF.write(str(0) + ",") ################## try with log2 values this time
      else:
        try:
          outF.write(str(round(log10(int(lfqI)))) + ",")
        except ValueError:
          print(inpItems)
          raise
    
    if inpLFQ[-1] == "_" or inpLFQ[-1] == "0": outF.write(str(0) + "\n")
    else: outF.write(str(round(log10(int(inpLFQ[-1])))) + "\n")
    
    
    rowCount += 1

  print("lfq parser 2x finished successfully")

def spectrum_parser():
  """remove 0 values from spectral counts and replace them with a random number between 0 and 1
  This is needed for ttseting later in R, as each measurement there has to have some sort of noise in it
  
  This function is a modification of the lfq_parser()
  """
  from tools import file_importer, file_outporter
  from random import random
  # from math import log10
  
  print("this is spectrum parser")
  
  relPath = "bob/processed/OST-24-05-2017_combined.csv"
  outPath = "bob/processed/OST-24-05-2017_combined_no0_spectrum.csv"
  inpF = file_importer(relPath)
  outF = file_outporter(outPath)  
  headerFlag = True
  rowCount = 0
  for inpLine in inpF:
    if headerFlag: 
      headerFlag = False
      spColCount = 0
      inpList = inpLine.split("\t")
      for headerI in inpList:
        if "Peptides ST-1|Peptides ST-2|Peptides ST-3" == headerI:
          break
        else: spColCount += 1
      outF.write("ID,Protein ID, Gene name,") # write header into output file
      for headI in inpList[spColCount].split("|"):
        outF.write(headI + ",")
      for headI in inpList[spColCount + 1].split("|")[:-1]:
        outF.write(headI + ",")
      outF.write(inpList[spColCount + 1].split("|")[-1] + "\n")
      rowCount += 1
      continue
    
    outF.write(str(rowCount) + ",")
    
    inpLine = inpLine.strip()
    inpItems = inpLine.split("\t")
    inpName = max(inpItems[0].split("|"), key=len) # get unique protein ID
    inpGeneName = max(inpItems[6].split("|"), key=len) # and gene name
    outF.write(inpName + "," + inpGeneName + ",")
    
    inpSP = inpItems[spColCount].split("|") + inpItems[spColCount + 1].split("|") # get lfq intensity scores
    # print inpSP
    for lfqI in inpSP[:-1]: # write lfq values
      if lfqI == "_" or lfqI == "0":
        outF.write(str(random()) + ",") ################## try with log10 values this time
      else:
        try:
          outF.write(str(lfqI) + ",")
        except ValueError:
          print(inpItems)
          raise
    
    if inpSP[-1] == "_" or inpSP[-1] == "0": outF.write(str(random()) + "\n")
    else: outF.write(inpSP[-1] + "\n")
    
    
    rowCount += 1

def entry_parser():
  """remove duplicate protein name and total peptide count cell entries from a proteinGroups.txt dataset
  
  columns:
  1) Protein IDs 
  2) Majority protein IDs  
  3) Peptide counts (all)  
  4) Peptide counts (razor+unique)  
  5) Peptide counts (unique)  
  6) Protein names  
  7) Gene names  
  8) Fasta headers  
  9) Number of proteins  
  10) Peptides  
  11) Razor + unique peptides  
  12) Unique peptides  
  13) Peptides ST  
  14) Peptides TI  
  15) Razor + unique peptides ST  
  16) Razor + unique peptides TI  
  17) Unique peptides ST  
  18) Unique peptides TI  
  19) Sequence coverage [%]  
  20) Unique + razor sequence coverage [%]  
  21) Unique sequence coverage [%]  
  22) Mol. weight [kDa]  
  23) Sequence length
  24) Sequence lengths  
  25) Q-value  
  26) Score  
  27) Identification type ST  
  28) Identification type TI  
  29) Sequence coverage ST [%]  
  30) Sequence coverage TI [%]  
  31) Intensity  
  32) Intensity ST  
  33) Intensity TI  
  34) iBAQ  
  35) iBAQ ST  
  36) iBAQ TI  
  37) LFQ intensity ST  
  38) LFQ intensity TI  
  39) MS/MS Count ST  
  40) MS/MS Count TI  
  41) MS/MS Count  
  42) Only identified by site  
  43) Reverse
  44) Potential contaminant  
  45) id  
  46) Peptide IDs  
  47) Peptide is razor  
  48) Mod. peptide IDs  
  49) Evidence IDs  
  50) MS/MS IDs  
  51) Best MS/MS  
  52) Deamidation (NQ) site IDs  
  53) Oxidation (M) site IDs  
  54) Deamidation (NQ) site positions  
  55) Oxidation (M) site positions
  
  """
  # from tools import file_importer, file_outporter
  from copy import copy
  from collections import defaultdict
  import os.path
  
  print("this is entry parser")
  
  # inPathL = ["bob/processed/proteinGroups - OST-1-09042017.txt","bob/processed/proteinGroups_OST2.txt","bob/processed/proteinGroups_OST3.txt"]
  inpathL = []
  inpF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "cav1ko", "txt_cav1ko-1-17082017", "proteinGroups.txt"),"r")
  # outPath = "bob/processed/OST-24-05-2017_combined.csv"
  fileCount = 1
  # outF = file_outporter(outPath)
  outF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "cav1ko", "processed", "cav1ko-1.csv"),"w")
  # newFlag = True
  
  finDict = defaultdict(list)
  cN = 0
  # for relPath in inPathL:
  outDict = {}
  # inpF = file_importer(relPath)
  headerFlag = True
  
  for inpLine in inpF:
    cN += 1
    if headerFlag:
      headerFlag = False
      headerLine = inpLine
      continue
    inpLine = inpLine.strip("\n\r")
    inpItem = inpLine.split("\t")
    geneL = inpItem[0].split(";")
    lenS = len(geneL[0])
    curGene = geneL[0]
    for geneI in geneL: # find gene name with the shortest length
      if len(geneI) < lenS:
        lenS = len(geneI)
        curGene = geneI
    if "__" in curGene: continue # get rid of contaminant lines
    try: # get rid of wonky lines introduced by excel
      int(curGene)
      continue
    except ValueError: 
      pass

    if curGene[-2] == "-":
      curGene = curGene[:-2]
    if curGene[-3] == "-":
      curGene = curGene[:-3]
    
    # remove ambiguities based on gene name from the entire entry:
    
    corrPos = geneL.index(curGene)
    corrLine = []
    targetCount = 46 # after the 45th item row in the list, peptide IDs and modification start to appear which are allowed to have multiple entries and do not need to be disambiguated
    currCount = 1
    pepFlag = True
    for inpE in inpItem:
      currCount += 1
      if currCount == targetCount:
        pepFlag = False
        # print inpE
      if ";" in inpE and pepFlag:
        try:
          corrLine.append(inpE.split(";")[corrPos])
        except IndexError:
          corrLine.append(inpE.split(";")[0])
      else:
        corrLine.append(inpE.rstrip("\n"))

      
    if inpItem[6] == "":
      # print "no protein name found. adding the uniprot ID."
      inpItem[6] = curGene
          
    """
    try:
      for inpN in inpItem[4:10]:
        inpItem[inpItem.index(inpN)] = int(inpN)
      countFlag = True
    except ValueError:
      print inpItem[4:10]
      countFlag = False
    if countFlag:
      if sum(inpItem[4:10]) == 0: continue # there are some unexpressed proteins in there
      
    """
    # print len(corrLine)
    if curGene in outDict: # handle duplicate protein entries and merge them together
      # print "%s is duplicate" % curGene
      if curGene == "Protein IDs": 
        """
        quickCount2 = 0
        for quickDictI in outDict[curGene]:
          print str(quickCount2) + " " + quickDictI
          quickCount2 += 1
        quickList = inpItem
        quickCount3 = 0
        for quickImp in quickList:
          print str(quickCount3) + " " + quickImp
          quickCount3 += 1         
        # print inpItem
        # print outDict[curGene]
        """
        continue
      combList = []
      
      """
      addL = []
      for i in outDict[curGene][3:]:
        addL.append(i)
      addL2 = []
      for j in corrLine[3:]:
        addL2.append(i)
      outL[3:] = map(add, addL, addL2) # admittedly this looks terrible
      """
      
      indexN = 0
      for cItem in corrLine:
        # print indexN
        # print "---"
        # print len(corrLine)
        if indexN < 18 or 30 <= indexN <= 43:
          try:
            currC = int(cItem)
            currC = currC + int(outDict[curGene][indexN]) # numbers like peptide counts or LFQ values are added up during merge
          except ValueError:
            currC = cItem
        
        elif 18 <= indexN <= 25 or 28 <= indexN <= 29: # sequence coverage and scores
          currC = max([float(cItem),float(outDict[curGene][indexN])])
        
        elif 26 <= indexN <= 27 or indexN == 44:
          """
          quickCount = 0
          for corrItem in corrLine:
            print str(quickCount) + " " + corrItem
            quickCount += 1
            
          import time
          
          print relPath
          print corrLine
          print outDict[curGene]
          print "++++++++++++++++++++++++"
          print indexN
          time.sleep(0.5)"""
          currC = cItem

          
        else:
          corrL = cItem.split(";")
          # print indexN
          # print corrLine
          # print outDict[curGene][indexN]
          dictL = outDict[curGene][indexN].split(";")
          mergeL = copy(dictL)
          for corrI in corrL:
            if corrI not in dictL:
              mergeL.append(corrI)
          
          currC = ";".join(mergeL)

        combList.append(currC)

        
        indexN +=1
      
      
      combList[-1] = "merged"    
      outDict[curGene] = combList 
      # print "merged:"
      # print combList
    else:
      corrLine.append("unique")
      outDict[curGene] = corrLine

    
  print(fileCount)
  

  #   if not newFlag: print fileCount, testKey, finDict[testKey]     
  # if newFlag:
  #   newFlag = False
  
  for outKey,outValue in list(outDict.items()): 
    if outKey in finDict: # add modified dicts together into single, unified dict
      # print fileCount, finDict[outKey]
      # print outValue
      outIndex = 0
      for outItem in outValue:
        finDict[outKey][outIndex].append(outItem)
        outIndex += 1
      # print finDict[outKey]

    else:  # or just add new entries
      if fileCount == 1:
        for outItem in outValue:
          finDict[outKey].append([outItem])
      
      else: # fill up entries that were not present in the previous cycle
        loopCount = 0
        while loopCount < fileCount - 1:
          for i in range(len(outValue)):
            if len(finDict[outKey]) == i:
              finDict[outKey].append([])
            else:
              finDict[outKey][i].append("")
          loopCount += 1
        outIndex = 0
        for outItem in outValue:
          # print finDict[outKey]
          finDict[outKey][outIndex].append(outItem)          
          outIndex += 1

  for testKey in finDict: # fill up entries in result dict which were not present in previous file
    if len(finDict[testKey][0]) < fileCount:
      for i in range(len(finDict[testKey])):
        finDict[testKey][i].append("")

  if len(inpathL) > 1: fileCount += 1 # this is needed if multiple files are parsed
  for finK, finV in list(finDict.items()):
    for finI in finV[-1]:
      if finI != "unique" and finI != "":
        print(finK, finV)

    
  
  outN = 0  
  # prepare header for file:
  headList = headerLine.strip("\n\r").split("\t")
  if fileCount > 1:
    for headerItem in headList[:-1]:
      headerI = headerItem.replace(",",".")
      headerCount = 1
      while headerCount < fileCount:
        outF.write(headerI + "-" + str(headerCount) + "|")
        headerCount += 1  
      outF.write(headerI + "-" + str(headerCount) + "\t")
      
    headerCount = 1
    while headerCount < fileCount:
      outF.write(headList[-1] + "-" + str(headerCount) + "|")
      headerCount += 1
    
    outF.write(headList[-1] + "-" + str(headerCount) + "\n")

  elif fileCount == 1:
    for headerItem in headList[:-1]:
      headerI = headerItem.replace(",",".")    
      outF.write(headerI + "\t")
    outF.write(headList[-1].replace(",",".") + "\n")
  
  else:
    print("number of input files should be at least one. Got less somehow")
    raise ValueError
    
  
  for outDK, outDV in list(finDict.items()): # write out assembled results to a file
    outN += 1
    if len(outDK) > 30: print("this line should not be displayed")
    # print outDV[1]
    # if outN == 100: break
    nameCount = 0
    for outI in outDV:
      # if nameCount == 0: print outI
      for outPiece in outI[:-1]:
        outU = outPiece.replace(",",".")
        if outU == "": outF.write("_|")
        else: outF.write(str(outU) + "|")
      if outI[-1] == "": # handle missing entries
        if nameCount == 6: outF.write(outDV[0][0] + "\t") # replace missing gene names with their uniprot ID
        else: outF.write("_\t")
      else: outF.write(str(outI[-1]).replace(",",".") + "\t")
      nameCount += 1
    outF.write("\n")
  

  print("unique proteins: ", outN)
  print("lines parsed: ", cN)
  # print headerLine
  inpF.close()
  outF.close()
  
def crapome_parser():
  """take crapome 1.1 tab delimited file as input and add it to the output of the lfq_parser method in proteingroups_parser.py.
  proteins not present in the database shall also be included."""
  import os.path
  
  # contTreshold = 30 # set this to the desired contamination score
  resD = {}
  
  # crapFile = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "cav1ko", "processed", "1503486016360_gp-1.txt"),"rU")
  crapFile = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "Crapome-OT1-upregulated.txt"),"rU")
  outF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "Crapome-OT1-upregulated_parsed.txt"),"w")
  
  headerFlag = True
  
  fileLength = 0
  for inpLine in crapFile: # parse crapome output file
    if headerFlag:
      headerFlag = False
      outF.write("protName\tscore\n")
      continue
    fileLength += 1
    lineList = inpLine.split("\t")
    if lineList[2] == "": continue
    elif len(lineList) > 2: contScore = int(lineList[2].split("/")[0])
    else: contScore = 0
    
    outF.write(lineList[1]+"\t"+ str(contScore) + "\n")
    
    # if contScore < contTreshold:
    resD[lineList[0]] = contScore
  
  # print "Contaminant treshold: " + str(contTreshold)
  
  print("lines parsed: " + str(fileLength))
  print("Number of results: " + str(len(resD)))
  
      
#   inpFile = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "cav1ko", "processed", "cav1ko-1_no0.csv"),"r")
#   outF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "cav1ko", "processed", "cav1ko-1_no0_crapome.csv"),"w")
#   
#   headerFlag = True
#   for inpLine in inpFile: # parse the input file for crapome and add crapome results to it
#     inpList = inpLine.rstrip("\n").split(",")
#     for inpI in inpList:
#       outF.write(inpI + ",")
#     
#     if headerFlag: 
#       outF.write("Crapome score")
#       headerFlag = False
#     elif inpList[2].upper() in resD: outF.write(str(resD[inpList[2].upper()]))
#     else: outF.write("0")
#     
#     outF.write("\n")
#   print "results written to file"

  
if __name__ == "__main__":
  main()