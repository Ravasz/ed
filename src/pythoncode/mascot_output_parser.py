'''
Created on 30 May 2019

@author: mate

take mascot output csv files from Kyriaki and assemble them to a single file while performing some analysis on them.

Use Pandas for this where possible
'''

##########################################################
# set variables here:
fileNames = ["/home/mate/code/ed/src/data/kyriaki/F013719_KN_WT-1.csv", "/home/mate/code/ed/src/data/kyriaki/F013712_KN_WT-2.csv", "/home/mate/code/ed/src/data/kyriaki/F013713_KN_WT-3.csv"]
fileNames.extend(["/home/mate/code/ed/src/data/kyriaki/F013716_KN_WAGO-1.csv", "/home/mate/code/ed/src/data/kyriaki/F013742_KN_WAGO-2.csv", "/home/mate/code/ed/src/data/kyriaki/F013743_KN_WAGO-3.csv"])

crapome_file = "152632541280_gp_ptpn22_r619w_FC_16052018.txt"

group1KeyWord = "WT"
group2KeyWord = "WAGO"

outputFileName = "/home/mate/code/ed/src/data/kyriaki/WAGO-IP.csv"

##########################################################


def main_function(fileNameL, outputFileS, g1K, g2K):
  """call functions in the correct order:
  - combine separate files into single dataframe
  - calculate fold changes and asign rank
  - calculaute p values
  - add crapome score - not done yet
  - write results to a file
  """
  
  import numpy as np
  from scipy.stats import ttest_ind
  
  
  # assemble full dataframe
  print("processing files", end = "")
  dataDF = dataframe_merger(fileNameL)
  print("Done.")
  
  wtList = []
  otherList = []
  for columnI in list(dataDF):
    if str(columnI).startswith("prot_score-" + g1K): wtList.append(columnI)
    elif str(columnI).startswith("prot_score-" + g2K): otherList.append(columnI)
  
  # calculate fold changes:
  dataDF["avg_" + g1K] = dataDF[wtList].mean(axis = 1, skipna = True).fillna(1.0)
  dataDF["avg_" + g2K] = dataDF[otherList].mean(axis = 1, skipna = True).fillna(1.0)
  dataDF["FC"] = dataDF["avg_" + g1K] / dataDF["avg_" + g2K]
  dataDF["FC"].replace(np.inf, 100, inplace = True)
  dataDF["FC"].replace(0, 0.01, inplace = True)
  
  # assign rank
  dataDF["rank"] = dataDF["FC"].rank(method = "min")
  # dataDF = dataDF.sort_values(by=["FC"], ascending = False)
  # dataDF["rank"] = range(1,len(dataDF)+1)
  
  # add P values
  dataDF["p_value"] = ttest_ind(dataDF[wtList],dataDF[otherList], axis = 1, nan_policy = "omit" , equal_var = False)[1]
  dataDF["p_value"] = dataDF["p_value"].fillna(1.0)
  
  # add CRAPOME score - not done yet
#   crapDF = crapome_merger()
#   print(crapDF[["frequency"]])
#   
#   dataDF = pd.merge(dataDF, crapDF[["frequency"]], how = "left", left_index = True, right_index = True)
  

  
  # write out results
  dataDF = dataDF.fillna(0.0)
  print(dataDF)
  dataDF.to_csv(outputFileS, index = False)


  
def crapome_merger():
  """"merge converted protein ID dataframe with crapome dataframe"""
  import pandas as pd
  crapomeDF = crapome_parser("/home/mate/code/ed/src/data/kyriaki/1559312746821_userCrapDB.txt")
  
  convertDF = crapomeDF[["refseq_unique_id"]].copy()
  convertDF.index = convertDF["refseq_unique_id"]
  
  uniprotDF = id_converter(convertDF)
  
  # print(crapomeDF)
  # print(uniprotDF)
  crapomeUniDF = pd.merge(crapomeDF, uniprotDF, how = "left", left_on='refseq_unique_id', right_on='refseq_id')  
  crapomeUniDF = crapomeUniDF.set_index("uniprot_id")
  crapomeUniDF = crapomeUniDF.drop_duplicates(keep = "first")
  crapomeUniDF = crapomeUniDF[crapomeUniDF.index.notnull()]
  print(crapomeUniDF)
  
  # print(crapomeUniDF)
  
  return crapomeUniDF



def id_converter(protDF, orgnID = "9606", inpDB = "refseqproteinaccession", outDB = "uniprotaccession"):
  """take in a list of uniprot entry names 
  and convert them to protein IDs using biodbnet.
  return a list of protein IDs. 
  bioDBnet often takes time to load so this might take several minutes to complete.
 
  """
  
  import urllib.request, json
  import sys
  import pandas as pd
  
  fullJason = []
  cutList = []
  countNum = 0
  
  protList = list(protDF["refseq_unique_id"])
  
  print("\nconnecting to biodbnet. This might take a while", end = "")
  for protI in protList:
    cutList.append(protI)
    countNum += 1
    if countNum == 190:
      urlStr = "http://biodbnet.abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&format=row&input=" + inpDB + "&inputValues=" + ",".join(cutList) + "&outputs=" + outDB + "&taxonId=" + orgnID  
      uParsed = urllib.request.urlopen(urlStr)  

      responseJson = uParsed.read()
      cutList = []
      countNum = 0
      # print(responseJson)
      print(".", end="")
      parsedJson = json.loads(responseJson.decode('utf-8'))
      fullJason.extend(parsedJson)
  
  if countNum > 0: # parse the last, smaller batch
      urlStr = "http://biodbnet.abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&format=row&input=" + inpDB + "&inputValues=" + ",".join(cutList) + "&outputs=" + outDB + "&taxonId=" + orgnID  
      uParsed = urllib.request.urlopen(urlStr)  

      responseJson = uParsed.read()
      cutList = []
      countNum = 0
      # print(responseJson)
      print(".", end="")
      parsedJson = json.loads(responseJson.decode('utf-8'))
      fullJason.extend(parsedJson)    
  
  print("Done.")

  if len(fullJason) == 0:
    print("No result returned. Site may be down or query may be incomplete")
    sys.exit(0)    
  
  resDF = pd.DataFrame(columns = ["refseq_id","uniprot_id"])
  notFoundCount = 0
  
  for inpAccD in fullJason:
    queryVL = inpAccD["UniProt Accession"].split("//")
    queryK = inpAccD["InputValue"]# .encode("ascii","ignore")
    foundFlag = False
    for queryI in queryVL:
      if queryI == "" or queryI == "-": continue
      if "-" in queryI: queryI = queryI.split("-")[0]
      resDF.loc[len(resDF)] = [queryK, queryI]
      foundFlag = True
      break
    
    if not foundFlag: 
      notFoundCount += 1
#       print("GI not found in this query:")
#       print(inpAccD)
  
  resDF.set_index("refseq_id")
  print("IDs not converted: " + str(notFoundCount))
  
  return resDF  
  
  
def crapome_parser(crapFile):
  """read in a crapome 1.1 database file and turn it into a dataframe"""
  
  def refseq_splitter(refSeries):
    """mini function: take in string from refseq cell in crapome input file and return a single, searchable refseq identifier"""
    
    refStr = str(refSeries)
    refL = refStr.split(";")
    minI = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    for itemI in refL: 
      if len(itemI) < len(minI): minI = itemI
    
    minI = minI.split(".")[0]
    
    return minI
  
  
  
  import pandas as pd
  pd.set_option('display.max_columns', 500)
  # pd.set_option('display.max_rows', 100)
  pd.set_option('display.width', 1000)
  
  # read in crapome file
  crapDF = pd.read_csv(crapFile, sep = "\t")
  crapDF[crapDF.columns[4:]] = crapDF[crapDF.columns[4:]].astype(int, errors = "ignore")
  crapDF[crapDF.columns[3]] = crapDF[crapDF.columns[3]].astype(float, errors = "ignore")
  
  # calculate frequencies
  crapDF["frequency"] = (crapDF[crapDF.columns[5:]] > 0).sum(1) / len(crapDF.columns[5:])
  
  # convert IDs to uniprot
  crapDF["refseq_unique_id"] = crapDF["REFSEQ_ID"].apply(refseq_splitter)
  # print(crapDF)  
  return crapDF
  



  
  
def dataframe_merger(fileList):
  """merge dataframes together from a list of filenames. each file is turned into a DF by calling csv_organizer, and merged into a the complete DF. 
  Returns the full dataframe"""
  
  import pandas as pd
  
  fullDF = pd.DataFrame()

  for fName in fileList:
    currDF = csv_organizer(fName)
    if len(fullDF) == 0: fullDF = currDF
    else: 
      fullDF = pd.merge(fullDF, currDF[currDF.columns[4:]], how = "outer", left_index=True, right_index=True) # add in new measurements with a full outer merge. Add new rows for new proteins, append existing for present proteins
      fullDF = fullDF.combine_first(currDF) # fill in descriptions by replacing nan values in the first DF with an non-nan value in the second DF
  
  return fullDF
    

def csv_organizer(fileName):
  """take out the protein data from the mascot result file. return a dataframe"""
  import pandas as pd
  pd.set_option('display.max_columns', 500)
  # pd.set_option('display.max_rows', 100)
  pd.set_option('display.width', 1000)
  
  fID = fileName.split("_")[-1].split(".")[0]
  
  with open(fileName,"r") as inFile:
    print(".", end = "")
    for inLine in inFile:
      if inLine.startswith("prot_hit_num"): # loop until the actual data starts and skip over metadata
        headerLRaw = inLine.rstrip("\n ").split(",")
        break
    
    headerL = headerLRaw[:4]
    for headerI in headerLRaw[4:]:
      headerL.append(headerI + "-" + fID)
    
    firstLine = next(inFile).rstrip("\n ").replace("\"","").split(",")[:-1]

    inDF = pd.DataFrame(firstLine).T # init dataframe
    inDF.columns = headerL
    inDF.index = inDF[inDF.columns[2]]
    
    for inLine in inFile: 
      inList = inLine.strip("\n ").replace("\"","").split(",")[:-1]
      
      strCount = 0
      for inI in inList: # handle commas in the gene description (which is acutally some fasta header)
        if inI == "":continue
        try:
          int(inI)
        except ValueError:
          strCount += 1
        
      if strCount > 2: # this means there was at least one comma in the header
        
        inList = inList[:3] + ["".join(inList[3:strCount + 2])] + inList[strCount + 2:]

      if len(inList) > 2: inDF.loc[inList[2]] = inList # add nonempty lines to DF
      
  inDF = inDF.apply(pd.to_numeric, errors = "ignore") # covert strings to numbers where appropriate
  
  # print(inDF)
  return inDF
      
# crapome_parser("/home/mate/code/ed/src/data/kyriaki/1559312746821_userCrapDB.txt")

main_function(fileNames, outputFileName, group1KeyWord, group2KeyWord)
