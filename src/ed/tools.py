'''
Created on 15 May 2015

@author: mate
bunch of functions useful for writing scripts with proteins. No action script.
'''
def main():
  print "you are in protein tools"
  print "you should really call functions from outside this script though"
  
  # testL = ["NP_055136"]
  
  # print prot_id_converter(testL, inpDB = "refseqproteinaccession", outDB = "genesymbol", orgnID="9606")  
  
  intact_parser()


def file_reader(fileStr,resType="dict"):
  """read in a file specified by fileStr (full path) and return its contents as a resType variable, 
  preferably a dict. Hmm... there is apparently a builtin csv reader. How about that.
  This needs further working on despite the simple problem to achieve the optimal solution"""
  
  import csv
  
  if resType == "dict":
    resD = {}
    with open(fileStr,"rU") as inpF:
      for inpLine in inpF:
        inpList = inpLine.strip("\n\r").split(",")
  


def ROutputFormatter():
  """take a terrible output file from R and format it it in a more nice way, 
  like remove leftover spaces and commas in it then add fold change and FDR score"""
  from math import log
  
  fdrN = 0.05
  def p_value_key(protItem):
    """mini function returning the last element of a list. just because I do not like unnamed functions"""
    return protItem[-1]
  
  protList = []
  headerFlag = True
  missCount = 1
  with open("bob/processed/OST-24-05-2017_combined_ttest.csv", "r") as inpF: # read and process the csv with protein names and p values
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
      if curLine[2] == "": # if no gene name is given, just add a placeholder
        curLine[2] = "Noname" + str(missCount)
        missCount += 1

      protList.append(curLine)
      
  protList.sort(key = p_value_key) # sort the whole list on p value (lowest to highest) 
  i = 0.0 
  m = float(len(protList))
  print "dataset length: ", int(m)
  with open("bob/processed/OST-24-05-2017_combined_ttest_formatted.csv","w") as outputF:
    outputF.write("#,UniprotID,geneName,OST1,OST2,OST3,WT1,WT2,WT3,pValue,FDR,log2FoldChange\n")
    for protListI in protList:
      i += 1
      critVal = (i/m)*fdrN # this is the benjamini-hochberg defined critical value
      protListI.append(critVal)
      try:
        FAvg = (protListI[2] + protListI[3] + protListI[4])/3.0 # OST
        SAvg = (protListI[5] + protListI[6] + protListI[7])/3.0 # OT1
      except TypeError:
        print curLine
        raise
      logFoldChange = log(SAvg/FAvg,2) # so positive numbers are more abundant in the OT1 cells, negatives number in the OST cells, at least for the OST IP mass spec file
      protListI.append(logFoldChange)
      
      for outI in protListI:
        outputF.write(str(outI))
        if outI is protListI[-1]:
          outputF.write("\n")
        else:
          outputF.write(",")
      
  print "formatting complete"
      
    

def protein_quantifier(protName = "Ptpn22"):
  """look up a protein by name in Bob's dataset, and return its abundance in ppm
  the ppm is calculated from the lfq values. Apparently, single lfq value, divided by its column total results in ppm.   
  """
  print "running protein quantifier"
  with open("bob/proteinGroups.txt", "r") as inpF:
    lfqTotal = [0.0,0.0,0.0] # this holds the total of each column
    headerFlag = True
    fullD = {} # complete list of names, abundances and even a rank
    for inpLine in inpF: 
      inpListed = inpLine.split("\r")
      for inpListedI in inpListed: # iterate through file
        if headerFlag:
          headerFlag = False
          continue
        inpL = inpListedI.split("\t")
        colCount = 0
        protNameList = inpL[6].split(";")
        if protName in protNameList: # find protein of interest
          print protName, "found as:", inpL[6]
        fullD[protNameList[0]] = []
        for inpI in inpL[83:86]: # do running sum
          if inpI != "NaN":
            lfqTotal[colCount] += float(inpI)
            fullD[protNameList[0]].append(float(inpI))
          colCount += 1
    for keyS, valueL in fullD.items(): # calculate mean ppm
      ppmSum = 0
      ppmCount = 0
      for valueI in valueL:
        ppmSum += ((valueI/lfqTotal[ppmCount])*1000000 )
        ppmCount += 1
      fullD[keyS].append(ppmSum/3)
    sortedFullD = sorted(fullD.items(), key=lambda (k, v): v[-1], reverse=True) # sort
    print lfqTotal
    for i in range(len(sortedFullD)): # add ranks
      sortedFullD[i][1].append(i+1)
    resD = {}
    for sortedI in sortedFullD:
      resD[sortedI[0]] = sortedI[1]
    
    if protName in resD:
      print protName
      print "average ppm in WT samples:", resD[protName][-2]
      print "abundance rank:", resD[protName][-1]
      
    else: 
      print protName, "does not match to any entry"

    with open("protein_abundances-ibaq.txt","w") as outF:
      for resI in sortedFullD:
        outF.write(resI[0] + "," + str(resI[1:]).strip("()[],") + "\n")
        
        

def file_importer(relPath, methodS = "r"):
  """open a file specified by the relPath variable, 
  even if its in another directory"""
  import os
  scriptDir = os.path.dirname(__file__) # absolute dir this script is in
  absFilePath = os.path.join(scriptDir, relPath)
  inpF = open(absFilePath, methodS)
  return inpF

def file_outporter(outPath):
  """create a file to write into even if its in another directory.
  deprecated. use file_importer and define the methodS variable as "w" instead."""
  import os
  scriptDir = os.path.dirname(__file__) # absolute dir this script is in
  absFilePath = os.path.join(scriptDir, outPath)
  outF = open(absFilePath,"w")
  return outF


def iso_e(protS):
  """return the isoelectric point of protS string protein sequence"""
  from Bio.SeqUtils.ProtParam import ProteinAnalysis
  protA = ProteinAnalysis(protS)
  return (protA.isoelectric_point())

def stop_watch():
  """this should not be called, its just a stopwatch 
  that can be copied into each script to measure how fast it runs"""
  import timeit
  start = timeit.default_timer()
  stop = timeit.default_timer()
  print stop - start 
  
def uniprot_dicter():
  """create dict of all known uniprot mouse proteins and return dict 
  with keys as uniprot IDs and values as protein sequences"""
  
  relPath = "datafiles/mouse_proteome.fasta"
  
  inpF = file_importer(relPath)
  geneD = {}
  
  for inpLine in inpF:
    if inpLine[0] == ">":
      geneKey = inpLine.split("|")[1]
    else:
      geneD[geneKey] = inpLine
  return geneD

def protein_groups_slicer():
  """open the large protein groups.txt file and take out parts of it into a new file for further analysis
  """
  with open("bob/proteinGroups.txt","r") as inpF, open("bob/bob_lfq.csv", "w") as outF:  
    for inpLine in inpF: # this unfortunately reads in the whole file as it is split by \r characters
      inpLine = inpLine.split("\r")
      cN = 0
      for inpL in inpLine: # read file line by line
        inpL = inpL.split("\t")
        commaN = 0
        for inpI in inpL[86:92]: # this loop is the important part. It writes some columns of the large proteinGorups.txt file into a new file
          if inpI == "NaN": inpI = 0
          if commaN == 5: outF.write(str(inpI)) # write commas between items, but not at the end of the line
          else: outF.write(str(inpI) + ",")
          commaN += 1
        outF.write("\n")
        cN += 1
    print str(cN) + " lines written to " + str(outF.name) 


def prot_id_converter(protList, orgnID = "10090", inpDB = "uniprotaccession", outDB = "genbankproteingi"):
  """take in a list of uniprot entry names 
  and convert them to protein IDs using biodbnet.
  return a list of protein IDs. 
  bioDBnet often takes time to load so this might take several minutes to complete.
  
  For orgnID use 10090 for mouse and 9606 for human.
  
  Human gene symbols are all caps: VAV1
  Mouse gene symbols are first caps: Vav1
  
  inpDB is input database type, use "uniprotaccession" for uniprot entry IDs or "genesymbol" for gene symbols
  outDB is output database type, use "geneid" for Gene ID or "genbankproteinaccession" or "genbankproteingi" or "refseqproteingi" for those
 
  """
  import urllib, json
  urlStr = "http://biodbnet.abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&format=row&input=" + inpDB + "&inputValues=" + ",".join(protList) + "&outputs=" + outDB + "&taxonId=" + orgnID  
  print "connecting to biodbnet. This might take a while..."
  uParsed = urllib.urlopen(urlStr)  
  print "connection successful"
  responseJson = uParsed.read()
  parsedJson = json.loads(responseJson)
  print parsedJson
  # parsedJson = [{u'Gene ID': u'54196', u'InputValue': u'Q8CCS6'}, {u'Gene ID': u'99982', u'InputValue': u'Q6ZQ88'}]
  # parsedJson = [{u'GenBank Protein Accession': u'BAC27741//Q8CCS6//EDL36322//AAH55866//NP_062275//XP_006519335//AAC00210////EDL36323', u'InputValue': u'Q8CCS6'}, {u'GenBank Protein Accession': u'AAH19417//XP_006539394//XP_006539393//NP_598633//AAH59885//CBY79415//CBY88367////XP_006539392//EDL29935//Q6ZQ88//BAC97980', u'InputValue': u'Q6ZQ88'}]
  # parsedJson = [{u'GenBank Protein GI': u'9506945//148704376//46396417//26328001//33585929//2351846////148704375//568988212', u'InputValue': u'NP_062275'}, {u'GenBank Protein GI': u'51315882//18044445//224994233////568932208//317440660//37589595//568932212//315003691//148697988//37360004//568932210', u'InputValue': u'NP_598633'}]
  # parsedJson = [{u'RefSeq Protein GI': u'//6005942', u'InputValue': u'P55072'}, {u'RefSeq Protein GI': u'530368795////46488944//767919614//578804849//31455611', u'InputValue': u'P43403'}, {u'RefSeq Protein GI': u'7108367//384551646//768003854//530425424//384551649', u'InputValue': u'P15498'}, {u'RefSeq Protein GI': u'767904317//112789546//767904319////112789548//767904315', u'InputValue': u'P06239'}, {u'RefSeq Protein GI': u'4502671//', u'InputValue': u'P07766'}, {u'RefSeq Protein GI': u'7108367//384551646//768003854//530425424//384551649', u'InputValue': u'P15498'}, {u'RefSeq Protein GI': u'767904317//112789546//767904319////112789548//767904315', u'InputValue': u'P06239'}, {u'RefSeq Protein GI': u'767910875//37595565//4557431////767910873', u'InputValue': u'P20963'}, {u'RefSeq Protein GI': u'4502671//', u'InputValue': u'P07766'}, {u'RefSeq Protein GI': u'530368795////46488944//767919614//578804849//31455611', u'InputValue': u'P43403'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'767910875//37595565//4557431////767910873', u'InputValue': u'P20963'}, {u'RefSeq Protein GI': u'767910875//37595565//4557431////767910873', u'InputValue': u'P20963'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'-', u'InputValue': u'P32577'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'//768033853//4507909', u'InputValue': u'P42768'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'767985679//767985664////767985659//767985662//767985670//20149528//767985683//578827539//767985657//767985674//767985677//767985681//767985668', u'InputValue': u'O43586'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'767985679//767985664////767985659//767985662//767985670//20149528//767985683//578827539//767985657//767985674//767985677//767985681//767985668', u'InputValue': u'O43586'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'//29725609//41327734//41327732//41327736', u'InputValue': u'P00533'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'744066863//16753212', u'InputValue': u'O75563'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'-', u'InputValue': u'EBI-9974954'}, {u'RefSeq Protein GI': u'744066863//16753212', u'InputValue': u'O75563'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}, {u'RefSeq Protein GI': u'-', u'InputValue': u'Q3UND0'}, {u'RefSeq Protein GI': u'-', u'InputValue': u'P97814'}, {u'RefSeq Protein GI': u'-', u'InputValue': u'Q99M15'}, {u'RefSeq Protein GI': u'530362357//301171669//301171662//767903702//767903706//224586929//767903708//767903704//767903700////815891380', u'InputValue': u'Q9Y2R2'}]
  # parsedJson = [{u'UniProt Accession': u'Q3TIM2//Q3TFH9//Q3TXN9//Q01853//Q6PI18//Q8BSR6//Q8CEG4', u'InputValue': u'Vcp'}, {u'UniProt Accession': u'P43404//P97455//Q80VV2//Q8CHJ3', u'InputValue': u'Zap70'}, {u'UniProt Accession': u'P27870//Q8BTV7', u'InputValue': u'Vav1'}, {u'UniProt Accession': u'Q91X65//Q62320//E9Q696//Q61794//P06240//Q61795', u'InputValue': u'Lck'}, {u'UniProt Accession': u'P22646', u'InputValue': u'Cd3e'}, {u'UniProt Accession': u'Q9D3G3//P29020//P24161', u'InputValue': u'Cd247'}, {u'UniProt Accession': u'Q03143//Q8VCW1//P41241//Q80WU4', u'InputValue': u'Csk'}, {u'UniProt Accession': u'Q64260//Q6GU11//P15116', u'InputValue': u'Cdh2'}, {u'UniProt Accession': u'P70424//Q61525//Q6ZPE0', u'InputValue': u'Erbb2'}, {u'UniProt Accession': u'E9QPE2//P05622', u'InputValue': u'Pdgfrb'}, {u'UniProt Accession': u'Q3UFB7', u'InputValue': u'Ntrk1'}, {u'UniProt Accession': u'A6H6U3//Q9Z2A0//Q9R1D8//Q9R215', u'InputValue': u'Pdpk1'}, {u'UniProt Accession': u'P16590//Q920Z3//Q9R264//P16882//Q61653//Q6DI66//Q80W86//Q8R1M5', u'InputValue': u'Ghr'}, {u'UniProt Accession': u'P70315', u'InputValue': u'Was'}, {u'UniProt Accession': u'Q60631//Q61240', u'InputValue': u'Grb2'}, {u'UniProt Accession': u'Q9EP98//Q01279', u'InputValue': u'Egfr'}, {u'UniProt Accession': u'Q8BK74//Q3UND0//Q9Z2K4', u'InputValue': u'Skap2'}, {u'UniProt Accession': u'Q4V9R4//P97814', u'InputValue': u'Pstpip1'}, {u'UniProt Accession': u'Q99M15//Q6GTF6//Q9Z189', u'InputValue': u'Pstpip2'}, {u'UniProt Accession': u'P08032//B2RWX6//P97502', u'InputValue': u'Spta1'}, {u'UniProt Accession': u'Q3U527//Q8CEA1//P22682', u'InputValue': u'Cbl'}]
  # parsedJson = [{u'KEGG Gene ID': u'mmu:16656', u'InputValue': u'A2A884'}, {u'KEGG Gene ID': u'mmu:11737', u'InputValue': u'O35381'}, {u'KEGG Gene ID': u'mmu:19763', u'InputValue': u'O35730'}]
  # parsedJson = [{u'InputValue': u'Q5SWU9', u'Gene Symbol': u'Acaca'}, {u'InputValue': u'Q8VDD5', u'Gene Symbol': u'Myh9'}, {u'InputValue': u'Q3T9S7', u'Gene Symbol': u'Pcx'}, {u'InputValue': u'B2RRX1', u'Gene Symbol': u'Actb'}, {u'InputValue': u'Q71LX8', u'Gene Symbol': u'Hsp90ab1'}, {u'InputValue': u'B2RRE2', u'Gene Symbol': u'Myo18a'}, {u'InputValue': u'Q3U2W2', u'Gene Symbol': u'Mybbp1a'}, {u'InputValue': u'Q3TII3', u'Gene Symbol': u'Eef1a1'}, {u'InputValue': u'P99024', u'Gene Symbol': u'Tubb5'}, {u'InputValue': u'E9QAS3', u'Gene Symbol': u'Ptpn22'}, {u'InputValue': u'Q99MR8', u'Gene Symbol': u'Mccc1'}, {u'InputValue': u'Q3THE2', u'Gene Symbol': u'Myl12b'}, {u'InputValue': u'D3YZ62', u'Gene Symbol': u'Myo5a'}, {u'InputValue': u'Q3UGC8', u'Gene Symbol': u'Pcca'}, {u'InputValue': u'Q6S385', u'Gene Symbol': u'Plec'}, {u'InputValue': u'B2RTP7', u'Gene Symbol': u'Krt2'}, {u'InputValue': u'B1AQ77', u'Gene Symbol': u'Krt15'}, {u'InputValue': u'D3Z6I8', u'Gene Symbol': u'Tpm3'}, {u'InputValue': u'B2RTM0', u'Gene Symbol': u'Hist2h4'}, {u'InputValue': u'Q8K0Z5', u'Gene Symbol': u'Tpm3'}, {u'InputValue': u'Q3TNH0', u'Gene Symbol': u'Tmpo'}, {u'InputValue': u'Q3TIG9', u'Gene Symbol': u'Myl6'}, {u'InputValue': u'D2KHZ9', u'Gene Symbol': u'GAPDH'}, {u'InputValue': u'Q6P5D8', u'Gene Symbol': u'Smchd1'}, {u'InputValue': u'Q4FZG4', u'Gene Symbol': u'Flna'}, {u'InputValue': u'F1DGF6', u'Gene Symbol': u'Prkcd'}, {u'InputValue': u'Q3TFG3', u'Gene Symbol': u'Eif4a1'}, {u'InputValue': u'B2RPX1', u'Gene Symbol': u'Iqcd'}, {u'InputValue': u'Q8BQ35', u'Gene Symbol': u'Sptbn1'}, {u'InputValue': u'E0CZ27', u'Gene Symbol': u'H3f3a'}, {u'InputValue': u'Q9CR57', u'Gene Symbol': u'Rpl14'}, {u'InputValue': u'Q0VG47', u'Gene Symbol': u'Hnrnpa3'}, {u'InputValue': u'Q8C553', u'Gene Symbol': u'Lmnb1'}, {u'InputValue': u'Q3T9U9', u'Gene Symbol': u'Rpl3'}, {u'InputValue': u'Q3KQJ4', u'Gene Symbol': u'Hspa8'}, {u'InputValue': u'Q3U7D2', u'Gene Symbol': u'Rpl15'}, {u'InputValue': u'A0PJE6', u'Gene Symbol': u'Pccb'}, {u'InputValue': u'Q68FG3', u'Gene Symbol': u'Spty2d1'}, {u'InputValue': u'Q0VB76', u'Gene Symbol': u'Gzmc'}, {u'InputValue': u'Q32P04', u'Gene Symbol': u'Krt5'}, {u'InputValue': u'D3Z6F5', u'Gene Symbol': u'Atp5a1'}, {u'InputValue': u'Q3U0I3', u'Gene Symbol': u'Cct3'}, {u'InputValue': u'Q3TJZ1', u'Gene Symbol': u'Eef2'}, {u'InputValue': u'Q3UI57', u'Gene Symbol': u'Mcm3'}]

  if "GenBank Protein Accession" in parsedJson[0]:
    return accession_wrangler(parsedJson)
  elif "GenBank Protein GI" in parsedJson[0]:
    return protein_gi_wrangler(parsedJson)
  elif "RefSeq Protein GI" in parsedJson[0]:
    return refseq_gi_wrangler(parsedJson)
  elif "UniProt Accession" in parsedJson[0]:
    return uniprot_wrangler(parsedJson)
  elif "KEGG Gene ID" in parsedJson[0]:
    return kegg_wrangler(parsedJson)
  elif "Gene Symbol" in parsedJson[0]:
    return gene_symbol_wrangler(parsedJson)
                         
  else: 
    print "was expecting Genbank protein accessions, genbank protein GIs, Uniprot accessions or RefSeq protein GIs, but got something else:"
    print parsedJson
    raise ValueError
 
def kegg_wrangler(inpAcc):
  """to be called by prot_id_converter
  Take in a loaded json file which has kegg gene IDs as results.
  Return a list of lists with where each sublist contains the results of a single query"""  
  print "processing KEGG Gene IDs Accessions"
  resL = []
  for inpAccD in inpAcc:
    queryL = inpAccD["KEGG Gene ID"].split("//")
    foundFlag = False
    curL = []
    for queryI in queryL:
      if queryI == "" or queryI == "-": continue
      curQuery = queryI.encode("ascii","ignore")
      curL.append(curQuery)
      foundFlag = True
    if not foundFlag: 
      print "GI not found in this query:"
      print inpAccD
    else:
      resL.append(curL)    

  return resL

def uniprot_wrangler(inpAcc):
  """to be called by prot_id_converter
  Take in a loaded json file which has uniprot accessions as results.
  Return a list of lists with where each sublist contains the results of a single query"""
  
  print "processing Uniprot Accessions"
  resL = []
  for inpAccD in inpAcc:
    queryL = inpAccD["UniProt Accession"].split("//")
    foundFlag = False
    curL = []
    for queryI in queryL:
      if queryI == "" or queryI == "-": continue
      curQuery = queryI.encode("ascii","ignore")
      curL.append(curQuery)
      foundFlag = True
    if not foundFlag: 
      print "GI not found in this query:"
      print inpAccD
    else:
      resL.append(curL)    

  return resL
  
def gene_symbol_wrangler(inpAcc):
  """to be called by prot_id_converter
  Take in a loaded json file which has gene symbols as results.
  Return a list with with the gene names from the query"""
  
  print "processing gene symbols"
  resL = []
  for inpAccD in inpAcc:
    queryI = inpAccD["Gene Symbol"]
    curQuery = queryI.encode("ascii","ignore")
    resL.append(curQuery)    

  return resL
  

def protein_gi_wrangler(inpAccessions):
  """to be called by prot_id_converter.
  Take in a loaded json entry from biodbnet db2db which has genbank protein accession ID results. 
  pick the largest GI number for each entry (this is really arbitrary though) and return a list of protein GIs. 
  These are unique identifiers (UIDs) that can be used by entrez.epost, entrez.efetch and such things.
  """
  print "processing GenBank Protein GIs"
  resL = []
  for inpAccD in inpAccessions:
    queryL = inpAccD["GenBank Protein GI"].split("//")
    maxQ = 1 
    for queryI in queryL:
      if queryI == "" or queryI == "-": continue
      curQuery = int(queryI.encode("ascii","ignore"))
      if curQuery > maxQ:
        maxQ = curQuery
    if maxQ == 1: 
      print "GI not found in this query:"
      print inpAccD
    else:
      resL.append(maxQ)    

  return resL
  
def refseq_gi_wrangler(inpAccessions):
  """to be called by prot_id_converter.
  Take in a loaded json entry from biodbnet db2db which has refseq protein accession ID results. 
  pick the smallest GI number for each entry (this should be a nice, full length sequence) and return a list of protein GIs. 
  These are unique identifiers (UIDs) that can be used by entrez.epost, entrez.efetch and such things.
  """
  print "processing GenBank Protein GIs"
  resL = []
  for inpAccD in inpAccessions:
    queryL = inpAccD["RefSeq Protein GI"].split("//")
    minQ = 999999999999999 # this should be big enough
    for queryI in queryL:
      if queryI == "" or queryI == "-": continue
      curQuery = int(queryI.encode("ascii","ignore"))
      if curQuery < minQ:
        minQ = curQuery
    if minQ == 999999999999999: 
      print "GI not found in this query:"
      print inpAccD
    else:
      resL.append(minQ)    

  return resL 
  
def accession_wrangler(inpAccessions):
  """to be called by prot_id_converter. 
  Take in a loaded json entry from biodbnet db2db which has genbank protein accession ID results. 
  Pick the refseq entry (NP) for each item and return a list of resulting refseq genbank protein accessions"""
  print "processing GenBank Protein Accessions"
  resL = []
  for inpAccD in inpAccessions:
    queryL = inpAccD["GenBank Protein Accession"].split("//")
    findFlag = False
    for queryI in queryL:
      if "NP" in queryI:
        resL.append(queryI.encode("ascii","ignore"))
        findFlag = True
        break
    if not findFlag:
      print "refseq entry not found in:"
      print inpAccD
  return resL
  
def prot_entrez_fetch(proteinList, retM="text", retT="fasta"):
  """take in a list of protein GI accessions and return their corresponding fasta sequences as a list using Entrez.
  Each returned list item is the fasta header + new line + sequence.
  Do not give it more than 200 GIs at once as requested by entrez.
  """
  from Bio import Entrez
  Entrez.email ="mate.ravasz@ed.ac.uk"
  for i in proteinList:
    try:
      int(i) # test if really a list of UIDs
    except ValueError:
      print "was expecting UIDs like \"12345678\", but got something else instead:"
      print i
      raise
    
  proteinList = map(str, proteinList)
  print "connecting to Entrez..."
  requestR = Entrez.epost("protein",id=",".join(proteinList)) # send all UIDs as a single query to entrez. 
  resultO = Entrez.read(requestR)
  webEnv = resultO["WebEnv"]
  queryKey = resultO["QueryKey"]
  handleO = Entrez.efetch(db="protein",retmode=retM, rettype=retT, webenv=webEnv, query_key=queryKey) # retrieve all results in batch
  print "connection successful"
  if retT == "fasta":
    return entrez_fasta_parser(handleO)
  elif retM == "text" and retT == "gp": # use "gp" for genpelt faltfile format, and "gb" for genbank flatfile for genes
    return handleO.read()
  else:
    print "this data format was not expected:"
    print "retmode: ", retM
    print "rettype: ", retT
    raise ValueError
    
def entrez_flatfile_parser(handleGB):
  """not written yet"""
  return handleGB.read()


def entrez_fasta_parser(handleFasta):
  """To be called by prot_entrez_fetch.
  Split and organise a single efetch query which is a string containing fasta sequences
  into a list of fasta sequences. Return the list."""
  fullList = handleFasta.read().split("\n") 
  resL = []
  seqFlag = False
  for fullLine in fullList:
    if fullLine == "":
      seqFlag = False
      continue
    elif fullLine[0] == ">":
      resL.append(fullLine + "\n")
      seqFlag = True
    elif seqFlag:
      resL[-1] += fullLine     
  return resL
  
def html_creator(titleS = "test HTML", bodyS = "test content", outpuFN = "test.html"):
  """create HTML file with titleS as title in header, bodyS as the body of the page and outpuFN as the filename"""
  with open(outpuFN, "w") as outF:
    outF.write(r"<!DOCTYPE html>" + "\n")
    outF.write(r"<html>" + "\n")
    outF.write(r"<head>" + "\n")
    outF.write(r"""<style>
p.monospace {
    font-family: "Courier New", Courier, monospace; line-height: 110%;
}

</style>""" + "\n")
    outF.write(r"<title>" + "\n")
    outF.write(titleS)
    outF.write(r"</title>" + "\n")
    outF.write(r"</head>" + "\n")
    outF.write(r"<body>" + "\n")
    outF.write(r"""<p class="monospace">""" + "\n")
    counterN = 0
    breakFlag = True
    needBreak = False
    for itemS in bodyS:
      outF.write(itemS)
      if itemS == "<": breakFlag = False
      if breakFlag: counterN += 1
      if itemS == ">": breakFlag = True
      
      if counterN % 100 == 0: needBreak = True        
      if needBreak and breakFlag:
        outF.write(r"<br>")
        needBreak = False
    outF.write(r"</p>" + "\n")
    outF.write(r"</body>" + "\n")
    outF.write(r"</html>" + "\n")
  

def go_term_advanced_lookup(protID):
  """like go_term_lookup but check first if GO annotation is already downloaded"""

  
  from os import listdir
  idList = [x[:x.find(".")] for x in listdir("/home/mate/workspace/katamari/src/ed/datafiles/go_terms")]
  
  if protID not in idList: 
    
    import urllib
  
    txtS = urllib.urlopen("http://www.ebi.ac.uk/QuickGO/GAnnotation?protein="+protID+"&db=UniProtKB&format=tsv")
    # print txtS.read()
    with open("/home/mate/workspace/katamari/src/ed/datafiles/go_terms/"+protID+".txt","w") as outF:
      for fileLine in txtS:
        outF.write(fileLine)
      
def go_term_lookup(protID):
  """take a uniprot id and rpint out any GO terms associated with it using quickGO
  
  To do:
  
  Take in a list of uniprot IDs, cut them into batches of 100, get their respective GO terms and return htem"""
  
  import urllib
  
  txtS = urllib.urlopen("http://www.ebi.ac.uk/QuickGO/GAnnotation?protein="+protID+"&db=UniProtKB&format=tsv")
  # print txtS.read()
  with open("/home/mate/workspace/katamari/src/ed/datafiles/go_terms/"+protID+".txt","w") as outF:
    for fileLine in txtS:
      outF.write(fileLine)

def biogrid_parser():
  '''
  Created on 23 Jun 2017
  
  @author: mate
  
  Take in an exported textfile from biogrid (biogrid TAB 2.0 format) and extract all interactor gene names. Return them in a list.
  
  '''
  
  inpFileName = "datafiles/BIOGRID-GENE-117604-3.4.149.tab2.txt"
  
  inpF = file_importer(inpFileName)
  
  interactorList = []
  
  headerFlag = True
  
  for inpLine in inpF:
    if headerFlag:
      headerFlag = False
      continue
    inpList = inpLine.rstrip("\n").split("\t")
    interactorList.append(inpList[8])
  
  inpF.close()
  return interactorList

def intact_parser():
  """open ptpn22.txt and extract prey protein uniprot accessions. 
  Convert those to refseq protein accessions.
  Return them as a list.
  
  added on 13.10.2016: 
  - check if interaction is actually between the protein of interest and a possible target -- done
  - check if interator name is acutally a valid uniprot ID -- done
  - remove isoform tags (-2, -2 etc) from uniprot accessions and just use the base name -- done
  
  added on 23-06-2017:
  - this is copied over from ppi_parser.py and modified here
  - instead of refseq GIs, return gene names
  - handle different organisms too using the taxids
  """
  
  baitStr = "PTPN22" # gene name of bait protein. Has to be entered all caps
  
  # relPath = "ptpn22_ppi_data/ptpn22.txt"
  relPath = "datafiles/ptpn22_interactions_intact.txt"
  inpF = file_importer(relPath, "r")
  headerFlag = True
  preyL = []
  nfCount = 0
  foundFlagA = False
  foundFlagB = False
  for inpLine in inpF:
    if headerFlag:
      headerFlag = False
      continue
    inpList = inpLine.split("\t")
    counterN = 0
    # print inpList
    while True: # find gene name in A interactor entry
      try:
        if inpList[4].split("|")[counterN].split(":")[1].split("(")[1][:-1] == "gene name": 
          foundFlagA = True
          break
      except IndexError:
        break
      
      counterN += 1
      
    if foundFlagA: 
      geneNameA = inpList[4].split("|")[counterN].split(":")[1].split("(")[0] # gene name A extracted here
      foundFlagA = False
    else: 
      geneNameA = "not found"
      nfCount += 1
    
    
    
    counterN = 0
    while True: # do the same for B interactor entry
      try:
        if inpList[5].split("|")[counterN].split(":")[1].split("(")[1][:-1] == "gene name": 
          foundFlagB = True
          break
      except IndexError:
        break      
      
      counterN += 1
      
    if foundFlagB:
      geneNameB = inpList[5].split("|")[counterN].split(":")[1].split("(")[0] # gene name B extracted here
      foundFlagB = False
    else: 
      geneNameB = "not found"
      nfCount += 1
    
    if geneNameA.upper() not in preyL:
      preyL.append(geneNameA.upper())
    if geneNameB.upper() not in preyL:
      preyL.append(geneNameB.upper())
    
    """
    if geneNameA.upper() == baitStr or geneNameB.upper() == baitStr: # check if bait protein is one of the interactors
      inpItem = inpList[1].split(":")[-1]
      inpSecItem = inpList[0].split(":")[-1]
      
      if "-" in inpItem: # remove - symbols
        inpItem = inpItem[:inpItem.index("-")]

      if "-" in inpSecItem:
        inpSecItem = inpSecItem[:inpSecItem.index("-")]        
    
    else: continue
    
    print inpItem, " ", inpSecItem
    
    if "\"" not in inpItem and inpItem != "not found" and len(inpItem)>3 and inpItem not in preyL:
      preyL.append(inpItem)
      
    if "\"" not in inpSecItem and inpSecItem != "not found" and len(inpSecItem)>3 and inpSecItem not in preyL:
      preyL.append(inpSecItem)      
      
    """
  inpF.close()
  
  try:
    preyL.remove("NOT FOUND")
  except ValueError:
    pass
  
  try:
    preyL.remove(baitStr)
  except ValueError:
    pass
      
  
  print preyL
  # idList = prot_id_converter(preyL, "", outDB="refseqproteingi") # convert uniprot ID to refseq accessions
  # return idList
      
  
if __name__ == "__main__":
  main()
  