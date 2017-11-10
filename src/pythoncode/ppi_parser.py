'''
Created on 25 Jun 2015

@author: mate
extract gene names from ptpn22.txt (retrieved from the IntAct PPI database 25.06.2015)
and convert them to GI accessions using biodbnet.
Then, get all their sequences from entrez and write these out in STDOut.
'''

def main():
  print("this is ppi_parser \n")
  # name_collector()
  # sh3_counter()
  protname_compare()
  # biogrid_parser()


def protname_compare():
  """quick script to compare biogrid interactor list to intact interactor list and return a list of all interactors present in at least one database"""
  import os.path
  
  intactList = intact_parser()

  intList = [x.lower() for x in intactList] # convert the list to lowercase

  
  inpF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "BIOGRID-GENE-113578-3.4.152.tab2-genes.txt"),"rU") 
  # inpF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "biorgid-ptpn22-interactors-08082017.txt"),"rU") # needs rU as it is from mac. 
  # This is stripped from the biogrid output file by loading it into excel and copying both columns with genes names in it to a single column in textedit. 
  # I could write one that parses the biogrid output file but I won't as I'm lazy
  
  fullList = []
  
  for inpLine in inpF:
    if inpLine.strip().lower() not in fullList:
      fullList.append(inpLine.strip().lower())
  
  for intactI in intList:
    if intactI.lower() not in fullList:
      fullList.append(intactI.lower())
  
  # print "---"
  
  titleL = []

  for outI in fullList:
    titleL.append(outI.upper())
    print(outI.upper())
  return titleL

def biogrid_parser():
  """open biogrid 2.0 tab format file and extract interactor protein gene IDs. Convert to refseq protein accessions. Return them a list.
  """
  
  import os.path
  from .tools import prot_id_converter
  
  print("processing biogrid file")
  
  inpF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "BIOGRID-GENE-117604-3.4.151.tab2.txt"),"r")
  
  idL = []
  
  headerFlag = True
  
  for inpLine in inpF:
    if headerFlag:
      headerFlag = False
      continue
    inpList = inpLine.split("\t")[1:3]
    for inpI in inpList:
      if inpI not in idL:
        idL.append(inpI)
  # print idL
  
  idList = prot_id_converter(idL, orgnID = "9606", inpDB = "geneid", outDB="refseqproteingi") # convert uniprot ID to refseq accessions or gene names
  # print idList
  return idList
  
def intact_parser(outDataType = "genesymbol"):
  """open ptpn22.txt and extract prey protein uniprot accessions. 
  Convert those to refseq protein accessions.
  Return them as a list.
  
  added on 13.10.2016: 
  - check if interaction is actually between the protein of interest and a possible target -- done
  - check if interator name is acutally a valid uniprot ID -- done
  - remove isoform tags (-2, -2 etc) from uniprot accessions and just use the base name -- done
  
  todo:
  - retrieve both gene name and full protein name
  """
  from .tools import prot_id_converter
  import os.path
  
  print("processing intact file")
  
  baitStr = "PTPN22" # gene name of bait protein. Has to be entered all caps
  
  # relPath = "ptpn22_ppi_data/ptpn22.txt"
  # relPath = "ptpn22_ppi_data/csk.txt"
  #relPath = "data/ptpn22-intact-mi-tab27-07-08-2017.txt"
  #inpF = file_importer(relPath, "r")
  # inpF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "ptpn22-intact-mi-tab27-07-08-2017.txt"),"r")
  inpF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "ptp4a1.txt"),"r")
  headerFlag = True
  preyL = []
  foundFlagA = False
  foundFlagB = False
  for inpLine in inpF:
    if headerFlag:
      headerFlag = False
      continue
    inpList = inpLine.split("\t")
    counterN = 0
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
    else: geneNameA = "not found"
    
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
      else: geneNameB = "not found"
    
    if geneNameA.upper() == baitStr or geneNameB.upper() == baitStr: # check if bait protein is one of the interactors
      inpItem = inpList[1].split(":")[-1]
      inpSecItem = inpList[0].split(":")[-1]
      
      if "-" in inpItem: # remove - symbols
        inpItem = inpItem[:inpItem.index("-")]

      if "-" in inpSecItem:
        inpSecItem = inpSecItem[:inpSecItem.index("-")]        
    
    else: continue
    
    # print inpItem, " ", inpSecItem
    
    if "\"" not in inpItem and inpItem != "not found" and len(inpItem)>3 and inpItem not in preyL:
      preyL.append(inpItem)
      
    if "\"" not in inpSecItem and inpSecItem != "not found" and len(inpSecItem)>3 and inpSecItem not in preyL:
      preyL.append(inpSecItem)      
      

  inpF.close()
  
  # idList = prot_id_converter(preyL, "", outDB="refseqproteingi") # convert uniprot ID to refseq accessions
  idList = prot_id_converter(preyL, "", outDB=outDataType) # convert uniprot ID to refseq accessions or gene names
  """for idI in idList:
    print idI"""
  return idList

def name_collector():
  """look up protein interactors, download their fasta sequences to extract their names. 
  Print the names to STDout"""
  from .tools import prot_entrez_fetch
  
  resD = {}
  intactIdList = intact_parser("refseqproteingi") 
  biogridIdList = biogrid_parser()
  
  combinedList = biogridIdList
  
  for intactI in intactIdList:
    if intactI not in combinedList:
      combinedList.append(intactI)
  
  fullL = prot_entrez_fetch(combinedList, retM="text", retT="gp").split("//\n")

  for fullI in fullL:
    longFlag = False
    shortFlag = False
    completeFlag = True
    longName = ""
    shortName = ""
    
    for fullIline in fullI.split("\n"):
      
      if longFlag and not completeFlag:
        if fullIline[:9] != "ACCESSION":
          longName += " "
          longName += fullIline.lstrip(" ").split("[")[0]
        longFlag = False
      
      if shortFlag:
        shortName = fullIline.split("\"")[1]
        shortFlag = False
        
      if "DEFINITION" in fullIline[:11]:
        completeFlag = False
        # print fullIline
        longName = fullIline[12:].split("[")[0]
        if len(fullIline[12:].split("[")) > 1:
          completeFlag = True
        longFlag = True
        
      if "     CDS" in fullIline[:11]:
        shortFlag = True
    
    if shortName.upper() not in resD and shortName != "":
      resD[shortName.upper()] = longName
      
  for j,k in list(resD.items()):
    print(j) #, ": ", k
  print(len(resD))
  
  """fastaL = prot_entrez_fetch(idList, retM="text", retT="fasta")
  for fastaItem in fastaL:
    nameS = fastaItem.split("\n")[0].split("|")[-1]
    print nameS"""


def sh3_counter():
  """look up a list of uniprot IDs, download their full genbank entries from the Entrez database 
  and count the number of SH3 domains the interactors have. Print the results to STDout"""
  from .tools import prot_entrez_fetch, prot_id_converter
  # from bobdata_parser import protein_name_collector
  fullPreyL = protname_compare()
  # fullPreyL = intact_parser() # to use for Ptpn22 interactome
  # fullPreyL = protein_name_collector() - use for Bob's data
  # fullPreyL = ["P20152", "Q8BFZ3", "P17182", "P17742", "P11499"]
  # print fullPreyL
  if len(fullPreyL) > 200: # chop up very large lists of uniprot IDs to batches of 100
    maxBatch = (len(fullPreyL)/200) + 1
    lenCount = 0
    idList = []
    seqL = []
    batchCount = 0
    preyL = []
    for listItem in fullPreyL:
      preyL.append(listItem)
      lenCount += 1
      if lenCount == 200:
        batchCount += 1
        print("processing batch number %d of Uniprot IDs..." % (batchCount, ))
        idList = prot_id_converter(preyL, "10090", outDB="refseqproteingi")
        seqL = seqL + prot_entrez_fetch(idList, retM="gb", retT="text").split("\n") # fetch the complete genbank entries from entrez using this function from tools.py
        lenCount = 0
        idList = []
        preyL = []
        
        print("this was batch number %d of %d" %(batchCount, maxBatch)) 
        print("")
    if lenCount != 0:
      batchCount += 1
      print("processing batch number %d of Uniprot IDs..." % (batchCount, ))
      idList = prot_id_converter(preyL, "10090", outDB="refseqproteingi")
      seqL = seqL + prot_entrez_fetch(idList, retM="gb", retT="text").split("\n") # fetch the complete genbank entries from entrez using this function from tools.py
      lenCount = 0
      idList = []
      preyL = []
      
      print("this was batch number %d of %d" %(batchCount, maxBatch)) 
      print("")
      
        
  else: 
    # idList = prot_id_converter(fullPreyL, "10090", outDB="refseqproteingi")
    idList = prot_id_converter(fullPreyL, "9606", inpDB = "genesymbol", outDB="refseqproteingi")
    # print len(idList)
    seqL = prot_entrez_fetch(idList, retM="text", retT="gp").split("\n") # fetch the complete genbank entries from entrez using this function from tools.py
    # print len(seqL)
  
  """
  dumpStr = ""
  for seqLine in seqL:
    seqStr = "***".join(seqLine) 
    dumpStr = dumpStr + "xxxxx" + seqStr # this might get huge
  
  outputDump = open("entrezdump.txt", "w")
  outputDump.write(dumpStr)
  """
  
  regionFlag = False
  sHFlag = True
  regionCount = 0
  sHCount = 0
  protCount = 0
  shProtCount = 0
  for flatLine in seqL:
    # print flatLine
    if flatLine[:10]  == "LOCUS     ": # this is usually the title of an entrez flatfile and contains the protien name 
      protCount+= 1
      newProtFlag = True
      # curProt = flatLine
    if flatLine[:11] == "     Region": # regions like this mark domains in the flatfile. Look here for SH3 domains
      regionCount += 1
      regionFlag = True
      continue
    if regionFlag:
      if flatLine[:11] == "           ":
        if "SH3" in flatLine and sHFlag:
          sHCount += 1
          sHFlag = False
          if newProtFlag:
            shProtCount += 1
            newProtFlag = False
          # print flatLine
          # print curProt
      else: 
        regionFlag = False
        sHFlag = True
  print("%d SH3 domains found in %d domains of %d proteins" % (sHCount, regionCount, protCount))
  print("%d proteins out of %d contain SH3 domains" % (shProtCount, protCount))

def intact_publications():
  """open orc6.txt and extract publication first authors and pubmed IDs. 
  print results to STDout"""
  
  import os.path
  with open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "orc6.txt"),"r") as inpF:
    #with open("orc6.txt","r") as inpF:
    headerFlag = True
    pubD = {}
    pubmedIDList = []
    for inpLine in inpF:
      if headerFlag:
        headerFlag = False
        continue
      inpList = inpLine.split("\t")
      if inpList[7] not in pubD:
        pubD[inpList[7]] = inpList[8] # 7th item is publication first author plus publication year. 8th item is article ID in often multiple databases
        pubIdList = inpList[8].split("|") # this par and the next few lines extract the pubmed IDs
        for listItem in pubIdList:
          if "pubmed" in listItem:
            pubmedIDList.append(listItem.split(":")[1])
      
  print(pubD)
  print(len(pubD))
  print(pubmedIDList)


if __name__ == "__main__":
  main()