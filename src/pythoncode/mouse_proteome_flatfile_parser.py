'''
Created on 11 Mar 2016

@author: mate
process the full mouse proteome with sh2_counter to find all SH2 domain containing proteins in the mouse, 
then cross-reference this list with the existing CTL proteomics dataset, RNA 24H T4 stim dataset, 
and protein 24H T4 datast to find proteins we actually detected.
'''

def main():
  """main function"""
  print("this is mouse_proteome_flatfile_parser")

def sh2_dicter():
  """read in the curated sh2 containing protein list and return it as a neat dict"""
  
  outD = {}
  
  with open ("../processed_data/sh2_containing_proteins_mouse_curated.csv","rU") as inpF:
    headerFlag = True
    for inpLine in inpF:
      if headerFlag:
        headerFlag = False
        continue
      inpList = inpLine.strip("\r\n").split(",")
      outD[inpList[0]] = inpList[1:]
  return outD

def sh2_counter():
  """open datafiles/entrez_flatfile_all_mouse_proteins.txt and find all proteins that have an SH2 domain in them. 
  Write out results into a new file."""
  
  with open("datafiles/entrez_flatfile_all_mouse_proteins.txt","r") as inpF:
    
    lineCount = 0
    resD = {}
    curGene = ""
    curDef = ""
    curAcc = ""
    curDom = ""
  
    newEntry = True
    defFlag = False
    cdsFlag = False
    foundFlag = False
  
    for inpLine in inpF:
      lineCount += 1
      # if lineCount == 1000000: break
      inpLine = inpLine.strip()
      if defFlag: # this bit is used to pull out the second line of the description if it is too long for a single line
        defFlag = False
        if inpLine[:9] == "ACCESSION":
          pass
        else:
          curDef = curDef + inpLine
          curDef = curDef.replace(",", "-")
          continue
      if cdsFlag: # pull out the gene name if the previous line was the CDS start
        curGene = inpLine[6:].strip("\"")
        cdsFlag = False
        continue
      if inpLine == "//":  # end of entry. collect data if needed here
        newEntry = True
        if foundFlag:
          foundFlag = False
          if curGene not in resD or resD[curGene][1][:9] == "PREDICTED":
            resD[curGene] = [curAcc, curDef, curDom]  
        continue
      if inpLine[:9] == "ACCESSION": # collect accession numbers
        curAcc = inpLine[12:]
        continue
      if inpLine[:10] == "DEFINITION": # collect protein name. This is a bit tricky. see handling defFlag above
        curDef = inpLine[12:]
        defFlag = True
        continue
      if inpLine[:3] == "CDS":
        cdsFlag = True  
        continue    
      if newEntry and inpLine[:12] == "/region_name":
        if "SH2" in inpLine[13:]:
          curDom = inpLine[13:].strip("\"").replace(",","-")
          newEntry = False
          foundFlag = True
    print("file parsed successfully.")
  
  with open("processed_data/sh2_containing_proteins_mouse.csv","w") as outF:
    
    outF.write("Gene_name,Accession,full_name,domain_found\n")
    
    for itemStr in resD:
      outF.write(itemStr)
      outF.write(",")
      for resItem in resD[itemStr][:-1]:
        outF.write(str(resItem))
        outF.write(",")
      outF.write(resD[itemStr][-1])
      outF.write("\n")
    print("%d outputs written to result file" %(len(resD)))


if __name__ == "__main__":
  main()