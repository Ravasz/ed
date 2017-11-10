'''
Created on 29 Feb 2016

@author: mate

module to analyize RNA seq dataset from Al.
This dataset is made using a 24 hour T4 stimulation of murine OT-1 T cells. 
PTPN22 KO and WT cells are compared. After 24H, RNA is isolated and sequenced. 
Dataset contains these RNAs, with ensembl ID, gene name, entrez id, fold change, 
average expression, raw P value, adjusted P, B.
Fold change is KO expression level divided by WT, log2. 
Hence, PTPN22 is around -4.8, meaning it is about 32 fold downregulated in KO cells over wt.
'''

def avg_dev(someD):
  """calculate number of entries, averages, and standard deviation for the input which is a dict of lists. """
  import numpy # calculate average fold changes and SD for each gene
  from math import isnan
  calD = {}
  # numpy.seterr(all="ignore")
  for keyS, valueLL in list(someD.items()):
    fCList = []
    for valueL in valueLL:
      fCList.append(valueL[8])
    # avgN = numpy.mean(fCList, dtype=numpy.float64) - wow, numpy mean does a very inaccurate floating point calculation even in 64 bits
    avgN = round(sum(fCList)/len(fCList),5)
    if len(fCList) == 1: stDev = 0
    else: stDev = round(numpy.std(fCList, ddof = 1),5)
    if isnan(stDev):
      calD[keyS] = [valueL[11], len(fCList), avgN, 0]
    else:
      calD[keyS] = [valueL[11], len(fCList), avgN, stDev]
  print("averages, standard deviations calculated")
  return calD

def log_converter(numberF):
  """convert log10 base of number to log2"""
  from math import log
  
  convF = log(10**numberF,2)
  
  return convF
    
def rna_parser():
  """process the tab delimited txt file (converted from mac excel) that contains the RNA seq results made by Al.
  return a dict with unique gene names as keys, and """
  from collections import defaultdict
  
  with open("../bob/002_WT-KO_genelist copy.txt","rU") as inpF:
    headerFlag = True
    parsedD = defaultdict(list)
    
    #FeatureID Symbol Description EntrezID Chrom ChromLoc PathwayID PathDesc logFC AveExpr t P.Value adj.P.Val B
    for inpLine in inpF:
      if headerFlag:
        headerFlag = False
        continue
      procList = []
      inpList = inpLine.split("\t")
      procList.extend(inpList[:3]) # add featureID, symbol, Description
      inpList[-1]=inpList[-1].rstrip()
      
      if inpList[3] == "NA": procList.append(0) # add entrezID
      else: procList.append(int(inpList[3]))
      
      try: procList.append(int(inpList[4])) # add chrom
      except ValueError: procList.append(inpList[4])
      
      try: procList.append(int(inpList[5])) # add chromloc
      except ValueError: procList.append(inpList[5])
      
      for curItem in inpList[6:8]: # PathwayID, PathDesc
        procList.append(curItem)
      
      try: procList.append(log_converter(float(inpList[8]))) # convert fold changes to log 2
      except ValueError: print(inpList[8])
      
      for curItem in inpList[9:]: # logFC AveExpr t P.Value adj.P.Val B
        try:
          procList.append(float(curItem))
        except ValueError:
          procList.append(curItem)
          
      #['ENSMUST00000147883', 'Ctsa', 'cathepsin A', 19025, 2, 164658372, '"04142, 04614"', '"Lysosome, Renin-angiotensin system"', 0.18, 2.1, 0.33, 0.75, 1.0, -8.1]
      if procList[1] != "NA" and procList[9] > 0: # filter out unexpressed genes
        parsedD[procList[1]].append(procList)
  print("file parsed successfully, %d genes found" % len(parsedD))
  
  extD = avg_dev(parsedD)
  with open("../bob/processed/001_KO-WT_genelist_export_Al_RNA-seq_unique.csv","w") as outF:
    outF.write("Gene name,number of unique mRNAs,Average fold change,standard deviation of fold change\n")
    for extI in extD:
      outF.write(extI + ",")
      for extS in extD[extI][:-1]:
        outF.write(str(extS) + ",")
      outF.write(str(extD[extI][-1]) + "\n")
  return extD
    

def prot_importer():
  """import processed and formatted proteomics data from the 24 hour bob dataset 
  and collect it into a dict with keys as gene names and values as a list of everything else
  """
  with open("../bob/processed/24h_bobdata_ed2_formatted.csv") as inpF:
    
    headerFlag = True
    protD = {}
    
    for inpLine in inpF:
      inpList = inpLine.split(",")
      if headerFlag:
        headerFlag = False
        continue     
       
      """if inpList[2] in protD:
        print protD[inpList[2]]"""
      protD[inpList[2]] = inpList

  print("found a total of %d proteins" % (len(protD)))
  return protD 

def fc_comparison(compareD):  
  """compare fold changes in the RNA and proteomics datasets.
  Print out results."""
  
  fCSigD = {}
  fCRNAD = {}
  fCProtD = {}
  fCnonSigD = {}
  
  for commonI in compareD:
    if abs(float(compareD[commonI][12])) > 1 and abs(float(compareD[commonI][15])) > 1:
      fCSigD[commonI] = compareD[commonI]
    elif abs(float(compareD[commonI][12])) > 1:
      fCProtD[commonI] = compareD[commonI]
    elif abs(float(compareD[commonI][15])) > 1:
      fCRNAD[commonI] = compareD[commonI]
    else:
      fCnonSigD[commonI] = compareD[commonI]
      
  print("double significant FC: ", len(fCSigD))
  print("RNA FC > 2", len(fCRNAD))
  print("protein FC > 2", len(fCProtD))
  print("FC < 2 in both datasets", len(fCnonSigD))

def p_value_comparison(compareD):  
  """compare p values in the RNA and proteomics datasets.
  Print out results."""
  
  fCSigD = {}
  fCRNAD = {}
  fCProtD = {}
  fCnonSigD = {}
  
  for commonI in compareD:
    if abs(float(compareD[commonI][10])) < 0.05 and abs(float(compareD[commonI][13])) < 0.05:
      fCSigD[commonI] = compareD[commonI]
    elif abs(float(compareD[commonI][10]))< 0.05:
      fCProtD[commonI] = compareD[commonI]
    elif abs(float(compareD[commonI][13])) < 0.05:
      fCRNAD[commonI] = compareD[commonI]
    else:
      fCnonSigD[commonI] = compareD[commonI]
      
  print("double significant p value: ", len(fCSigD))
  print("RNA p value < 0.05", len(fCRNAD))
  print("protein p value < 0.05", len(fCProtD))
  print("p value > 0.05 in both datasets", len(fCnonSigD))

def string_analyzer(targetD):
  """take stringDB output of protein interactors and map it against the RNA and protein dataset to find matches"""
  resL = [] # read in string file and put gene names into a list
  with open("../datafiles/tabdelimited_PDK1_network_string.txt","r") as inpF:
    for inpLine in inpF:
      inpList = inpLine.split("\t")
      if inpList[0] not in resL:
        resL.append(inpList[0])
  
  with open("../datafiles/tabdelimited_PDK1_network_string_matched.txt","w") as outF:
    hitCount = 0
    for baitI in resL:
      if baitI in targetD:
        hitCount += 1
        outF.write(str(targetD[baitI][2]) + "," + str(targetD[baitI][12]).rstrip("\n") + "," + str(targetD[baitI][15]) + "\n")
    
  print("out of %d proteins, %d was found in the dataset" %(len(resL),hitCount))


    

def main():
  from ed.mouse_proteome_flatfile_parser import sh2_dicter as sh2d
  
  rnaD = rna_parser() # build RNA dictionary
  protDict = prot_importer() # build protein dictionary
  shD = sh2d() # build dict of sh2 containing proteins from mouse proteome flatfile parser
  
  # group genes in these 3 dicts based on whether they are present in one dataset, the other or both
  commonD = {}
  protOnlyD = {}
  rnaOnlyD = {}
  
  # group of dicts distributing the SH2 domain containing proteins
  sh2inprotD = {}
  sh2inrnaD = {}
  sh2CommonD = {}
  
  for protItem in protDict:
    if protItem in rnaD:
      commonD[protItem] = protDict[protItem] + rnaD[protItem]
    else:
      protOnlyD[protItem] = protDict[protItem]
    
    if protItem in shD:
      sh2inprotD[protItem] = protDict[protItem] + shD[protItem]
      
  for rnaI in rnaD:
    if rnaI not in commonD:
      rnaOnlyD[rnaI] = rnaD[rnaI]
    if rnaI in shD:
      sh2inrnaD[rnaI] = rnaD[rnaI] + shD[rnaI]
  
  for comI in commonD:
    if comI in shD:
      sh2CommonD[comI] = commonD[comI] + shD[comI]
  
  print("number of genes in proteomics data in total: ", len(protDict))
  print("number of genes found in RNA seq: ", len(rnaD))
  print("number of genes found in both datasets: ", len(commonD))
  print("number of genes found in proteomics dataset only: ", len(protOnlyD))
  print("number of genes found in RNA seq only: ", len(rnaOnlyD))
  print("-----")
  print("number of proteins with SH2 domains", len(sh2inprotD))
  print("number of mRNAs encoding SH2 domains", len(sh2inrnaD))
  print("number of genes with both protein and RNA having SH2 domains", len(sh2CommonD))
  print("-----")
  """
  with open("../bob/processed/24h_bobgenes_in_common.csv","w") as outF:
    outF.write("oldID,UniprotID,GeneName,uniquePeptides,LFQKO1,LFQKO2,LFQKO3,LFQWT1,LFQWT2,LFQWT3,pValue,FDR,l2FC,raw P value in RNA seq,number of unique mRNAs,Average log 2 fold change RNA seq,standard deviation of fold change in RNA seq\n")
    for commonI in commonD:
      for cI in commonD[commonI][:-1]:
        outF.write(str(cI).rstrip("\n") + ",")
      outF.write(str(commonD[commonI][-1]) + "\n")"""
  
  with open("../bob/processed/24h_bobgenes_in_common_with_sh2.csv","w") as outF:
    outF.write("oldID,UniprotID,GeneName,uniquePeptides,LFQKO1,LFQKO2,LFQKO3,LFQWT1,LFQWT2,LFQWT3,pValue,FDR,l2FC,raw P value in RNA seq,number of unique mRNAs,Average log 2 fold change RNA seq,standard deviation of fold change in RNA seq,Accession,full_name,domain_found\n")
    for commI in sh2CommonD:
      for cI in sh2CommonD[commI][:-1]:
        outF.write(str(cI).rstrip("\n") + ",")
      outF.write(str(sh2CommonD[commI][-1]) + "\n")
    

  fc_comparison(commonD)
  p_value_comparison(commonD)
  # string_analyzer(commonD)


  

if __name__ == "__main__":
  main()