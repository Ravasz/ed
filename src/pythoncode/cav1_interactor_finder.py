'''
Created on 1 Apr 2019

@author: mate

Take in intact interaction file, biogrid interaction file, and the group 1 vs group 2 venn diagram datafile and prepare a 3 circle venn
with group1, group2 and interactors.
'''
def venn_maker(procFlag = True):
  """draw venn diagrams of 3 groups"""
  from tools import intact_parser, prot_id_converter
   
  intactFile = "/home/mate/code/ed/src/data/cav1ko/cav1-intact-interactome.txt"
  biogridFile = "/home/mate/code/ed/src/data/cav1ko/BIOGRID-GENE-107305-3.5.170.tab2.txt"
  groupFile = open("/home/mate/code/ed/src/data/cav1ko/processed/vennlists_cav1ko_vs_wt.txt","r")
   
  intact_parser(intactFile, biogridFile, "/home/mate/code/ed/src/data/cav1ko/processed/cav1_interactor.txt")
  
  bindFile = open("/home/mate/code/ed/src/data/cav1ko/processed/cav1_interactor.txt","r")
  bindL = []
  for bindLine in bindFile: bindL.append(bindLine.strip())
  
  bindSet = set(bindL)
  
  if procFlag:
  
    group1L = []
    group2L = []
    for groupLine in groupFile:
      groupL = groupLine.rstrip("\n").split("\t")
      group1L.append(groupL[0])
      if groupL[1] != "": group2L.append(groupL[1])
      
    
    group1Names = prot_id_converter(group1L, "10090", "uniprotaccession", "genesymbol")
    group2Names = prot_id_converter(group2L, "10090", "uniprotaccession", "genesymbol")
    
    group1Up = [x.upper() for x in group1Names]
    group2Up = [x.upper() for x in group2Names]
    
    group1UpSet = set(group1Up)
    group2UpSet = set(group2Up)
    
    if "CAV1" in group1UpSet: print("group1")
    if "CAV1" in group2UpSet: print("group2")
    
    return(bindSet, group1UpSet, group2UpSet)
    
  else:

    return(bindSet, None, None)

  


  
  return(bindSet, group1UpSet, group2UpSet)
  

def venn_plotter(group1, group2, group3):
  """prepare venn diagram images and save as file"""
  import matplotlib_venn
  import matplotlib.pyplot as plt
  import os
  
  outFolder = "/home/mate/code/ed/src/data/cav1ko/processed/"
  outFigName = "Cav1_venn"
  
  plt.figure()
  matplotlib_venn.venn3([group1, group2, group3])
  i = 1
  while os.path.exists(os.path.join(outFolder, (outFigName + "-" + str(i) + ".png"))): i += 1
  plt.savefig(os.path.join(outFolder, (outFigName + "-"  + str(i) + ".png"))) # write out to new file
  plt.close()
  
  bindGroupFSet = group1.intersection(group2.difference(group3))
  bindGroupSSet = group1.intersection(group3.difference(group2))
  bindBothSet = group1.intersection(group2.intersection(group3))
  group1only = group2.difference(group3)
  group2only = group3.difference(group2)
      
  print("\nonly group1:")
  for setI in bindGroupFSet:
    print(setI)
    
  print("\nonly group2:")
  for setI in bindGroupSSet:
    print(setI)
    
  print("\ncommon")
  print(len(bindBothSet))
  for setI in bindBothSet:
    print(setI)
    
  print("\nGroup1:")
  print(len(group1only))
#   for setI in group1only:
#     print(setI)

  print("\nGroup2:")
  print(len(group2only))
#   for setI in group2only:
#     print(setI)
    
  print(len(group2.union(group3)))

  with open(os.path.join(outFolder, "fullList_2.txt"), "w") as outF:
    for setI in group2.union(group3):
      outF.write(setI + "\n")
      
# 
#   print("\nWT only")
#   for setI in group2.difference(group3):
#     print(setI)
    
  # seq_writer(list(group2.difference(group3)))

def seq_writer(groupList):
  
  
  from tools import prot_entrez_fetch, prot_id_converter
  
  wtList = [x[0].upper() + x[1:].lower() for x in groupList]
  print(wtList)
  wtAccList = prot_id_converter(wtList, "10090", "genesymbol", "refseqproteingi")
  print(wtAccList)
  protSeqList = prot_entrez_fetch(wtAccList)
  
  with open("/home/mate/code/ed/src/data/cav1ko/processed/exosome_sequences_18-04-2019-1.txt","w") as outF:
    for seqItem in protSeqList:
      outF.write(seqItem)
      outF.write("\n")
    

def adhesome_parser():
  """open and parse the list of proteins in adhesomes from http://www.adhesome.org/downloads.html"""
  
  import pandas as pd
  pd.set_option('display.max_columns', 500)
  pd.set_option('display.width', 1000)
  
  # inFile = open("/home/mate/code/ed/src/data/cav1ko/components.csv","r")
  
  # inDF = pd.read_csv("/home/mate/code/ed/src/data/cav1ko/components.csv")
  # return set(inDF["Official Symbol"])
  inDF = pd.read_csv("/home/mate/code/ed/src/data/cav1ko/cholesterol_binders.csv", header = None)
  # inDF = pd.read_csv("/home/mate/code/ed/src/data/cav1ko/processed/null_names_2.txt", header = None)
  
  return set(inDF[0])
  
  # return set(list(inDF["Official Symbol"]))
  
#   print(inDF.shape)
#   print(inDF)
#   print(inDF.columns)
#   print(inDF[:2])
#   print(inDF.columns[0])
#   print(inDF["Official Symbol"])

def venn_lister():
  """create set of gene names from analyzed proteingroups file based on preset conditions"""
  import pandas as pd
  # import numpy as np
  pd.set_option('display.max_columns', 500)
  pd.set_option('display.width', 1000)  
  
  #inDF = pd.read_csv("/home/mate/code/ed/src/data/cav1ko/processed/proteinGroups_EV_matched_samples_24-11-2018_exosome_cav1ko_paper_11-12-2019_combined_2019-12-17-14.csv")
  
  inDF = pd.read_csv("/home/mate/code/ed/src/data/cav1ko/processed/proteinGroups_EV_matched_samples_24-11-2018_exosome_cav1ko_paper_08-01-2020_combined_2020-01-08-1.csv")


  
  # inDF = pd.read_csv("/home/mate/code/ed/src/data/cav1ko/processed/proteinGroups_EV_matched_samples_24-11-2018_exosome_cav1ko_vs_wt_all_datasets_24-11-2018_combined_2019-04-17-5.csv")

  inDF["Gene names"] = inDF["Gene names"].str.upper()
  
#   inDF["Venn groups"] = np.where(inDF["Avg1"] * 1.5 <= inDF["Avg2"], "group2", "")
#   inDF["Venn groups"] = np.where(inDF["Avg2"] * 1.5 <= inDF["Avg1"], "group1", inDF["Venn groups"])

#   inDF["Venn groups"] = np.where(inDF["Avg1"] * 1.5 <= inDF["Avg2"], "group2", "")
#   inDF["Venn groups"] = np.where(inDF["Avg2"] * 1.5 <= inDF["Avg1"], "group1", inDF["Venn groups"])
  
  venn1DF = inDF[inDF["Avg1"] > 1000000] # prots present in group 1
  
  venn2DF = inDF[inDF["Avg2"] > 1000000] # prots present in group 2
  
  print(inDF[inDF["Gene names"] == "CAV1"])
  print(inDF[inDF["Gene names"] == "NPC1"])
  print(inDF[inDF["Gene names"] == "NPC2"])
  
  # print(list(inDF["Gene names"]))
  
  # return(set(list(venn1DF["Majority protein IDs"])),set(list(venn2DF["Majority protein IDs"])))
  return(set(list(venn1DF["Gene names"])),set(list(venn2DF["Gene names"])))

  
def runner():
  """run the functions needed"""
  vennO = venn_lister()

  # venn_plotter(venn_maker(False)[0],vennO[0],vennO[1])
  venn_plotter(adhesome_parser(),vennO[0],vennO[1])
  # venn_maker(False)
  
def row_finder():
  """find rows based on the gene name and save the rows to a new single file"""
  
  
  
if __name__ == "__main__":
  runner()