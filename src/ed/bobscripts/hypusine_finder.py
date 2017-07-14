'''
Created on 27 Jun 2016

@author: mate
used to parse hypusine peptide list files from maxquant output

'''
def dohh_hyp_finder():
  """take list of DHS and DOHH interactors from the paper 
  "Protein-protein-interaction Network Organization of the Hypusine Modification System" 
  and check if any of them are in the new hypusinylation dataset I made from Bob's 24H T4 data."""
  with open("/home/mate/workspace/katamari/src/root/ed/datafiles/DOHH_interactors_published.csv","r") as dohF:
    with open("/home/mate/workspace/katamari/src/root/ed/datafiles/hypusine-24H-T4-proteomics-27062016.txt","rU") as hypF:
      qL = []
      for inpLine in dohF:
        inpL = inpLine.rstrip().split(" ")
        for inpI in inpL:
          if inpI not in qL:
            qL.append(inpI)
      
      hL = []
      for hypLine in hypF:
        hypL = hypLine.rstrip().split("\t")
        hypLI = hypL[5].split(";")
        for hypI in hypLI:
          if hypI not in hL:
            hL.append(hypI)
        
  
  print qL
  print hL
  
  for hI in hL:
    if hI in qL:
      print hI


def hyp_comp():
  """take list of hypusine containing proteins for two different datasets (24H TCR and Jens' TCR data pxd002928) and compare them to see if they overlap greatly. 
  They should reasonably overlap if the hypusine finding went well."""
  with open("/home/mate/workspace/katamari/src/root/ed/datafiles/HypusineSites_pxd002928.txt","rU") as bobF:
    with open("/home/mate/workspace/katamari/src/root/ed/datafiles/hypusine-24H-T4-proteomics-27062016.txt","rU") as jensF:
      bobList = []
      for bobLine in bobF:
        bobL = bobLine.rstrip().split("\t")
        bobLI = bobL[5].split(";")        
        bobList.append(bobLI[0])
      
      jensList = []
      for jensLine in jensF:
        jensL = jensLine.rstrip().split("\t")
        jensLI = jensL[5].split(";")      
        jensList.append(jensLI[0])      
        
      for bobI in bobList:
        if bobI in jensList:
          print bobI
        



hyp_comp()