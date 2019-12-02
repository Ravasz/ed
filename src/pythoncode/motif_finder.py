'''
Created on 11 Apr 2019

@author: mate

Use regex to find simple motifs in protein sequences
'''


def starter():
  """run all the scripts"""
  from tools import uniprot_dicter, go_term_advanced_lookup
  
  mouseD = uniprot_dicter()
  
  with open("/home/mate/code/ed/src/data/cav1ko/processed/proteinGroups_EV_matched_samples_24-11-2018_exosome_Cav1ko_vs_wt_all_datasets_24-11-2018_combined_2019-04-10-3.csv","r") as inpF:
    with open("/home/mate/code/ed/src/data/cav1ko/processed/proteinGroups_EV_matched_samples_24-11-2018_exosome_Cav1ko_vs_wt_all_datasets_24-11-2018_combined_2019-04-10-3_cholesterol_2.csv","w") as outF:

      headerS = next(inpF)  
      outF.write(headerS.rstrip("\n") + ", Protein sequence, Cholesterol binding\n")
      
      lineCount = 0
      for inpLine in inpF:
        lineCount += 1
        print(".", end = "")
        if lineCount == 100:
          lineCount = 0
          print("\n")
        inpL = inpLine.rstrip().split(",")
        protID = inpL[0]
        
        go_term_advanced_lookup(protID)
        
        protSeq = mouseD[protID]
        inpL.append(protSeq)
        # protSeq = "MSGGKYVDSEGHLYTVPIREQGNIYKPNNKAMADELSEKQVYDAHTKEIDLVNRDPKHLNDDVVKIDFEDVIAEPEGTHSFDGIWKASFTTFTVTKYWFYRLLSALFGIPMALIWGIYFAILSFLHIWAVVPCIKSFLIEIQCISRVYSIYVHTVCDPLFEAVGKIFSNVRINLQKEI"
        inpL.append(str(motif_finder(protSeq)))
        
        for outI in inpL[:-1]:
          outF.write(outI + ",")
          
        outF.write(inpL[-1] + "\n")

        
  


def motif_finder(protSeq):
  """do the regex on a prtein sequence string
  regular expressions are hard coded below"""
  import re
    
  searchS = r"[LV].{0,5}Y.{0,5}[KR]" # CRAC motif (L/V)-X1−5-(Y)-X1−5-(K/R)
  searchSBack = r"[KR].{0,5}[YF].{0,5}[LV]" # CARC motif (K/R)-X1−5-(Y/F)-X1−5-(L/V)
  
  resO = re.search(searchS,protSeq)
  resOBack = re.search(searchSBack,protSeq)
  
  if resO: return resO.group()
  elif resOBack: return resOBack.group()
  else: return 0
  
#   print(resO)
#   print(resO.span())
#   print(resO.string)
#   print(resO.group())


if __name__ == "__main__":
  starter()