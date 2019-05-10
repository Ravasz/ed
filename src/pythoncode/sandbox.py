'''
Created on 13 May 2015

@author: mate

This is a sandbox script for testing any sort of code snippets. 
Can be freely edited or deleted anytime.
'''


def main():  
  print("you are in the sandbox")
  pandas_tester()
  
  
  
  
def test_bio():
  from Bio.Seq import Seq
  from Bio.Alphabet import IUPAC
  new_seq = Seq ("GATCAGAAG", IUPAC.unambiguous_dna)
  print(new_seq)    
  
  import time
  print((time.strftime("%d-%m-%Y")))
  
  print(int("150"))
  a = [1, 3, 5]
  a[:] = [x + 2 for x in a]
  print(a)
  
  import Bio
  print(Bio.__version__)
  


def pandas_tester():
  import pandas as pd
  
  df = pd.DataFrame({'targetCol': [True, False, False, True, True], 'secondCol': ["one", "two", "three", "four", "five"], 'numCol': [4, 8, 2, 16, 32]}, index=[10,20,30,40,50])  
  print(df)
  
  df["newCol"] = df["numCol"]
  del df["numCol"]
  
  print(df)
   


def floater():
  from decimal import Decimal
  a = 0.1
  print(a)
  print(Decimal(a))


def finder_thing():
  """just find a line in a file..."""
  with open("bob/processed/24h_bobdata_ed2.csv", "r") as inpF: # read and process the csv with protein names and p values
    for inpLine in inpF:
      if "Ptpn22" in inpLine:
        print(inpLine)


  

def p_value_explorer():
  """find all entries in bob's data that have the ominous 0.4226 for p value and print them out in STDout"""
  countD = {}
  with open("bob/bobdata2_ed2.csv", "r") as inpF:
    headerFlag = True
    for inpLine in inpF:
      if headerFlag:
        headerFlag = False
        continue
      inpList = inpLine.split("\" \"")
      inpP = inpList[10].strip(" \"\n")
      if inpP == "0.4226497308": print(inpList)
      if inpP not in countD:
        countD[inpP] = 1
      else:
        countD[inpP] += 1
    
  v=list(countD.values())
  k=list(countD.keys())
  print(k[v.index(max(v))], max(v))


def box():
  print("box module reporting")


if __name__ == "__main__":
  main()