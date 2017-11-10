'''
Created on 29 Jun 2016

@author: mate

take gene name lists outputted by mass_spec_result_parser and compare gene nems to one another. Easy.
'''

with open("/home/mate/workspace/katamari/src/ed/datafiles/ost-mass-spec-gene-names-29062016.csv","r") as inpF:
  ostL = []
  for inpLine in inpF:
    ostL.append(inpLine.rstrip())
  print(ostL)

with open("/home/mate/workspace/katamari/src/ed/datafiles/ot1-mass-spec-gene-names-29062016.csv","r") as inpF:
  otL = []
  for inpLine in inpF:
    otL.append(inpLine.rstrip())
  print(otL)
  for oI in otL:
    print(oI)
  print("=====")
  
for ostI in otL:
  if ostI not in ostL:
    print(ostI)



