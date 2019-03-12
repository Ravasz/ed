'''
Created on 12 Mar 2019

@author: mate
quick script to convert mrna sequence ids to hgnc ids for Caludio's project
'''

from tools import prot_id_converter

with open("/home/mate/code/ed/src/data/tissue_mrna_names_processed.csv","w") as outF:
  with open("/home/mate/code/ed/src/data/tissue_mrna_names.csv","r") as inF:
    nameL = []
    next(inF)
    nameCount = 1
    for inLine in inF:
      nameCount += 1
      nameL.append(inLine.rstrip())
      if nameCount % 100 == 0:
        resL = prot_id_converter(nameL, inpDB = "ensemblgeneid", outDB = "hgncid", orgnID="9606")  
        # print(resL)
        for resI in resL:
          for keyS, valueS in resI.items():
            if keyS == "InputValue":
              outF.write(valueS.rstrip() + ",")
            else:
              if "//" in valueS:
                valL = valueS.split("//")
                for j in valL[:-1]:
                  outF.write(j.rstrip() + ";")
                outF.write(valL[-1].rstrip() + "\n")
              else:
                outF.write(valueS.rstrip() + "\n")

                  
        nameL = []
      
    resL = prot_id_converter(nameL, inpDB = "ensemblgeneid", outDB = "hgncid", orgnID="9606")  
    # print(resL)
    for resI in resL:
      for keyS, valueS in resI.items():
        if keyS == "InputValue":
          outF.write(valueS.rstrip() + ",")
        else:
          if "//" in valueS:
            valL = valueS.split("//")
            for j in valL[:-1]:
              outF.write(j.rstrip() + ";")
            outF.write(valL[-1].rstrip() + "\n")
          else:
            outF.write(valueS.rstrip() + "\n")
    nameL = []
