'''
Created on 27 May 2016

@author: mate

Take the output file from abs_quant.py and the folder with all the GO terms and fuse them into a single file.
'''

# with open("/home/mate/workspace/katamari/src/ed/bob/processed/abs_quant_24H_T4-14-24-27-05-20162.txt","r") as inpF:
#   with open("/home/mate/workspace/katamari/src/ed/bob/processed/abs_quant_24H_T4-14-24-27-05-20164.txt","w") as outF:
with open("/home/mate/code/ed/src/data/processed/julia_time_course/julia_t0-15-50-08-08-2019_2.txt","r") as inpF:
  with open("/home/mate/code/ed/src/data/processed/julia_time_course/julia_t0-15-50-08-08-2019_3.txt","w") as outF:
  
    headerFlag = True
    
    for inpLine in inpF:
      inpL = inpLine.split(",")
      if headerFlag:
        headerFlag = False
        outF.write(inpLine.rstrip("\n") + "," + "GO Function," + "GO Process," + "GO Component\n")
        continue
      # this block here opens the GO term file associated with the query at hand, and collects out the GO terms into 3 lists
      with open("/home/mate/workspace/katamari/src/ed/datafiles/go_terms/"+inpL[1]+".txt","r") as goF:
        goHead = True
        goFun = []
        goProc = []
        goComp = []
        for goLine in goF:
          if goHead:
            goHead = False
            continue

          goL = goLine.split("\t")
          currGo = goL[7].replace(",","-")
          # print(goLine)
          if currGo == "molecular_function" or currGo == "biological_process" or currGo == "cellular_component":
            continue
          
          if goL[11] == "Function" or goL[11] == "molecular_function":
            if currGo not in goFun:
              goFun.append(currGo)
          elif goL[11] == "Process" or goL[11] == "biological_process":
            if currGo not in goProc:
              goProc.append(currGo)
          elif goL[11] == "Component" or goL[11] == "cellular_component":
            if currGo not in goComp:
              goComp.append(currGo)        
          else:
            print("got something odd:")
            print(goLine)
            raise ValueError
 
      funS = ""
      for funI in goFun:
        funS += funI
        if funI is not goFun[-1]:
          funS += ";"
      
      procS = ""
      for procI in goProc:
        procS += procI
        if procI is not goProc[-1]:
          procS += ";"
          
      compS = ""
      for compI in goComp:
        compS += compI
        if compI is not goComp[-1]:
          compS += ";"                  
      
      outF.write(inpLine.rstrip("\n") + "," + funS + "," + procS + "," + compS + "\n")
      

        
    
