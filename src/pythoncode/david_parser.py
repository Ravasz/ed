'''
Created on 18 Jan 2016

@author: mate

take files outputted by DAVID and do stuff with em:
- cluster_miner will process clusters and output the first 3 terms of significant clusters into a new file
- chart miner will take a functional annotation chart and output all terms with a bonferroni p value < 0.01
'''


def main():
  print("this is DAVID parser")
  cluster_miner()

def cluster_miner():
  """take a functional annotation clustering output file 
  and extract the first 3 terms from each cluster with an enrichment score greater than 2 (-log(0.01))"""
  
  import os.path
  
  print("this is cluster_miner")
  
  inpF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "cav1ko", "processed", "t2t_5979A7E9CBA41503507470597.txt"),"rU")
  outF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "cav1ko", "processed", "david_clusters_t2t_5979A7E9CBA41503507470597_members.txt"),"w")
  
  # with open("bob/processed/DAVID_FAC_24hbobprots_down_18012016_clusters.txt", "w") as outF:
    # with open("bob/processed/DAVID_FAC_24hbobprots_down_18012016.txt", "r") as inpF:
  clusterFlag = False
  entryCount = 1
  headerFlag = True
  startFlag = True
  protL = []
  outF.write("Enrichment score\t1st term\t2nd term\t3rd term\tnumber of proteins\n")
  for inpLine in inpF:
    inpList = inpLine.split("\t")
    if inpList[0] == "\n": 
      clusterFlag = False
      continue # skip empty lines
    
    if len(inpList) > 4:
      currProtL = inpList[5].rstrip("\t").split(",") # count member genes here
      for currProtI in currProtL:
        if currProtI not in protL:
          protL.append(currProtI)
        
    if clusterFlag: # take out first 3 terms from each cluster
      if entryCount < 4:
        if headerFlag:
          headerFlag = False
          continue
        if "~" in inpList[1]: # remove go term ID, just display the name
          outF.write(inpList[1][inpList[1].index("~")+1:] + "\t")
          # print repr(inpList[1][inpList[1].index("~")+1:])
        else:
          outF.write(inpList[1] + "\t")
          # print inpList[1]
        entryCount += 1
      else:
        clusterFlag = False

          
        
        
    if inpList[0][:18] == "Annotation Cluster": # start a new cluster here
      if not startFlag: outF.write(str(len(protL)) + "\t") 
      protL = []
      EnrScore = inpList[1][17:].strip()
      if float(EnrScore) < 2: break
      if startFlag:
        outF.write(str(EnrScore) + "\t")
        startFlag = False
      else:
        outF.write("\n" + str(EnrScore) + "\t")
      # print ""
      # print EnrScore
      headerFlag = True
      clusterFlag = True
      entryCount = 1
  outF.write("\n")    
  inpF.close()
  outF.close()
  print("cluster_miner completed")
          
def chart_miner():
  """take a functional annotation chart output from DAVID and extract all terms 
  with a bonferroni p value above 0.01"""  
  
  import os.path
  
  print("this is chart_miner")
  inpF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "cav1ko", "processed", "chart_5979A7E9CBA41503507854395.txt"),"rU")
  outF = open(os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "cav1ko", "processed", "david_chart_5979A7E9CBA41503507854395.txt"),"w")
  """
  with open("bob/processed/DAVID_FACh_24hbobprots_up_18012016_terms.txt", "w") as outF:
    with open("bob/processed/DAVID_FACh_24hbobportos_up_18012016.txt", "r") as inpF:  
  """
  headerFlag = True
  for inpLine in inpF:
    if headerFlag:
      headerFlag = False
      outF.write("P value (bonferroni)\tterm\n")
      continue
    inpList = inpLine.split("\t")
    if float(inpList[-3]) <= 0.01:
      if "~" in inpList[1]:
        outF.write(str(float(inpList[-3])) + "\t" + inpList[1][inpList[1].index("~")+1:] + "\n")
      else:
        outF.write(str(float(inpList[-3])) + "\t" + inpList[1] + "\n")
   
    
  print("chart miner completed")

if __name__ == "__main__":
  main()