'''
Created on 18 May 2016

@author: mate

compare protein abundance estimates calculated using the proteomic ruler method in Bob's dataset 
to the ones published in the Hukelmann dataset that were made with the same method

So it appears just comparing the calculated protein amounts per cell won't work as Jens loaded smaller amounts than Bob. 

Also, I don't have all the information from Bob's 24h dataset, I should get that still from Sara.

Bob's CTL dataset: 
16 fractions per sample, each fraction with 1.2ug of peptides, so each sample has 19.2ug peptides.

Jens' dataset:
12 fraction per sample, each fraction with 1ug of peptides, so each sample has 12ug peptides.  
'''

# open hukelmann data
with open("/home/mate/workspace/katamari/src/root/ed/datafiles/ni.3314-f1.txt", "rU") as inpHukF:
  # open Bob's data
  with open("/home/mate/workspace/katamari/src/root/ed/bob/processed/abs_quant_CTL-21-55-24-05-20162.txt","r") as inpBobF:
    # abs_quant_24H_TCR-17-21-19-05-20162.txt
    # abs_quant_CTL-17-22-19-05-20162.txt
    
    bobD = {} # contains all lines in the protein copy number estimation dataset with keys as gene names, 
    # and values as a list of everything
    commonD = {}
    correctD = {}
    matchCount = 0
    
    for bobLine in inpBobF:
      bobL = bobLine.split(",")
      bobL[-1] = bobL[-1].rstrip("\n")
      bobD[bobL[3]] = bobL
    
    countN = 0
    for inpHukLine in inpHukF:

      countN += 1
      # if countN == 20: break
      inpHukL = inpHukLine.split("\t")
      inpHukL[-1] = inpHukL[-1].rstrip("\n")
      
      """
      if len(inpHukL) != 72: 
        print inpHukL
        print len(inpHukL)
      # test if splitting of lines works as intended - it does
      """
      
      hukI = inpHukL[1].split(";")
          
      for geneI in hukI:
        if geneI in bobD and geneI != "Gene names":
          if geneI == "Hist1h4a": 
            print inpHukL
            print inpHukL[8]
            print int(bobD[geneI][-2])*0.5
            print int(bobD[geneI][-2])*2
          matchCount += 1       
          commonD[geneI] = bobD[geneI] + inpHukL           
          if inpHukL[8] != "" and float(bobD[geneI][-1])*0.5 < float(inpHukL[8]) < float(bobD[geneI][-1])*2: 
            # this is the important bit where the comparison is made
            correctD[geneI] = commonD[geneI]
            # print commonD[geneI]
          break
          
    print len(bobD)
    print countN
    print matchCount
    print "number of proteins within the acceptable range:", len(correctD)
    
    outF = open("/home/mate/workspace/katamari/src/root/ed/bob/processed/abs_quant_CTL-21-55-24-05-20163.txt","w")
    
    for outDV in sorted(commonD.items(), key=lambda x:int(x[1][0])): # sort the dict based on the ID they have
      for outI in outDV[1][:-1]:
        outF.write(str(outI) + ",")
      outF.write(str(outDV[1][-1]) + "\n")
    
    