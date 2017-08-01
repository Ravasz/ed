'''
Created on 7 Aug 2015

@author: mate

open peptides.txt and find peptides in there corresponding to the gene name specified in the queryS variable.
Then find the protein refseq GI using db2db
Then, download the full protein sequence from Entrez.
Then, mark the found peptides in the protein sequence and output the resulting marked sequence in a html file as protein name.html
Then, use the phosphorylation site dataset in the datafiles folder which is the complete phosphosite.org phosphosite map and mark phosphosites on the HTML
'''


def main():
  from ed.tools import html_creator, prot_id_converter, prot_entrez_fetch
  from copy import deepcopy
  print "this is peptide parser"
  
  queryS = "Eif5a" # this is the search term that will be worked on. It should be a protein name like "Ptpn22"
  
  
  print "working on: " + queryS
  with open("../bob/peptides.txt","r") as inpF: 
    pepL, uniId = peptide_finder(targetFile=inpF, targetS = queryS) # find peptides from peptides.txt for the protein name
  print pepL
  print "peptides found"
  targetL = [queryS]
  idList = prot_id_converter(targetL, "10090", inpDB = "genesymbol",outDB="refseqproteingi")
  seqL = prot_entrez_fetch(idList, retM="gb", retT="fasta")
  for seqItem in seqL:
    seqS = seqItem.split("\n")[1]
    print seqS
  print "protein sequence found"
  annotS = deepcopy(seqS)
  pStartL = []
  pEndL = []
  for pepItem in pepL: # locate peptides in full protein sequence and store positions for starts and ends. merge overlapping peptides.
    pepStart = seqS.index(pepItem)
    pepEnd = pepStart + len(pepItem)
    
    startCount = 0 # handle starts
    for startItem in pStartL:
      if startItem <= pepStart:
        startCount += 1
    endCount = 0
    for endItem in pEndL:
      if endItem <= pepStart:
        endCount += 1
    if startCount == endCount and pepStart not in pEndL: # start new peptide
      pStartL.append(pepStart)
    elif startCount == endCount and pepStart in pEndL: # start a new peptide at the end of another peptide
      pEndL.remove(pepStart)
    
    overlapCount = 0
    for startItem in pStartL[:]: # handle ends
      if pepStart<startItem<=pepEnd:
        pStartL.remove(startItem) # remove extra starts
        overlapCount += 1
        
    for endItem in pEndL[:]:
      if pepStart<=endItem<=pepEnd: # remove extra ends
        pEndL.remove(endItem)
        overlapCount -= 1
    
    if pepEnd not in pEndL and overlapCount <= 0: # add end
      curStart = 500000
      for pSI in pStartL:
        if curStart > pSI > pepEnd:
          curStart = pSI
      curEnd = 500000
      for pEI in pEndL:
        if curEnd > pEI > pepEnd:
          curEnd = pEI
      if curStart <= curEnd:  # check if next tag is start or end. if start, add end. if end, do nothing
        pEndL.append(pepEnd)
     
  print uniId
  phL = []
  with open("../datafiles/Phosphorylation_site_dataset") as phInp: # now for the phosphosite data
    for phLine in phInp:
      phList = phLine.split("\t")
      try:
        if uniId == phList[1]: 
          phL.append(int(phList[4][1:-2]) - 1)
      except IndexError:
        continue
  print phL
  
  fullL = pStartL + pEndL
  for phItem in phL:
    if phItem not in fullL: fullL.append(phItem)
  fullL.sort()
  
  offsetN = 0
  for posI in fullL: # from the resulting intervals, create emphasis in html file, mark phosphosites in red
    if posI in pStartL:
      annotS = annotS[:posI+offsetN] + "<mark>" + annotS[posI+offsetN:]
      offsetN += 6
    elif posI in pEndL:
      annotS = annotS[:posI+offsetN] + "</mark>" + annotS[posI+offsetN:]
      offsetN += 7
    if posI in phL:
      annotS = annotS[:posI+offsetN] + r"""<strong style="color: red;">""" + annotS[posI+offsetN] + r"""</strong>""" + annotS[posI+offsetN + 1 :]
      offsetN += 37
        
  print annotS
  html_creator(queryS + " peptides", annotS, queryS + ".html")
  print "found peptides marked in the file: ",
  print queryS + ".html"  


  
  
  
def peptide_finder(targetFile, targetS = "Ptpn22"):
  """Find all peptides which are matched to the targetS gene name in peptides.txt
  return peptide sequences as a list of strings and uniprot ID as string. The two together work as a tuple.
  """
  headerFlag = True
  peptideL = []
  uniprotId = ""
  for lineS in targetFile:
    if headerFlag:
      headerFlag = False
      continue
    lineL = lineS.split("\t")
    geneNames = lineL[37].split(";")
    if geneNames[0] == "": continue
    if geneNames[0] == targetS:
      if lineL[0] not in peptideL:
        peptideL.append(lineL[0])
        uniprotId = lineL[34]
        if "-" in uniprotId:
          uniprotId = uniprotId[:uniprotId.index("-")]
  return peptideL, uniprotId
  
if __name__ == "__main__":
  main()