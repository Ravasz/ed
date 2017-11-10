'''
Created on 14 May 2015

@author: mate

explore evidence.txt from Bob's mass spec data 
and look for different modified peptides, gene names and such.

evidence.txt seems to be a list of all peptides
with their modifications and identification scores listed

take all peptides, the protein they are mapped to, and the fraction of isoelectric focusing they were found in, 
and using the protein sequence information found in mouse_proteome.fasta, 
calculate if the isoelectric point of the protein matches the fraction where each peptide was found.

After running the program, it seems that 380624 matching, and 1058049 unmatching peptides were found.
'''  

def main():
  from .tools import file_importer, uniprot_dicter, iso_e
  from collections import defaultdict
  print("this is evidence parser")
  inpEvF = file_importer("bob/evidence.txt")
  seqD = uniprot_dicter()
  resD = defaultdict(list)
  for keyS, valueS in list(seqD.items()): # create dict with keys as all mouse genes and values as a list, with the protein sequence as first item in the list
    resD[keyS].append(valueS.strip())
  geneL = []
  goodC = 0
  badC = 0
  for inpLine in inpEvF:
    inpList = inpLine.split("\t")
    inpL = inpList[14].split(";") # this is the leading protein name column. It holds uniprot IDs separated by ; 
    try:
      inpFrac = int(inpList[20]) # this is the "fraction" column. it holds the isoelectric fractionation number
    except ValueError:
      print(inpList[20])
      continue
    for inpI in inpL:
      if inpI[-2] == "-":
        inpI = inpI[:-2]
      if inpI not in geneL: # make a list of gene names
        geneL.append(inpI)
      if inpI[:3] == "REV" or inpI[:3] == "CON": continue # there are some weird entries. CON is contaminants methinks.
      try:
        if inpFrac - 2 < iso_e(resD[inpI][0]) < inpFrac + 2 : # check if the isoelectric point of the protein to which this peptide belongs to matches the fraction it was retrieved from
          goodC += 1
        else:
          badC += 1
          
        """resD[inpI].append()
        resD[inpI].append(iso_e(resD[inpI][0])) # add isoelectric point to the genes found"""
      except IndexError:
        print("%s has no sequence" % inpI)
  print(str(len(geneL)) + " genes found in evidence.txt")
  print(str(goodC) + " matching, and " + str(badC) + " unmatching peptides found") 
  inpEvF.close()  


  
if __name__ == "__main__":
  main()
  