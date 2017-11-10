'''
Created on 18 May 2015

@author: mate
parse the complete mouse proteome from the original uniprot database named UP000000589_10090.fasta (29.04.2015).
This file was downloaded from:
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/

Parse into a slightly different, but more desirable format and output to mouse_proteome.fasta. 

The result should look like this:

>tr|A0A075B5I0|A0A075B5I0_MOUSE Uncharacterized protein OS=Mus musculus GN=Ighv1-13 PE=4 SV=1
SLTIIVSSVLRCQASELHVGCVGGCVYSKCGLAFELLIVGNTTIFRMNSSNPFQSLLQTLLHPLNFIVIKIVARSLAGQLHLSPRYYQLWSRLLKLDL
>tr|A0A075B5I2|A0A075B5I2_MOUSE Protein Trbv4 (Fragment) OS=Mus musculus GN=Trbv4 PE=4 SV=1
MGCRLLSCVAFCLLGIGPLETAVFQTPNYHVTQVGNEVSFNCKQTLGHDTMYWYKQDSKKLLKIMFSYNNKQLIVNETVPRRFSPQSSDKAHLNLRIKSVEPEDSAVYLCASS
>tr|A0A075B5J4|A0A075B5J4_MOUSE Protein Trbc2 (Fragment) OS=Mus musculus GN=Trbc2 PE=4 SV=1
XDLRNVTPPKVSLFEPSKAEIANKQKATLVCLARGFFPDHVELSWWVNGKEVHSGVSTDPQAYKESNYSYCLSSRLRVSATFWHNPRNHFRCQVQFHGLSEEDKWPEGSPKPVTQNISAEAWGRADCGITSASYHQGVLSATILYEILLGKATLYAVLVSGLVLMAMVKKKNS
>tr|A0A075B5J5|A0A075B5J5_MOUSE Protein Trbv31 OS=Mus musculus GN=Trbv31 PE=4 SV=1
MLYSLLAFLLGMFLVSAQTIHQWPVAEIKAVGSPLSLGCTIKGKSSPNLYWYWQATGGTLQQLFYSITVGQVESVVQLNLSASRPKDDQFILSTEKLLLSHSGFYLCAWSL

this script was used only once and is not needed now

'''
def main():
  
  from .tools import file_importer, file_outporter
  
  relPath = "datafiles/UP000000589_10090_additional.fasta"
  relOut = "datafiles/mouse_proteome.fasta"
  
  inpF = file_importer(relPath)
  outF = file_outporter(relOut)
  # 13628909

  seqS = ""
  for inpLine in inpF:
    inpS = inpLine.strip()
    if inpS[0] == ">":
      if seqS == "":
        outF.write(inpS + "\n")
        continue
      outF.write(seqS + "\n")
      outF.write(inpS + "\n")
      seqS = ""
    else: 
      seqS = seqS + inpS
  
    
    
  inpF.close()

if __name__ == "__main__":
  main()

    