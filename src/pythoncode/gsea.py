'''
Created on 11 Jan 2019

@author: mate

needs own virtual environment with gsea python package installed (gsea-p3)


1) convert mass spec analysis output file to the format of gene name, description, log2 lfq1, log2 lfq2..... log2 lfqn - done
2) create cls file as outlined here: http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats - file wasn't needed, only list. done.
3) create gene set file from vesiclepedia top 100 proteins and make it into a gmt file outlined in the link above
4) place them all in the code outlined on line 28 in this link: https://gseapy.readthedocs.io/en/master/gseapy_example.html
5) pray
'''
def gsea_data_maker(g1Names, g2Names, inpF):
  """create data file for GSEA from mass spec analysis output file 
  (made by proteingroups_analyzer.py and file_combiner within it.)"""
  
  import pandas as pd
  from math import log2
  
  pd.set_option("display.expand_frame_repr", False) # this prevents the splitting of dataframes to multiple rows upon printing
  pd.set_option("display.max_columns", 50)
  
  inpFile = open(inpF,"r")
  inpDF = pd.read_csv(inpFile)
  inpFile.close()
  
  # print(inpDF)
  
  colNameList = []
  orderedColNameList = []
  corrColL = []
  dataDF = pd.DataFrame()
  
  for colName in inpDF.columns:
    if colName.startswith("LFQ intensity"): colNameList.append(colName) # collect out all the needed column names. store them in a list
  
  groupCount = 1  
  for groupI in g1Names: # arrange list into list of lists by group
    for colI in colNameList: 
      if groupI in colI.upper(): 
        if len(orderedColNameList) < groupCount: orderedColNameList.append([])
        
        orderedColNameList[groupCount - 1].append(colI)
    groupCount += 1
    
  for groupI in g2Names: # arrange list into list of lists by group
    for colI in colNameList: 
      if groupI in colI: 
        if len(orderedColNameList) < groupCount: orderedColNameList.append([])
          
        orderedColNameList[groupCount - 1].append(colI)
    groupCount += 1
  
  for groupGroup in orderedColNameList:
    corrColL.extend(sorted(groupGroup,key = lambda b: int(b.split("-")[1].strip()))) # organize subgroups on number at the end and add them to a final list
    
  dataDF["NAME"] = inpDF["Gene names"].str.upper()
  dataDF["DESCRIPTION"] = inpDF["Protein names"]
  for colI in corrColL:
    dataDF[colI.split(" ")[-1]] = inpDF[colI].apply(lambda x: round(log2(int(round(x,0))),4)) # convert LFQ values to log scale
  
  dataDF.set_index("NAME", inplace = True)
  
  return dataDF


def gsea_cls_maker(g1Names, g2Names, testDF):
  """creates a class vector for GSEA analysis. Outlined here:
  https://gseapy.readthedocs.io/en/master/gseapy_example.html#4.-GSEA-Example"""
  
  clsL = []
  
  for g1I in g1Names: 
    for colI in testDF.columns: 
      if g1I in colI.upper(): clsL.append("EXO")

  for g2I in g2Names: 
    for colI in testDF.columns: 
      if g2I in colI: clsL.append("CTRL")  
      
  return clsL
  

def enricher_library_preparer_vesiclepedia():
  """append the exosome data from exosome datasets from vesiclepedia to the KEGG2016 collection."""
  
  exoFileL = [open("/home/mate/code/ed/src/data/cav1ko/vesiclepedia_results_1127.csv","r"),open("/home/mate/code/ed/src/data/cav1ko/vesiclepedia_results_1007.csv","r"), open("/home/mate/code/ed/src/data/cav1ko/vesiclepedia_results_192.csv","r")]
  nameL = ["Vesiclepedia_1127","Vesiclepedia_1007", "Vesiclepedia_192"]

  libFile = open("/home/mate/code/ed/src/data/cav1ko/KEGG_2016_EV.gmt","r")
  
  modFile = open("/home/mate/code/ed/src/data/cav1ko/KEGG_2016_EV_Vesi.gmt","w")
  
  for libLine in libFile:
    libL = libLine.split("\t")
    modFile.write(libL[0])
    modFile.write("\tNA\t")
    libS = "\t".join(libL[1:])
    modFile.write(libS)
  
  libFile.close()
  
  for i in range(len(exoFileL)):
    exoFile = exoFileL[i]
    next(exoFile) # skip header
    exoL = []
    for exoLine in exoFile:
      exoL.append(exoLine.split(",")[-1].strip()) # add upp exosome terms into a list
    
    exoFile.close()
    
    modFile.write(nameL[i] + "\tNA\t")
    for exoI in exoL[:-1]:
      modFile.write(exoI.upper())
      modFile.write("\t")
      
    modFile.write(exoL[-1].upper())
    modFile.write("\n")
    
    
  modFile.close()
  
def enricher_library_preparer():
  """append the exosome data from EV_TOP_100.txt - from vesiclepedia - to the KEGG2016 collection."""
  
  exoFile = open("/home/mate/code/ed/src/data/cav1ko/EV_TOP_100.txt","r")
  next(exoFile) # skip header
  exoL = []
  for exoLine in exoFile:
    exoL.append(exoLine.split("\t")[0].upper().strip()) # add upp exosome terms into a list
  
  exoFile.close()
  
  libFile = open("/home/mate/code/ed/src/data/cav1ko/KEGG_2016","r")
  
  modFile = open("/home/mate/code/ed/src/data/cav1ko/KEGG_2016_EV.gmt","w")
  
  for libLine in libFile:
    libL = libLine.split("\t")
    modFile.write(libL[0])
    modFile.write("\t")
    libS = "\t".join(libL[1:])
    modFile.write(libS)

    
  
  modFile.write("top 100 exosome proteins\t")
  for exoI in exoL[:-1]:
    modFile.write(exoI)
    modFile.write("\t")
    
  modFile.write(exoL[-1])
  modFile.write("\n")
  
  libFile.close()
  modFile.close()

  

def gsea_calculator(dataDF, geneSetO, clsList, idL):
  import gseapy as gp  
  from gseapy.plot import gseaplot, heatmap

  # run gsea
  # enrichr libraries are supported by gsea module. Just provide the name
  
  gs_res = gp.gsea(data = dataDF, # or data='./P53_resampling_data.txt'
                   gene_sets = geneSetO, # enrichr library names
                   cls = clsList, #'./data/P53.cls', # cls=class_vector
                   # set permutation_type to phenotype if samples >=15
                   # permutation_type='phenotype',
                   permutation_num=500, # reduce number to speed up test
                   max_size = 10000, # set max size of groups. larger ones are excluded
                   outdir="/home/mate/code/ed/src/data/cav1ko/processed/log/",  # do not write output to disk
                   no_plot=False, # Skip plotting
                   method="ratio_of_classes", # 'signal_to_noise',
                   processes=4,
                   format='png')

  
    
  for idI in idL:
  
    print(gs_res.res2d.loc[idI])
    gseaplot(gs_res.ranking, term=idI, ofname = "/home/mate/code/ed/src/data/cav1ko/processed/gseatest-" + idI +".png", **gs_res.results[idI])

    # plotting heatmap
    genes = gs_res.res2d.genes[idI].split(";")
    heatmap(df = gs_res.heatmat.loc[genes], z_score=0, ofname = "/home/mate/code/ed/src/data/cav1ko/processed/gseatestheatmap-" + idI + ".png", title=idI, figsize=(18,6))


  
  
def gsea_maker():  
  import gseapy as gp
  print(gp.__version__)
  names = gp.get_library_name()
  print(names[:100])
  
  phenoA, phenoB, class_vector =  gp.parser.gsea_cls_parser("/home/mate/code/ed/src/data/gsea/P53.cls")
  
  print(phenoA, phenoB, class_vector)
  print(len(class_vector))
  
  
  # run enrichr
  # if you are only intrested in dataframe that enrichr returned, please set no_plot=True
  
  # list, dataframe, series inputs are supported
  enr = gp.enrichr(gene_list="/home/mate/code/ed/src/data/gsea/gene_list.txt",
                   # or gene_list=glist
                   description='test_name',
                   gene_sets='KEGG_2016',
                   # or gene_sets='KEGG_2016,KEGG_2013',
                   # gene_sets=['KEGG_2016','KEGG_2013'],
                   outdir='test/enrichr_kegg',
                   cutoff=0.5, # test dataset, use lower value of range(0,1)
                   no_plot=False
                  )
  
  
  print(enr.results)

def single_gsea(dataFrame, gene_sets, outdir="ssGSEA_", sample_norm_method='rank', min_size=15, max_size=2000,
                permutation_num=0, weighted_score_type=0.25, scale=True, ascending=False, processes=4,
                figsize=(7,6), graph_num=20, no_plot=False, seed=None, verbose=False):
  """single sample gene set enrichment analysis. see package documentation for details, or here: https://gseapy.readthedocs.io/en/master/_modules/gseapy/gsea.html#ssgsea"""
  import gseapy as gp  
  from gseapy.plot import gseaplot, heatmap
  
  resO = gp.ssgsea(dataFrame, gene_sets = gene_sets, outdir = outdir, permutation_num=0)
  print(resO.res2d)
  print(resO.res2d.loc["Vesiclepedia_1007"])
  # gseaplot(resO.ranking, term="Vesiclepedia_1007", ofname = "/home/mate/code/ed/src/data/cav1ko/processed/ssgseatest-" + "Vesiclepedia_1007" +".png", **resO.results["Vesiclepedia_1007"])
  
  print(resO.resultsOnSamples)
  heatmap(resO.res2d, ofname = "/home/mate/code/ed/src/data/cav1ko/processed/ssgseaheatmap2.png", figsize = (20,40), fontN = 9)
  
def gmt_file_converter():
  """remove _homo sapiens_hsaxxx tags from group names to get shorter descriptions"""  
  
  with open("/home/mate/code/ed/src/data/cav1ko/KEGG_2016_EV_Vesi.gmt","r") as inF:
    with open("/home/mate/code/ed/src/data/cav1ko/KEGG_2016_EV_Vesi_short.gmt","w") as outF:
      for inLine in inF:
        inL = inLine.split("_Homo sapiens_")
        if len(inL) > 1:
          outF.write(inL[0])
          outF.write(inL[1][8:])
        else: outF.write(inLine)
        

def main_gsea_function():  
  """call all the other GSEA functions and run the whole pipeline to perform the analysis"""
  
  group1Names =  ["WT","CAV1KO"] # ["WT","CAV1KO","IL7"]
  group2Names = ["IL7"]
  inFile = "/home/mate/code/ed/src/data/cav1ko/processed/proteinGroups_EV_matched_samples_24-11-2018_exosome_Cav1ko_vs_wt_all_datasets_26-11-2018_combined_2019-01-09-1.csv"
  
  IDList = ["top 100 exosome proteins","Vesiclepedia_1007","Vesiclepedia_192"]
    
  datasetDF = gsea_data_maker(g1Names = group1Names, g2Names = group2Names, inpF = inFile)
  clsList = gsea_cls_maker(g1Names = group1Names, g2Names = group2Names, testDF = datasetDF)
  
  print(datasetDF)
  
  print(clsList)
  
  gmt_file_converter()
  
  geneSetName = "/home/mate/code/ed/src/data/cav1ko/KEGG_2016_EV_Vesi_short.gmt"
  
  single_gsea(dataFrame = datasetDF, gene_sets = geneSetName, outdir = "/home/mate/code/ed/src/data/cav1ko/processed/gsea/")
    
  # gsea_calculator(dataDF = datasetDF, geneSetO = geneSetName, clsList = clsList, idL = IDList)
  
  
  
main_gsea_function()

# enricher_library_preparer_vesiclepedia()