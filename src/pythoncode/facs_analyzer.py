'''
Created on 24 Oct 2019

@author: mate

take in .csv file with gated flow cytometry data, and convert it to a pandas dataframe object. Then analyze samples and output graphs comparing the effect of Pstpip1 KO in different Ptpn22 and CD2 backgrounds.
To do this, take in one channel at a time, normalize based on the unstimulated sample, then plot stimulation curve pairs of WT VS Pstpip1 KO pairs with Ptpn22 presence or absence, and with CD2 presence or absence.
'''

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

pd.set_option("display.expand_frame_repr", False) # this prevents the splitting of dataframes to multiple rows upon printing
pd.set_option("display.max_columns", 50)

inpDF = pd.read_csv("/home/mate/code/ed/src/data/facs/CD2_KO_VS_WT.csv")

rawDF = inpDF.rename({"cells/Single Cells/live cells/CD2+ | Geometric Mean (B1-A)" : "CD2+",
                      "cells/Single Cells/live cells/CD2+ | Geometric Mean (B2-A)" : "CD69+",
                      "cells/Single Cells/live cells/CD2+ | Geometric Mean (B3-A)": "CD8+",
                      "cells/Single Cells/live cells/CD2+ | Geometric Mean (B4-A)": "CD71+",
                      "cells/Single Cells/live cells/CD2+ | Geometric Mean (V2-A)": "LD+",
                      "cells/Single Cells/live cells/CD2+ | Geometric Mean (V1-A)": "CD25+",
                      "cells/Single Cells/live cells/CD2+ | Geometric Mean (R1-A)": "CD44+",
                      "cells/Single Cells/live cells/CD2+ | Geometric Mean (R2-A)": "TCR+",
                      "cells/Single Cells/live cells/CD2+ | Geometric Mean (FSC-A)": "FSC+",
                      "cells/Single Cells/live cells/CD2+ | Geometric Mean (SSC-A)": "SSC+",
                      "cells/Single Cells/live cells/CD2- | Geometric Mean (B1-A)" : "CD2-",
                      "cells/Single Cells/live cells/CD2- | Geometric Mean (B2-A)" : "CD69-",
                      "cells/Single Cells/live cells/CD2- | Geometric Mean (B3-A)": "CD8-",
                      "cells/Single Cells/live cells/CD2- | Geometric Mean (B4-A)": "CD71-",
                      "cells/Single Cells/live cells/CD2- | Geometric Mean (V2-A)": "LD-",
                      "cells/Single Cells/live cells/CD2- | Geometric Mean (V1-A)": "CD25-",
                      "cells/Single Cells/live cells/CD2- | Geometric Mean (R1-A)": "CD44-",
                      "cells/Single Cells/live cells/CD2- | Geometric Mean (R2-A)": "TCR-",
                      "cells/Single Cells/live cells/CD2- | Geometric Mean (FSC-A)": "FSC-",
                      "cells/Single Cells/live cells/CD2- | Geometric Mean (SSC-A)": "SSC-",
                      "Sample:": "Sample"
                      }, axis = 1)

print(rawDF.head())

wtDF = rawDF.loc[(rawDF["Sample"].str.endswith("0001.fcs") & (rawDF["Sample"].str.contains("H") == False)) | rawDF["Sample"].str.endswith("A1.0002.fcs")] # select WT samples (first plate)
wtDF = wtDF.reindex([wtDF.columns[0]] + sorted(wtDF.columns[1:], key = lambda x: x[-1]), axis=1)
wtDF = wtDF.reindex([wtDF.columns[0]] + list(wtDF.columns[1:(int(len(wtDF.columns)/2) + 1)]) + list(wtDF.columns[int(len(wtDF.columns)/2) + 1:])[::-1], axis=1)

# print(wtDF)


rwDF = rawDF.loc[((rawDF["Sample"].str.endswith("0002.fcs") & (rawDF["Sample"].str.contains("A1") == False)) & (rawDF["Sample"].str.contains("H") == False)) | rawDF["Sample"].str.endswith("A1.0003.fcs")] # select R619W samples (second plate)
rwDF = rwDF.reindex([rwDF.columns[0]] + sorted(rwDF.columns[1:], key = lambda x: x[-1]), axis=1)
rwDF = rwDF.reindex([rwDF.columns[0]] + list(rwDF.columns[1:(int(len(rwDF.columns)/2) + 1)]) + list(rwDF.columns[int(len(rwDF.columns)/2) + 1:])[::-1], axis=1)


pepkoDF = rawDF.loc[((rawDF["Sample"].str.endswith("0003.fcs") & (rawDF["Sample"].str.contains("A1") == False)) & (rawDF["Sample"].str.contains("H") == False)) | rawDF["Sample"].str.endswith("A1.0004.fcs")] # select PEPKO samples (third plate)
pepkoDF = pepkoDF.reindex([pepkoDF.columns[0]] + sorted(pepkoDF.columns[1:], key = lambda x: x[-1]), axis=1)
pepkoDF = pepkoDF.reindex([pepkoDF.columns[0]] + list(pepkoDF.columns[1:(int(len(pepkoDF.columns)/2) + 1)]) + list(pepkoDF.columns[int(len(wtDF.columns)/2) + 1:])[::-1], axis=1)

print(pepkoDF)

facsL = [wtDF, rwDF, pepkoDF]
cd2L = []

for facsDF in facsL:
  cd2L.append(facsDF.loc[facsDF["Sample"].str.contains("B|D|F|G")])


cd2WTL = []
cd2KOL = []
skL = []

for cd2I in cd2L:
  cd2SKWT = cd2I.loc[cd2I["Sample"].str.contains("B"),:][cd2I.columns[1:]].iloc[:-1,:].append(cd2I.loc[cd2I["Sample"].str.contains("F"),:][cd2I.columns[1:]].iloc[5:,:]).append(cd2I.loc[cd2I["Sample"].str.contains("B"),:][cd2I.columns[1:]].iloc[-1,:]) 
  # Pstpip1 WT
  cd2SKWTNorm = cd2SKWT.subtract(cd2SKWT.iloc[-1,:]).iloc[:-1,:] # normalize to unstimulated
  cd2SKWTNorm[cd2SKWTNorm < 0] = 1
  # cd2SKWTNorm /= cd2SKWTNorm.max() # scale
  print(cd2SKWTNorm)
  print()
#   sns.heatmap(cd2SKWTNorm[["FSC+","SSC+","FSC-","SSC-"]].T, fmt = ".0f", cmap="YlGnBu", annot=True, xticklabels = False, vmin = 0, vmax = 50000)
#   plt.show()
#   sns.heatmap(cd2SKWTNorm.drop(columns = ["FSC+","FSC-","SSC+","SSC-", "LD+", "LD-", "CD2-"]).T, fmt = ".0f", cmap="YlGnBu", annot=True, xticklabels = False, vmin = 0, vmax = 20000)
#   plt.show()
  
#   
  
  cd2SKWTNorm_cd2RatioDF =  pd.DataFrame()
  for colN in cd2SKWTNorm.columns[:int(len(cd2SKWTNorm.columns)/2) + 1]:
    colN = colN[:-1]
    cd2SKWTNorm_cd2RatioDF[colN] = cd2SKWTNorm[colN + "+"] - cd2SKWTNorm[colN + "-"]
  
  sns.heatmap(cd2SKWTNorm_cd2RatioDF.drop(columns = ["FSC","SSC", "LD"]).T, fmt = ".0f", cmap="Oranges", annot=True, xticklabels = False, vmin = 0, vmax = 8000)
  plt.show()    
    
  cd2WTL.append(cd2SKWTNorm_cd2RatioDF.T)
#   print(cd2SKWTNorm_cd2RatioDF)
#   print()
  
  # Pstpip1 KO
  cd2SKKO = cd2I.loc[cd2I["Sample"].str.contains("D"),:][cd2I.columns[1:]].iloc[:-1,:].append(cd2I.loc[cd2I["Sample"].str.contains("G"),:][cd2I.columns[1:]].iloc[5:,:]).append(cd2I.loc[cd2I["Sample"].str.contains("D"),:][cd2I.columns[1:]].iloc[-1,:]) 
  cd2SKKONorm = cd2SKKO.subtract(cd2SKKO.iloc[-1,:]).iloc[:-1,:] # normalize to unstimulated
  cd2SKKONorm[cd2SKKONorm < 0] = 1
  # cd2SKKONorm /= cd2SKKONorm.max()
  
  print(cd2SKKONorm)
  print("=============================")
#   sns.heatmap(cd2SKKONorm[["FSC+","SSC+","FSC-","SSC-"]].T, fmt = ".0f", cmap="YlGnBu", annot=True, xticklabels = False, vmin = 0, vmax = 50000)
#   plt.show()
#   sns.heatmap(cd2SKKONorm.drop(columns = ["FSC+","FSC-","SSC+","SSC-", "LD+", "LD-", "CD2-"]).T, fmt = ".0f", cmap="YlGnBu", annot=True, xticklabels = False, vmin = 0, vmax = 20000)
#   plt.show()

  
  cd2SKKONorm_cd2RatioDF =  pd.DataFrame() # effect of CD2 presence in Pstpip1 WT cells
  for colN in cd2SKKONorm.columns[:int(len(cd2SKKONorm.columns)/2) + 1]:
    colN = colN[:-1]
    cd2SKKONorm_cd2RatioDF[colN] = cd2SKKONorm[colN + "+"] - cd2SKKONorm[colN + "-"]

  sns.heatmap(cd2SKKONorm_cd2RatioDF.drop(columns = ["FSC","SSC", "LD"]).T, fmt = ".0f", cmap="Oranges", annot=True, xticklabels = False, vmin = 0, vmax = 8000)
  plt.show()  
    
  cd2KOL.append(cd2SKKONorm_cd2RatioDF.T)
#   print(cd2SKWTNorm_cd2RatioDF)
#   print("====================")
  
  cd2SKKONorm.index = cd2SKWTNorm.index
  
  skL.append(cd2SKWTNorm.subtract(cd2SKKONorm).T)
  

for skI in skL:
  print(skI)
#   sns.heatmap(skI.T.drop(columns = ["FSC+","FSC-","SSC+","SSC-", "LD+", "LD-", "CD2-"]).T, fmt = ".0f", cmap=sns.color_palette("RdBu_r",200), vmin = -3000, vmax = 5000, center = 0, annot=True, xticklabels = False)
#   plt.show()
  
  
  


  

