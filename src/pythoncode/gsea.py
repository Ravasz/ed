'''
Created on 11 Jan 2019

@author: mate

needs own virtual environment with gsea python package installed (gsea-p3)


1) convert mass spec analysis output file to the format of gene name, description, log2 lfq1, log2 lfq2..... log2 lfqn
2) create cls file as outlined here: http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
3) create gene set file from vesiclepedia top 100 proteins and make it into a gmt file outlined in hte link above
4) place them all in the code outlined on line 28 in this link: https://gseapy.readthedocs.io/en/master/gseapy_example.html
5) pray
'''
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

pd.set_option("display.expand_frame_repr", False) # this prevents the splitting of dataframes to multiple rows upon printing
pd.set_option("display.max_columns", 50)

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
                 no_plot=True
                )


print(enr.results)