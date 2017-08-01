#instal SPADE:
R.Version()
install.packages("devtools")
devtools::install_github("nolanlab/Rclusterpp")
source("http://bioconductor.org/biocLite.R")
devtools::install_github("nolanlab/spade")
devtools::install_github("RGLab/flowCore",ref='trunk')
library(spade)

#download sample dataset and run sample analysis:
setwd("/home/mate/workspace/katamari/src/root/ed/Rcode")

markers <- c("Cd(110,111,112,114)","Cell_length","Dy(163.929)-Dual","Er(165.930)-Dual","Er(166.932)-Dual","Er(167.932)-Dual","Er(169.935)-Dual","Eu(150.919)-Dual","Eu(152.921)-Dual","Gd(155.922)-Dual","Gd(157.924)-Dual","Gd(159.927)-Dual","Ho(164.930)-Dual","In(114.903)-Dual","Ir(190.960)-Dual","La(138.906)-Dual","Lu(174.940)-Dual","Nd(141.907)-Dual","Nd(143.910)-Dual","Nd(144.912)-Dual","Nd(145.913)-Dual","Nd(147.916)-Dual","Nd(149.920)-Dual","Pr(140.907)-Dual","Sm(146.914)-Dual","Sm(151.919)-Dual","Sm(153.922)-Dual","Tb(158.925)-Dual","Tm(168.934)-Dual","Yb(170.936)-Dual","Yb(171.936)-Dual","Yb(173.938)-Dual","Yb(175.942)-Dual")
PANELS <- list(list(panel_files=c("Bendall_et_al_Science_2011_Marrow_1_SurfacePanel_Live_CD44pos_Singlets.fcs"), median_cols=NULL,reference_files=c("Bendall_et_al_Science_2011_Marrow_1_SurfacePanel_Live_CD44pos_Singlets.fcs"),fold_cols=c()))

SPADE.driver("Bendall_et_al_Science_2011_Marrow_1_SurfacePanel_Live_CD44pos_Singlets.fcs", out_dir="output", cluster_cols=markers, panels=PANELS, transforms=flowCore::arcsinhTransform(a=0, b=0.2), layout=SPADE.layout.arch, downsampling_target_percent=0.1, downsampling_target_number=NULL, downsampling_target_pctile=NULL, downsampling_exclude_pctile=0.01, k=200, clustering_samples=50000)

#generate PDFs from analysis:
layout <- read.table(file.path("output","layout.table"))
mst <- read.graph(file.path("output","mst.gml"),format="gml")
SPADE.plot.trees(mst,"output",file_pattern="*fcs*Rsave",layout=as.matrix(layout),out_dir=file.path("output","pdf"),size_scale_factor=1.2)




#2nd test with different layout (kamada.kawai, much better):

#download sample dataset and run sample analysis:
setwd("/home/mate/workspace/katamari/src/root/ed/Rcode")

markers <- c("Cd(110,111,112,114)","Cell_length","Dy(163.929)-Dual","Er(165.930)-Dual","Er(166.932)-Dual","Er(167.932)-Dual","Er(169.935)-Dual","Eu(150.919)-Dual","Eu(152.921)-Dual","Gd(155.922)-Dual","Gd(157.924)-Dual","Gd(159.927)-Dual","Ho(164.930)-Dual","In(114.903)-Dual","Ir(190.960)-Dual","La(138.906)-Dual","Lu(174.940)-Dual","Nd(141.907)-Dual","Nd(143.910)-Dual","Nd(144.912)-Dual","Nd(145.913)-Dual","Nd(147.916)-Dual","Nd(149.920)-Dual","Pr(140.907)-Dual","Sm(146.914)-Dual","Sm(151.919)-Dual","Sm(153.922)-Dual","Tb(158.925)-Dual","Tm(168.934)-Dual","Yb(170.936)-Dual","Yb(171.936)-Dual","Yb(173.938)-Dual","Yb(175.942)-Dual")
PANELS <- list(list(panel_files=c("Bendall_et_al_Science_2011_Marrow_1_SurfacePanel_Live_CD44pos_Singlets.fcs"), median_cols=NULL,reference_files=c("Bendall_et_al_Science_2011_Marrow_1_SurfacePanel_Live_CD44pos_Singlets.fcs"),fold_cols=c()))

SPADE.driver("Bendall_et_al_Science_2011_Marrow_1_SurfacePanel_Live_CD44pos_Singlets.fcs", out_dir="output", cluster_cols=markers, panels=PANELS, transforms=flowCore::arcsinhTransform(a=0, b=0.2), layout=igraph::layout.kamada.kawai, downsampling_target_percent=0.1, downsampling_target_number=NULL, downsampling_target_pctile=NULL, downsampling_exclude_pctile=0.01, k=200, clustering_samples=50000)

#generate PDFs from analysis:
layout <- read.table(file.path("output","layout.table"))
mst <- read.graph(file.path("output","mst.gml"),format="gml")
SPADE.plot.trees(mst,"output",file_pattern="*fcs*Rsave",layout=as.matrix(layout),out_dir=file.path("output","pdf"),size_scale_factor=1.2)
