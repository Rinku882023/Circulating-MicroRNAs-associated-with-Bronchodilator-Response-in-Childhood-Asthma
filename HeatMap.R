library(ComplexHeatmap)
library(circlize)

x=c("Q","hsa.miR.143.3p","hsa.miR.378i","hsa.miR.26b.5p","hsa.miR.223.3p","hsa.miR.885.5p","hsa.miR.378a.3p","hsa.miR.500b.3p","hsa.miR.1290","hsa.miR.1246","hsa.miR.200b.3p","hsa.miR.877.5p")
heat=CRA_norm[x]
heat=heat[order(heat$Q),]
h=t(heat[2:12])
colnames(h)=heat$Q
row.names(h)=c("hsa-miR-143-3p" , "hsa-miR-378i"   , "hsa-miR-26b-5p" , "hsa-miR-223-3p" , "hsa-miR-885-5p" , "hsa-miR-378a-3p","hsa-miR-500b-3p", "hsa-miR-1290" ,   "hsa-miR-1246" ,   "hsa-miR-200b-3p", "hsa-miR-877-5p")

heatmap=Heatmap(h,clustering_distance_rows = "pearson",clustering_distance_columns = "pearson",col = colorRamp2(c(0,1,5,10,15,20), c("yellow", "green","blue","orange", "red","darkred")),column_split =  paste0(colnames(h)),row_split = 6,border=TRUE, row_dend_width = unit(20, "mm"), column_dend_height = unit(30, "mm"),show_column_names = FALSE,row_names_gp = gpar(fontface = "bold"),heatmap_legend_param = list(
    legend_direction = "vertical",title="log Norm Expression",title_position = "lefttop-rot",legend_height = unit(4, "cm")))
draw(heatmap, heatmap_legend_side = "left")
