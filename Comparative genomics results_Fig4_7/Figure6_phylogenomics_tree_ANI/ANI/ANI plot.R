setwd("/Users/sihanbu/Library/CloudStorage/OneDrive-UniversityofConnecticut/Research/whole\ genome\ assembly/2026Mar04_Biolog_Analysis_new_data/selected\ genomes\ comparision/01_PANGENOME/New\ with\ 24\ genomes/ANI")

ANI<- read.delim("ANIb_percentage_identity.txt",stringsAsFactors = FALSE)
ANI<- as.data.frame(ANI)
ANI<-column_to_rownames(ANI, var = "key")

name_map <- c(
  "Ho04_635" = "Ho04-635",
  "b1_3" = "b1-3",
  "Hv98_241" = "Hv98-241",
  "Hv98_571" = "Hv98-571",
  "AV1708_29120" = "1708-29120",
  "AV183026" = "183026",
  "Hv98_231" = "Hv98-231",
  "Hv98_232" = "Hv98-232",
  "Hv98_225" = "Hv98-225")

# Rename row names
idx <- rownames(ANI) %in% names(name_map)
rownames(ANI)[idx] <- name_map[rownames(ANI)[idx]]

# Rename column names
idx <- colnames(ANI) %in% names(name_map)
colnames(ANI)[idx] <- name_map[colnames(ANI)[idx]]

my_levels <- c("AV0110","Hv98-225", "AV040", 
               "FDAARGOS_632", "1708-29120",  "b1-3", "KMM358", "KMM366", "C4","KMM289",
               "Hm22","Hv98-571","Hv98-231", "Hv98-232","Hm21" ,"Ho04-635",
               "Hv98-241",  "FC951","D3","ANYA_18196","183026",
               "ZfB1","ANYA_13997", "ANYA_14430")

 
ANI <- ANI[my_levels, my_levels]


pheatmap(
  ANI,
  color = colorRampPalette(c("white", "#D55E00"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey90",
  fontsize_row = 10,
  fontsize_col = 8)


#### claude assignment
 
library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  seq(min(ANI), max(ANI), length.out = 100),
  colorRampPalette(c("white", "#D55E00"))(100)
)

aniheatmap<-Heatmap(
  as.matrix(ANI),
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = gpar(col = "grey90", lwd = 0.5),   # border_color equivalent
  row_names_side = "left",                      # y labels on LEFT
  column_names_side = "top",                    # x labels on TOP
  row_names_gp = gpar(fontsize = 10),           # fontsize_row
  column_names_gp = gpar(fontsize = 10),         # fontsize_col                
  heatmap_legend_param = list(title = "ANI"))



#export as a svg file
svglite("ANI_heatmap.svg",width = 6,height = 5)
print(aniheatmap)      
dev.off()