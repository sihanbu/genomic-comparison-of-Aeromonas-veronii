setwd("/Users/sihanbu/Library/CloudStorage/OneDrive-UniversityofConnecticut/Research/whole\ genome\ assembly/2026Mar04_Biolog_Analysis_new_data/selected\ genomes\ comparision/AMR_and_VF")

VF <- read.delim("vfdb_all.tab", check.names = FALSE) 
group <- read.csv("Aeromonas_VFs_grouping.csv", header = T)

AMR <- read.delim("ncbi_amr_all.tab", check.names = FALSE) 
AMR$genome_name <- sub("\\.fasta$", "", basename(AMR$`#FILE`))
#write.csv(AMR, "AMR_info.csv", row.names = F)

#remove Hv96-22.fasta
VF<-filter(VF, `#FILE`  != "Hv96-22.fasta")
AMR<-filter(AMR, `#FILE`  != "Hv96-22.fasta")       
     

#create a new col as genome name 
VF$genome_name <- sub("\\.fasta$", "", basename(VF$`#FILE`))
summary(VF$genome_name)
#write.csv(VF, "VF_info.csv", row.names = F)


##pick the %IDENTITY>=90 to make a heatmap
VF<- VF[VF$`%IDENTITY` >= 90, ]


#62 genes 
pa <- VF %>%
  select(genome_name, GENE) %>%
  distinct() %>%                    # one row per genome-gene pair
  mutate(present = 1) %>%
  pivot_wider(names_from = GENE,
              values_from = present,
              values_fill = 0)   

#pa_1 = 62 genes
pa <-as.data.frame(pa)
pa_1<-column_to_rownames(pa, var = "genome_name")
pa_1<-t(pa_1)
pa_1<-as.data.frame(pa_1)
pa_1<-rownames_to_column(pa_1, var="Genes")  # 62 genes

# 64 genes-2 = 62 genes 
## removed one for each, need to add them by in svg
# two flgC - Lateral flagella/Polar flagella - Motility
# two fliG - Lateral flagella/Polar flagella - Motility
merge<-merge(pa_1, group, by = "Genes")
merge<-merge[-c(48,51),]
merge<-as.data.frame(merge)
rownames(merge) <- NULL

#core, accessory gene type, create gene_type column
merge$gene_type <- ifelse(
  rowSums(merge[, -c(1, 27:28)] == 1) == 25,
  "core",
  "accessory")

accessory<-filter(merge, merge$gene_type =="accessory")
summary(accessory)

accessory$Group<-as.factor(accessory$Group)


merge<-column_to_rownames(merge, var="Genes")  # 62 genes


my_levels <- c("AV0110","Hv98-225", "AV040", 
               "FDAARGOS_632", "1708-29120",  "b1-3", "KMM358", "KMM366", "C4","KMM289",
               "Hm22","Hv98-571","Hv98-231", "Hv98-232","Hm21" ,"Ho04-635",
               "Hv98-241",  "FC951","D3","ANYA_18196","183026",
               "ZfB1","ANYA_13997", "ANYA_14430")

 
mat <- merge%>%
  select( -Group, -Virulence.factors) %>%
  as.matrix()
mode(mat) <- "numeric"
mat<- mat[, my_levels]


## removed one for each, need to add them by in svg
# two flgC - Lateral flagella/Polar flagella - Motility
# two fliG - Lateral flagella/Polar flagella - Motility

group_row <- merge %>%
  select(Virulence.factors, Group) %>%
  rename(virulence_group = Virulence.factors)  


my_color <- c("grey90", "#b36692")      
ph_breaks <- c(-0.5, 0.5, 1.5)  


group_levels <- sort(unique(group_row$virulence_group))
group_levels_2 <- sort(unique(group_row$Group))


group_colors <- list(
  virulence_group = setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(length(group_levels)),
    nm = group_levels),
  Group = setNames(
    colorRampPalette(brewer.pal(3, "Set3"))(length(group_levels_2)),
    group_levels_2))

png("vf_heatmap_vf_group.png", width = 3500, height = 4500, res = 400)
p<-pheatmap(mat,
         #设置颜色（离散的0/1色块）
         color = my_color,
         breaks = ph_breaks,
         legend_breaks = c(0, 1),
         legend_labels = c("Absent", "Present"),
         labels_row = merge$Genes,
         clustering_distance_rows = "binary",
         clustering_distance_cols = "binary",
         treeheight_row = 12,
         treeheight_col = 10,
         annotation_row = group_row,
         annotation_colors = group_colors,
         annotation_names_row = FALSE,
         #设置行/列的字体
         fontsize_row = 11,
         fontsize_col = 12,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         na_col = "grey80")
dev.off()

svglite("vf_heatmap_vf_group.svg",width = 8,height = 11)
print(p)      
dev.off()


####  AMR brief analysis

AMR <- read.delim("ncbi_amr_all.tab", check.names = FALSE) 
AMR$genome_name <- sub("\\.fasta$", "", basename(AMR$`#FILE`))
#write.csv(AMR, "AMR_info.csv", row.names = F)

library(dplyr)

amr_core <- AMR %>%
  group_by(RESISTANCE) %>%
  summarise(
    n_genomes = n_distinct(genome_name),
    .groups = "drop"
  ) %>%
  filter(n_genomes == 25)

amr_core




