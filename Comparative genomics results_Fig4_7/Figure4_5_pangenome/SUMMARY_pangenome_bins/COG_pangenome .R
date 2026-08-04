setwd("/Users/sihanbu/Library/CloudStorage/OneDrive-UniversityofConnecticut/Research/whole\ genome\ assembly/2026Mar04_Biolog_Analysis_new_data/selected\ genomes\ comparision/01_PANGENOME/New\ with\ 24\ genomes/SUMMARY_pangenome_bins")

cog_pan_sum<- read.delim("All-GENOMES_COG_gene_clusters_summary.txt",stringsAsFactors = FALSE)

#how many gene clusters have COG annotation? 
#5017
summary(unique(cog_pan$gene_cluster_id))
cog_pan  %>%
  summarise(n_gene_clusters_with_COG = n_distinct(gene_cluster_id))

#use this info in txt, use cog_pan_sum regardless of cog annotation
summary_bin<- cog_pan_sum %>%
  group_by(bin_name) %>%
  summarise(n = n(), .groups = "drop")


summary(cog_pan)


summary <- cog_pan %>%
  group_by(bin_name,genome_name,COG24_CATEGORY) %>%
  summarise(n = n(), .groups = "drop")

#most abundant regardless bin and genome name
summary_category<-cog_pan %>%
  group_by(COG24_CATEGORY) %>%
  summarise(n = n(), .groups = "drop")

factor(summary$genome_name)
#change the names 
# Ho04_635 to  Ho04-635
# b1_3  to b1-3 
# Hv98_241  to Hv98-241 
# Hv98_571  to  Hv98-571
summary <- summary %>%
  mutate(
    genome_name = dplyr::recode(
      genome_name,
      "Ho04_635" = "Ho04-635",
      "b1_3" = "b1-3",
      "Hv98_241" = "Hv98-241",
      "Hv98_571" = "Hv98-571",
      "AV1708_29120" = "1708-29120",
      "AV183026" = "183026",
      "Hv96_22" = "Hv96-22",
      "Hv98_571" = "Hv98-571",
      "Hv98_241"="Hv98-241",
      "Hv98_231" ="Hv98-231",
      "Hv98_232" = "Hv98-232",
      "Hv98_225" = "Hv98-225"))

 

my_levels <- c("Hv98-571", "Hm22", "KMM366", "KMM358",
               "Hv98-231", "Hv98-232", "Hv98-241", "Hm21", "Ho04-635",
               "b1-3", "1708-29120", "ANYA_14430", "ANYA_13997", "FC951",
               "AV040", "KMM289", "D3", "Hv98-225", "ANYA_18196",
               "FDAARGOS_632", "AV0110", "C4", "ZfB1", "183026")


top25 <- summary %>%
  group_by(COG24_CATEGORY) %>%
  summarise(total = sum(n), .groups = "drop")


#Make heatmap

##top 25 in core gene cluster
top25_core <- summary %>%
  filter(bin_name == "Core") %>%
  group_by(COG24_CATEGORY) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 25) %>%
  pull(COG24_CATEGORY)


#make heatmaps on core and unique 
core <- summary %>%
  filter(bin_name == "Core",COG24_CATEGORY %in% top25_core) %>%
  select(-bin_name) %>%
  pivot_wider(names_from = COG24_CATEGORY,values_from = n,values_fill = 0) 
  
core <- core[match(my_levels, core$genome_name), ]
core<-column_to_rownames(core, var = "genome_name")
core_mat <- t(as.matrix(core))



rownames(core_mat) <- str_wrap(rownames(core_mat), width = 35)

p1<-pheatmap(
  core_mat,
  color = colorRampPalette(c("white", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey90",
  fontsize_row = 10,
  fontsize_col = 8,
  angle_col = "45",
  main = "Core gene clusters of pan-genome")


#export as a svg file
svglite("Core_gene_clusters_heatmap.svg",width = 10,height = 10)
print(p1)      
dev.off()


top25_unique <- summary %>%
  filter(bin_name == "Unique") %>%
  group_by(COG24_CATEGORY) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 25) %>%
  pull(COG24_CATEGORY)

unique <- summary %>%
  filter(bin_name == "Unique",COG24_CATEGORY %in% top25_unique) %>%
  select(-bin_name) %>%
  pivot_wider(names_from = COG24_CATEGORY,values_from = n,values_fill = 0) 

unique <- unique[match(my_levels, unique$genome_name), ]
unique<-as.data.frame(unique)
# 4 of them doesn't have unique gene clusters annotation - same as the circular chart
unique$genome_name[is.na(unique$genome_name)] <-
  c("Hv98-571", "Hm22", "Hv98-231", "Hv98-232")


unique<-column_to_rownames(unique, var = "genome_name")
unique_mat <- t(as.matrix(unique))


rownames(unique_mat) <- str_wrap(rownames(unique_mat), width = 35)

p2<-pheatmap(
  unique_mat,
  color = colorRampPalette(c("white", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey90",
  fontsize_row = 10,
  fontsize_col = 8,
  angle_col = "45",
  main = "Unique gene clusters of pan-genome")

#export as a svg file
svglite("Unique_gene_clusters_heatmap.svg",width = 10,height = 8)
print(p2)      
dev.off()

###Extra analysis
top25_acc <- summary %>%
  filter(bin_name == "Accessory") %>%
  group_by(COG24_CATEGORY) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 25) %>%
  pull(COG24_CATEGORY)

##
per_cluster <- cog_pan_sum %>%
  distinct(gene_cluster_id, .keep_all = TRUE) %>%
  mutate(bin_name = tidyr::replace_na(bin_name, "Accessory"))   # 1 blank cluster, present in 3 genomes

per_cluster %>%
  count(bin_name) %>%
  mutate(pct = round(100 * n / sum(n), 2))


per_gene<- cog_pan_sum %>%
  distinct(unique_id, .keep_all = TRUE) %>%
  mutate(bin_name = tidyr::replace_na(bin_name, "Accessory"))   # 1 blank cluster, present in 3 genomes

per_gene %>%
  count(bin_name) %>%
  mutate(pct = round(100 * n / sum(n), 2))

 