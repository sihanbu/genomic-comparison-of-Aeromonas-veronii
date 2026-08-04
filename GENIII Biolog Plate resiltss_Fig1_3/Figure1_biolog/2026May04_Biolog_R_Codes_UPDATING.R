require(MESS)
require(tibble)
require(rstatix)
require(RColorBrewer)
require(tibble)
require(dplyr)
require(stringr)
require(tidyr)
require(MESS)
require(rstatix)
require(ggplot2)
require(RColorBrewer)
require(ggh4x)
require(patchwork)

setwd("/Users/sihanbu/Library/CloudStorage/OneDrive-UniversityofConnecticut/Research/whole\ genome\ assembly/2026Mar04_Biolog_Analysis_new_data/R\ analysis")


###########################################################
########   Run3: AUC statistics, 72h  #####################
###########################################################

run3_biolog_org<-read.csv("Run3_Feature_Table.csv",check.names = FALSE,fileEncoding = "UTF-8")

run3_biolog<-run3_biolog_org[,-1]
run3_biolog<-column_to_rownames(run3_biolog, var="SUBSTRATE")
run3_biolog<-t(run3_biolog)

run3_biolog<-as.data.frame(run3_biolog)
run3_biolog<-rownames_to_column(run3_biolog, var="SampleID")

#Add 3 columns: strain, replicate, timepoint
run3_biolog <- run3_biolog %>%
  mutate(
    Strain     = str_sub(SampleID, 1, 3),
    Replicate  = str_extract(SampleID, "P\\d+") |> str_remove("P") |> as.integer(),
    Timepoint  = str_extract(SampleID, "\\d+$") |> as.integer()
  ) %>%
  select(Strain, Replicate, Timepoint, everything())

run3_biolog_72<- run3_biolog %>%
  mutate(Strain = paste0("KMM", Strain)) %>%
  filter(Timepoint != 96 & Timepoint != 120)


#calculate AUC per replicate across time points, per substrate (columns)
run3_biolog_AUC_72 <-run3_biolog_72 %>%
  pivot_longer(
    cols = -c(Strain, Replicate, Timepoint, SampleID),
    names_to = "Substrate",
    values_to = "OD"
  ) %>%
  group_by(Strain, Replicate, Substrate) %>%
  arrange(Timepoint, .by_group = TRUE) %>%
  summarise(AUC = auc(Timepoint, OD), .groups = "drop")


#stats on comparing 3 strains using one way ANOVA
#overall p-values
summary(run3_biolog_AUC_72)
run3_biolog_AUC_72$Strain<-as.factor(run3_biolog_AUC_72$Strain)
run3_biolog_AUC_72$Replicate<-as.factor(run3_biolog_AUC_72$Replicate)
run3_biolog_AUC_72$Substrate<-as.factor(run3_biolog_AUC_72$Substrate)

Run3_anova_AUC_72<- run3_biolog_AUC_72 %>%
  group_by(Substrate) %>%
  anova_test(AUC ~ Strain)  %>%
  as.data.frame() %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

Run3_anova_AUC_72_sig<- Run3_anova_AUC_72 %>%
  as.data.frame() %>%
  filter(p_adj < 0.05)

#pairwise comparsion 
Run3_tukey_pairwise_72 <- run3_biolog_AUC_72 %>%
  group_by(Substrate) %>%
  tukey_hsd(AUC ~ Strain)  

#subset the significant ones from pairwise comparsion 
sig_sub_72 <- Run3_tukey_pairwise_72 %>%
  filter(p.adj < 0.05) %>%
  pull(Substrate) %>%
  unique()

Run3_sig_all_72 <- Run3_tukey_pairwise_72 %>%
  filter(Substrate %in% sig_sub_72)


 



###########################################################
########   Run2: AUC statistics, 72h  #####################
###########################################################

run2_biolog_org<-read.csv("Run2_Feature_Table.csv",check.names = FALSE,fileEncoding = "UTF-8")

run2_biolog<-run2_biolog_org[,-1]
run2_biolog<-column_to_rownames(run2_biolog, var="SUBSTRATE")
run2_biolog<-t(run2_biolog)

run2_biolog<-as.data.frame(run2_biolog)
run2_biolog<-rownames_to_column(run2_biolog, var="SampleID")

#Add 3 columns: strain, replicate, timepoint
run2_biolog <- run2_biolog %>%
  mutate(
    Strain     = str_sub(SampleID, 1, 3),
    Replicate  = str_extract(SampleID, "P\\d+") |> str_remove("P") |> as.integer(),
    Timepoint  = str_extract(SampleID, "\\d+$") |> as.integer()
  ) %>%
  select(Strain, Replicate, Timepoint, everything())

run2_biolog_72<- run2_biolog %>%
  mutate(Strain = paste0("KMM", Strain)) %>%
  filter(Timepoint != 96 & Timepoint != 120)


#calculate AUC per replicate across time points, per substrate (columns)
run2_biolog_AUC_72 <-run2_biolog_72 %>%
  pivot_longer(
    cols = -c(Strain, Replicate, Timepoint, SampleID),
    names_to = "Substrate",
    values_to = "OD"
  ) %>%
  group_by(Strain, Replicate, Substrate) %>%
  arrange(Timepoint, .by_group = TRUE) %>%
  summarise(AUC = auc(Timepoint, OD), .groups = "drop")


#stats on comparing 3 strains using one way ANOVA
#overall p-values
summary(run2_biolog_AUC_72)
run2_biolog_AUC_72$Strain<-as.factor(run2_biolog_AUC_72$Strain)
run2_biolog_AUC_72$Replicate<-as.factor(run2_biolog_AUC_72$Replicate)
run2_biolog_AUC_72$Substrate<-as.factor(run2_biolog_AUC_72$Substrate)

Run2_anova_AUC_72<- run2_biolog_AUC_72 %>%
  group_by(Substrate) %>%
  anova_test(AUC ~ Strain)  %>%
  as.data.frame() %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

Run2_anova_AUC_72_sig<- Run2_anova_AUC_72 %>%
  as.data.frame() %>%
  filter(p_adj < 0.05)

#pairwise comparsion 
Run2_tukey_pairwise_72 <- run2_biolog_AUC_72 %>%
  group_by(Substrate) %>%
  tukey_hsd(AUC ~ Strain)  

#subset the significant ones from pairwise comparsion 
sig_sub_72 <- Run2_tukey_pairwise_72 %>%
  filter(p.adj < 0.05) %>%
  pull(Substrate) %>%
  unique()

Run2_sig_all_72 <- Run2_tukey_pairwise_72 %>%
  filter(Substrate %in% sig_sub_72)


##############################################
########   Run4: AUC statistics, 72h ############# 
##############################################
run4_biolog_org<-read.csv("Run4_Feature_Table.csv",check.names = FALSE,fileEncoding = "UTF-8")

run4_biolog<-run4_biolog_org[,-1]
run4_biolog<-column_to_rownames(run4_biolog, var="SUBSTRATE")
run4_biolog<-t(run4_biolog)

run4_biolog<-as.data.frame(run4_biolog)
run4_biolog<-rownames_to_column(run4_biolog, var="SampleID")

#Add 3 columns: strain, replicate, timepoint
run4_biolog <- run4_biolog %>%
  mutate(
    Strain     = str_sub(SampleID, 1, 3),
    Replicate  = str_extract(SampleID, "P\\d+") |> str_remove("P") |> as.integer(),
    Timepoint  = str_extract(SampleID, "\\d+$") |> as.integer()
  ) %>%
  select(Strain, Replicate, Timepoint, everything())

run4_biolog_72<- run4_biolog %>%
  mutate(Strain = paste0("KMM", Strain)) %>%
  filter(Timepoint != 96 & Timepoint != 120)



#calculate AUC per replicate across time points, per substrate (columns)
run4_biolog_72_AUC <-run4_biolog_72  %>%
  pivot_longer(
    cols = -c(Strain, Replicate, Timepoint, SampleID),
    names_to = "Substrate",
    values_to = "OD"
  ) %>%
  group_by(Strain, Replicate, Substrate) %>%
  arrange(Timepoint, .by_group = TRUE) %>%
  summarise(AUC = auc(Timepoint, OD), .groups = "drop")

#stats on comparing 3 strains using one way ANOVA
#overall p-values
summary(run4_biolog_72_AUC)
run4_biolog_72_AUC$Strain<-as.factor(run4_biolog_72_AUC$Strain)
run4_biolog_72_AUC$Replicate<-as.factor(run4_biolog_72_AUC$Replicate)
run4_biolog_72_AUC$Substrate<-as.factor(run4_biolog_72_AUC$Substrate)

Run4_anova_AUC<- run4_biolog_72_AUC %>%
  group_by(Substrate) %>%
  anova_test(AUC ~ Strain)  %>%
  as.data.frame() %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

Run4_anova_AUC_sig<- Run4_anova_AUC %>%
  as.data.frame() %>%
  filter(p_adj < 0.05)

#pairwise comparsion 
Run4_tukey_pairwise <- run4_biolog_72_AUC %>%
  group_by(Substrate) %>%
  tukey_hsd(AUC ~ Strain)  


#subset the significant ones from pairwise comparsion 
sig_sub <- Run4_tukey_pairwise %>%
  filter(p.adj < 0.05) %>%
  pull(Substrate) %>%
  unique()

Run4_sig_all <- Run4_tukey_pairwise %>%
  filter(Substrate %in% sig_sub)





######################################################################
########  Find consistent substrates in 3 tukey results  ############# 
######################################################################

consistent_substrate_72 <- Reduce(
  intersect,
  list( unique(Run2_sig_all_72$Substrate),
    unique(Run3_sig_all_72$Substrate),
    unique(Run4_sig_all$Substrate)))

Run4_consistent_72 <- Run4_sig_all %>%
  filter(Substrate %in% consistent_substrate_72)


#get the grouping info
grouping<-read.csv("tabula-Origanum heracleoticum paper_STable1_GENIII_grouping_KEGG.csv", header = T,fileEncoding = "UTF-8")


#avergae the AUC for heatmap 
auc_mean <- run4_biolog_72_AUC %>%
  group_by(Strain, Substrate) %>%
  summarise(mean_AUC = mean(AUC))


auc_heatmap<-Run4_consistent_72  %>%
  left_join(
    auc_mean %>%
      select(Substrate, Strain, mean_AUC) %>%
      rename(group1 = Strain, mean_AUC_group1 = mean_AUC),
    by = c("Substrate", "group1")) %>%
  left_join(
    auc_mean %>%
      select(Substrate, Strain, mean_AUC) %>%
      rename(group2 = Strain, mean_AUC_group2 = mean_AUC),
    by = c("Substrate", "group2")) %>%
  mutate(diff_meanAUC = mean_AUC_group1 - mean_AUC_group2) %>%
  mutate(compare = paste(group1, group2, sep = "-"))  %>%
  mutate(sign.mark = ifelse(p.adj.signif == "ns", "", p.adj.signif))

auc_heatmap<-left_join(auc_heatmap, grouping, by="Substrate")
summary(auc_heatmap)

auc_heatmap$compare<-as.factor(auc_heatmap$compare)
auc_heatmap$Group<-as.factor(auc_heatmap$Group)
auc_heatmap$Substrate<-as.factor(auc_heatmap$Substrate)
auc_heatmap$diff_meanAUC <-as.numeric(auc_heatmap$ diff_meanAUC)
auc_heatmap$Updated_Carbon_subgroup<-as.factor(auc_heatmap$Updated_Carbon_subgroup)


#visualize
heatmap <- ggplot(auc_heatmap, aes(compare, Substrate)) +
  geom_tile(aes(fill = diff_meanAUC)) +
  geom_text(aes(label = sign.mark), color = "black", vjust = 0.75, size = 6) +
  scale_fill_gradientn(colours = rev(brewer.pal(11, "RdYlBu")), name = "ΔAUC") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))+
  ggh4x::facet_nested(
    Group + Updated_Carbon_subgroup ~ .,
    scales = "free_y",
    space  = "free_y",
    switch = "y",
    strip = ggh4x::strip_nested(
      text_y = ggh4x::elem_list_text(
        element_text(face = "bold", size = 10),  # outer strip (Group)
        element_text(face = "bold", size = 10)   # inner strip (grey boxes)
      ),
      background_y = ggh4x::elem_list_rect(
        fill  = c("white", "grey85"),   # outer Group = white; inner subgroup = grey
        color = c("white", "grey60")),by_layer_y = TRUE)) +
  theme_bw() +
  theme(axis.title = element_blank(),  
        axis.text.x = element_text(size = 10,face = "bold"),  
        axis.text.y = element_text(size = 10,face = "bold"),   
        legend.title = element_text(size = 10,face = "bold"),   
        legend.text = element_text(size = 10,face = "bold"),   
        legend.key.width = unit(0.3, "cm"),       
        legend.key.height = unit(0.4, "cm"),     
        strip.background.y = element_rect(color="white", fill = "white"),  
        strip.placement = "outside",   
        strip.switch.pad.grid = unit(1, "cm"),
        strip.text.x = element_text(size = 9,face = "bold"),   
        strip.text.y.left = element_text(size = 9, angle = 90, hjust = 0.5,face = "bold"),  
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        axis.ticks = element_line(linewidth = 0.6),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'))


heatmap

ggsave(filename = "Consistent_Sig_Pairwise_heatmap_use_Run4AUC.png",
       plot = heatmap,
       width = 8, height = 6,
       dpi =800,
       device = "png")



################################
###### Growth curve  ########

regroup<-data.frame(Substrate = c(
  "Sucrose",
  "N-Acetyl-D-Galactosamine",
  "N-Acetyl-β-D-Mannosamine"  ,
  "D-Mannose"  ,
  "D-Fructose-6-PO4"   ,
  "D-Fructose" ,
  "Propionic Acid",
  "p-Hydroxy-Phenylacetic Acid" , 
  "Formic Acid",
  "Tetrazolium Blue" ,
  "Rifamycin SV",
  "Niaproof 4" ),
  relabel_group = c( "Carbon source | Carbohydrate",
                     "Carbon source | Carbohydrate",
                     "Carbon source | Carbohydrate",
                     "Carbon source | Carbohydrate",
                     "Carbon source | Carbohydrate",
                     "Carbon source | Carbohydrate",
                    "Carbon source | Carboxylic acid",
                    "Carbon source | Carboxylic acid",
                    "Carbon source | Carboxylic acid",
                    "Chemical sensitivity",
                    "Chemical sensitivity",
                    "Chemical sensitivity")) 

order<-c(
  "Sucrose",
  "N-Acetyl-D-Galactosamine",
  "N-Acetyl-β-D-Mannosamine"  ,
  "D-Mannose"  ,
  "D-Fructose-6-PO4"   ,
  "D-Fructose" ,
  "Propionic Acid",
  "p-Hydroxy-Phenylacetic Acid" , 
  "Formic Acid",
  "Tetrazolium Blue" ,
  "Rifamycin SV",
  "Niaproof 4" )
  
  
run4_biolog_filtered <- run4_biolog %>%
  select(1:3, any_of(consistent_substrate_72))

run4_biolog_filtered <- filter(run4_biolog_filtered, Timepoint != 96 &Timepoint != 120 )

 
long_df <- run4_biolog_filtered %>%
  pivot_longer(cols = 4:ncol(.),names_to = "Substrate",values_to = "Value")

long_df <-left_join(long_df, regroup, by = "Substrate")


long_df$relabel_group<-as.factor(long_df$relabel_group)
long_df$Timepoint<-as.factor(long_df$Timepoint)

long_df <- long_df%>%
  mutate(Strain = paste0("KMM", Strain)) %>%
  mutate(Substrate = factor(Substrate, levels = order))
 

growth_curve<-ggplot(long_df, aes(x = Timepoint, y = Value, color = Strain)) +
  stat_summary(aes(group = Strain), 
               geom = "errorbar",fun.data = mean_se, width = 0.4, size = 0.6) +
  stat_summary(aes(group = Strain),geom = "line", fun = mean, linewidth = 1) +
  stat_summary(aes(fill = Strain), geom = "point", 
               fun = mean, size = 2, shape = 21, stroke = 0.8)+
  facet_wrap(~ relabel_group+Substrate, scales = "free_y") +
  labs(x = "Incubation time (h)", y = "Corrected OD value") +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    strip.placement = "outside",
    legend.position = "bottom" ,
    strip.text.x = element_text(face = "bold", size = 9),  # Substrate
    strip.text.y = element_text(face = "bold", size = 9),   # relabel_group
    strip.background.y = element_rect(fill = "grey90", color = NA),
    axis.title.y = element_text(margin = margin(r = 6)))


 
combined_plot <- ( heatmap+growth_curve ) +
  plot_layout(widths = c( 2,4.8)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20, face = "bold"))

ggsave("Run4_Combined_Figure.png",
       combined_plot,
       width = 18, height = 9,
       dpi = 500)








#######################################################
########### KEGG results map to biolog  #######
############################################


#find the compound names of the significant taxa 
#download KEGG module and co database 
co_module_full<-read.table("https://rest.kegg.jp/link/compound/module", quote="", sep="\t")
colnames(co_module_full)<- c("Module","Compound")
co_module_full$Compound <- sub("cpd:", "", co_module_full$Compound)
co_module_full$Module     <- sub("md:", "", co_module_full$Module)

grouping<-read.csv("tabula-Origanum heracleoticum paper_STable1_GENIII_grouping_KEGG.csv", header = T,fileEncoding = "UTF-8")

co_mo_biolog<-merge(grouping, co_module_full, by.x = "KEGG", by.y = "Compound")   

compound_modules <- co_mo_biolog %>%
  group_by(KEGG) %>%
  summarise(
    Modules = paste(unique(Module), collapse = ", "))

#Supplementary table3
merge<-full_join(grouping, compound_modules, by="KEGG")
write.csv(merge, "tabula-Origanum heracleoticum paper_STable1_GENIII_grouping_KEGG_Module.csv", row.names = F)

##################
#pick significant taxa 
######################

target_substrates <- c(
  "Sucrose",
  "N-Acetyl-D-Galactosamine",
  "N-Acetyl-β-D-Mannosamine"  ,
  "D-Mannose"  ,
  "D-Fructose-6-PO4"   ,
  "D-Fructose" ,
  "Propionic Acid",
  "p-Hydroxy-Phenylacetic Acid" , 
  "Formic Acid",
  "Tetrazolium Blue" ,
  "Rifamycin SV",
  "Niaproof 4" )

picked <- merge %>%
  dplyr::filter(Substrate %in% target_substrates)


anvio_kegg_module <- read.table(
  "/Users/sihanbu/Library/CloudStorage/OneDrive-UniversityofConnecticut/Research/whole\ genome\ assembly/2026Mar04_Biolog_Analysis_new_data/Anvio\ KEGG\ analysis/289_358_366_metabolism_modules_custom_in\ supp\ table1.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "")

#mapping 
#format the picked 
picked_long <- picked %>%
  separate_rows(Modules, sep = ",\\s*")

mapping<-merge(picked_long,anvio_kegg_module, by.x="Modules", by.y = "module")

kegg_completeness_wide <- mapping %>%
  filter(db_name %in% c("KMM289", "KMM358", "KMM366")) %>%
  distinct() %>%
  pivot_wider(
    names_from = db_name,
    values_from = c(pathwise_module_completeness,
                    stepwise_module_completeness),
    names_glue = "{.value}_{db_name}")

write.csv(kegg_completeness_wide, "kegg_enrich_completeness_biolog.csv", row.names = F)


