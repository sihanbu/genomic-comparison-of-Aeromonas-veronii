setwd("/Users/sihanbu/Library/CloudStorage/OneDrive-UniversityofConnecticut/Research/whole\ genome\ assembly/2026Mar04_Biolog_Analysis_new_data/Anvio\ KEGG\ analysis")
metabolism_modules <- read.delim("289_358_366_metabolism_modules_custom_in supp table1.txt",header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)

compare_all <- metabolism_modules %>%
  filter(db_name %in% c("KMM289", "KMM358", "KMM366")) %>%
  select(module, module_name, db_name,
         pathwise_module_completeness,
         stepwise_module_completeness) %>%
  pivot_wider(
    names_from = db_name,
    values_from = c(pathwise_module_completeness,
                    stepwise_module_completeness),
    names_glue = "{db_name}_{.value}")


compare_all_ordered <- compare_all %>%
  
  
  
  # Remove completely blank rows (all 3 pathwise are NA)
  
  filter(
    
    !(is.na(KMM289_pathwise_module_completeness) &
        
        is.na(KMM358_pathwise_module_completeness) &
        
        is.na(KMM366_pathwise_module_completeness))
    
  ) %>%
  
  
  
  mutate(
    
    # Any NA in the 3 pathwise columns
    
    has_na_pathwise = if_any(
      
      c(KMM289_pathwise_module_completeness,
        
        KMM358_pathwise_module_completeness,
        
        KMM366_pathwise_module_completeness),
      
      is.na
      
    ),
    
    
    
    # Difference ignoring NA
    
    pathwise_difference = pmax(
      
      KMM289_pathwise_module_completeness,
      
      KMM358_pathwise_module_completeness,
      
      KMM366_pathwise_module_completeness,
      
      na.rm = TRUE
      
    ) -
      
      pmin(
        
        KMM289_pathwise_module_completeness,
        
        KMM358_pathwise_module_completeness,
        
        KMM366_pathwise_module_completeness,
        
        na.rm = TRUE
        
      ),
    
    
    
    # TRUE only if all 3 are present and identical
    
    same_pathwise = !has_na_pathwise & (
      
      KMM289_pathwise_module_completeness ==
        
        KMM358_pathwise_module_completeness &
        
        KMM358_pathwise_module_completeness ==
        
        KMM366_pathwise_module_completeness
      
    )
    
  ) %>%
  
  
  
  # Order:
  
  # 1. Partial NA rows first
  
  # 2. Different rows next
  
  # 3. Same rows last
  
  arrange(
    
    desc(has_na_pathwise),
    
    same_pathwise,
    
    desc(pathwise_difference)
    
  )
  


group_info<-read.delim("289_358_366_metabolism_modules.txt",header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)

group_info<-group_info[,c(1,4,5,6,7)]

merge<-merge(compare_all_ordered, group_info, by="module",all.x = TRUE, sort = FALSE)

write.csv(merge,"Supp_Table1.csv", row.names = F)

