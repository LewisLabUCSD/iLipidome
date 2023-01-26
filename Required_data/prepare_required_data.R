library(tidyverse)
library(igraph)

file <- dirname(rstudioapi::getSourceEditorContext()$path)

#---------------Prepare lipid substructure data-----------------

network_node <- read.csv(file.path(file,'network_node.csv'))
network_edge <- read.csv(file.path(file,'network_edge.csv'))

graph <- graph_from_data_frame(network_edge[c('S1','P1')], directed = T, vertices = network_node$Abbreviation)

find_path_from_G3P <- function(network_node, network_edge){
  lipid <- character()
  G3P_start_path <- list()
  lipid_list <- unique(c(network_edge$S1,network_edge$P1))
  num_class <- 1
  num_path <- 1
  
  while(!is.na(lipid_list[num_class])){
    class <- network_node %>% filter(Abbreviation==lipid_list[num_class]) %>% .$Class

    if(class=='Sphingolipid'){
      shortest_path <- tryCatch({all_simple_paths(graph,'Serine+Palmitoyl CoA',lipid_list[num_class], mode='out')},
                                error=function(e){NULL})
    }
    else if(class=='Deoxysphingolipid'){
      shortest_path <- tryCatch({all_simple_paths(graph,'Alanine+Palmitoyl CoA',lipid_list[num_class], mode='out')},
                                error=function(e){NULL})
    }
    else{
      shortest_path <- tryCatch({all_simple_paths(graph,'G3P',lipid_list[num_class], mode='out')},
                                error=function(e){NULL})
    }
    
    if(!is.null(shortest_path)){
      if(length(shortest_path)!=0){
        
        for (path_num in 1:length(shortest_path)) {
          lipid[num_path] <- lipid_list[num_class]
          G3P_start_path[[num_path]] <- shortest_path[[path_num]] %>% 
            attributes() %>% .$names
          num_path <- num_path+1
          
        }
      }
      else{
        lipid[num_path] <- lipid_list[num_class]
        G3P_start_path[[num_path]] <- ''
        num_path <- num_path+1
      }
      
    }
    else{
      lipid[num_path] <- lipid_list[num_class]
      G3P_start_path[[num_path]] <- ''
      num_path <- num_path+1
    }
    num_class <- num_class+1
  }
  names(G3P_start_path) <- lipid
  
  
  lipid_substructure <- plyr::ldply(G3P_start_path, rbind) %>% unique()
  
  lipid_substructure[is.na(lipid_substructure)] <- ''
  
  colnames(lipid_substructure) <- c('Lipid', str_c('Unit', 1:(ncol(lipid_substructure)-1)))
  
  
  for(a in 1:nrow(lipid_substructure)){
    FA_num <- network_node %>% filter(Abbreviation==lipid_substructure$Lipid[a]) %>% .$FA
    if(FA_num!=0){
      lipid_substructure[-1][lipid_substructure[-1]==lipid_substructure$Lipid[a]] <- str_c(lipid_substructure$Lipid[a], '_FA', as.character(FA_num))
    }
  }
  
  return(lipid_substructure)
  
}

lipid_substructure <- find_path_from_G3P(network_node, network_edge)
#write.csv(lipid_substructure, file.path(file, 'lipid_substructure.csv'), row.names = F)


#---------------Prepare lipidmaps reference data-----------------
lipidmaps <- read.delim(file.path(file,'LIPIDMAPS.txt'))

lipidmaps_sum <- lipidmaps %>% mutate(sub_class_id=str_extract(sub_class,'\\[\\w+\\]')) %>% 
  filter(sub_class_id %in% network_node$LIPIDMAPS_subclass_id) %>% 
  filter(!is.na(sub_class_id)) %>% 
  mutate(FA_sum=str_extract(abbrev, '\\d+:\\d+')) %>% 
  mutate(FA_split=str_extract_all(abbrev_chains, '\\d+:\\d+') %>% 
           map(.f = function(x){str_c(x, collapse = '_')}) %>% unlist) %>% 
  .[c('name','sub_class_id','FA_sum','FA_split')] %>% filter(FA_split!='')

#write.csv(lipidmaps_sum, file.path(file, 'lipidmaps_sum.csv'), row.names = F)

LM_mapping <- lipidmaps_sum %>% filter(sub_class_id=='[GP1201]')
CL <- character()
possible_MLCL_sum <- character()
possible_MLCL_split <- character()
possible_PG_sum <- character()
possible_PG_split <- character()
add <- 1
for(num1 in 1:nrow(LM_mapping)){
  CL_FA_sum <- LM_mapping$FA_sum[num1] %>% str_extract_all('\\d+') %>% 
    unlist() %>% as.integer()
  CL_FA_split <-LM_mapping$FA_split[num1] %>% str_split('_') %>% 
    unlist()
  for(num2 in CL_FA_split){
    FA_split <- str_extract_all(num2, '\\d+') %>% unlist() %>% 
      as.integer()
    
    MLCL_split <- which(CL_FA_split==num2)[1]
    MLCL_split <- CL_FA_split[-MLCL_split] %>% str_c(collapse = '_')
    
    MLCL_FA_split <-MLCL_split %>% str_split('_') %>% 
      unlist()
    
    
    for(num3 in MLCL_FA_split){
      CL[add] <- LM_mapping$name[num1]
      
      possible_MLCL_sum[add] <- str_c(CL_FA_sum[1]-FA_split[1], ':', 
                                      CL_FA_sum[2]-FA_split[2], ';0')
      
      possible_MLCL_split[add] <- MLCL_split
      
      MLCL_FA_sum <- possible_MLCL_sum[add] %>% str_extract_all('\\d+') %>% 
        unlist() %>% as.integer()
      
      FA_split2 <- str_extract_all(num3, '\\d+') %>% unlist() %>% 
        as.integer()
      
      possible_PG_sum[add] <- str_c(MLCL_FA_sum[1]-FA_split2[1], ':', 
                                    MLCL_FA_sum[2]-FA_split2[2], ';0')
      
      PG_split <- which(MLCL_FA_split==num3)[1]
      PG_split <- MLCL_FA_split[-PG_split] %>% str_c(collapse = '_')
      
      possible_PG_split[add] <- PG_split
      
      add <- add+1
    }
    
  }
}

MLCL_mapping <- LM_mapping %>% 
  left_join(unique(data.frame(name=CL,
                              possible_MLCL_sum=possible_MLCL_sum,
                              possible_MLCL_split=possible_MLCL_split,
                              possible_PG_sum=possible_PG_sum,
                              possible_PG_split=possible_PG_split)),
            by='name') %>% unique()

#write.csv(MLCL_mapping, file.path(file, 'MLCL_mapping.csv'), row.names = F)

