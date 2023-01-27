library(tidyverse)
library(dplyr)
library(igraph)
library(visNetwork)
library(xlsx)
library(data.table)
library(MKmisc)
library(gplots)
library(gtools)

build_char_table <- function(raw_data, network_node){
  
  class <- raw_data$feature %>% str_extract('[A-Za-z]+( O-)*')
  class_included <- class %in% network_node$Abbreviation
  
  
  totallength <- raw_data$feature %>% str_extract_all('\\d+:') %>% 
    map_int(~str_sub(.x, end = -2) %>% as.integer() %>% sum)
  
  totaldb <- raw_data$feature %>% str_extract_all('\\d+;') %>% 
    map_int(~str_sub(.x, end = -2) %>% as.integer() %>% sum)
  
  totaloh <- raw_data$feature %>% str_extract_all(';\\d+') %>% 
    map_int(~str_sub(.x, start = 2) %>% as.integer() %>% sum)
  
  FA_sum <- str_c(totallength, ':', totaldb, ';', totaloh)
  
  FA_split <- raw_data$feature %>% str_extract_all('\\d+:\\d+;\\d+') %>% 
    map_chr(~str_c(.x, collapse = '_'))
  
  char_table <- data.frame(feature=raw_data$feature, class=class, totallength=totallength,
                           totaldb=totaldb, totaloh=totaloh, FA_sum=FA_sum, FA_split=FA_split) %>% 
    .[class_included,] %>% left_join(network_node[c('Abbreviation','FA')], by=c('class'='Abbreviation'))
  
  colnames(char_table)[8] <- 'FA_num'
  FA_exact <- map_int(str_split(char_table$FA_split, '_'), ~length(.x))==char_table$FA_num
  char_table$FA_split[!FA_exact] <- ''
  
  each_FA <- str_extract_all(char_table$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})
  
  char_table <- char_table %>% mutate(each_FA=each_FA)
  
  raw_data <- raw_data[class_included,]
  
  return(list(raw_data, char_table))
}


unprocessed_data_test <- function(exp_data, char_table, method='t.test',
                                    significant=p_value, ctrl_group, 
                                    exp_group){
  
  ctrl_group <- ctrl_group+1
  exp_group <- exp_group+1
  
  #FA expression
  
  char_sel <- char_table %>% filter(feature %in% exp_data$feature)
  
  
  each_FA <- str_extract_all(char_sel$FA_split, '\\d+:\\d+;\\d+')
  
  FA_mapping <- unique(unlist(each_FA)) %>% sort()
  
  Lipid_FA_freq_list <- each_FA  %>% 
    map(.f = function(x){factor(sort(x[x!='0:0;0']), levels = FA_mapping) %>% 
        table() %>% as.integer()})
  
  Lipid_FA_freq_matrix <- matrix(unlist(Lipid_FA_freq_list), ncol = length(Lipid_FA_freq_list))
  
  exp_matrix <- exp_data[-1] %>% as.matrix()
  
  
  FA_exp <- Lipid_FA_freq_matrix %*% exp_matrix %>% as.data.frame()
  
  FA_exp <- FA_exp %>% mutate(feature=FA_mapping, type='FA') %>% 
    dplyr::select(feature, type, everything())
  
  #lipid class expression
  
  exp_data[-1][is.na(exp_data[-1])] <- 0
  
  lipid_class_exp <- exp_data %>% left_join(char_table[c('feature', 'class')], by='feature') %>% 
    dplyr::select(-1)
  
  lipid_class_exp <- lipid_class_exp[!is.na(lipid_class_exp[['class']]),]
  if(nrow(lipid_class_exp)==0){
    lipid_class_exp <- data.frame()
  }else{
    lipid_class_exp <- lipid_class_exp %>% 
      aggregate(as.formula(str_c('. ~ class')), ., sum)
  }
  
  colnames(lipid_class_exp)[1] <- 'feature'
  
  lipid_class_exp <- lipid_class_exp %>% mutate(type='class') %>% 
    dplyr::select(feature, type, everything())
  
  lipid_species_exp <- exp_data %>% mutate(type='species') %>% 
    dplyr::select(feature, type, everything())
  
  all_exp_data <- rbind(lipid_species_exp, lipid_class_exp, FA_exp)
  
  #statistics test
  lipid <- character()
  type <- character()
  
  mean_ctrl <- numeric()
  mean_exp <- numeric()
  mean_all <- numeric()
  
  FC <- numeric()
  
  sd_ctrl <- numeric()
  sd_exp <- numeric()
  p_value <- numeric()
  statistics <- numeric()
  
  ctrl_group <- ctrl_group+1
  exp_group <- exp_group+1
  if(method=='mod.t.test'){
    
    g2 <- rep("group 2", length(c(ctrl_group, exp_group)))
    g2[(exp_group-2)] <- "group 1"
    g2 <- factor(g2)
    print(g2)
    mod_t_test <- tryCatch({mod.t.test(as.matrix(all_exp_data[-c(1:2)]), group = g2)},
                           error=function(e){return(NULL)})
    
    if(!is.null(mod_t_test)){
      p_value <- mod_t_test$p.value
      statistics <- mod_t_test$t
    }else{
      p_value <- NA
      statistics <- NA
    }
    
  }
  for (a  in 1:nrow(all_exp_data)) {
    lipid[a] <- all_exp_data$feature[a]
    type[a] <- all_exp_data$type[a]
    
    mean_ctrl[a] <- all_exp_data[a, ctrl_group] %>% unlist() %>% mean(na.rm=T)
    mean_exp[a] <- all_exp_data[a, exp_group] %>% unlist() %>% mean(na.rm=T)
    mean_all[a] <- all_exp_data[a, c(ctrl_group, exp_group)] %>% unlist() %>% mean(na.rm=T)
    
    FC[a] <- mean_exp[a]/mean_ctrl[a]
    sd_ctrl[a] <- all_exp_data[a, ctrl_group] %>% unlist() %>% sd(na.rm=T)
    sd_exp[a] <- all_exp_data[a, exp_group] %>% unlist() %>% sd(na.rm=T)
    
    if(method=='t.test'){
      p_value[a] <- tryCatch({t.test(unlist(all_exp_data[a, exp_group]), unlist(all_exp_data[a, ctrl_group]), var.equal=T)$p.value},
                             error=function(e){return(NA)})
      statistics[a] <- tryCatch({t.test(unlist(all_exp_data[a, exp_group]), unlist(all_exp_data[a, ctrl_group]), var.equal=T)$statistic},
                                error=function(e){return(NA)})
    }else if(method=='wilcox.test'){
      p_value[a] <- tryCatch({wilcox.test(unlist(all_exp_data[a, exp_group]), unlist(all_exp_data[a, ctrl_group]), var.equal=T)$p.value},
                             error=function(e){return(NA)})
      statistics[a] <- tryCatch({wilcox.test(unlist(all_exp_data[a, exp_group]), unlist(all_exp_data[a, ctrl_group]), var.equal=T)$statistic},
                                error=function(e){return(NA)})
    }
    
  }
  
  result <- data.frame(lipid=lipid, type=type, mean_all=mean_all, mean_ctrl=mean_ctrl, mean_exp=mean_exp, sd_ctrl=sd_ctrl,
                       sd_exp=sd_exp, FC=FC, log2FC=log2(FC), statistics=statistics, p_value=p_value, mlog10p=-log10(p_value))
  
  result <- result %>% group_by(type) %>% 
    mutate(adj_p_value = stats::p.adjust(p_value, method = 'fdr')) %>% 
    mutate(mlog10padj=-log10(adj_p_value))
  
  
  if(significant=='p_value'){
    result <- result %>% mutate(sig=ifelse(p_value<0.05, 'yes', 'no'))
  }else{
    result <- result %>% mutate(sig=ifelse(adj_p_value<0.05, 'yes', 'no'))
  }
  
  return(list(all_exp_data, result))
}

build_FA_net <- function(FA_network, unprocessed_data_result){
  
  non_processed_data_result <- unprocessed_data_result
  all_FA <- non_processed_data_result[[2]] %>% 
    filter(type=='FA') %>% .$lipid %>% unique()

  FA_below_16 <- min(which(FA_network$S1[1:7] %in%all_FA))
  
  if(is.infinite(FA_below_16)){
    FA_network <- FA_network[-c(1:7),]
    endo_start <- data.frame(S1=str_c('endo: 2:0;0--14:0;0'),
                             P1='16:0;0',
                             S1_detail='endo: 2:0;0--14:0;0',
                             P1_detail='16:0;0',pathway='Non_essential_FA_synthesis')
    FA_network <- FA_network %>% 
      rbind(endo_start)
    exo <- ''
  }
  else{
    exo <- FA_network[FA_below_16:7,]
    exo$S1 <- str_c('exo: ',exo$S1)
    exo$P1 <- str_c('exo: ',exo$P1)
    exo$S1_detail <- str_c('exo: ',exo$S1_detail)
    exo$P1_detail <- str_c('exo: ',exo$P1_detail)
    endo <- data.frame(lapply(exo, function(x){gsub('exo','endo',x)}))
    endo_start <- data.frame(S1=str_c('endo: 2:0;0--',2*(FA_below_16-1),':0;0'),
                             P1=endo[1,1],
                             S1_detail=str_c('endo: 2:0;0--',2*(FA_below_16-1),':0;0'),
                             P1_detail=endo[1,1],pathway='Non_essential_FA_synthesis')
    
    FA_network <- FA_network[-c(1:7),] %>% 
      rbind(exo) %>% rbind(endo) %>% rbind(endo_start)
    FA_network <- data.frame(lapply(FA_network, function(x){gsub('exo: 16:0;0|endo: 16:0;0','16:0;0',x)}))
  }
 
  return(FA_network)
  
}

FA_sub_transform <- function(FA_network, unprocessed_data_result,
                             unmapped_FA=NULL){
  
  non_processed_data_result <- unprocessed_data_result
  all_FA <- non_processed_data_result[[2]] %>% 
    filter(type=='FA') %>% .$lipid %>% unique()
  
  endo_start <- FA_network %>% filter(str_detect(S1, 'endo: 2:0;0'))
  exo <- FA_network %>% filter(str_detect(S1, 'exo: ')) %>% .[1,]

  graph <- graph_from_data_frame(FA_network[c('S1','P1')], directed = T,
                                 vertices = unique(c(FA_network$S1,FA_network$P1)))
  
  
  De_novo <- FA_network %>% filter(pathway=='Non_essential_FA_synthesis') %>% unlist()
  Omega6 <- FA_network %>% filter(pathway=='Omega_6_FA_synthesis') %>% unlist()
  Omega3 <- FA_network %>% filter(pathway=='Omega_3_FA_synthesis') %>% unlist()
  
  lipid <- character()
  G3P_start_path <- list()
  lipid_list <- unique(c(FA_network$S1,FA_network$P1))
  num_class <- 1
  num_path <- 1
  
  while(!is.na(lipid_list[num_class])){
    if(str_detect(lipid_list[num_class], 'exo|w6-18:2;0|w3-18:3;0')){
      lipid[num_path] <- lipid_list[num_class]
      G3P_start_path[[num_path]] <- lipid_list[num_class]
      num_path <- num_path+1
      num_class <- num_class+1
      next
    }
    if(lipid_list[num_class] %in% De_novo){
      start=c(exo[1,1],endo_start[1,1])
    }
    else if(lipid_list[num_class] %in% Omega6){
      start='w6-18:2;0'
    }
    else if(lipid_list[num_class] %in% Omega3){
      start='w3-18:3;0'
    }
    for(starts in start){
      
      
      shortest_path <- tryCatch({all_simple_paths(graph,starts,
                                                  lipid_list[num_class],
                                                  mode='out')},
                                error=function(e){NULL})
      
      
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
      
    }
    num_class <- num_class+1
  }
  
  names(G3P_start_path) <- lipid
  FA_substructure <- plyr::ldply(G3P_start_path, rbind) %>% unique()
  
  FA_substructure[is.na(FA_substructure)] <- ''
  
  colnames(FA_substructure) <- c('FA', str_c('Unit', 1:(ncol(FA_substructure)-1)))
  
  
  FA_substructure <- FA_substructure %>%  filter(Unit1!='')
  FA_substructure$FA <- str_extract(FA_substructure$FA, '\\d+:.+')
  
  special_FA <- all_FA[all_FA!='']
  special_FA <- special_FA[!special_FA %in% FA_substructure$FA]
  
  FA_substructure <- FA_substructure %>% bind_rows(data.frame(FA=special_FA, Unit1=special_FA))
  
  FA_substructure[is.na(FA_substructure)]=''
  
  if(!is.null(unmapped_FA)){
    rm_FA <- which(apply(FA_substructure, MARGIN = 1, FUN = function(x){last(x[x!=''])}) %in% unmapped_FA)
    
    FA_substructure <- FA_substructure[-rm_FA,]
  }
  
  return(FA_substructure)
  
}

FA_sub_extract <- function(char_table, FA_substructure,
                           unprocessed_data_result,
                           exact_FA='no', exo_lipid=NULL){
  
  if(class(unprocessed_data_result)=='list'){
    unprocessed_data_result <- unprocessed_data_result[[2]] %>% filter(type=='FA')
  }
  non_processed_data_result <- unprocessed_data_result %>% filter(lipid %in% FA_substructure$FA)
  if(exact_FA=='no'){
    FA_t_test <- FA_substructure[-1] %>% t() %>% as.data.frame() %>% 
      apply(MARGIN = 2, FUN = function(x){as.data.frame(x) %>% filter(x!='') %>% 
          mutate(Lipid=str_extract(x, '\\d+:\\d+;\\d+')) %>% #modify
          left_join(non_processed_data_result,by=c('Lipid'='lipid'))})
  }
  else{
    FA_t_test <- FA_substructure[-1] %>% t() %>% as.data.frame() %>% 
      apply(MARGIN = 2, FUN = function(x){as.data.frame(x) %>% filter(x!='') %>% 
          mutate(Lipid=x) %>% #modify
          mutate(Lipid=str_replace_all(Lipid, 'endo: ','')) %>% 
          mutate(Lipid=str_replace_all(Lipid, 'exo: ','')) %>% 
          left_join(non_processed_data_result,by=c('Lipid'='lipid'))})
  }
  
  if(!is.null(exo_lipid)){
    exo_lipid_with_neighbor <- unlist(filter(FA_network, S1==exo_lipid | P1==exo_lipid)[c('S1','P1')]) %>% unique()
  }
  else{
    exo_lipid_with_neighbor <- ''
  }
  FA_sub_stop <- list()
  FA_name <- character()
  for(num in 1:length(FA_t_test)){
    FA_sub <- FA_t_test[[num]]
    FA_name[num] <- last(FA_sub$x)
    
    if(is.na(last(FA_sub$log2FC))){
      FA_sub_stop[[num]] <- ''
    }
    else if(last(FA_sub$x)=='16:0;0'){
      FA_sub_stop[[num]] <- FA_sub$x
    }
    
    else if(exo_lipid_with_neighbor!='' && last(FA_sub$x) %in% exo_lipid_with_neighbor){
      
      FA_sub_stop[[num]] <- last(FA_sub$x)
    }
    else{
      if(last(FA_sub$log2FC>0)){
        stop_point <- which(FA_sub$log2FC<0)
        if(length(stop_point)==0){
          FA_sub_stop[[num]] <- FA_sub$x
        }
        else if(max(stop_point)==(nrow(FA_sub)-1)){
          FA_sub_stop[[num]] <- ''
        }
        else{
          stop_loc <- max(stop_point)+1
          FA_sub_stop[[num]] <- FA_sub$x[stop_loc:length(FA_sub$x)]
        }
        
      }
      else{
        stop_point <- which(FA_sub$log2FC>0)
        if(length(stop_point)==0){
          FA_sub_stop[[num]] <- FA_sub$x
        }
        else if(max(stop_point)==(nrow(FA_sub)-1)){
          FA_sub_stop[[num]] <- ''
        }
        else{
          stop_loc <- max(stop_point)+1
          FA_sub_stop[[num]] <- FA_sub$x[stop_loc:length(FA_sub$x)]
        }
      }
    }
    names(FA_sub_stop) <- FA_name #modify
  }
  FA_sub_stop <- plyr::ldply(FA_sub_stop, rbind)
  FA_sub_stop[is.na(FA_sub_stop)] <- ''
  
  colnames(FA_sub_stop) <- c('FA', str_c('Unit', 1:(ncol(FA_sub_stop)-1)))
  FA_sub_stop <- unique(FA_sub_stop) %>% filter(Unit1!='')
  
  exist_FA <- FA_substructure %>% filter(FA %in% non_processed_data_result$lipid) %>% 
    apply(MARGIN = 1, FUN = function(x){last(x[x!=''])}) %>% unique()
  FA_sub_exist_FA <- FA_sub_stop %>% apply(MARGIN = 1, FUN = function(x){last(x[x!=''])}) %>% 
    unique()
  lost_FA <- exist_FA[!exist_FA%in%FA_sub_exist_FA]
  
  if(length(lost_FA)!=0){
    FA_sub_stop <- FA_sub_stop %>% add_row(FA = lost_FA, Unit1 =lost_FA)
    FA_sub_stop[is.na(FA_sub_stop)] <- ''
  }
  if(exact_FA=='no'){
    FA_sub_stop <- FA_sub_stop %>% 
      mutate(FA=str_extract(FA, '\\d+:\\d+;\\d+'))
  }
  else{
    FA_sub_stop <- FA_sub_stop %>% 
      mutate(FA=ifelse(str_detect(FA, '(exo)|(endo)'), str_extract(FA, '\\d+:\\d+;\\d+'), FA))
  }
  
  FA_substructure_transform <- function(char_table, FA_substructure){
    
    #Filter lipid and FA with substructure
    char <- char_table %>% filter(FA_split!='')
    
    lipid_name <- list()
    lipid_sub_list <- list()
    
    add <- 1
    
    for (var1 in 1:nrow(char)) {
      char[var1, 'feature']
      FA <- str_split(char[var1, 'FA_split'], '_') %>% unlist()
      FA_sub <- FA_sub_trans(char[var1, 'feature'], FA, FA_substructure)
      
      
      lipid_name[[add]] <- FA_sub[[1]]
      lipid_sub_list[[add]] <- FA_sub[[2]]
      add <- add+1
      
    }
    lipid_name <- unlist(lipid_name, recursive = F)
    lipid_sub_list <- unlist(lipid_sub_list, recursive = F)
    
    return(list(lipid_name, lipid_sub_list))
  }
  FA_sub_trans <- function(lipid, FA_list, FA_substructure){
    add <- 1
    #FA substructure: all possible path method
    all_FA_sub <- list()
    for(FA_num in 1:length(FA_list)){
      FA_sub_all <- FA_substructure %>% filter(FA==FA_list[FA_num])
      
      if(nrow(FA_sub_all)==0){
        FA_sub_all <- data.frame(FA=FA_list[FA_num],Unit1=FA_list[FA_num])
      }
      all_FA_sub[[FA_num]] <- FA_sub_all
    }
    
    lipid_name <- character()
    FA_sub_list <- list()
    
    for(var1 in 1:nrow(all_FA_sub[[1]])){
      if(length(all_FA_sub)>1){
        for(var2 in 1:nrow(all_FA_sub[[2]])){
          if(length(all_FA_sub)>2){
            for(var3 in 1:nrow(all_FA_sub[[3]])){
              if(length(all_FA_sub)>3){
                for(var4 in 1:nrow(all_FA_sub[[4]])){
                  FA1_sub <- all_FA_sub[[1]][var1,] %>% dplyr::select(-FA) %>% unlist()
                  FA1_sub <- FA1_sub[FA1_sub!='']
                  
                  FA2_sub <- all_FA_sub[[2]][var2,] %>% dplyr::select(-FA) %>% unlist()
                  FA2_sub <- FA2_sub[FA2_sub!='']
                  
                  FA3_sub <- all_FA_sub[[3]][var3,] %>% dplyr::select(-FA) %>% unlist()
                  FA3_sub <- FA3_sub[FA3_sub!='']
                  
                  FA4_sub <- all_FA_sub[[4]][var4,] %>% dplyr::select(-FA) %>% unlist()
                  FA4_sub <- FA4_sub[FA4_sub!='']
                  
                  FA_sub <- c(FA1_sub, FA2_sub, FA3_sub, FA4_sub)
                  
                  names(FA_sub) <- c(rep('FA1', length(FA1_sub)), 
                                     rep('FA2', length(FA2_sub)),
                                     rep('FA3', length(FA3_sub)),
                                     rep('FA4', length(FA4_sub)))
                  
                  lipid_name[add] <- lipid
                  FA_sub_list[[add]] <- FA_sub
                  add <- add+1
                }
              }
              else{
                FA1_sub <- all_FA_sub[[1]][var1,] %>% dplyr::select(-FA) %>% unlist()
                FA1_sub <- FA1_sub[FA1_sub!='']
                
                FA2_sub <- all_FA_sub[[2]][var2,] %>% dplyr::select(-FA) %>% unlist()
                FA2_sub <- FA2_sub[FA2_sub!='']
                
                FA3_sub <- all_FA_sub[[3]][var3,] %>% dplyr::select(-FA) %>% unlist()
                FA3_sub <- FA3_sub[FA3_sub!='']
                
                
                FA_sub <- c(FA1_sub, FA2_sub, FA3_sub)
                
                names(FA_sub) <- c(rep('FA1', length(FA1_sub)), 
                                   rep('FA2', length(FA2_sub)),
                                   rep('FA3', length(FA3_sub)))
                
                lipid_name[add] <- lipid
                FA_sub_list[[add]] <- FA_sub
                add <- add+1
              }
              
              
            }  
          }
          else{
            FA1_sub <- all_FA_sub[[1]][var1,] %>% dplyr::select(-FA) %>% unlist()
            FA1_sub <- FA1_sub[FA1_sub!='']
            FA2_sub <- all_FA_sub[[2]][var2,] %>% dplyr::select(-FA) %>% unlist()
            FA2_sub <- FA2_sub[FA2_sub!='']
            
            FA_sub <- c(FA1_sub, FA2_sub)
            
            names(FA_sub) <- c(rep('FA1', length(FA1_sub)), 
                               rep('FA2', length(FA2_sub)))
            
            lipid_name[add] <- lipid
            FA_sub_list[[add]] <- FA_sub
            add <- add+1
          }
          
        }
      }
      else{
        FA_sub <- all_FA_sub[[1]][var1,] %>% dplyr::select(-FA) %>% unlist()
        #20:2;0: 18:1;0,9, 18:2;0,6,9, 20:2;0,8,11
        FA_sub <- FA_sub[FA_sub!='']
        
        names(FA_sub) <- c(rep('FA1', length(FA_sub)))
        
        lipid_name[add] <- lipid
        FA_sub_list[[add]] <- FA_sub
        add <- add+1
      }
    }
    return(list(lipid_name, FA_sub_list))
  }
  
  FA_sub_stop <- FA_substructure_transform(char_table, FA_sub_stop)
  
  return(FA_sub_stop)
}


lipid_sub_matrix <- function(exp_data, sub_data,
                                      sub_type='species'){
  sub_mat <- list()
  
  if(sub_type=='FA'){
    sub_data[[2]] <- sub_data[[2]] %>% map(.f = function(x){x[str_detect(names(x), 'FA')]})
  }
  else if(sub_type=='Class'){
    sub_data[[2]] <- sub_data[[2]] %>% map(.f = function(x){x <- x[!str_detect(names(x), 'FA')]
    x <- unique(str_replace(x, '_.+', ''))})
  }
  else{
    sub_data[[2]] <- sub_data[[2]] %>% map(.f = function(x){x <- x[!str_detect(names(x), 'FA')]})
  }
  
  all_sub <- sub_data[[2]] %>% unlist() %>% unique() %>% sort()
  
  for (num in 1:nrow(exp_data)) {
    
    loc <- which(sub_data[[1]]==exp_data[num,1])
    if(length(loc)!=0){
      
      #Substructure weight: all possible path divided by path number
      sub_mat[[num]] <-  as.integer(table(factor(unlist(sub_data[[2]][loc]), levels = all_sub)))/length(loc)
      
    }else{
      sub_mat[[num]] <- integer(length(all_sub))
    }
  }
  
  sub_mat <- sub_mat %>% as.data.frame()
  
  colnames(sub_mat) <- exp_data$feature
  
  rownames(sub_mat) <- sub_data[[2]] %>% unlist() %>% unique() %>% sort()
  
  sub_mat <- sub_mat[rowSums(sub_mat)!=0,]
  
  sub_mat <- as.matrix(sub_mat)
  
  exp_mat <- exp_data[-1] %>% as.matrix()
  
  sub_exp <- sub_mat %*% exp_mat 
  
  return(list(sub_mat,exp_mat, sub_exp))
}


t_test <- function(data, ctrl, exp, method='t.test', significant='p_value'){
  lipid <- character()
  mean_ctrl <- numeric()
  mean_exp <- numeric()
  mean_all <- numeric()
  
  FC <- numeric()
  
  sd_ctrl <- numeric()
  sd_exp <- numeric()
  p_value <- numeric()
  statistics <- numeric()
  
  if(method=='mod.t.test'){
    
    g2 <- rep("group 2", length(c(ctrl, exp)))
    g2[exp] <- "group 1"
    g2 <- factor(g2)
    data_m <- as.matrix(data)
    mod_t_test <- tryCatch({mod.t.test(data_m, group = g2)},
                           error=function(e){return(NULL)})
    
    if(!is.null(mod_t_test)){
      p_value <- mod_t_test$p.value
      statistics <- mod_t_test$t
    }else{
      p_value <- NA
      statistics <- NA
    }
    
  }
  
  
  for (a  in 1:nrow(data)) {
    lipid[a] <- rownames(data)[a]
    mean_ctrl[a] <- data[a, ctrl] %>% unlist() %>% mean(na.rm=T)
    mean_exp[a] <- data[a, exp] %>% unlist() %>% mean(na.rm=T)
    mean_all[a] <- data[a, c(ctrl, exp)] %>% unlist() %>% mean(na.rm=T)
    FC[a] <- mean_exp[a]/mean_ctrl[a]
    sd_ctrl[a] <- data[a, ctrl] %>% unlist() %>% sd(na.rm=T)
    sd_exp[a] <- data[a, exp] %>% unlist() %>% sd(na.rm=T)
    if(method=='t.test'){
      p_value[a] <- tryCatch({t.test(unlist(data[a, exp]), unlist(data[a, ctrl]), var.equal=T)$p.value},
                             error=function(e){return(NA)})
      statistics[a] <- tryCatch({t.test(unlist(data[a, exp]), unlist(data[a, ctrl]), var.equal=T)$statistic},
                                error=function(e){return(NA)})
    }
    else if(method=='wilcox.test'){
      p_value[a] <- tryCatch({wilcox.test(unlist(data[a, exp]), unlist(data[a, ctrl]), var.equal=T)$p.value},
                             error=function(e){return(NA)})
      statistics[a] <- tryCatch({wilcox.test(unlist(data[a, exp]), unlist(data[a, ctrl]), var.equal=T)$statistic},
                                error=function(e){return(NA)})
    }
    
  }
  adj_p_value <- p.adjust(p_value, method = 'fdr', n = length(p_value))
  result <- data.frame(lipid=lipid, mean_all=mean_all, mean_ctrl=mean_ctrl, mean_exp=mean_exp, sd_ctrl=sd_ctrl,
                       sd_exp=sd_exp, FC=FC, log2FC=log2(FC), statistics=statistics, p_value=p_value, mlog10p=-log10(p_value),
                       adj_p_value=adj_p_value, mlog10padj=-log10(adj_p_value))
  
  if(significant=='p_value'){
    result <- result %>% mutate(sig=ifelse(p_value<0.05, 'yes', 'no'))
  }
  else{
    result <- result %>% mutate(sig=ifelse(adj_p_value<0.05, 'yes', 'no'))
  }
  
  result$lipid <- str_replace(result$lipid, '_FA\\d+','')
  
  return(result)
}

draw_network <- function(network_data, DE_data, if_species=F, significant='p_value',
                         path_scoring_result,  reaction_scoring_result,
                         top_n=3, path_type='both'){
  
  if(path_type=='both'){
    active=T
    suppressed=T
  }
  else if(path_type=='active'){
    active=T
    suppressed=F
  }
  else{
    active=F
    suppressed=T
  }
  FA_node <- data.frame(lipid=unique(c(network_data$S1,network_data$P1)), stringsAsFactors = F) %>% 
      left_join(DE_data, by=c('lipid'))

  
  replace_inf <- max(abs(FA_node$log2FC[!is.infinite(FA_node$log2FC)]), na.rm = T)+0.5
  
  if(significant=='p_value'){
    FA_node <- FA_node %>% mutate(sig=ifelse(p_value<0.05, 'yes','no'))
  }
  else{
    FA_node <- FA_node %>% mutate(sig=ifelse(adj_p_value<0.05, 'yes','no')) %>% 
      mutate(mlog10p=-log10(adj_p_value))
  }
  
  FA_node <- FA_node %>% 
    mutate(in_ref=ifelse(!is.na(FC), 'yes','no')) %>% 
    mutate(sig=ifelse(is.na(sig),'na',sig)) %>% 
    mutate(log2FC=ifelse(is.infinite(log2FC) & log2FC>0, replace_inf, log2FC)) %>% 
    mutate(log2FC=ifelse(is.infinite(log2FC) & log2FC<0, -replace_inf, log2FC)) %>% 
    mutate(mlog10p=ifelse(is.na(mlog10p), 0, mlog10p)) %>% 
    mutate(mean_all=ifelse(is.na(mean_all), 0, mean_all))
  
  FC_color <- cut(FA_node$log2FC, breaks = seq(-replace_inf-0.1, replace_inf+0.1, length.out=100), labels = bluered(99)) %>% 
    as.character()
  
  
  
  FC_color[is.na(FC_color)] <- 'black'
  
  
  node <- data.frame(id=FA_node$lipid,
                     label=FA_node$lipid,
                     color.background=FC_color,
                     color.border=recode(FA_node$sig,
                                         "yes"='purple',
                                         "no"='black',
                                         'na'='transparent'),
                     borderWidth=recode(FA_node$sig,
                                        "yes"=2.5,
                                        "no"=1,
                                        'na'=0),
                     #value=FA_node$mlog10p,
                     value=(FA_node$mlog10p),
                     font.size = 35)
  edge <- data.frame(from=network_data$S1,to=network_data$P1,
                     #label = paste("Edge", 1:8), 
                     color='gray',
                     arrows=c('to'),
                     length=100)
  
  path_color <- function(network_edge, path_scoring_result,
                         reaction_scoring_result, top_n=3,
                         active=T, suppressed=T){
    
    if(active==T){
      #top paths
      
      sig_active_path <- path_scoring_result %>% filter(Type=='Active') %>% 
        .[!duplicated(.$rep_sub_path),] %>% filter(Significant=='yes')
      if(nrow(sig_active_path)==0){
        top_n_new <- nrow(sig_active_path)
        topN_active_path <- data.frame()
      }
      else if(nrow(sig_active_path)<top_n){
        top_n_new <- nrow(sig_active_path)
        topN_active_path <- sig_active_path[1:top_n_new,] %>% 
          mutate(rank=1:top_n_new)
      }
      else{
        top_n_new <- top_n
        topN_active_path <- sig_active_path[1:top_n_new,] %>% 
          mutate(rank=1:top_n_new)
      }

      
      #top edges
      
      sig_active_edge <- reaction_scoring_result %>% filter(Mode=='Increase') %>% 
        filter(p_value<0.05)
      
      if(nrow(sig_active_edge)==0){
        topN_active_edge <- data.frame()
      }
      else if(nrow(sig_active_edge)<top_n){
        top_n_new <- nrow(sig_active_edge)
        topN_active_edge <- sig_active_edge[1:top_n_new,] %>% 
          mutate(rank=1:top_n_new)
      }
      else{
        top_n_new <- top_n
        topN_active_edge <- sig_active_edge[1:top_n_new,] %>% 
          mutate(rank=1:top_n_new)
      }
    }
    else{
      topN_active_path <- data.frame()
      topN_active_edge <- data.frame()
    }
    
    if(suppressed==T){
      #top paths
      
      sig_suppressed_path <- path_scoring_result %>% filter(Type=='Suppressed') %>% 
        arrange(cal_score) %>%.[!duplicated(.$rep_sub_path),] %>% 
        filter(Significant=='yes')
      if(nrow(sig_suppressed_path)==0){
        topN_suppressed_path <- data.frame()
      }
      else if(nrow(sig_suppressed_path)<top_n){
        top_n_new <- nrow(sig_suppressed_path)
        topN_suppressed_path <- sig_suppressed_path[1:top_n_new,] %>% 
          mutate(rank=c(rep(5, top_n_new)+1:top_n_new))
      }
      else{
        top_n_new <- top_n
        topN_suppressed_path <- sig_suppressed_path[1:top_n_new,] %>% 
          mutate(rank=c(rep(5, top_n_new)+1:top_n_new))
      }

      
      #top edges
      
      sig_suppressed_edge <- reaction_scoring_result %>% filter(Mode=='Decrease') %>% 
        arrange(perturbation_score) %>% filter(p_value<0.05)
      if(nrow(sig_suppressed_edge)==0){
        topN_suppressed_edge <- data.frame()
      }
      else if(nrow(sig_suppressed_edge)<top_n){
        top_n_new <- nrow(sig_suppressed_edge)
        topN_suppressed_edge <- sig_suppressed_edge[1:top_n_new,] %>% 
          mutate(rank=c(rep(5, top_n_new)+1:top_n_new))
      }
      else{
        top_n_new <- top_n
        topN_suppressed_edge <- sig_suppressed_edge[1:top_n_new,] %>% 
          mutate(rank=c(rep(5, top_n_new)+1:top_n_new))
      }
    }
    else{
      topN_suppressed_path <- data.frame()
      topN_suppressed_edge <- data.frame()
    }
      
        
    topN_path <- rbind(topN_active_path, topN_suppressed_path)
   
    topN_edge <- rbind(topN_active_edge, topN_suppressed_edge)
    
    path <- data.frame(rank=1:10, edge_color=bluered(100)[c(seq(100,56, length.out=5),seq(1,45, length.out=5))],
                       width=rep(seq(12,4, length.out=5),2))
    
    #
    edge <- data.frame(rank=1:10, perturb_type=c(rep('Increase',5),rep('Decrease',5)),
                       fontsize=rep(rep(40,5),2),
                       fontcolor=c(rep('red',5),rep('blue',5)))
    if(nrow(topN_path)!=0){
      path_score <- topN_path %>% left_join(path, by='rank')
    }
    else{
      path_score <- data.frame()
    }
    if(nrow(topN_edge)!=0){
      perturbation_score <- topN_edge %>% left_join(edge, by='rank')
    }
    else{
      perturbation_score <- data.frame()
    }

    edge_color <- character()
    width <- numeric()
    perturb_type <- character()
    fontsize <- numeric()
    fontcolor <- character()
    
    for (edge_num in 1:nrow(network_edge)) {
      from <- network_edge$from[edge_num]
      to <- network_edge$to[edge_num]
      
      edge_in_path <- str_detect(path_score$path, str_c(from,' --> ',to))
      
      if(sum(edge_in_path)==0){
        edge_color[edge_num] <- 'black'
        width[edge_num] <- 1
      }
      else{
        names(edge_in_path) <- path_score$rank
        path_rank <- names(edge_in_path)[which.max(edge_in_path)]
        edge_color[edge_num] <- filter(path_score, rank==path_rank)$edge_color
        width[edge_num] <- filter(path_score, rank==path_rank)$width
      }
      
      edge_in_perturb <- perturbation_score$edge_name==str_c(from,' --> ',to)
      
      if(sum(edge_in_perturb)==0){
        perturb_type[edge_num] <- ''
        fontsize[edge_num] <- 0
        fontcolor[edge_num] <- 'black'
      }
      else{
        
        names(edge_in_perturb) <- perturbation_score$rank
        perturb_rank <- names(edge_in_perturb)[which.max(edge_in_perturb)]
        perturb_type[edge_num] <- filter(perturbation_score, rank==perturb_rank)$perturb_type
        fontsize[edge_num] <- filter(perturbation_score, rank==perturb_rank)$fontsize
        fontcolor[edge_num] <- filter(perturbation_score, rank==perturb_rank)$fontcolor
        
      }
    }
    
    network_edge <- network_edge %>% mutate(color=edge_color, 
                                            width=width,
                                            label=perturb_type,
                                            font.size=fontsize,
                                            font.color=fontcolor)
    return(network_edge)
  }
  
  if(if_species==T){
    
    if(active==T && suppressed==T){
      topN_rep_path <- c(path_scoring_result %>% filter(Type=='Active') %>% 
                           .[!duplicated(.$rep_sub_path),] %>% .[1:top_n,] %>% .$rep_sub_path, 
                         path_scoring_result %>% filter(Type=='Suppressed') %>% 
                           arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>%
                           .[1:top_n,] %>% .$rep_sub_path)
    }
    else if(active==T && suppressed==F){
      topN_rep_path <- path_scoring_result %>% filter(Type=='Active') %>% 
        .[!duplicated(.$rep_sub_path),] %>% .[1:top_n,] %>% .$rep_sub_path
    }
    else if(active==F && suppressed==T){
      topN_rep_path <-path_scoring_result %>% filter(Type=='Suppressed') %>% 
        arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>%
        .[1:top_n,] %>% .$rep_sub_path
    }

    topN_rep_path <- path_scoring_result %>% filter(rep_sub_path%in%topN_rep_path, Significant=='yes') %>% 
      .$path %>% str_split(' --> ') %>% unlist() %>% unique()
    
    edge_in_topN_path <- reaction_scoring_result$edge_name %>% str_split(' --> ') %>% 
      map_lgl(.f = function(x){x[1] %in% topN_rep_path && x[2] %in% topN_rep_path})
    edge_in_topN_path <- reaction_scoring_result[edge_in_topN_path,]$edge_name
    
    topN_net_edge <- network_data %>% filter(S1 %in% topN_rep_path, P1 %in% topN_rep_path)
    
    topN_net_edge <- topN_net_edge %>% mutate(color='gray', arrows='to',length=100)
    
    colnames(topN_net_edge)[1:2] <- c('from','to')
    
    edge <- topN_net_edge
    
    topN_net_node <- node %>% 
      filter(id %in% unlist(topN_net_edge))
    

    reaction_scoring_result <- reaction_scoring_result%>% 
      filter(edge_name %in% edge_in_topN_path)

  }
  
  edge <- path_color(edge, path_scoring_result,
                     reaction_scoring_result, top_n,
                     active, suppressed)
  
  if(if_species==T){
    node <- topN_net_node %>% filter(id %in%c(edge$from, edge$to))
  }
  return(list(node, edge))
  
}

path_scoring <- function(network, sub_t,  calibrate=T, data_type='FA'){
  options(scipen = 999)
  
  sub_t <- sub_t %>% dplyr::select(-sig)
  sub_t <- sub_t %>% 
    mutate(zscore=qnorm(1-sub_t$p_value/2)) %>% 
    mutate(zscore=ifelse(is.infinite(zscore), 8, zscore)) %>% 
    mutate(zscore=ifelse(statistics>0, zscore, -zscore))
  
  Network_w <- data.frame(node=unique(unlist(network[c('S1','P1')]))) %>% 
    left_join(sub_t, by=c('node'='lipid')) %>% 
    mutate(weight=zscore) %>% 
    dplyr::select(node, weight)
  
  graph <- graph_from_data_frame(network[c('S1','P1')], directed = T, vertices = unique(c(network$S1,network$P1)))
  
  all_node <- network[c('S1','P1')] %>%  unlist() %>% unique()
  
  combo <- permutations(length(all_node), 2, all_node)
  path <- character()
  path_full <- character()
  score <- numeric()
  type <- character()
  from <- character()
  to <- character()
  
  num_path <- 1
  for (num in 1:nrow(combo)) {
    
    if((combo[num,1] %in% sub_t$lipid) && (combo[num,2] %in% sub_t$lipid)){
      
      all_path <- all_simple_paths(graph, combo[num,1], combo[num,2], mode='out')
      
      
      if(length(all_path)!=0){
        for (path_num in 1:length(all_path)) {
          
          path_node <- all_path[[path_num]] %>% attributes() %>% .$names
          #c(2:0;0, 4:0;0, 6:0;0, 8:0;0, 10:0;0, 12:0;0)
          
          path[num_path] <- str_c(path_node, collapse = ' --> ')
          #c(2:0;0_4:0;0_6:0;0_8:0;0_10:0;0_12:0;0)
          
          from[num_path] <- path_node[1]
          to[num_path] <- tail(path_node,1)
          
          S1_P1_w <- Network_w %>% filter(node %in% path_node) %>% .$weight
          score[num_path] <- sum(S1_P1_w)/sqrt(length(path_node))
          
          num_path <- num_path+1
        }
      }
    }
    
  }
  path_score <- data.frame(path=path, from=from, to=to, score=score) %>% 
    filter(!is.na(score)) %>% 
    mutate(Type=ifelse(score>0, 'Active', 'Suppressed')) %>% 
    arrange(score)
  
  
  if(calibrate==T){
    all_z <- sub_t$zscore
    max_length <- path_score$path %>% str_split(' --> ') %>% 
      map(length) %>% unlist() %>% max
    
    path <- numeric()
    mean_each_length <- numeric()
    sd_each_length <- numeric()
    
    for(path_length in 2:max_length){
      path[path_length-1] <- path_length
      
      cal <- replicate(10000, sample(all_z, path_length, replace = F)) %>% 
        colSums() 
      cal <- cal/sqrt(path_length)
      
      mean_each_length[path_length-1] <- mean(cal)
      sd_each_length[path_length-1] <- sd(cal)
    }
    
    cal_score <- numeric()
    for(each_path in 1:nrow(path_score)){
      path_length <- path_score$path[each_path] %>% str_split(' --> ') %>% 
        unlist() %>% length()
      cal_score[each_path] <- (path_score$score[each_path]-mean_each_length[path_length-1])/sd_each_length[path_length-1]
    }
    path_score <- path_score %>% mutate(cal_score=cal_score) %>% 
      dplyr::select(path, path, from, to, score, cal_score, everything())
    
    path_score <- path_score %>% 
      mutate(Significant=ifelse(abs(cal_score)>1.96, 'yes', 'no')) %>% 
      arrange(desc(cal_score)) %>% 
      mutate(Type=ifelse(cal_score>0, 'Active','Suppressed'))
  }
  else{
    path_score <- path_score %>% 
      mutate(Significant=ifelse(abs(score)>1.96, 'yes', 'no')) %>% 
      arrange(desc(score))
  }
  
  if(data_type=='FA'){
    rep_sub <- function(path_score, type){
      
      path_score <- path_score %>% filter(Type==type)
      
      if(nrow(path_score)==0){
        return(NULL)
      }
      score_list <- path_score$path
      score_list_name <- path_score$path
      
      rep_sub_path_list <- list()
      rep_sub_path_num <- 1
      rep_sub_path_name <- character()
      #path_name <- character()
      
      
      for(num in 1:length(score_list)){
        sub_path <- unlist(str_split(score_list[num], ' --> '))
        sub_path <- str_c(sub_path[-length(sub_path)], sub_path[-1], sep = '_')
        overlap_prop <- rep_sub_path_list %>% map_dbl(.f = function(x){sum(sub_path %in% x)/length(sub_path)})
        
        if(sum(overlap_prop>0.5)==0){
          #path_name[num] <- score_list_name[num]
          
          rep_sub_path_list[[rep_sub_path_num]] <- sub_path
          names(rep_sub_path_list)[rep_sub_path_num] <- score_list_name[num]
          
          rep_sub_path_name[num] <- score_list_name[num]
          rep_sub_path_num <- rep_sub_path_num+1
        }
        else{
          rep_sub_path_name[num] <- names(rep_sub_path_list)[which.max(overlap_prop)] 
          #path_name[num] <- score_list_name[num]
        }
      }
      return(rep_sub_path_name)
    }
    
    rep_sub_path <- c(rep_sub(path_score ,'Active'), rev(rep_sub(path_score[nrow(path_score):1,], 'Suppressed')))
    path_score <- path_score %>% mutate(rep_sub_path=rep_sub_path) 
  }
  else if(data_type=='Class'){
    rep_sub <- function(path_score, type){
      
      path_score <- path_score %>% filter(Type==type)
      if(nrow(path_score)==0){
        return(NULL)
      }
      score_list <- path_score$path
      score_list_name <- path_score$path
      
      rep_sub_path_list <- list()
      rep_sub_path_num <- 1
      rep_sub_path_name <- character()
      #path_name <- character()
      
      
      for(num in 1:length(score_list)){
        sub_path <- unlist(str_split(score_list[num], ' --> '))
        
        sub_path <- str_c(sub_path[-length(sub_path)], sub_path[-1], sep = '_')
        
        overlap_prop <- rep_sub_path_list %>% map_dbl(.f = function(x){sum(sub_path %in% x)/length(sub_path)})
      
        if(sum(overlap_prop>=0.5)==0){
          #path_name[num] <- score_list_name[num]
          
          rep_sub_path_list[[rep_sub_path_num]] <- sub_path
          names(rep_sub_path_list)[rep_sub_path_num] <- score_list_name[num]
          
          rep_sub_path_name[num] <- score_list_name[num]
          rep_sub_path_num <- rep_sub_path_num+1
        }
        else{
          rep_sub_path_name[num] <- names(rep_sub_path_list)[which.max(overlap_prop)] 
          #path_name[num] <- score_list_name[num]
        }
      }
      return(rep_sub_path_name)
    }
    
    rep_sub_path <- c(rep_sub(path_score ,'Active'), rev(rep_sub(path_score[nrow(path_score):1,], 'Suppressed')))
    path_score <- path_score %>% mutate(rep_sub_path=rep_sub_path) 
  }
  else{
    if(sum(str_detect(path_score$path, '\\d+:\\d+;0_\\d+:\\d+;0'))!=0){
      
      rep_sub_path <- path_score$path%>% 
        str_extract_all('\\d+:\\d+;0_\\d+:\\d+;0') %>% 
        map(.f = function(x){k=table(x) %>% sort(decreasing = T);
        if(length(k)==0){''}
        else if(sum(k!=1)==0){sort(names(k), decreasing = T)[1]}
        else if(sum(k!=1)!=0){names(k)[1]}}) %>% 
        unlist
    }
    else{
      
      rep_sub_path <- path_score$path%>% 
        str_extract_all('\\d+:\\d+;\\d+') %>% 
        map(.f = function(x){k=table(x) %>% sort(decreasing = T);
        if(length(k)==0){''}
        else if(sum(k!=1)==0){sort(names(k), decreasing = T)[1]}
        else if(sum(k!=1)!=0){sort(names(k[k==max(k)]), decreasing = T)[1]}}) %>% 
        unlist
    }

    
    
    path_score <- path_score%>% 
      mutate(rep_sub_path=rep_sub_path)
    
  }
  return(path_score)
}


reaction_scoring <- function(network, sub_exp, sub_t, ctrl=1:7, 
                             exp=8:13, stat='p', Species='human'){
  
  sub_exp <- sub_exp %>% apply(MARGIN = 1, FUN = function(x){y=x; 
  y[y==0] <- min(y[y!=0])*0.5; y },simplify = T) 
  sub_exp <- as.data.frame(t(sub_exp))
  
  
  edge_name <- character(nrow(network))
  
  
  edge_type <- character(nrow(network))
  node1_log2FC <- numeric(nrow(network))
  node2_log2FC <- numeric(nrow(network))
  p_value <- numeric(nrow(network))
  statistics <- numeric(nrow(network))
  FC_ctrl <- numeric(nrow(network))
  FC_exp <- numeric(nrow(network))
  FC_exp_ctrl <- numeric(nrow(network))
  
  num <- 1
  for(edge in 1:nrow(network)){
    if((network[edge,1] %in% rownames(sub_exp)) && (network[edge,2] %in% rownames(sub_exp))){
      
      edge_name[num] <- str_c(network[edge,1], ' --> ', network[edge,2])
      
      a <- sub_t %>% filter(lipid==network[edge,1]) %>% .$log2FC
      b <- sub_t %>% filter(lipid==network[edge,2]) %>% .$log2FC
      
      node1_log2FC[num] <- a
      node2_log2FC[num] <- b
      
      if(a<0 && b>0){
        edge_type[num] <- 'Increase'
      }else if(a>0 && b<0){
        edge_type[num] <- 'Decrease'
      }else{
        edge_type[num] <- 'Non-change'
      }
      
      from <- which(rownames(sub_exp)==network[edge,1])
      to <- which(rownames(sub_exp)==network[edge,2])
      
      edge_FC <- unlist(sub_exp[to,])/unlist(sub_exp[from,])
      
      FC_ctrl[num] <- mean(edge_FC[ctrl])
      FC_exp[num] <- mean(edge_FC[exp])
      FC_exp_ctrl[num] <- FC_exp[num]/FC_ctrl[num]
      
      p_value[num] <- tryCatch({t.test(edge_FC[exp],edge_FC[ctrl], var.equal=T)$p.value},
                               error=function(e){return(NA)})
      statistics[num] <- tryCatch({t.test(edge_FC[exp],edge_FC[ctrl], var.equal=T)$statistic},
                                  error=function(e){return(NA)})
      
      num <- num+1
    }else{}
  }
  
  
  result <- data.frame(edge_name=edge_name, FC_ctrl=FC_ctrl, FC_exp=FC_exp,
                       FC_exp_ctrl=FC_exp_ctrl,statistics=statistics,
                       p_value=p_value,edge_type=edge_type, 
                       node1_log2FC=node1_log2FC,node2_log2FC=node2_log2FC) %>% 
    filter(edge_name!='') %>% 
    mutate(adj_p_value=p.adjust(.$p_value,method = 'fdr'),
           log2FC_exp_ctrl=log2(FC_exp_ctrl),
           mlog10p=-log10(p_value),mlog10padj=-log10(adj_p_value))
  
  result <- result %>% 
    mutate(perturbation_score=log2FC_exp_ctrl*mlog10p) %>% 
    mutate(Mode=ifelse(perturbation_score>0, 'Increase', 'Decrease'))
  
  result <- result %>% dplyr::select(edge_name, FC_ctrl, FC_exp,
                                     FC_exp_ctrl,statistics,
                                     p_value,mlog10p,
                                     adj_p_value,mlog10padj,
                                     perturbation_score,Mode,edge_type,
                                     node1_log2FC,node2_log2FC) %>% 
    mutate(perturbation_score=ifelse(is.na(perturbation_score), 0,perturbation_score)) %>% 
    arrange(desc(perturbation_score))
  
  
  edge_lipid1 <- str_split(result$edge_name, ' --> ') %>% map_chr(~.[1])
  edge_lipid2 <- str_split(result$edge_name, ' --> ') %>% map_chr(~.[2])
  
  
  FA_change <- map2_chr(edge_lipid1, edge_lipid2, .f = function(x,y){a=str_extract_all(x, '\\d+:\\d+;\\d+') %>% unlist();
  if(length(a)==0){a='0:0;0'};
  a1=str_extract_all(a,'\\d+:') %>% str_sub(end = -2) %>% as.integer() %>% sum;
  a2=str_extract_all(a,'\\d+;') %>% str_sub(end = -2) %>% as.integer() %>% sum;
  a3=str_extract_all(a,'\\d+$') %>%  as.integer() %>% sum;
  b=str_extract_all(y, '\\d+:\\d+;\\d+') %>% unlist();
  if(length(b)==0){b='0:0;0'};
  b1=str_extract_all(b,'\\d+:') %>% str_sub(end = -2) %>% as.integer() %>% sum;
  b2=str_extract_all(b,'\\d+;') %>% str_sub(end = -2) %>% as.integer() %>% sum;
  b3=str_extract_all(b,'\\d+$') %>%  as.integer() %>% sum;
  return(str_c(' (',abs(a1-b1),':',abs(a2-b2),';',abs(a3-b3),')'))})
  
  FA_change[FA_change==' (0:0;0)'] <- ''
  
  result <- result %>% mutate(FA_change=str_replace(FA_change,';0',''))
  
  
  lipid_reaction <- str_split(result$edge_name, ' --> ') %>% 
    map(~str_replace_all(.x, 'endo: ','') %>% 
          str_replace_all('exo: ','') %>% 
          str_replace_all('2:0;0--','')) %>% 
    map_chr(.f = function(x){a <- 
      if(str_detect(x[1], '_')){str_extract(x[1],'[A-Za-z -]+')}
    else{x[1]};
    b <- if(str_detect(x[2], '_')){str_extract(x[2],'[A-Za-z -]+')}
    else{x[2]};
    str_c(a, b, sep  = '_')})
  
  genes <- data.frame(Reaction=lipid_reaction) %>% 
    left_join(filter(reaction_gene_mapping, species==Species), by='Reaction') %>% 
    .$gene
  
  result <- result %>% mutate(genes=genes)
  
  
  return(result)
}



species_sub_transform <- function(char, lipid_substructure, network_node){

  char <- char %>% filter(class %in% lipid_substructure$Lipid)
  lipid_name <- character()
  lipid_sub_list <- list()
  
  add <- 1
  add2 <- 1
  
  FA_pool <- char$FA_split[char$FA_split!=''] %>% str_split('_') %>% unlist() %>% unique()
  
  for (var1 in 1:nrow(char)) {
    #PC, PE
    #print(var1)
    
    LMID <- filter(network_node, Abbreviation==char$class[var1])$LIPIDMAPS_subclass_id
    
    #lipid substructure: all possible path method
    lipid_sub_all <- lipid_substructure %>% filter(Lipid%in%char$class[var1]) 
    #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
    #PE: G3P, LPA_FA1, PA_FA2, DAG_FA2, PE_FA2
    if(char[var1, 'FA_num']==2){
      
      if(char[var1, 'FA_split']==''){
        
        FA2 <- char$FA_sum[var1] #c('34:1;0')
        
        #chain_db_oh <- str_extract_all(char[var1,'FA_sum'], '\\d+') %>% unlist()
        #c('34','1','0')
        if_sphingolipid <- filter(network_node, Abbreviation==char$class[var1])$Class
        
        if(if_sphingolipid %in% c('Sphingolipid','Deoxysphingolipid')){
          for(var2 in 1:nrow(lipid_sub_all)){
            lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
            #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub!='']
            
            lipid_sub[str_detect(lipid_sub, 'FA2')] <- str_c(lipid_sub[str_detect(lipid_sub, 'FA2')], FA2, sep='_')
            
            lipid_sub_FA1 <- lipid_sub
            
            if(!char$class[var1] %in% c('dhCer', 'dhSM', 'doxdhCer','doxdhSM')){
              
              length_db_oh <- str_extract_all(FA2, '\\d+') %>% unlist %>% as.integer()
              
              lipid_sub_FA1[str_detect(lipid_sub_FA1,'^dhCer')] <- str_c('dhCer_FA2_',length_db_oh[1], ':',
                                                                         length_db_oh[2]-1, ';2')
              lipid_sub_FA1[str_detect(lipid_sub_FA1,'doxdhCer')] <- str_c('doxdhCer_FA2_',length_db_oh[1], ':',
                                                                           length_db_oh[2]-1, ';1')
            }
            if(str_detect(char$class[var1], 'dox')){
              lipid_sub_FA1[str_detect(lipid_sub_FA1, 'ketoSPB')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'ketoSPB')], '18:0;0' , sep='_')
              lipid_sub_FA1[str_detect(lipid_sub_FA1, 'dhSPB')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'dhSPB')], '18:0;1' , sep='_')
            }
            else{
              lipid_sub_FA1[str_detect(lipid_sub_FA1, 'ketoSPB')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'ketoSPB')], '18:0;1' , sep='_')
              lipid_sub_FA1[str_detect(lipid_sub_FA1, 'dhSPB')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'dhSPB')], '18:0;2' , sep='_')
            }
            
            
            names(lipid_sub_FA1) <- rep('Lipid',length(lipid_sub_FA1))
            
            lipid_sub_list[[add]] <- lipid_sub_FA1
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add+1
            
          }  
        }
        else{
          LM_mapping <- lipidmaps_sum %>% filter(sub_class_id==LMID, FA_sum==str_extract(FA2, '\\d+:\\d+'))
          
          if(nrow(LM_mapping)!=0){
            possible_FA <- LM_mapping$FA_split %>% str_split('_') %>% unlist() %>% unique()
            
            possible_FA <- str_c(possible_FA, ';0')
            
            if(length(FA_pool)!=0){
              possible_FA <- possible_FA[possible_FA%in%FA_pool]
            }
            
          }
          else{
            for(var2 in 1:nrow(lipid_sub_all)){
              lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
              #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
              
              lipid_sub <- lipid_sub[lipid_sub!='']
              
              lipid_sub[str_detect(lipid_sub, 'FA2')] <- str_c(lipid_sub[str_detect(lipid_sub, 'FA2')], FA2, sep='_')
              
              names(lipid_sub) <- rep('Lipid',length(lipid_sub))
              
              lipid_sub_list[[add]] <- lipid_sub
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add+1
            } 
            next
          }
          for(var2 in 1:nrow(lipid_sub_all)){
            lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
            #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub!='']
            
            lipid_sub[str_detect(lipid_sub, 'FA2')] <- str_c(lipid_sub[str_detect(lipid_sub, 'FA2')], FA2, sep='_')
            
            last_FA1_lipid <- last(lipid_sub[str_detect(lipid_sub, 'FA1')]) %>% 
              str_replace('_FA1','')
            
            possible_FA1_lipid <- char %>% filter(class==last_FA1_lipid, FA_split %in%possible_FA)
            
            if(nrow(possible_FA1_lipid)!=0){
              possible_FA1 <- possible_FA1_lipid$FA_split
            }
            else{
              possible_FA1 <- possible_FA
            }
            for(var3 in 1:length(possible_FA1)){
              
              lipid_sub_FA1 <- lipid_sub
              
              lipid_sub_FA1[str_detect(lipid_sub_FA1, 'FA1')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'FA1')],
                                                                       possible_FA1[var3], sep='_')
              
              names(lipid_sub_FA1) <- rep('Lipid',length(lipid_sub_FA1))
              
              lipid_sub_list[[add]] <- lipid_sub_FA1
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add+1
              
            }
          }
          
        }
        
        #FA_pool_chain_db_oh <- str_extract_all(FA_pool, '\\d+')
        #c('16','0','0'), c('18','1','0')...
        
        #possible_FA <- FA_pool_chain_db_oh %>% map(.f = function(x){str_c(as.integer(chain_db_oh)-as.integer(x), collapse = ':')}) %>% 
        #  unlist()
        #c('18:1:0', '16:0:0')
        
        #str_sub(possible_FA, -2,-2) <- ';'
        #c('18:1;0', '16:0;0')
      }
      else if(char[var1, 'FA_split']!=''){
        
        FA2 <- unlist(char[var1, 'each_FA']) #c('16:0;0', '18:0;0')
        
        
        for(var2 in 1:nrow(lipid_sub_all)){
          lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
          #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub!='']
          
          lipid_sub[str_detect(lipid_sub, 'FA2')] <- str_c(lipid_sub[str_detect(lipid_sub, 'FA2')], str_c(FA2, collapse = '_') , sep='_')
          
          if_sphingolipid <- filter(network_node, Abbreviation==char$class[var1])$Class
          
          if(if_sphingolipid %in% c('Sphingolipid','Deoxysphingolipid')){
            lipid_sub_FA1 <- lipid_sub
            
            if(str_detect(char$class[var1], 'dox')){
              lipid_sub_FA1[str_detect(lipid_sub_FA1,'doxdhCer')] <- str_c('doxdhCer_FA2_18:0;1_',FA2[2])
              lipid_sub_FA1[str_detect(lipid_sub_FA1, 'ketoSPB')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'ketoSPB')], '18:0;0' , sep='_')
              lipid_sub_FA1[str_detect(lipid_sub_FA1, 'dhSPB')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'dhSPB')], '18:0;1' , sep='_')
            }
            else{
              lipid_sub_FA1[str_detect(lipid_sub_FA1,'dhCer')] <- str_c('dhCer_FA2_18:0;2_',FA2[2])
              lipid_sub_FA1[str_detect(lipid_sub_FA1, 'ketoSPB')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'ketoSPB')], '18:0;1' , sep='_')
              lipid_sub_FA1[str_detect(lipid_sub_FA1, 'dhSPB')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'dhSPB')], '18:0;2' , sep='_')
              
            }
            
            
            names(lipid_sub_FA1) <- rep('Lipid',length(lipid_sub_FA1))
            
            lipid_sub_list[[add]] <- lipid_sub_FA1
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add+1
            
            FA2 <- FA2[2]
          }
          else{
            
            if(filter(network_node, Abbreviation==char$class[var1])$Class=='Ether lipid'){
              
              lipid_sub_FA1 <- lipid_sub
              
              lipid_sub_FA1[str_detect(lipid_sub_FA1, 'FA1')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'FA1')], FA2[1] , sep='_')
              
              names(lipid_sub_FA1) <- rep('Lipid',length(lipid_sub_FA1))
              
              lipid_sub_list[[add]] <- lipid_sub_FA1
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add+1
              next
            }
            else{
              #first FA
              possible_FA1_lipid <- char %>% filter(class=='LPA', FA_split %in% FA2)
              
              if(nrow(possible_FA1_lipid)!=0){
                FA1 <- possible_FA1_lipid$FA_split
              }
              else{
                FA1 <- FA2
              }
              for(var3 in 1:length(FA1)){
                
                
                lipid_sub_FA1 <- lipid_sub
                
                lipid_sub_FA1[str_detect(lipid_sub_FA1, 'FA1')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'FA1')], FA1[var3], sep='_')
                
                names(lipid_sub_FA1) <- rep('Lipid',length(lipid_sub_FA1))
                
                lipid_sub_list[[add]] <- lipid_sub_FA1
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add+1
              }
            }
          }
        }
        
      }
    }
    else if(char[var1, 'FA_num']==1){
      
      FA1 <- unlist(char[var1, 'each_FA']) #c('16:0;0')
      
      if_sphingolipid <- filter(network_node, Abbreviation==char$class[var1])$Class
      
      
      for(var2 in 1:nrow(lipid_sub_all)){
        
        lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
        #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
        
        lipid_sub <- lipid_sub[lipid_sub!='']
        lipid_sub[str_detect(lipid_sub, 'FA1')] <- str_c(lipid_sub[str_detect(lipid_sub, 'FA1')], FA1, sep='_')
        
        if(if_sphingolipid %in% c('Sphingolipid','Deoxysphingolipid')){
          
          if(!str_detect(char$class[var1], 'dox')){
            lipid_sub <- str_replace(lipid_sub, 'ketoSPB_FA1_.+','ketoSPB_FA1_18:0;1')
            lipid_sub <- str_replace(lipid_sub, 'dhSPB_FA1_.+','dhSPB_FA1_18:0;2')
          }
          else{
            lipid_sub <- str_replace(lipid_sub, 'doxketoSPB_FA1_.+','doxketoSPB_FA1_18:0;0')
            lipid_sub <- str_replace(lipid_sub, 'doxdhSPB_FA1_.+','doxdhSPB_FA1_18:0;1')
          }
        }
        
        #LPE: LPA_FA1_16:0;0
        
        
        if(sum(str_detect(lipid_sub, 'FA2'))!=0){
          
          last_FA2_lipid <- last(lipid_sub[str_detect(lipid_sub, 'FA2')]) %>% 
            str_replace('_FA2','')
          #PC for LPC
          if(filter(network_node, Abbreviation==last_FA2_lipid)$Class=='Ether lipid'){
            possible_FA2_lipid <- char %>% filter(class== last_FA2_lipid) %>% 
              filter(str_detect(FA_split, str_c(FA1,'_')))
          }
          else{
            possible_FA2_lipid <- char %>% filter(class== last_FA2_lipid) %>% 
              filter(str_detect(FA_split, FA1))
          }
          
          
          if(nrow(possible_FA2_lipid)!=0){
            for(var3 in 1:nrow(possible_FA2_lipid)){
              
              lipid_sub_FA1 <- lipid_sub
              
              possible_FA2 <- possible_FA2_lipid$FA_split[var3]
              
              
              lipid_sub_FA1[str_detect(lipid_sub_FA1, 'FA2')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'FA2')],possible_FA2, sep='_')
              
              if(if_sphingolipid %in% c('Sphingolipid','Deoxysphingolipid')){
                lipid_sub_FA1 <- str_replace(lipid_sub_FA1, '^dhCer_FA2_18:1;2','dhCer_FA2_18:0;2')
                lipid_sub_FA1 <- str_replace(lipid_sub_FA1, 'doxdhCer_FA2_18:1;1','dhCer_FA2_18:0;1')
              }
              
              names(lipid_sub_FA1) <- rep('Lipid',length(lipid_sub_FA1))
              
              lipid_sub_list[[add]] <- lipid_sub_FA1
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add+1
              
            }
            
          }
          else if(nrow(possible_FA2_lipid)==0){
            LMID <- filter(network_node, Abbreviation==last_FA2_lipid)$LIPIDMAPS_subclass_id
            
            LM_mapping <- lipidmaps_sum %>% filter(sub_class_id==LMID) %>% 
              filter(str_detect(FA_split, str_extract(FA1, '\\d+:\\d+')))
            
            if(nrow(LM_mapping)!=0){
              possible_FA2_sum <- LM_mapping$FA_sum %>%  unique()
              possible_FA2_sum <- str_c(possible_FA2_sum, ';0')
            }
            else{
              lipid_sub_FA1 <- lipid_sub
              
              names(lipid_sub_FA1) <- rep('Lipid',length(lipid_sub_FA1))
              
              lipid_sub_list[[add]] <- lipid_sub_FA1
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add+1
              
              next
            }
            
            possible_FA2_lipid <- char %>% filter(FA_split=='') %>% 
              filter(class==last_FA2_lipid, FA_sum %in% possible_FA2_sum)
            
            if(nrow(possible_FA2_lipid)!=0){
              for(var3 in 1:nrow(possible_FA2_lipid)){
                
                lipid_sub_FA1 <- lipid_sub
                
                possible_FA2 <- possible_FA2_lipid$FA_sum[var3]
                
                lipid_sub_FA1[str_detect(lipid_sub_FA1, 'FA2')] <- str_c(lipid_sub_FA1[str_detect(lipid_sub_FA1, 'FA2')],possible_FA2, sep='_')
                
                if(if_sphingolipid %in% c('Sphingolipid','Deoxysphingolipid')){
                  
                  length_db_oh <- str_extract_all(possible_FA2, '\\d+') %>% as.integer()
                  
                  lipid_sub_FA1[str_detect(lipid_sub_FA1,'^dhCer')] <- str_c('dhCer_FA2_',length_db_oh[1], ':',
                                                                             length_db_oh[2]-1, ';2', sep = '_')
                  lipid_sub_FA1[str_detect(lipid_sub_FA1,'doxdhCer')] <- str_c('doxdhCer_FA2_',length_db_oh[1], ':',
                                                                               length_db_oh[2]-1, ';1', sep = '_')
                }
                
                
                names(lipid_sub_FA1) <- rep('Lipid',length(lipid_sub_FA1))
                
                lipid_sub_list[[add]] <- lipid_sub_FA1
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add+1
                
              }
            }
            else{
              lipid_sub_FA1 <- lipid_sub
              
              names(lipid_sub_FA1) <- rep('Lipid',length(lipid_sub_FA1))
              
              lipid_sub_list[[add]] <- lipid_sub_FA1
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add+1
            }
            
          }
        }
        else{
          
          names(lipid_sub) <- rep('Lipid',length(lipid_sub))
          
          lipid_sub_list[[add]] <- lipid_sub
          
          lipid_name[add] <- char$feature[var1]
          
          add <- add+1
        }
      }
      
    }
    else if(char$class[var1]=='TAG'){
      
      if(char[var1, 'FA_split']!=''){
        FA3 <- unlist(char[var1, 'each_FA'])
        #c('16:0;0', '18:0;0','18:0;0')
        
        #combo1 <- str_c(FA3[1], FA3[2],sep = '_')
        #combo2 <- str_c(FA3[2], FA3[3],sep = '_')
        #combo3 <- str_c(FA3[1], FA3[1],sep = '_')
        #all_combo <- c(combo1, combo2,combo3)
        
        possible_FA2_lipid <- char %>% filter(class== 'DAG') %>% 
          filter(sum(!unlist(str_split(FA_split,'_')) %in% FA3)==0) %>% 
          filter(sum(FA3 %in% unlist(str_split(FA_split,'_')))>=2)
        
        
        for(var2 in 1:nrow(lipid_sub_all)){
          lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
          #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub!='']
          
          lipid_sub[lipid_sub=='TAG_FA3'] <- str_c('TAG_FA3',str_c(FA3,collapse = '_'), sep='_')
          
          if(nrow(possible_FA2_lipid)!=0){
            
            for(var3 in 1:nrow(possible_FA2_lipid)){
              
              lipid_sub_FA3 <- lipid_sub
              
              possible_FA2 <- possible_FA2_lipid$FA_split[var3]
              
              lipid_sub_FA3[str_detect(lipid_sub_FA3, 'FA2')] <- str_c(lipid_sub_FA3[str_detect(lipid_sub_FA3, 'FA2')], 
                                                                       possible_FA2, sep='_')
              possible_FA1 <- unlist(str_split(possible_FA2))
              
              possible_FA1_lipid <- char %>% filter(class=='LPA', FA_split %in%possible_FA1)
              
              if(nrow(possible_FA1_lipid)!=0){
                possible_FA1 <- possible_FA1_lipid$FA_split
              }
              for(var4 in 1:length(possible_FA1)){
                lipid_sub_FA3_2 <- lipid_sub_FA3
                
                lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, 'FA1')] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, 'FA1')], 
                                                                             possible_FA1[var4], sep='_')
                names(lipid_sub_FA3_2) <- rep('Lipid',length(lipid_sub_FA3_2))
                
                lipid_sub_list[[add]] <- lipid_sub_FA3_2
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add+1
              }
            }
          }
          
          else if(nrow(possible_FA2_lipid)!=0){
            
            lipid_sub_FA3 <- lipid_sub
            
            names(lipid_sub_FA3) <- rep('Lipid',length(lipid_sub_FA3))
            
            lipid_sub_list[[add]] <- lipid_sub_FA3
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add+1
          }
        }
        
      }
      else if(char[var1, 'FA_split']==''){
        
        FA3 <- char$FA_sum[var1] #c('52:2;0')
        LM_mapping <- lipidmaps_sum %>% filter(sub_class_id==LMID, FA_sum==str_extract(FA3, '\\d+:\\d+'))
        
        if(nrow(LM_mapping)!=0){
          possible_DAG <- LM_mapping$FA_split %>% str_split('_') %>% 
            map(.f = function(x){a <- str_c(x[1],x[2],sep = '_');
            b <- str_c(x[1],x[3],sep = '_'); c <- str_c(x[2],x[3],sep = '_');
            return(c(a,b,c))}) %>% unlist() %>% unique()
          possible_DAG <- str_replace_all(possible_DAG, '_',';0_')
          possible_DAG <- str_c(possible_DAG, ';0')
          
          if(length(FA_pool)!=0){
            
            FA_pool_filter <- str_split(possible_DAG, '_') %>% 
              map(.f = function(x){sum(!x %in%FA_pool)==0}) %>% unlist()
            possible_DAG <- possible_DAG[FA_pool_filter]
          }
          
        }
        if(nrow(LM_mapping)==0 || length(possible_DAG)==0){
          for(var2 in 1:nrow(lipid_sub_all)){
            lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
            #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub!='']
            
            lipid_sub[lipid_sub=='TAG_FA3'] <- str_c('TAG_FA3', FA3, sep='_')
            
            names(lipid_sub) <- rep('Lipid',length(lipid_sub))
            
            lipid_sub_list[[add]] <- lipid_sub
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add+1
            
          }
          next
        }
        
        possible_FA2_lipid <- char %>% filter(class=='DAG', FA_split %in%possible_DAG)
        
        for(var2 in 1:nrow(lipid_sub_all)){
          lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
          #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub!='']
          
          lipid_sub[lipid_sub=='TAG_FA3'] <- str_c('TAG_FA3', FA3, sep='_')
          
          
          if(nrow(possible_FA2_lipid)!=0){
            
            for(var3 in 1:nrow(possible_FA2_lipid)){
              lipid_sub_FA3 <- lipid_sub
              possible_FA2 <- possible_FA2_lipid$FA_split[var3]
              lipid_sub_FA3[str_detect(lipid_sub_FA3, 'FA2')] <- str_c(lipid_sub_FA3[str_detect(lipid_sub_FA3, 'FA2')],possible_FA2, sep='_')
              
              possible_FA1 <- unlist(str_split(possible_FA2,'_'))
              possible_FA1_lipid <- char %>% filter(class=='LPA', FA_split %in%possible_FA1)
              
              if(nrow(possible_FA1_lipid)!=0){
                possible_FA1 <- possible_FA1_lipid$FA_split
              }
              
              for(var4 in 1:length(possible_FA1)){
                lipid_sub_FA3_2 <- lipid_sub_FA3
                
                
                lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, 'FA1')] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, 'FA1')], 
                                                                             possible_FA1[var4], sep='_')
                names(lipid_sub_FA3_2) <- rep('Lipid',length(lipid_sub_FA3_2))
                
                lipid_sub_list[[add]] <- lipid_sub_FA3_2
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add+1
              }              
            }  
          }
          else if(nrow(possible_FA2_lipid)==0){
            
            possible_FA_sum <- possible_DAG %>% str_extract_all('\\d+') %>% 
              map(.f = function(x){a <- as.integer(x[1])+as.integer(x[4]);
              b <- as.integer(x[2])+as.integer(x[5]); c <- str_c(a,b,sep = ':');
              return(c)}) %>% unlist()
            
            possible_FA_sum <- str_c(possible_FA_sum, ';0')
            
            
            possible_FA2_lipid_2 <- char %>% filter(FA_split=='') %>% 
              filter(class=='DAG', FA_sum %in%possible_FA_sum)
            
            
            if(nrow(possible_FA2_lipid_2)!=0){
              for(var3 in 1:nrow(possible_FA2_lipid_2)){
                
                lipid_sub_FA3 <- lipid_sub
                possible_FA2 <- possible_FA2_lipid_2$FA_sum[var3]
                lipid_sub_FA3[str_detect(lipid_sub_FA3, 'FA2')] <- str_c(lipid_sub_FA3[str_detect(lipid_sub_FA3, 'FA2')],possible_FA2, sep='_')
                
                possible_FA1 <- possible_DAG[which(possible_FA_sum==possible_FA2)] %>% 
                  str_split('_') %>% unlist()
                possible_FA1_lipid <- char %>% filter(class=='LPA', FA_split %in%possible_FA1)
                
                if(nrow(possible_FA1_lipid)!=0){
                  possible_FA1 <- possible_FA1_lipid$FA_split
                }
                
                for(var4 in 1:length(possible_FA1)){
                  lipid_sub_FA3_2 <- lipid_sub_FA3
                  
                  lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, 'FA1')] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, 'FA1')], 
                                                                               possible_FA1[var4], sep='_')
                  names(lipid_sub_FA3_2) <- rep('Lipid',length(lipid_sub_FA3_2))
                  
                  lipid_sub_list[[add]] <- lipid_sub_FA3_2
                  
                  lipid_name[add] <- char$feature[var1]
                  
                  add <- add+1
                }  
                
              }  
            }
            else if(nrow(possible_FA2_lipid_2)==0){
              lipid_sub_FA3 <- lipid_sub
              
              names(lipid_sub_FA3) <- rep('Lipid',length(lipid_sub_FA3))
              
              lipid_sub_list[[add]] <- lipid_sub_FA3
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add+1  
            }
          }
        }
      }
    }
    else if(char$class[var1]=='CL'){
      
      if(char[var1, 'FA_split']!=''){
        FA4 <- unlist(char[var1, 'each_FA'])
        
        possible_FA2_lipid <- char %>% filter(class== 'PG') %>% 
          filter(sum(!unlist(str_split(FA_split,'_')) %in% FA4)==0) %>% 
          filter(sum(FA4 %in% unlist(str_split(FA_split,'_')))>=2)
        
        for(var2 in 1:nrow(lipid_sub_all)){
          lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
          #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub!='']
          
          lipid_sub[lipid_sub=='CL_FA4'] <- str_c('CL_FA4',str_c(FA4,collapse = '_'), sep='_')
          
          if(nrow(possible_FA2_lipid)!=0){
            
            for(var3 in 1:nrow(possible_FA2_lipid)){
              
              lipid_sub_FA4 <- lipid_sub
              
              possible_FA2 <- possible_FA2_lipid$FA_split[var3]
              
              lipid_sub_FA4[str_detect(lipid_sub_FA4, 'FA2')] <- str_c(lipid_sub_FA4[str_detect(lipid_sub_FA4, 'FA2')], 
                                                                       possible_FA2, sep='_')
              
              possible_FA1 <- unlist(str_split(possible_FA2, '_'))
              
              possible_FA1_lipid <- char %>% filter(class=='LPA', FA_split %in%possible_FA1)
              
              if(nrow(possible_FA1_lipid)!=0){
                possible_FA1 <- possible_FA1_lipid$FA_split
              }
              for(var4 in 1:length(possible_FA1)){
                lipid_sub_FA4_2 <- lipid_sub_FA4
                
                lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, 'FA1')] <- str_c(lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, 'FA1')], 
                                                                             possible_FA1[var4], sep='_')
                names(lipid_sub_FA4_2) <- rep('Lipid',length(lipid_sub_FA4_2))
                
                lipid_sub_list[[add]] <- lipid_sub_FA4_2
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add+1
              }              
            }
          }
          
          else if(nrow(possible_FA2_lipid)==0){
            
            lipid_sub_FA4 <- lipid_sub
            
            names(lipid_sub_FA4) <- rep('Lipid',length(lipid_sub_FA4))
            
            lipid_sub_list[[add]] <- lipid_sub_FA4
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add+1
          }
        }
        
      }
      else if(char[var1, 'FA_split']==''){
        
        FA4 <- char$FA_sum[var1] #c('52:2;0')
        
        LM_mapping <- lipidmaps_sum %>% filter(sub_class_id==LMID, FA_sum==str_extract(FA4, '\\d+:\\d+'))
        
        if(nrow(LM_mapping)!=0){
          possible_PG <- MLCL_mapping %>% filter(FA_sum==str_extract(FA4, '\\d+:\\d+')) %>% 
            .$possible_PG_split %>% unique()
          
          possible_PG <- str_replace_all(possible_PG, '_',';0_')
          possible_PG <- str_c(possible_PG, ';0')
          
          if(length(FA_pool)!=0){
            
            FA_pool_filter <- str_split(possible_PG, '_') %>% 
              map(.f = function(x){sum(!x %in%FA_pool)==0}) %>% unlist()
            possible_PG <- possible_PG[FA_pool_filter]
          }
          
        }
        if(nrow(LM_mapping)==0 || length(possible_PG)==0){
          for(var2 in 1:nrow(lipid_sub_all)){
            lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
            #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub!='']
            
            lipid_sub[lipid_sub=='CL_FA4'] <- str_c('CL_FA4', FA4, sep='_')
            
            names(lipid_sub) <- rep('Lipid',length(lipid_sub))
            
            lipid_sub_list[[add]] <- lipid_sub
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add+1
            
          }
          next
        }
        
        possible_FA2_lipid <- char %>% filter(class=='PG', FA_split %in%possible_PG)
        
        if(nrow(possible_FA2_lipid)!=0){
          for(var2 in 1:nrow(lipid_sub_all)){
            lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
            #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub!='']
            
            lipid_sub[lipid_sub=='CL_FA4'] <- str_c('CL_FA4', FA4, sep='_')
            for(var3 in 1:nrow(possible_FA2_lipid)){
              lipid_sub_FA4 <- lipid_sub
              
              possible_FA2 <- possible_FA2_lipid$FA_split[var3]
              
              lipid_sub_FA4[str_detect(lipid_sub_FA4, 'FA2')] <- str_c(lipid_sub_FA4[str_detect(lipid_sub_FA4, 'FA2')],possible_FA2, sep='_')
              
              possible_FA1 <- unlist(str_split(possible_FA2, '_'))
              
              possible_FA1_lipid <- char %>% filter(class=='LPA', FA_split %in%possible_FA1)
              
              if(nrow(possible_FA1_lipid)!=0){
                possible_FA1 <- possible_FA1_lipid$FA_split
              }
              for(var4 in 1:length(possible_FA1)){
                lipid_sub_FA4_2 <- lipid_sub_FA4
                
                lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, 'FA1')] <- str_c(lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, 'FA1')], 
                                                                             possible_FA1[var4], sep='_')
                names(lipid_sub_FA4_2) <- rep('Lipid',length(lipid_sub_FA4_2))
                
                lipid_sub_list[[add]] <- lipid_sub_FA4_2
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add+1
              }
            }  
            
          }  
        }
        else if(nrow(possible_FA2_lipid)==0){
          
          possible_PG_sum <- possible_PG %>% str_extract_all('\\d+') %>% 
            map(.f = function(x){a <- as.integer(x[1])+as.integer(x[4]);
            b <- as.integer(x[2])+as.integer(x[5]); c <- str_c(a,b,sep = ':');
            return(c)}) %>% unlist()
          
          possible_PG_sum <- str_c(possible_PG_sum, ';0')
          
          
          possible_FA2_lipid_2 <- char %>% filter(FA_split=='') %>% 
            filter(class=='PG', FA_sum %in%possible_PG_sum)
          
          for(var2 in 1:nrow(lipid_sub_all)){
            lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
            #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub!='']
            
            lipid_sub[lipid_sub=='CL_FA4'] <- str_c('CL_FA4', FA4, sep='_')
            
            if(nrow(possible_FA2_lipid_2)!=0){
              for(var3 in 1:nrow(possible_FA2_lipid_2)){
                
                lipid_sub_FA4 <- lipid_sub
                
                possible_FA2 <- possible_FA2_lipid_2$FA_sum[var3]
                
                lipid_sub_FA4[str_detect(lipid_sub_FA4, 'FA2')] <- str_c(lipid_sub_FA4[str_detect(lipid_sub_FA4, 'FA2')],possible_FA2, sep='_')
                
                possible_FA1 <- possible_PG[which(possible_PG_sum==possible_FA2)] %>% 
                  str_split('_') %>% unlist()
                
                possible_FA1_lipid <- char %>% filter(class=='LPA', FA_split %in%possible_FA1)
                
                if(nrow(possible_FA1_lipid)!=0){
                  possible_FA1 <- possible_FA1_lipid$FA_split
                }
                for(var4 in 1:length(possible_FA1)){
                  lipid_sub_FA4_2 <- lipid_sub_FA4
                  
                  lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, 'FA1')] <- str_c(lipid_sub_FA4_2[str_detect(lipid_sub_FA4_2, 'FA1')], 
                                                                               possible_FA1[var4], sep='_')
                  names(lipid_sub_FA4_2) <- rep('Lipid',length(lipid_sub_FA4_2))
                  
                  lipid_sub_list[[add]] <- lipid_sub_FA4_2
                  
                  lipid_name[add] <- char$feature[var1]
                  
                  add <- add+1
                }
                
              }  
            }
            else if(nrow(possible_FA2_lipid_2)==0){
              lipid_sub_FA4 <- lipid_sub
              
              names(lipid_sub_FA4) <- rep('Lipid',length(lipid_sub_FA4))
              
              lipid_sub_list[[add]] <- lipid_sub_FA4
              
              lipid_name[add] <- char$feature[var1]
              
              add <- add+1  
            }
          }  
          
        }
      }
    }
    else if(char$class[var1]=='MLCL'){
      
      if(char[var1, 'FA_split']!=''){
        FA3 <- unlist(char[var1, 'each_FA'])
        #c('16:0;0', '18:0;0','18:0;0')
        
        
        possible_FA4_lipid <- char %>% filter(class== 'CL') %>% 
          filter(sum(!unlist(str_split(FA_split,'_')) %in% FA3)==0) %>% 
          filter(sum(!FA3 %in% unlist(str_split(FA_split,'_')))==0)
        
        for(var2 in 1:nrow(lipid_sub_all)){
          
          lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
          #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub!='']
          
          lipid_sub[lipid_sub=='MLCL_FA3'] <- str_c('MLCL_FA3',str_c(FA3,collapse = '_'), sep='_')
          
          if(nrow(possible_FA4_lipid)!=0){
            
            for(var3 in 1:nrow(possible_FA4_lipid)){
              
              lipid_sub_FA3 <- lipid_sub
              
              possible_FA4 <- possible_FA4_lipid$FA_split[var3]
              
              lipid_sub_FA3[lipid_sub=='CL_FA4'] <- str_c('CL_FA4',possible_FA4, sep='_')
              
              FA4 <- str_split(possible_FA4, '_') %>% unlist()
              
              possible_FA2_lipid <- char %>% filter(class== 'PG') %>% 
                filter(sum(!unlist(str_split(FA_split,'_')) %in% FA4)==0) %>% 
                filter(sum(FA4 %in% unlist(str_split(FA_split,'_')))>=2)
              
              if(nrow(possible_FA2_lipid)!=0){
                
                for(var4 in 1:nrow(possible_FA2_lipid)){
                  
                  lipid_sub_FA3_2 <- lipid_sub_FA3
                  
                  possible_FA2 <- possible_FA2_lipid$FA_split[var4]
                  
                  lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, 'FA2')] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, 'FA2')], 
                                                                               possible_FA2, sep='_')
                  
                  possible_FA1 <- unlist(str_split(possible_FA2, '_'))
                  
                  possible_FA1_lipid <- char %>% filter(class=='LPA', FA_split %in%possible_FA1)
                  
                  if(nrow(possible_FA1_lipid)!=0){
                    possible_FA1 <- possible_FA1_lipid$FA_split
                  }
                  for(var4 in 1:length(possible_FA1)){
                    lipid_sub_FA3_3 <- lipid_sub_FA3_2
                    
                    
                    lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, 'FA1')] <- str_c(lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, 'FA1')], 
                                                                                 possible_FA1[var4], sep='_')
                    names(lipid_sub_FA3_3) <- rep('Lipid',length(lipid_sub_FA3_3))
                    
                    lipid_sub_list[[add]] <- lipid_sub_FA3_3
                    
                    lipid_name[add] <- char$feature[var1]
                    
                    add <- add+1
                  }
                  
                }
              }
              
              else if(nrow(possible_FA2_lipid)==0){
                
                names(lipid_sub_FA3) <- rep('Lipid',length(lipid_sub_FA3))
                
                lipid_sub_list[[add]] <- lipid_sub_FA3
                
                lipid_name[add] <- char$feature[var1]
                
                add <- add+1
              }
            }
          }
          else if(nrow(possible_FA4_lipid)==0){
            
            
            lipid_sub_FA3 <- lipid_sub
            
            names(lipid_sub_FA3) <- rep('Lipid',length(lipid_sub_FA3))
            
            lipid_sub_list[[add]] <- lipid_sub_FA3
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add+1
          }
        }
        
      }
      else if(char[var1, 'FA_split']==''){
        
        FA3 <- char$FA_sum[var1] #c('52:2;0')
        
        LM_mapping <- MLCL_mapping %>% filter(possible_MLCL_sum==FA3) %>% .[3:5] %>% 
          unique()
        
        if(nrow(LM_mapping)!=0){
          
          possible_CL <- LM_mapping$FA_split
          possible_CL_sum <- LM_mapping$FA_sum
          
          possible_CL <- str_replace_all(possible_CL, '_',';0_')
          possible_CL <- str_c(possible_CL, ';0')
          possible_CL_sum <- str_c(possible_CL_sum, ';0')
          
          if(length(FA_pool)!=0){
            
            FA_pool_filter <- str_split(possible_CL, '_') %>% 
              map(.f = function(x){sum(!x %in%FA_pool)==0}) %>% unlist()
            
            possible_CL <- possible_CL[FA_pool_filter]
            possible_CL_sum <- possible_CL_sum[FA_pool_filter]
            
          }
          
        }
        
        if(nrow(LM_mapping)==0 || length(possible_CL)==0){
          for(var2 in 1:nrow(lipid_sub_all)){
            lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
            #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
            
            lipid_sub <- lipid_sub[lipid_sub!='']
            
            lipid_sub[lipid_sub=='MLCL_FA3'] <- str_c('MLCL_FA3', FA3, sep='_')
            
            names(lipid_sub) <- rep('Lipid',length(lipid_sub))
            
            lipid_sub_list[[add]] <- lipid_sub
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add+1
            
          }
          next
        }
        
        possible_FA4_lipid <- char %>% filter(class=='CL', FA_sum %in%possible_CL_sum)
        
        
        for(var2 in 1:nrow(lipid_sub_all)){
          lipid_sub <- lipid_sub_all[var2,] %>% dplyr::select(-Lipid) %>% unlist()
          #PE: G3P, LPA_FA1, PA_FA2, PS_FA2, PE_FA2
          
          lipid_sub <- lipid_sub[lipid_sub!='']
          
          lipid_sub[lipid_sub=='MLCL_FA3'] <- str_c('MLCL_FA3', FA3, sep='_')
          
          
          if(nrow(possible_FA4_lipid)!=0){
            
            for(var3 in 1:nrow(possible_FA4_lipid)){
              
              lipid_sub_FA3 <- lipid_sub
              FA4 <- possible_FA4_lipid$FA_sum[var3]
              lipid_sub_FA3[lipid_sub_FA3=='CL_FA4'] <- str_c('CL_FA4', FA4, sep='_')
              
              possible_FA4 <- possible_CL[which(possible_CL_sum==FA4)]
              
              FA4s <- str_split(possible_FA4, '_')
              
              
              for (var4 in 1:length(FA4s)){
                FA4 <- FA4s[[var4]]
                
                possible_FA2_lipid <- char %>% filter(class== 'PG') %>% 
                  filter(sum(!unlist(str_split(FA_split,'_')) %in% FA4)==0) %>% 
                  filter(sum(FA4 %in% unlist(str_split(FA_split,'_')))>=2)
                
                
                if(nrow(possible_FA2_lipid)!=0){
                  
                  for(var5 in 1:nrow(possible_FA2_lipid)){
                    
                    lipid_sub_FA3_2 <- lipid_sub_FA3
                    
                    possible_FA2 <- possible_FA2_lipid$FA_split[var5]
                    
                    lipid_sub_FA4[str_detect(lipid_sub_FA3_2, 'FA2')] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, 'FA2')], 
                                                                               possible_FA2, sep='_')
                    
                    possible_FA1 <- unlist(str_split(possible_FA2, '_'))
                    
                    possible_FA1_lipid <- char %>% filter(class=='LPA', FA_split %in%possible_FA1)
                    
                    
                    if(nrow(possible_FA1_lipid)!=0){
                      possible_FA1 <- possible_FA1_lipid$FA_split
                    }
                    for(var6 in 1:length(possible_FA1)){
                      lipid_sub_FA3_3 <- lipid_sub_FA3_2
                      
                      lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, 'FA1')] <- str_c(lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, 'FA1')], 
                                                                                   possible_FA1[var6], sep='_')
                      names(lipid_sub_FA3_3) <- rep('Lipid',length(lipid_sub_FA3_3))
                      
                      lipid_sub_list[[add]] <- lipid_sub_FA3_3
                      
                      lipid_name[add] <- char$feature[var1]
                      
                      add <- add+1
                    }
                    
                  }
                }
                
                else if(nrow(possible_FA2_lipid)==0){
                  FA4 <- str_c(FA4, collapse = '_') %>% 
                    str_replace_all(pattern = ';0','')
                  
                  possible_FA_sum_2 <- MLCL_mapping %>% filter(FA_split==FA4) %>% 
                    .$possible_PG_sum %>% unique()
                  
                  possible_FA2_lipid_2 <- char %>% filter(FA_split=='') %>% 
                    filter(class=='PG', FA_sum %in%possible_FA_sum_2)
                  
                  
                  if(nrow(possible_FA2_lipid_2)!=0){
                    for(var5 in 1:nrow(possible_FA2_lipid_2)){
                      
                      lipid_sub_FA3_2 <- lipid_sub_FA3
                      
                      possible_FA2 <- possible_FA2_lipid_2$FA_sum[var5]
                      lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, 'FA2')] <- str_c(lipid_sub_FA3_2[str_detect(lipid_sub_FA3_2, 'FA2')],possible_FA2, sep='_')
                      
                      
                      LM_mapping_2 <- lipidmaps_sum %>% filter(sub_class_id=='[GP0401]', FA_sum==str_extract(possible_FA2, '\\d+:\\d+'))
                      
                      
                      if(nrow(LM_mapping_2)!=0){
                        possible_FA1 <- LM_mapping_2$FA_split %>% str_split('_') %>% unlist() %>% unique()
                        
                        possible_FA1 <- str_c(possible_FA1, ';0')
                        
                        if(length(FA_pool)!=0){
                          possible_FA1 <- possible_FA1[possible_FA1%in%FA_pool]
                        }
                        
                      }
                      else{
                        names(lipid_sub_FA3_2) <- rep('Lipid',length(lipid_sub_FA3_2))
                        
                        lipid_sub_list[[add]] <- lipid_sub_FA3_2
                        
                        lipid_name[add] <- char$feature[var1]
                        
                        add <- add+1
                        next
                      }
                      
                      
                      
                      possible_FA1_lipid <- char %>% filter(class=='LPA', FA_split %in%possible_FA1)
                      
                      if(nrow(possible_FA1_lipid)!=0){
                        possible_FA1 <- possible_FA1_lipid$FA_split
                      }
                      for(var6 in 1:length(possible_FA1)){
                        lipid_sub_FA3_3 <- lipid_sub_FA3_2
                        
                        lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, 'FA1')] <- str_c(lipid_sub_FA3_3[str_detect(lipid_sub_FA3_3, 'FA1')], 
                                                                                     possible_FA1[var6], sep='_')
                        names(lipid_sub_FA3_3) <- rep('Lipid',length(lipid_sub_FA3_3))
                        
                        lipid_sub_list[[add]] <- lipid_sub_FA3_3
                        
                        lipid_name[add] <- char$feature[var1]
                        
                        add <- add+1
                      }
                      
                    }  
                  }
                  else if(nrow(possible_FA2_lipid)==0){
                    names(lipid_sub_FA3) <- rep('Lipid',length(lipid_sub_FA3))
                    
                    lipid_sub_list[[add]] <- lipid_sub_FA3
                    
                    lipid_name[add] <- char$feature[var1]
                    
                    add <- add+1
                  }
                }
                
              }
            }  
          }
          else if(nrow(possible_FA4_lipid)==0){
            
            names(lipid_sub) <- rep('Lipid',length(lipid_sub))
            
            lipid_sub_list[[add]] <- lipid_sub
            
            lipid_name[add] <- char$feature[var1]
            
            add <- add+1  
            
            
            
          }
        }
        
      }
      
    }
    
  }
  sub_species <- list(lipid_name, lipid_sub_list)
  species_substructure <- sub_species[[2]] %>% map(.f = function(x){names(x) <- 1:length(x);return(x)}) %>% 
    plyr::ldply(rbind) %>% mutate(Lipid=sub_species[[1]]) %>% 
    dplyr::select(Lipid, everything())
  
  species_substructure[is.na(species_substructure)] <- ''
  
  colnames(species_substructure) <- c('Lipid', str_c('Unit', 1:(ncol(species_substructure)-1)))
  
  
  return(species_substructure)
}



species_sub_extract <- function(lipid_substructure, unprocessed_data_result,
                                type='species', pct_limit=0.3,
                                exo_lipid=NULL){
  
  if(class(unprocessed_data_result)=='list'){
    unprocessed_data_result <- unprocessed_data_result[[2]] %>% filter(type==type)
  }
  
  non_pro <- unprocessed_data_result
  #non_pro <- non_processed_data_result %>% filter(lipid %in% lipid_substructure$Lipid)
  lipid_substructure <- lipid_substructure %>% filter(Lipid %in% non_pro$lipid)
  
  lipid_test <- lipid_substructure %>% mutate(num=1:nrow(.)) %>% 
    gather(-Lipid,-num, key='Unit', value='sub') %>% 
    mutate(Unit=factor(Unit, levels = unique(Unit))) %>% 
    mutate(sub=str_replace(sub, '_FA\\d','')) %>% 
    left_join(non_pro[c('lipid','log2FC')], by=c('sub'='lipid')) %>% .[-4] %>% 
    spread(key='Unit',value = 'log2FC') %>% arrange(num) %>% 
    dplyr::select(-num)
  
  lipid_sub_trans <- list(nrow(lipid_test))
  lipid_name <- character(nrow(lipid_test))
  for(num in 1:nrow(lipid_test)){
    #print(num)
    lipid_sub <- lipid_substructure[-1][num,] %>% unlist() %>% 
      str_replace('_FA\\d','')
    lipid_sub_FC <- lipid_test[-1][num,] %>% unlist()
    
    lipid_sub_FC <- lipid_sub_FC[lipid_sub!='']
    lipid_sub <- lipid_sub[lipid_sub!='']
    
    lipid_name[num] <- last(lipid_sub)
    
    Not_NA_pct <- sum(!is.na(lipid_sub_FC))/length(lipid_sub_FC)
    if(Not_NA_pct<pct_limit){
      lipid_sub_trans[[num]] <- ''
    }
    else{
      if(last(lipid_sub_FC)>0){
        stop_point <- which(lipid_sub_FC<0)
        if(length(stop_point)==0){
          lipid_sub_trans[[num]] <- lipid_sub
        }
        else if(max(stop_point)==(length(lipid_sub)-1)){
          lipid_sub_trans[[num]] <- ''
        }
        else{
          stop_loc <- max(which(lipid_sub_FC<0))+1
          lipid_sub_trans[[num]] <- lipid_sub[stop_loc:length(lipid_sub)]
        }
      }
      else{
        stop_point <- which(lipid_sub_FC>0)
        if(length(stop_point)==0){
          lipid_sub_trans[[num]] <- lipid_sub
        }
        else if(max(stop_point)==(length(lipid_sub)-1)){
          lipid_sub_trans[[num]] <- ''
        }
        else{
          stop_loc <- max(which(lipid_sub_FC>0))+1
          lipid_sub_trans[[num]] <- lipid_sub[stop_loc:length(lipid_sub)]
        }
      }
    }
  }
  names(lipid_sub_trans) <- lipid_name
  
  lipid_sub_trans <- plyr::ldply(lipid_sub_trans, rbind)
  lipid_sub_trans[is.na(lipid_sub_trans)] <- ''
  
  
  colnames(lipid_sub_trans) <- c('Lipid', str_c('Unit', 1:(ncol(lipid_sub_trans)-1)))
  lipid_sub_trans <- unique(lipid_sub_trans) %>% filter(Unit1!='')
  
  
  if(type=='class'){
    exist_lipid <- non_pro %>% filter(type=='class') %>% .$lipid
    lost_lipid <- exist_lipid[!exist_lipid%in%lipid_sub_trans$Lipid]
    lipid_sub_trans <- lipid_sub_trans %>% add_row(Lipid = lost_lipid, 
                                                   Unit1 =lost_lipid)
    lipid_sub_trans[is.na(lipid_sub_trans)] <- ''
  }
  else if(type=='species'){
    exist_lipid <- non_pro %>% filter(type=='species') %>% .$lipid
    lost_lipid <- exist_lipid[!exist_lipid%in%lipid_sub_trans$Lipid]
    lipid_sub_trans <- lipid_sub_trans %>% add_row(Lipid = lost_lipid, 
                                                   Unit1 =lost_lipid)
    lipid_sub_trans[is.na(lipid_sub_trans)] <- ''
  }
  
  species_sub_stop <- list(lipid_sub_trans[[1]], apply(lipid_sub_trans[-1], MARGIN = 1, FUN = function(x){x[x!='']}))
  if(!is.null(exo_lipid)){
    exo_lipid_with_neighbor <- c(which(map_chr(species_sub_stop[[2]], ~last(.)) %in% exo_lipid), 
                                 which(map_chr(species_sub_stop[[2]], ~last(., 2)[1]) %in% exo_lipid)) %>% unique()
    
    species_sub_stop[[2]][exo_lipid_with_neighbor] <- 
      species_sub_stop[[2]][exo_lipid_with_neighbor] %>% map(.f = function(x){y <- last(x); names(y) <- 'Unit1';y})
  }
  return(species_sub_stop)
}


add_rev_rection <- function(network_edge, species_net){
  rev_reaction <- network_edge %>% filter(reverse==1) %>% 
    apply(MARGIN = 1, FUN = function(x){c(x[1],x[2])})
  
  species_net_pair <- map2(species_net$S1, species_net$P1,
                           .f = function(x,y){c(str_split(x,'_')[[1]][1], str_split(y, '_')[[1]][1])})
  
  add_reaction <- species_net_pair %>% map_lgl(.f = function(x){max(colSums(matrix(rev_reaction %in% x, nrow=2)))==2})
  
  species_net_w_rev <- data.frame(S1=species_net$P1[add_reaction], P1=species_net$S1[add_reaction]) %>% 
    rbind(species_net) %>% unique()
  
  return(species_net_w_rev)
}

build_species_net <- function(species_substructure){
  species_network <- split(species_substructure[-1], seq(nrow(species_substructure)))
  species_network <- species_network %>% map(.f = function(x){x[x!='']})
  
  species_network <- species_network %>% map(.f = function(x){head(rep(x, each=2)[-1],-1)})
  species_network <-  matrix(unlist(species_network), ncol=2, byrow = T,dimnames = list(NULL, c('S1','P1'))) %>% 
    as.data.frame() %>%  unique()
  
  species_network <- species_network %>% mutate(S1=str_replace(S1, '_FA\\d',''),
                                                P1=str_replace(P1, '_FA\\d',''))
  return(species_network)
  
}

