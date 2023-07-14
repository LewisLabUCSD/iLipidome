library(ggsci)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(ggtext)
file <- dirname(rstudioapi::getSourceEditorContext()$path)

#-------------------Data upload-------------------

source(file.path(dirname(file),'Required_function/required_function.R'))
load(file.path(dirname(file),'Required_data/required_data.RData'))


exp_raw <- read.csv(file.path(file, 'lipidome_data/exp_DR.csv'))

char_raw <- read.csv(file.path(file, 'lipidome_data/char_DR.csv'))



data_process <- function(exp_data, exclude_var_missing=T, 
                         missing_pct_limit=50,
                         replace_zero=T, zero2what='min', xmin=0.5,
                         replace_NA=T, NA2what='min', ymin=0.5,
                         pct_transform=T,
                         data_transform=T,trans_type='log',
                         centering=F,
                         scaling=F){
  require(tidyverse)
  require(stringr)
  exp_data2 <- exp_data[-1]
  #replace 0 with: min, specfic num
  if(replace_zero==T){
    if(is.numeric(zero2what)){
      exp_data2[exp_data2==0] <- zero2what
    }else if(zero2what=='min'){
      for(a in 1:nrow(exp_data2)){
        num <- unlist(exp_data2[a,])[unlist(exp_data2[a,])!=0]
        exp_data2[a,][exp_data2[a,]==0] <- xmin*min(num, na.rm = T)
      }
    }else if(zero2what=='NA'){
      exp_data2[exp_data2==0] <- NA
    }
    exp_data <- cbind(exp_data[1],exp_data2)
  }
  
  
  
  exp_data2 <- exp_data[-1]
  
  #exclude_var_missing
  if(exclude_var_missing==T){
    missing_pct <- apply(exp_data2, 1, function(x){ sum(is.na(x)) / length(x) })
    maintain_var <- missing_pct*100 < missing_pct_limit
    exp_data <- exp_data[maintain_var,]
  }
  if(nrow(exp_data)==0){
    #stop('no species remains')
    return(NULL)
  }
  
  
  exp_data2 <- exp_data[-1]
  
  
  #replace NA with: min, mean, median, specfic num
  if(replace_NA==T){
    if(NA2what=='min'){
      for(a in 1:nrow(exp_data2)){
        exp_data2[a,][is.na(exp_data2[a,])] <- ymin*min(unlist(exp_data2[a,]), na.rm = T)
      }
    }else if(NA2what=='mean'){
      for(a in 1:nrow(exp_data2)){
        exp_data2[a,][is.na(exp_data2[a,])] <- mean(unlist(exp_data2[a,]), na.rm = T)
      }
    }else if(NA2what=='median'){
      for(a in 1:nrow(exp_data2)){
        exp_data2[a,][is.na(exp_data2[a,])] <- median(unlist(exp_data2[a,]), na.rm = T)
      }
    }else if(is.numeric(NA2what)){
      exp_data2[is.na(exp_data2)] <- NA2what
    }
    exp_data <- cbind(exp_data[1],exp_data2)
  }
  
  if(pct_transform==T){
    exp_data[-1] <- map2(exp_data[-1], colSums(exp_data[-1], na.rm = T), ~.x/.y*100) %>% 
      as.data.frame()
  }
  
  #data_transform
  if(data_transform==T){
    if(trans_type=='log'){
      exp_data[-1] <- log10(exp_data[-1])
    }
  }
  
  #centering
  if(centering==T){
    exp_data[-1] <- exp_data[-1] %>% t() %>% scale(scale = F) %>% t() %>% as.data.frame()
  }
  
  #scaling
  if(scaling==T){
    exp_data[-1] <- exp_data[-1] %>% t() %>% scale() %>% t() %>% as.data.frame()
  }
  
  return(exp_data)
}



#-----------------------------lipid class analysis--------------------------------

exp_pro_class <- data_process(exp_raw[-c(1:43),], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                              replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp_pro_class <- remove_rownames(exp_pro_class)

char_pro_class <- char_raw %>% filter(feature %in% exp_pro_class$feature)


each_FA <- str_extract_all(char_pro_class$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_pro_class <- char_pro_class %>% mutate(each_FA=each_FA)

#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_pro_class, char_pro_class, 't.test', 'adj_p_value', 1:4,  5:8)



#Lipid class substructure analysis

#Extract class substructures using fold change


class_t <- no_sub_t[[2]] %>% filter(type=='class')

class_sub_stop <- species_sub_extract(filter(lipid_substructure, Lipid %in% class_t$lipid),
                                      class_t, 'class', pct_limit = 0.01)


#Transform lipid exp into substructure exp

Species2Char <- function(exp_data, lipid_char_table, char_var){
  require(tidyverse)
  require(stringr)
  
  exp_data[-1][is.na(exp_data[-1])] <- 0
  
  if(str_detect(char_var, 'FA_')){
    max_chain_num <- str_split(lipid_char_table[[char_var]],  pattern = ',') %>% 
      map_dbl(length) %>% max()
    
    suppressWarnings(
      transform_table <- exp_data %>% left_join(lipid_char_table[c('feature', char_var)], by='feature') %>% 
        separate(eval(parse(text = char_var)), letters[1:max_chain_num]) %>% 
        gather(letters[1:max_chain_num], key='key', value = 'value') %>% 
        filter(!is.na(value))
    )
    if(nrow(transform_table)==0){
      transform_table <- data.frame()
    }else{
      transform_table <- transform_table %>% 
        dplyr::select(-feature,-key) %>% 
        aggregate(. ~ value, ., sum)
      transform_table[[1]] <- as.character(transform_table[[1]])
      colnames(transform_table)[1] <- char_var
    }
    
  }else{
    transform_table <- exp_data %>% left_join(lipid_char_table[c('feature', char_var)], by='feature') %>% 
      dplyr::select(-1)
    transform_table <- transform_table[!is.na(transform_table[[char_var]]),]
    if(nrow(transform_table)==0){
      transform_table <- data.frame()
    }else{
      transform_table <- transform_table %>% 
        aggregate(as.formula(str_c('. ~ ',char_var)), ., sum)
    }
  }
  
  #if(pct_transform==T  && nrow(transform_table)!=0){}
  # transform_table[-1] <- map2(transform_table[-1], colSums(transform_table[-1], na.rm = T), ~.x/.y*100) %>% 
  #   as.data.frame()
  
  transform_table[-1][transform_table[-1]==0] <- NA
  
  return(transform_table)
}


char_exp <- Species2Char(exp_pro_class, char_pro_class, 'class')

char_exp[-1][is.na(char_exp[-1])] <- 0

class_sub_exp <- lipid_sub_matrix(char_exp, class_sub_stop, 'class')


#Differential expression analysis for substructures

class_sub_exp_t <- t_test(class_sub_exp[[3]], 1:4, 5:8, 't.test', 'adj_p_value')

#Class biosynthetic network data transformation

class_network <- network_edge[c('S1','P1')] %>% filter(S1 %in% class_sub_exp_t$lipid,
                                                       P1 %in% class_sub_exp_t$lipid)


#Essential pathway analysis for class substructures

set.seed(1)

path_score_class <-  path_scoring(class_network, class_sub_exp_t,
                                  calibrate = T, data_type = 'Class')



#Essential edges (reactions) analysis for class substructures



reaction_score_class <- reaction_scoring(network_edge,
                                         class_sub_exp[[3]],
                                         class_sub_exp_t,
                                         ctrl=1:4, exp=5:8,
                                         Species = 'mouse')


#Lipid class biosynthetic network construction

class_network_data <- draw_network(network_data = class_network,
                                   DE_data = class_sub_exp_t,
                                   if_species = F,significant = 'adj_p_value',
                                   path_scoring_result = path_score_class,
                                   reaction_scoring_result = reaction_score_class,
                                   top_n = 3, path_type = 'both')

visNetwork(class_network_data[[1]], mutate(class_network_data[[2]], label=''))


SF6a_data <- class_network_data[[2]]%>% left_join(class_sub_exp_t, by=c('from'='lipid')) %>% 
  left_join(class_sub_exp_t, by=c('to'='lipid')) %>% 
  .[c(1:9,16,21,22,29,34,35)]



#write.xlsx(SF6a_data, file.path(file, 'source_data/SF6a.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)

#-----------------------------FA analysis--------------------------------


exp_pro_FA <- data_process(exp_raw[1:43,], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                           replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = F, data_transform = F)

exp_pro_FA <- remove_rownames(exp_pro_FA)

char_pro_FA <- char_raw %>% filter(feature %in% exp_pro_FA$feature)


each_FA <- str_extract_all(char_pro_FA$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_pro_FA <- char_pro_FA %>% mutate(each_FA=each_FA)


#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_pro_FA, char_pro_FA, 't.test', 'adj_p_value', 1:4,  5:8)

#FA substructure analysis
#FA biosynthetic network data transformation

FA_network_DR <- build_FA_net(FA_network, no_sub_t)

#Decompose lipids into FA substructures


FA_substructure <- FA_sub_transform(FA_network_DR, no_sub_t)


#Extract FA substructure using fold changes


FA_t <- no_sub_t[[2]] %>% filter(type=='FA')

FA_sub_stop <- FA_sub_extract(char_pro_FA, FA_substructure, FA_t,
                              exact_FA='no', exo_lipid=NULL)

#Transform lipid FA exp data into substructure exp data
FA_sub_exp <- lipid_sub_matrix(exp_pro_FA, FA_sub_stop, 'FA')


#Differential expression analysis for FA substructures
FA_sub_exp_t <- t_test(FA_sub_exp[[3]], 1:4, 5:8, 't.test', 'adj_value')


#Essential pathway analysis for FA substructures

set.seed(1)
path_data_DR <- path_scoring(FA_network_DR, FA_sub_exp_t,
                             calibrate = T, data_type ='FA')


#Essential edges (reactions) analysis for FA substructures

reaction_data_DR <- reaction_scoring(network = FA_network_DR, 
                                     sub_exp = FA_sub_exp[[3]],
                                     sub_t = FA_sub_exp_t, 
                                     ctrl = 1:4, exp = 5:8, 
                                     Species = 'mouse')

#FA biosynthetic network construction
FA_network_data <- draw_network(network_data = FA_network_DR,
                                DE_data = FA_sub_exp_t,
                                if_species = F, significant = 'adj_p_value',
                                path_scoring_result = path_data_DR,
                                reaction_scoring_result = reaction_data_DR,
                                top_n = 3, path_type = 'both')




visNetwork(FA_network_data[[1]],mutate(FA_network_data[[2]], label='')) %>% 
  visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                  physics = F, smooth = TRUE, randomSeed =3) 

SF7a_data <- FA_network_data[[2]]%>% left_join(FA_sub_exp_t, by=c('from'='lipid')) %>% 
  left_join(FA_sub_exp_t, by=c('to'='lipid')) %>% 
  .[c(1:9,16,21,22,29,34,35)]



#write.xlsx(SF7a_data, file.path(file, 'source_data/SF7a.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)


#-----------------------------FA residue analysis--------------------------------


exp_pro_FAresidue <- data_process(exp_raw[-c(1:43),], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                                  replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = F, data_transform = F)

exp_pro_FAresidue <- remove_rownames(exp_pro_FAresidue)

char_pro_FAresidue <- char_raw %>% filter(feature %in% exp_pro_FAresidue$feature)


each_FA <- str_extract_all(char_pro_FAresidue$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_pro_FAresidue <- char_pro_FAresidue %>% mutate(each_FA=each_FA)

#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_pro_FAresidue, char_pro_FAresidue, 't.test', 'adj_p_value', 1:4,  5:8)

#FA substructure analysis
#FA biosynthetic network data transformation

FA_network_DR <- build_FA_net(FA_network, no_sub_t)

#Decompose lipids into FA substructures


FA_substructure <- FA_sub_transform(FA_network_DR, no_sub_t)


#Extract FA substructure using fold changes


FA_t <- no_sub_t[[2]] %>% filter(type=='FA')

FA_sub_stop <- FA_sub_extract(char_pro_FAresidue, FA_substructure, FA_t,
                              exact_FA='no', exo_lipid=NULL)

#Transform lipid FA exp data into substructure exp data
FA_sub_exp <- lipid_sub_matrix(exp_pro_FAresidue, FA_sub_stop, 'FA')


#Differential expression analysis for FA substructures
FA_sub_exp_t <- t_test(FA_sub_exp[[3]], 1:4, 5:8, 't.test', 'adj_value')


#Essential pathway analysis for FA substructures

set.seed(1)
path_data_DR <- path_scoring(FA_network_DR, FA_sub_exp_t,
                             calibrate = T, data_type ='FA')



#Essential edges (reactions) analysis for FA substructures

reaction_data_DR <- reaction_scoring(network = FA_network_DR, 
                                     sub_exp = FA_sub_exp[[3]],
                                     sub_t = FA_sub_exp_t, 
                                     ctrl = 1:4, exp = 5:8, 
                                     Species = 'mouse')

#FA biosynthetic network construction
FA_network_data <- draw_network(network_data = FA_network_DR,
                                DE_data = FA_sub_exp_t,
                                if_species = F, significant = 'adj_p_value',
                                path_scoring_result = path_data_DR,
                                reaction_scoring_result = reaction_data_DR,
                                top_n = 3, path_type = 'both')




visNetwork(FA_network_data[[1]], mutate(FA_network_data[[2]], label='')) %>% 
  visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                  physics = F, smooth = TRUE, randomSeed =3) 

SF7b_data <- FA_network_data[[2]]%>% left_join(FA_sub_exp_t, by=c('from'='lipid')) %>% 
  left_join(FA_sub_exp_t, by=c('to'='lipid')) %>% 
  .[c(1:9,16,21,22,29,34,35)]



#write.xlsx(SF7b_data, file.path(file, 'source_data/SF7b.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)


#-----------------------------Correlation between FA and FA residue analysis-----------------------------------------



FA_data <- unprocessed_data_test(exp_pro_FA, char_pro_FA, 't.test', 'adj_p_value', 1:4,  5:8)

FA_residue_data <- unprocessed_data_test(exp_pro_FAresidue, char_pro_FAresidue, 't.test', 'adj_p_value', 1:4,  5:8)


options(scipen = 4)
SF7d_data <- filter(FA_data[[1]], type=='FA')[-2] %>% 
  gather(-feature, key='key',value='value') %>% 
  left_join(filter(FA_residue_data[[1]], type=='FA')[-2] %>% 
              gather(-feature, key='key',value='value'),
            by=c('feature','key')) %>% 
  `colnames<-`(c('feature','sample_name','Fatty acyls','FA residues'))

SF7d_data %>% ggplot(aes(x=log10(`Fatty acyls`), y=log10(`FA residues`)))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  stat_cor(method = "pearson")+
  labs(x='log10 (Fatty acyls expression)', y='log10 (FA residues expression)')+
  theme_bw()



#write.xlsx(SF7d_data, file.path(file, 'source_data/SF7d.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)


#-------------------Save data-------------------
#save.image(file.path(file, 'Dietary_restriction.RData'))
