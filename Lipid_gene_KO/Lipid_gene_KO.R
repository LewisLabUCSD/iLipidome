library(ggsci)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(ggtext)
file <- dirname(rstudioapi::getSourceEditorContext()$path)

#-------------------Data upload-------------------

source(file.path(dirname(file),'Required_function/required_function.R'))
load(file.path(dirname(file),'Required_data/required_data.RData'))


raw_data1 <- read.csv(file.path(file, 'lipidome_data/raw_data_deletion_t.csv'))
raw_data2 <- read.csv(file.path(file, 'lipidome_data/raw_data_trap_t.csv'))

raw_data1$feature[str_detect(raw_data1$feature, '(Cer)|(SM)') & str_detect(raw_data1$feature, ':0')] <- 
  str_c('dh',raw_data1$feature[str_detect(raw_data1$feature, '(Cer)|(SM)') & str_detect(raw_data1$feature, ':0')])

raw_data2$feature[str_detect(raw_data2$feature, '(Cer)|(SM)') & str_detect(raw_data2$feature, ':0')] <- 
  str_c('dh',raw_data2$feature[str_detect(raw_data2$feature, '(Cer)|(SM)') & str_detect(raw_data2$feature, ':0')])


exp_raw_1 <- build_char_table(raw_data1, network_node)[[1]]

exp_raw_1$feature <- str_replace(exp_raw_1$feature, ' O-', '>O- ') %>% 
  str_replace(' ','_') %>% str_replace('>O',' O') %>% str_replace('\\/','_')

char_raw_1 <- build_char_table(raw_data1, network_node)[[2]]
char_raw_1$feature <- exp_raw_1$feature

exp_raw_2 <- build_char_table(raw_data2, network_node)[[1]]
exp_raw_2$feature <- str_replace(exp_raw_2$feature, ' O-', '>O- ') %>% 
  str_replace(' ','_') %>% str_replace('>O',' O') %>% str_replace('\\/','_')

char_raw_2 <- build_char_table(raw_data2, network_node)[[2]]
char_raw_2$feature <- exp_raw_2$feature

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



#-------------------CERS2-------------------

exp_CERS2 <- data_process(exp_raw_1[1:10], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                          replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp_CERS2 <- remove_rownames(exp_CERS2)

char_CERS2 <- char_raw_1 %>% filter(feature %in% exp_CERS2$feature)


each_FA <- str_extract_all(char_CERS2$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_CERS2 <- char_CERS2 %>% mutate(each_FA=each_FA)

#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_CERS2, char_CERS2, 't.test', 'adj_p_value', 1:3,  4:9)

#Lipid species substructure analysis
#Decompose lipids into species substructures

species_substructure <- species_sub_transform(char_CERS2, lipid_substructure, 
                                              network_node)

#Extract species substructures using fold change

species_t <- no_sub_t[[2]] %>% filter(type=='species')


species_sub_stop <- species_sub_extract(species_substructure, species_t, 'species', pct_limit = 0.01)

#Transform lipid exp into substructure exp

species_sub_exp <- lipid_sub_matrix(exp_CERS2, species_sub_stop, 'species')

#Differential expression analysis for substructures

species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:3, 4:9, 
                            't.test', 'adj_p_value')

#Species biosynthetic network data transformation


species_network <- build_species_net(species_substructure)

SL_network <- species_substructure[str_detect(species_substructure$Lipid, '(Cer)|(SM)'),] %>% 
  unlist() %>% unique() %>% str_replace('FA\\d+_','')

SL_network <- species_network %>% filter(S1 %in% SL_network, P1 %in% SL_network)

#Essential pathway analysis for species substructures

set.seed(1)
path_score_species <-  path_scoring(SL_network, species_sub_exp_t,
                                    calibrate = T, data_type = 'species')


path_data_CERS2 <- rbind(path_score_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_species %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(path_data_CERS2, file.path(file, 'source_data/path_data_CERS2.csv'), row.names = F)

path_data_CERS2 %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  #scale_y_continuous(limits = c(-7,7))+
  scale_fill_manual(values = bluered(100)[c(1,10,20,30,40,100,90,80,70,60)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')


#Essential edges (reactions) analysis for species substructures

species_net_w_rev <- add_rev_rection(network_edge, SL_network)

reaction_score_species <- reaction_scoring(species_net_w_rev,
                                           species_sub_exp[[3]],
                                           species_sub_exp_t,
                                           ctrl=1:3, exp=4:9,
                                           Species = 'human')


top5_rep_path <- c(path_score_species %>% filter(Type=='Active') %>% 
                     arrange(desc(cal_score)) %>% 
                     .[!duplicated(.$rep_sub_path),] %>% .[1:3,] %>% .$rep_sub_path
                   ,path_score_species %>% filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path)


top5_rep_path <- path_score_species %>% filter(rep_sub_path%in%top5_rep_path, Significant=='yes') %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top5_path <- reaction_score_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top5_rep_path && x[2] %in% top5_rep_path})
edge_in_top5_path <- reaction_score_species[edge_in_top5_path,]$edge_name



reaction_data_CERS2 <- reaction_score_species%>% 
  filter(edge_name %in%edge_in_top5_path) %>% 
  filter(p_value<0.05) %>% 
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color, FA_change))

#write.csv(reaction_data_CERS2[-c(20:22)], file.path(file, 'source_data/reaction_data_CERS2.csv'), row.names = F)

reaction_data_CERS2 %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase','Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(reaction_data_CERS2$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', legend.box = "vertical",
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = rev(pal_lancet()(2)))+
  scale_color_manual(values = c('gold','white'))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 significant edges',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))


#-------------------FADS3-------------------

exp_FADS3 <- data_process(exp_raw_1[c(1:4,11:13)], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                          replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp_FADS3 <- remove_rownames(exp_FADS3)

char_FADS3 <- char_raw_1 %>% filter(feature %in% exp_FADS3$feature)


each_FA <- str_extract_all(char_FADS3$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_FADS3 <- char_FADS3 %>% mutate(each_FA=each_FA)

#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_FADS3, char_FADS3, 't.test', 'adj_p_value', 1:3,  4:6)


#Lipid species substructure analysis
#Decompose lipids into species substructures

species_substructure <- species_sub_transform(char_FADS3, lipid_substructure, 
                                              network_node)

#Extract species substructures using fold change

species_t <- no_sub_t[[2]] %>% filter(type=='species')

species_sub_stop <- species_sub_extract(species_substructure, species_t, 'species', pct_limit = 0.01)


#Transform lipid exp into substructure exp

species_sub_exp <- lipid_sub_matrix(exp_FADS3, species_sub_stop, 'species')

#Differential expression analysis for substructures
species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:3, 4:6, 't.test', 'adj_p_value')

#Species biosynthetic network data transformation
species_network <- build_species_net(species_substructure)
SL_network <- species_substructure[str_detect(species_substructure$Lipid, '(Cer)|(SM)'),] %>% 
  unlist() %>% unique() %>% str_replace('FA\\d+_','')

SL_network <- species_network %>% filter(S1 %in% SL_network, P1 %in% SL_network)

#Essential pathway analysis for species substructures

set.seed(1)
path_score_species <-  path_scoring(SL_network, species_sub_exp_t,
                                    calibrate = T, data_type = 'species')


path_data_FADS3 <- rbind(path_score_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_species %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(path_data_FADS3, file.path(file, 'source_data/path_data_FADS3.csv'), row.names = F)

path_data_FADS3 %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  #scale_y_continuous(limits = c(-7,7))+
  scale_fill_manual(values = bluered(100)[c(1,10,20,30,40,100,90,80,70,60)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')


#-------------------SGMS1-------------------

exp_SGMS1 <- data_process(exp_raw_2[c(1:4,29:31)], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                          replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp_SGMS1 <- remove_rownames(exp_SGMS1)

char_SGMS1 <- char_raw_2 %>% filter(feature %in% exp_SGMS1$feature)


each_FA <- str_extract_all(char_SGMS1$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_SGMS1 <- char_SGMS1 %>% mutate(each_FA=each_FA)



#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_SGMS1, char_SGMS1, 't.test', 'adj_p_value', 1:3,  4:6)

#Lipid species substructure analysis
#Decompose lipids into species substructures

species_substructure <- species_sub_transform(char_SGMS1, lipid_substructure, 
                                              network_node)

#Extract species substructures using fold change


species_t <- no_sub_t[[2]] %>% filter(type=='species')

species_sub_stop <- species_sub_extract(species_substructure, species_t, 'species', pct_limit = 0.01)


#Transform lipid exp into substructure exp
species_sub_exp <- lipid_sub_matrix(exp_SGMS1, species_sub_stop, 'species')


#Differential expression analysis for substructures
species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:3, 4:6, 't.test', 'adj_p_value')

#Species biosynthetic network data transformation
species_network <- build_species_net(species_substructure)
SL_network <- species_substructure[str_detect(species_substructure$Lipid, '(Cer)|(SM)'),] %>% 
  unlist() %>% unique() %>% str_replace('FA\\d+_','')

SL_network <- species_network %>% filter(S1 %in% SL_network, P1 %in% SL_network)

#Essential pathway analysis for species substructures


set.seed(1)
path_score_species <-  path_scoring(SL_network, species_sub_exp_t,
                                    calibrate = T, data_type = 'species')


path_data_SGMS1 <- rbind(path_score_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_species %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(path_data_SGMS1, file.path(file, 'source_data/path_data_SGMS1.csv'), row.names = F)

path_data_SGMS1 %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  #scale_y_continuous(limits = c(-7,7))+
  scale_fill_manual(values = bluered(100)[c(1,10,20,30,40,100,90,80,70,60)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')


#Essential edges (reactions) analysis for species substructures


species_net_w_rev <- add_rev_rection(network_edge, SL_network)


reaction_score_species <- reaction_scoring(species_net_w_rev,
                                           species_sub_exp[[3]],
                                           species_sub_exp_t,
                                           ctrl=1:3, exp=4:6,
                                           Species = 'human')


top5_rep_path <- c(path_score_species %>% filter(Type=='Active') %>% 
                     arrange(desc(cal_score)) %>% 
                     .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path
                   ,path_score_species %>% filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path)


top5_rep_path <- path_score_species %>% filter(rep_sub_path%in%top5_rep_path, Significant=='yes') %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top5_path <- reaction_score_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top5_rep_path && x[2] %in% top5_rep_path})
edge_in_top5_path <- reaction_score_species[edge_in_top5_path,]$edge_name


reaction_data_SGMS1 <- reaction_score_species%>% 
  filter(edge_name %in%edge_in_top5_path) %>% 
  filter(p_value<0.05) %>% 
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color, FA_change))

#write.csv(reaction_data_SGMS1[-c(20:22)], file.path(file, 'source_data/reaction_data_SGMS1.csv'), row.names = F)

reaction_data_SGMS1 %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase','Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(reaction_data_SGMS1$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', legend.box = "vertical",
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = rev(pal_lancet()(2)))+
  scale_color_manual(values = c('gold','white'))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 significant edges',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))

#-------------------MPDU1-------------------

exp_MPDU1 <- data_process(exp_raw_2[c(1:4,20:22)], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                          replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp_MPDU1 <- remove_rownames(exp_MPDU1)

char_MPDU1 <- char_raw_2 %>% filter(feature %in% exp_MPDU1$feature)


each_FA <- str_extract_all(char_MPDU1$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_MPDU1 <- char_MPDU1 %>% mutate(each_FA=each_FA)



#Analysis for unprocessed data



no_sub_t <- unprocessed_data_test(exp_MPDU1, char_MPDU1, 't.test', 'adj_p_value', 1:3,  4:6)

#Lipid species substructure analysis
#Decompose lipids into species substructures
species_substructure <- species_sub_transform(char_MPDU1, lipid_substructure, 
                                              network_node)

#Extract species substructures using fold change

species_t <- no_sub_t[[2]] %>% filter(type=='species')

species_sub_stop <- species_sub_extract(species_substructure, species_t, 'species', pct_limit = 0.01)


#Transform lipid exp into substructure exp
species_sub_exp <- lipid_sub_matrix(exp_MPDU1, species_sub_stop, 'species')


#Differential expression analysis for substructures
species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:3, 4:6, 't.test', 'adj_p_value')

#Species biosynthetic network data transformation
species_network <- build_species_net(species_substructure)

SL_network <- species_substructure[str_detect(species_substructure$Lipid, '(Cer)|(SM)'),] %>% 
  unlist() %>% unique() %>% str_replace('FA\\d+_','')

SL_network <- species_network %>% filter(S1 %in% SL_network, P1 %in% SL_network)

#Essential pathway analysis for species substructures

set.seed(1)
path_score_species <-  path_scoring(SL_network, species_sub_exp_t,
                                    calibrate = T, data_type = 'species')


path_data_MPDU1 <- rbind(path_score_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_species %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(path_data_MPDU1, file.path(file, 'source_data/path_data_MPDU1.csv'), row.names = F)

path_data_MPDU1 %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  #scale_y_continuous(limits = c(-7,7))+
  scale_fill_manual(values = bluered(100)[c(1,10,20,30,40,100,90,80,70,60)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')

#-------------------ORMDL2-------------------

exp_ORMDL2 <- data_process(exp_raw_2[c(1:4,23:28)], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                          replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp_ORMDL2 <- remove_rownames(exp_ORMDL2)

char_ORMDL2 <- char_raw_2 %>% filter(feature %in% exp_ORMDL2$feature)


each_FA <- str_extract_all(char_ORMDL2$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_ORMDL2 <- char_ORMDL2 %>% mutate(each_FA=each_FA)



#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_ORMDL2, char_ORMDL2, 't.test', 'adj_p_value', 1:3,  4:9)

#Lipid species substructure analysis
#Decompose lipids into species substructures
species_substructure <- species_sub_transform(char_ORMDL2, lipid_substructure, 
                                              network_node)

#Extract species substructures using fold change
species_t <- no_sub_t[[2]] %>% filter(type=='species')

species_sub_stop <- species_sub_extract(species_substructure, species_t, 'species', pct_limit = 0.01)



#Transform lipid exp into substructure exp

species_sub_exp <- lipid_sub_matrix(exp_ORMDL2, species_sub_stop, 'species')

#Differential expression analysis for substructures
species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:3, 4:9, 't.test', 'adj_p_value')

#Species biosynthetic network data transformation

species_network <- build_species_net(species_substructure)

SL_network <- species_substructure[str_detect(species_substructure$Lipid, '(Cer)|(SM)'),] %>% 
  unlist() %>% unique() %>% str_replace('FA\\d+_','')

SL_network <- species_network %>% filter(S1 %in% SL_network, P1 %in% SL_network)

#Essential pathway analysis for species substructures

set.seed(1)
path_score_species <-  path_scoring(SL_network, species_sub_exp_t,
                                    calibrate = T, data_type = 'species')


path_data_ORMDL2 <- rbind(path_score_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_species %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(path_data_ORMDL2, file.path(file, 'source_data/path_data_ORMDL2.csv'), row.names = F)

path_data_ORMDL2 %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  #scale_y_continuous(limits = c(-7,7))+
  scale_fill_manual(values = bluered(100)[c(1,10,20,30,40,100,90,80,70,60)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')

#-------------------GNPAT-------------------

exp_GNPAT <- data_process(exp_raw_1[c(1:4,14:16)], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                          replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp_GNPAT <- remove_rownames(exp_GNPAT)

char_GNPAT <- char_raw_1 %>% filter(feature %in% exp_GNPAT$feature)


each_FA <- str_extract_all(char_GNPAT$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_GNPAT <- char_GNPAT %>% mutate(each_FA=each_FA)



#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_GNPAT, char_GNPAT, 't.test', 'adj_p_value', 1:3,  4:6)

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

char_exp_GNPAT <- Species2Char(exp_GNPAT, char_GNPAT, 'class')

char_exp_GNPAT[-1][is.na(char_exp_GNPAT[-1])] <- 0

class_sub_exp <- lipid_sub_matrix(char_exp_GNPAT, class_sub_stop, 'Class')

#Differential expression analysis for substructures

class_sub_exp_t <- t_test(class_sub_exp[[3]], 1:3, 4:6, 't.test', 'adj_p_value')

#Class biosynthetic network data transformation


class_network <- network_edge[c('S1','P1')] %>% filter(S1 %in% class_sub_exp_t$lipid,
                                                       P1 %in% class_sub_exp_t$lipid)

GPL_GL_network <- network_node %>% filter(Class %in%c('Glycerophospholipid', 'Ether lipid')) %>% 
  .$Abbreviation

GPL_GL_network <- class_network %>% filter(S1 %in% GPL_GL_network, P1 %in% GPL_GL_network)

#Essential pathway analysis for class substructures


set.seed(1)

path_score_class <-  path_scoring(GPL_GL_network, class_sub_exp_t,
                                    calibrate = T, data_type = 'Class')



path_data_GNPAT <- rbind(path_score_class %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                     path_score_class %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(path_data_GNPAT, file.path(file, 'source_data/path_data_GNPAT.csv'), row.names = F)

path_data_GNPAT %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  #scale_y_continuous(limits = c(-7,7))+
  scale_fill_manual(values = bluered(100)[c(1,100,90,80,70,60)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')


#Essential edges (reactions) analysis for species substructures

reaction_score_class <- reaction_scoring(network_edge,
                                           class_sub_exp[[3]],
                                           class_sub_exp_t,
                                           ctrl=1:3, exp=4:6,
                                           Species = 'human')

reaction_data_GNPAT <- reaction_score_class%>% 
  filter(p_value<0.05) %>% 
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color, FA_change))

#write.csv(reaction_data_GNPAT[-c(20:22)], file.path(file, 'source_data/reaction_data_GNPAT.csv'), row.names = F)

reaction_data_GNPAT %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase','Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(reaction_data_GNPAT$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', legend.box = "vertical",
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = rev(pal_lancet()(2)))+
  scale_color_manual(values = c('gold','white'))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 significant edges',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))


#-------------------CEPT1-------------------

exp_CEPT1 <- data_process(exp_raw_2[c(1:4,8:10)], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                           replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)


exp_CEPT1 <- remove_rownames(exp_CEPT1)

char_CEPT1 <- char_raw_2 %>% filter(feature %in% exp_CEPT1$feature)


each_FA <- str_extract_all(char_CEPT1$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_CEPT1 <- char_CEPT1 %>% mutate(each_FA=each_FA)



#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_CEPT1, char_CEPT1, 't.test', 'adj_p_value', 1:3,  4:6)

#Lipid class substructure analysis

#Extract class substructures using fold change


class_t <- no_sub_t[[2]] %>% filter(type=='class')

class_sub_stop <- species_sub_extract(filter(lipid_substructure, Lipid %in% class_t$lipid),
                                      class_t, 'class', pct_limit = 0.01)


#Transform lipid exp into substructure exp



char_exp_CEPT1 <- Species2Char(exp_CEPT1, char_CEPT1, 'class')

char_exp_CEPT1[-1][is.na(char_exp_CEPT1[-1])] <- 0

class_sub_exp <- lipid_sub_matrix(char_exp_CEPT1, class_sub_stop, 'class')


#Differential expression analysis for substructures

class_sub_exp_t <- t_test(class_sub_exp[[3]], 1:3, 4:6, 't.test', 'adj_p_value')

#Class biosynthetic network data transformation


class_network <- network_edge[c('S1','P1')] %>% filter(S1 %in% class_sub_exp_t$lipid,
                                                       P1 %in% class_sub_exp_t$lipid)

GPL_GL_network <- network_node %>% filter(Class %in%c('Glycerophospholipid', 'Ether lipid')) %>% 
  .$Abbreviation

GPL_GL_network <- class_network %>% filter(S1 %in% GPL_GL_network, P1 %in% GPL_GL_network)

#Essential pathway analysis for class substructures

set.seed(1)

path_score_class <-  path_scoring(GPL_GL_network, class_sub_exp_t,
                                    calibrate = T, data_type = 'Class')



path_data_CEPT1 <- rbind(path_score_class %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                     path_score_class %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(path_data_CEPT1, file.path(file, 'source_data/path_data_CEPT1.csv'), row.names = F)

path_data_CEPT1 %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  #scale_y_continuous(limits = c(-7,7))+
  scale_fill_manual(values = bluered(100)[c(1,10,20,30,100,90,80)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')


#Essential edges (reactions) analysis for class substructures



reaction_score_class <- reaction_scoring(network_edge,
                                         class_sub_exp[[3]],
                                         class_sub_exp_t,
                                         ctrl=1:3, exp=4:6,
                                         Species = 'human')


reaction_data_CEPT1 <- reaction_score_class%>% 
  filter(p_value<0.05) %>% 
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color, FA_change))

#write.csv(reaction_data_CEPT1[-c(20:22)], file.path(file, 'source_data/reaction_data_CEPT1.csv'), row.names = F)

reaction_data_CEPT1 %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase','Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(reaction_data_CEPT1$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', legend.box = "vertical",
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = rev(pal_lancet()(2)))+
  scale_color_manual(values = c('gold','white'))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 significant edges',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))


#-------------------ACOT7-------------------

exp_ACOT7 <- data_process(exp_raw_2[c(1:4,4:6)], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                           replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp_ACOT7 <- remove_rownames(exp_ACOT7)

char_ACOT7 <- char_raw_2 %>% filter(feature %in% exp_ACOT7$feature)


each_FA <- str_extract_all(char_ACOT7$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_ACOT7 <- char_ACOT7 %>% mutate(each_FA=each_FA)



#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_ACOT7, char_ACOT7, 't.test', 'adj_p_value', 1:3,  4:6)

#FA substructure analysis
#FA biosynthetic network data transformation

FA_network_ACOT7 <- build_FA_net(FA_network, no_sub_t)

#Decompose lipids into FA substructures
#18:2 and 20:4 are majorly omega-6 FAs, so we only kept omega-6 forms of them


FA_substructure <- FA_sub_transform(FA_network_ACOT7, no_sub_t,  c('w9-18:2;0','w3-20:4;0'))


#Extract FA substructure using fold changes


FA_t <- no_sub_t[[2]] %>% filter(type=='FA')

FA_sub_stop <- FA_sub_extract(char_ACOT7, FA_substructure, FA_t,
                              exact_FA='no', exo_lipid=NULL)

#Transform lipid FA exp data into substructure exp data
FA_sub_exp <- lipid_sub_matrix(exp_ACOT7, FA_sub_stop, 'FA')


#Differential expression analysis for FA substructures
FA_sub_exp_t <- t_test(FA_sub_exp[[3]], 1:3, 4:6, 't.test', 'adj_value')


#Essential pathway analysis for FA substructures


set.seed(1)
path_score_FA <- path_scoring(FA_network_ACOT7, FA_sub_exp_t,
                              calibrate = T, data_type ='FA')


path_data_ACOT7 <- rbind(path_score_FA %>% filter(Type=='Active') %>% 
                    .[!duplicated(.$rep_sub_path),] %>% .[1:5,], 
                  path_score_FA %>% filter(Type=='Suppressed') %>% 
                    arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>%.[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(path_data_ACOT7, file.path(file, 'source_data/path_data_ACOT7.csv'), row.names = F)

path_data_ACOT7 %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(1.96,-1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(limits = c(-4.5,4.5))+
  scale_fill_manual(values = bluered(100)[c(100,90,80,70, 1,10,20)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 1),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Significant representative pathways')

#-------------------DECR2-------------------

exp_DECR2 <- data_process(exp_raw_2[c(1:4,11:13)], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                          replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp_DECR2 <- remove_rownames(exp_DECR2)

char_DECR2 <- char_raw_2 %>% filter(feature %in% exp_DECR2$feature)


each_FA <- str_extract_all(char_DECR2$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_DECR2 <- char_DECR2 %>% mutate(each_FA=each_FA)



#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_DECR2, char_DECR2, 't.test', 'adj_p_value', 1:3,  4:6)

#FA substructure analysis
#FA biosynthetic network data transformation

FA_network_DECR2 <- build_FA_net(FA_network, no_sub_t)

#Decompose lipids into FA substructures
#18:2 and 20:4 are majorly omega-6 FAs, so we only kept omega-6 forms of them
FA_substructure <- FA_sub_transform(FA_network_DECR2, no_sub_t,  c('w9-18:2;0','w3-20:4;0'))



#Extract FA substructure using fold changes


FA_t <- no_sub_t[[2]] %>% filter(type=='FA')

FA_sub_stop <- FA_sub_extract(char_DECR2, FA_substructure, FA_t,
                              exact_FA='no', exo_lipid=NULL)


#Transform lipid FA exp data into substructure exp data



FA_sub_exp <- lipid_sub_matrix(exp_DECR2, FA_sub_stop, 'FA')
#Differential expression analysis for FA substructures
FA_sub_exp_t <- t_test(FA_sub_exp[[3]], 1:3, 4:6, 't.test', 'adj_value')


#Essential pathway analysis for FA substructures



set.seed(1)
path_score_FA <- path_scoring(FA_network_DECR2, FA_sub_exp_t,
                              calibrate = T, data_type ='FA')


path_data_DECR2 <- rbind(path_score_FA %>% filter(Type=='Active') %>% 
                    .[!duplicated(.$rep_sub_path),] %>% .[1:5,], 
                  path_score_FA %>% filter(Type=='Suppressed') %>% 
                    arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>%.[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(path_data_DECR2, file.path(file, 'source_data/path_data_DECR2.csv'), row.names = F)

path_data_DECR2 %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(1.96,-1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(limits = c(-4.5,4.5))+
  scale_fill_manual(values = bluered(100)[c(100,90,1,10,20,30)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 1),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Significant representative pathways')

#-------------------ELOVL5-------------------

exp_ELOVL5 <- data_process(exp_raw_2[c(1:4,14:16)], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                          replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp_ELOVL5 <- remove_rownames(exp_ELOVL5)

char_ELOVL5 <- char_raw_2 %>% filter(feature %in% exp_ELOVL5$feature)


each_FA <- str_extract_all(char_ELOVL5$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_ELOVL5 <- char_ELOVL5 %>% mutate(each_FA=each_FA)



#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_ELOVL5, char_ELOVL5, 't.test', 'adj_p_value', 1:3,  4:6)

#FA substructure analysis
#FA biosynthetic network data transformation
FA_network_ELOVL5 <- build_FA_net(FA_network, no_sub_t)


#Decompose lipids into FA substructures
#18:2 and 20:4 are majorly omega-6 FAs, so we only kept omega-6 forms of them
FA_substructure <- FA_sub_transform(FA_network_ELOVL5, no_sub_t,  c('w9-18:2;0','w3-20:4;0'))



#Extract FA substructure using fold changes

FA_t <- no_sub_t[[2]] %>% filter(type=='FA')

FA_sub_stop <- FA_sub_extract(char_ELOVL5, FA_substructure, FA_t,
                              exact_FA='no', exo_lipid=NULL)

#Transform lipid FA exp data into substructure exp data

FA_sub_exp <- lipid_sub_matrix(exp_ELOVL5, FA_sub_stop, 'FA')

#Differential expression analysis for FA substructures

FA_sub_exp_t <- t_test(FA_sub_exp[[3]], 1:3, 4:6, 't.test', 'adj_value')

#Essential pathway analysis for FA substructures


set.seed(1)
path_score_FA <- path_scoring(FA_network_ELOVL5, FA_sub_exp_t,
                              calibrate = T, data_type ='FA')

path_data_ELOVL5 <- rbind(path_score_FA %>% filter(Type=='Active') %>% 
                    .[!duplicated(.$rep_sub_path),] %>% .[1:5,], 
                  path_score_FA %>% filter(Type=='Suppressed') %>% 
                    arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>%.[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(path_data_ELOVL5, file.path(file, 'source_data/path_data_ELOVL5.csv'), row.names = F)

path_data_ELOVL5 %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(1.96,-1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(limits = c(-4.5,4.5))+
  scale_fill_manual(values = bluered(100)[c(100,90,80,60,70,1,10,20,30)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 1),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Significant representative pathways')



#Essential edges (reactions) analysis for FA substructures

reaction_score_FA <- reaction_scoring(FA_network, 
                                      FA_sub_exp[[3]], FA_sub_exp_t, 
                                      ctrl = 1:3, exp = 4:6, 
                                      Species = 'human')


reaction_data_ELOVL5 <-reaction_score_FA %>%
  filter(!is.na(perturbation_score)) %>% 
  filter(p_value<0.05) %>% 
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color))

#write.csv(reaction_data_ELOVL5[-c(20:22)], file.path(file, 'source_data/reaction_data_ELOVL5.csv'), row.names = F)


reaction_data_ELOVL5 %>%
  mutate(Mode=factor(Mode, levels = c("Increase","Decrease"))) %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(reaction_data_ELOVL5$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', 
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = rev(pal_lancet()(2)))+
  scale_color_manual(values = c('gold','white'))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 significant edges',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))


#-------------------HSD17B12-------------------

exp_HSD17B12 <- data_process(exp_raw_2[c(1:4,17:19)], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                           replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp_HSD17B12 <- remove_rownames(exp_HSD17B12)

char_HSD17B12 <- char_raw_2 %>% filter(feature %in% exp_HSD17B12$feature)


each_FA <- str_extract_all(char_HSD17B12$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_HSD17B12 <- char_HSD17B12 %>% mutate(each_FA=each_FA)



#Analysis for unprocessed data

no_sub_t <- unprocessed_data_test(exp_HSD17B12, char_HSD17B12, 't.test', 'adj_p_value', 1:3,  4:6)

#FA substructure analysis
#FA biosynthetic network data transformation

FA_network_HSD17B12 <- build_FA_net(FA_network, no_sub_t)

#Decompose lipids into FA substructures
#18:2 and 20:4 are majorly omega-6 FAs, so we only kept omega-6 forms of them

FA_substructure <- FA_sub_transform(FA_network_HSD17B12, no_sub_t,  c('w9-18:2;0','w3-20:4;0'))


FA_t <- no_sub_t[[2]] %>% filter(type=='FA')

#Extract FA substructure using fold changes



FA_sub_stop <- FA_sub_extract(char_HSD17B12, FA_substructure, FA_t,
                              exact_FA='no', exo_lipid=NULL)


#Transform lipid FA exp data into substructure exp data
FA_sub_exp <- lipid_sub_matrix(exp_HSD17B12, FA_sub_stop, 'FA')


#Differential expression analysis for FA substructures

FA_sub_exp_t <- t_test(FA_sub_exp[[3]], 1:3, 4:6, 't.test', 'adj_value')

#Essential pathway analysis for FA substructures


set.seed(1)
path_score_FA <- path_scoring(FA_network_HSD17B12, FA_sub_exp_t,
                              calibrate = T, data_type ='FA')


path_data_HSD17B12 <- rbind(path_score_FA %>% filter(Type=='Active') %>% 
                    .[!duplicated(.$rep_sub_path),] %>% .[1:5,], 
                  path_score_FA %>% filter(Type=='Suppressed') %>% 
                    arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>%.[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(path_data_HSD17B12, file.path(file, 'source_data/path_data_HSD17B12.csv'), row.names = F)

path_data_HSD17B12 %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(1.96,-1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(limits = c(-4.5,4.5))+
  scale_fill_manual(values = bluered(100)[c(1,10,20,30)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 1),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Significant representative pathways')


#Essential edges (reactions) analysis for FA substructures

reaction_score_FA <- reaction_scoring(FA_network, 
                                      FA_sub_exp[[3]], FA_sub_exp_t, 
                                      ctrl = 1:3, exp = 4:6, 
                                      Species = 'human')


reaction_data_HSD17B12 <-reaction_score_FA %>%
  filter(!is.na(perturbation_score)) %>% 
  filter(p_value<0.05) %>% 
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color))

#write.csv(reaction_data_HSD17B12[-c(20:22)], file.path(file, 'source_data/reaction_data_HSD17B12.csv'), row.names = F)


reaction_data_HSD17B12 %>%
  mutate(Mode=factor(Mode, levels = c("Increase","Decrease"))) %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(reaction_data_HSD17B12$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', 
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = rev(pal_lancet()(2)))+
  scale_color_manual(values = c('gold','white'))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 significant edges',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))

#-------------------Save data-------------------
#save.image(file.path(file, 'Lipid_gene_KO.RData'))
