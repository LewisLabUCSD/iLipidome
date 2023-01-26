library(ggsci)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(ggtext)


file <- dirname(rstudioapi::getSourceEditorContext()$path)

#-------------------Data upload-------------------

source(file.path(dirname(file),'Required_function/required_function.R'))
load(file.path(dirname(file),'Required_data/required_data.RData'))


exp_raw <- read.csv(file.path(file,'lipidome_data/exp_LPCAT1_raw.csv'))
char_raw <- read.csv(file.path(file,'lipidome_data/char_LPCAT1_raw.csv'))

char_raw <- char_raw %>% 
  mutate(FA_sum=str_c(totallength,':',totaldb,';0'),
         FA_num=ifelse(totallength>24,2,1)) %>% 
  mutate(FA_split=ifelse(FA_num==1, FA_sum,'')) %>% 
  mutate(feature=str_c(class,'_',FA_sum))

exp_raw$feature <- char_raw$feature

each_FA <- str_extract_all(char_raw$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char_unique <- char_raw %>% mutate(each_FA=each_FA) %>% unique() %>% 
  remove_rownames()

lipid <- character()
sum_exp <- list()
add <- 1
for(num in unique(exp_raw$feature)){
  lipid[add] <- num
  sum_exp[[add]] <- exp_raw %>% filter(feature==num) %>% .[-1] %>% colSums()
  add <- add+1
}

exp_unique <- Reduce(rbind, sum_exp) %>% as.data.frame()

exp_unique <- exp_unique %>% mutate(feature=lipid) %>% 
  dplyr::select(feature, everything()) %>% remove_rownames()



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


exp <- rbind(data_process(exp_unique[-c(81:84),], exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                    replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F),
             exp_unique[c(81:84),])

char <- char_unique


#-------------------Analysis for unprocessed data-------------------
#Remove Free fatty acids
lipid_exp <- exp[-c(81:84),]
lipid_char <- char[-c(81:84),]

no_sub_t <- unprocessed_data_test(lipid_exp, lipid_char, 't.test', 'adj_p_value', 1:5,  6:10)
#-------------------Conventional lipidomics analysis-------------------

SF3a_data <- DE_lipid_data <- no_sub_t[[2]] %>% filter(type=='species') %>%  
  mutate(log2FC=ifelse(is.infinite(log2FC),10*sign(log2FC),log2FC)) %>% 
  mutate(Significance=ifelse(log2FC>0, 'Increase','Decrease')) %>% 
  mutate(Significance=ifelse(sig=='yes', Significance, 'No change'))

#write.csv(SF3a_data, file.path(file, 'source_data/SF3a.csv'), row.names = F)


SF3a <- SF3a_data %>% 
  mutate(label=ifelse(lipid %in% c('PC_28:0;0', 'PC_30:0;0', 'PC_32:0;0'), lipid, '')) %>% 
  ggplot(aes(x=log2FC, y=mlog10padj, col=Significance)) +
  geom_point() + 
  scale_color_manual(values=c("blue", "red", "gray")) +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='dashed')+
  geom_text_repel(aes(label=label))+
  theme_classic()+
  theme(legend.position = 'right')+
  labs(y='-log10(padj)')

#-------------------Lipid species substructure analysis-------------------

#Decompose lipids into substructures

species_substructure <- species_sub_transform(lipid_char, lipid_substructure, 
                                              network_node)


#Original data provide exact FA for PS, allowing us to correct PS-LPS connections

LPS_modify <- function(lipid_char, lipid_original_name, species_substructure){
  LPS_new_sub <- list()
  num <- 1
  PS <- lipid_original_name %>% filter(str_detect(feature, '^PS')) %>% 
    .$original_name %>% str_extract_all('\\d+.\\d+') %>% 
    map_chr(.f = function(x){str_c(sort(x)[1],';0_',sort(x)[2],';0')})
  
  LPS <- lipid_char %>% filter(class=='LPS') %>% .$feature
  for(a in LPS){
    LPS_FA <- lipid_char %>% filter(feature==a) %>% .$FA_split
    LPS_sub <- species_substructure %>% filter(Lipid==a)
    
    PS_sub <- species_substructure %>% filter(Lipid==a) %>%
      apply(MARGIN = 1, FUN = function(x){tail(x[x!=''],2)[1]})
    exist_PS <- PS[str_detect(PS, LPS_FA)]
    if(length(exist_PS)==0){
      LPS_new_sub[[num]] <- LPS_sub
      num <- num+1
    }
    else{
      exist_PS_FA_sum <- str_split(exist_PS, '_') %>% 
        map_chr(.f = function(x){a=unlist(str_extract_all(x[1],'\\d+')) %>% as.integer();
        b=unlist(str_extract_all(x[2],'\\d+')) %>% as.integer();
        str_c(a[1]+b[1],':',a[2]+b[2],';',a[3]+b[3])})
      exist_PS_FA_sum <- str_c(exist_PS_FA_sum, collapse = '|')
      exist_PS_FA_sum <- str_detect(PS_sub, exist_PS_FA_sum)
      if(length(exist_PS)!=0){
        LPS_new_sub[[num]] <- LPS_sub[which(exist_PS_FA_sum),]
        num <- num+1
      }
      else{
        LPS_new_sub[[num]] <- LPS_sub
        num <- num+1
      }
    }
  }
  return(LPS_new_sub)
  
}
lipid_original_name <- read.csv(file.path(file,'lipidome_data/lipid_original_name.csv'))

LPS_m <- LPS_modify(lipid_char, lipid_original_name, species_substructure)
LPS_m <- rbindlist(LPS_m)

species_substructure <- species_substructure %>% filter(!str_detect(Lipid,'LPS')) %>% 
  rbind(LPS_m)


#Extract species substructures using fold changes

species_t <- no_sub_t[[2]] %>% filter(type=='species')

species_sub_stop <- species_sub_extract(species_substructure, species_t,
                                        'species', pct_limit = 0.1)


#Transform lipid exp into substructure exp

species_sub_exp <- lipid_sub_matrix(lipid_exp, species_sub_stop, 'species')


#Differential expression analysis for substructures

species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:5, 6:10,
                            't.test', 'adj_p_value')




#Species biosynthetic network data transformation

species_network <- build_species_net(species_substructure)


#Essential pathway analysis for species substructures

set.seed(1)
path_score_species <-  path_scoring(species_network, species_sub_exp_t,
                                    calibrate = T, data_type = 'species')


F5c_data <- path_score_species %>% 
  arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% 
  .[1:5,] %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


F5c_data$path[1:5] <- F5c_data$path[1:5] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], tail(x,2)[1],tail(x,2)[2],sep=' --> ')})


#write.csv(F5c_data, file.path(file, 'source_data/F5c.csv'), row.names = F)




F5c <- F5c_data %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  scale_fill_manual(values = bluered(100)[c(1,10,20,30,40)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 suppressed pathways')


#Essential edges (reactions) analysis for species substructures

species_net_w_rev <- add_rev_rection(network_edge, species_network)


reaction_score_species <- reaction_scoring(species_net_w_rev,
                                           species_sub_exp[[3]],
                                           species_sub_exp_t,
                                           ctrl=1:5, exp=6:10,
                                           Species = 'human')


top5_rep_path <- path_score_species %>% filter(Type=='Suppressed') %>% 
  arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path

top5_rep_path <- path_score_species %>% filter(rep_sub_path%in%top5_rep_path, Significant=='yes') %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top5_path <- reaction_score_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top5_rep_path && x[2] %in% top5_rep_path})
edge_in_top5_path <- reaction_score_species[edge_in_top5_path,]$edge_name



F5d_data <- reaction_score_species%>% 
  filter(edge_name %in%edge_in_top5_path) %>% 
  filter(str_detect(edge_name, 'L.+ --> P.+|P.+ --> L.+')) %>% 
  filter(!str_detect(edge_name, 'LPA')) %>% 
  .[c((nrow(.)-4):nrow(.)),] %>% 
  filter(p_value<0.05) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color, FA_change))


#write.csv(F5d_data[-c(20:22)], file.path(file, 'source_data/F5d.csv'), row.names = F)


F5d <- F5d_data %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(F5d_data$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', legend.box = "vertical",
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = pal_lancet()(2))+
  scale_color_manual(values = c('gold','white'))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 suppressed edges (Lands cycle)',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))



#Construct species biosynthetic network


species_network_data <- draw_network(species_net_w_rev, species_sub_exp_t,
                                     if_species = T,significant = 'adj_p_value',
                                     path_scoring_result = path_score_species,
                                     reaction_scoring_result = filter(reaction_score_species, !str_detect(edge_name, 'LPA')),
                                     top_n = 5, path_type = 'suppressed')


PC_LPC_path <- path_score_species %>% filter(Type=='Suppressed') %>% 
  arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[2:4,] %>% .$rep_sub_path

PC_LPC_path <- path_score_species %>% filter(rep_sub_path%in%PC_LPC_path, Significant=='yes') %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

PC_LPC_network_data <- list()
PC_LPC_network_data[[2]] <- species_network_data[[2]] %>% filter(from  %in% PC_LPC_path, to %in% PC_LPC_path)
PC_LPC_network_data[[1]] <- species_network_data[[1]] %>% filter(id  %in% PC_LPC_path)

F5b <- visNetwork(PC_LPC_network_data[[1]], PC_LPC_network_data[[2]])


species_sub_t_net <- species_sub_exp_t

maxinf <- ceiling(max(species_sub_t_net$log2FC[is.finite(species_sub_t_net$log2FC)]))
mininf <- floor(min(species_sub_t_net$log2FC[is.finite(species_sub_t_net$log2FC)]))
species_sub_t_net$log2FC[species_sub_t_net$log2FC>0 & is.infinite(species_sub_t_net$log2FC)] <- maxinf
species_sub_t_net$log2FC[species_sub_t_net$log2FC<0 & is.infinite(species_sub_t_net$log2FC)] <- mininf


F5b_data <- PC_LPC_network_data[[2]]%>% left_join(species_sub_t_net, by=c('from'='lipid')) %>% 
  left_join(species_sub_t_net, by=c('to'='lipid')) %>% 
  mutate(from=str_replace_all(from, ';0','')) %>% 
  mutate(to=str_replace_all(to, ';0','')) %>% 
  .[c(1:9,16,21,22,29,34,35)]


#write.xlsx(F5b_data, file.path(file, 'source_data/F5b.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)



#BIOPAN's algorithm uses free FA data to build the connections involving FA transfer
biopan_network_edge <- PC_LPC_network_data[[2]] %>% filter(from %in% lipid_char$feature, to %in% lipid_char$feature)

biopan_network_edge <- biopan_network_edge %>% left_join(char[c('feature','class')], by=c('from'='feature')) %>% 
  left_join(char[c('feature','class')], by=c('to'='feature')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.x'='Abbreviation')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.y'='Abbreviation')) %>% 
  mutate(FA=str_c(abs(as.integer(str_extract(from, '\\d+'))-
                        as.integer(str_extract(to, '\\d+'))),':',
                  abs(as.integer(str_extract_all(from, '\\d+',simplify = T)[1,2])-
                        as.integer(str_extract_all(to, '\\d+',simplify = T)[1,2])))) %>% 
  filter(FA.x==FA.y | FA %in% str_extract(char$feature[char$class=='FFA'], '\\d+:\\d+'))

biopan_network_node <- PC_LPC_network_data[[1]]  %>% filter(id %in% char$feature)

SF5b <- visNetwork(biopan_network_node, biopan_network_edge)

SF5b_data1 <- biopan_network_node %>% 
  left_join(species_sub_t_net, by=c('id'='lipid')) %>% 
  mutate(id=str_replace_all(id, ';0','')) %>% 
  mutate(label=str_replace_all(label, ';0','')) %>% 
  .[c(1:7,14,19,20)]

SF5b_data2 <- biopan_network_edge[1:9] %>% left_join(species_sub_t_net, by=c('from'='lipid')) %>% 
  left_join(species_sub_t_net, by=c('to'='lipid')) %>% 
  mutate(from=str_replace_all(from, ';0','')) %>% 
  mutate(to=str_replace_all(to, ';0','')) %>% 
  .[c(1:9,16,21,22,29,34,35)]

#write.xlsx(SF5b_data1, file.path(file, 'source_data/SF5b1.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)
#write.xlsx(SF5b_data2, file.path(file, 'source_data/SF5b2.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)

#-------------------Species substructures improve statistical power-------------------


SF3b_data <- species_sub_exp_t %>% 
  left_join(species_t, by=c('lipid')) %>% 
  .[c('lipid', 'mlog10padj.y', 'mlog10padj.x')] %>% 
  `colnames<-`(c('lipid','Lipid species','Substructure')) %>% 
  filter(`Lipid species`>-log10(0.05)| Substructure>-log10(0.05)) 

#write.csv(SF3b_data, file.path(file, 'source_data/SF3b.csv'), row.names = F, na = '')


SF3b <- SF3b_data %>% 
  filter(!is.na(`Lipid species`)) %>% 
  ggpaired(cond1 = "Lipid species", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = '-log10(padj)')+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed')+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.24, label.y = 7.3)+
  theme(legend.position = 'none')



#-------------------Free fatty acid analysis-------------------
options(scipen = -1)
F5e_data <- exp[81,] %>% 
  gather(-feature, key='label',value='value') %>% 
  mutate(group=str_extract(label,'[A-Za-z]+')) %>% 
  mutate(group=ifelse(group=='Ctrl','shNT','shLPCAT1'))

#write.csv(F5e_data, file.path(file, 'source_data/F5e.csv'), row.names = F, na='')


F5e <- F5e_data%>% 
  ggboxplot(x = "group", y = "value",
            color = "group", palette = 'jco',
            add = "jitter")+ 
  stat_compare_means(aes(group=group),method = "t.test",label.y = 440, 
                     label.x = 1.18)+
  theme_classic() +
  #scale_y_continuous(limits = c(0,25))+
  labs(x='',title='', y='16:0 FA abundance')
#-------------------Raw lipid species analysis-------------------


#Essential pathway analysis for raw data

raw_species_net <- species_network

raw_species_net <- raw_species_net%>% filter(S1%in%species_t$lipid,
                                             P1%in%species_t$lipid)

set.seed(1)
path_score_raw_species <-  path_scoring(raw_species_net, species_t,
                                        calibrate = T, data_type = 'species')


SF3c_data <- path_score_raw_species %>% 
  arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% 
  .[1:5,] %>% 
  mutate(path=str_replace_all(path, ';0',''))


#write.csv(SF3c_data, file.path(file, 'source_data/SF3c.csv'), row.names = F)


SF3c <- SF3c_data %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score), fill='gray')+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 suppressed pathways')




#Essential edges (reactions) analysis for raw data

raw_species_net_w_rev <- add_rev_rection(network_edge, raw_species_net)

reaction_score_raw_species <- reaction_scoring(raw_species_net_w_rev,
                                               column_to_rownames(lipid_exp, var = 'feature'),
                                               species_t,ctrl=1:5, exp=6:10,
                                               Species = 'human')



top5_rep_path_raw <- path_score_raw_species %>% filter(Type=='Suppressed') %>% 
  arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path

top5_rep_path_raw <- path_score_raw_species %>% filter(rep_sub_path%in%top5_rep_path_raw) %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top5_path_raw <- reaction_score_raw_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top5_rep_path_raw && x[2] %in% top5_rep_path_raw})
edge_in_top5_path_raw <- reaction_score_raw_species[edge_in_top5_path_raw,]$edge_name



SF3d_data <- reaction_score_raw_species%>% 
  filter(edge_name %in%edge_in_top5_path_raw) %>% 
  filter(str_detect(edge_name, 'L.+ --> P.+|P.+ --> L.+')) %>% 
  filter(!str_detect(edge_name, 'LPA')) %>% 
  .[c((nrow(.)-4):nrow(.)),] %>% 
  filter(p_value<0.05) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color, FA_change))

#write.csv(SF3d_data[-c(20:22)], file.path(file, 'source_data/SF3d.csv'), row.names = F)

SF3d <- SF3d_data %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(SF3d_data$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', legend.box = "vertical",
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = pal_lancet()(2))+
  scale_color_manual(values = c('gold','white'))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 suppressed edges (Lands cycle)',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))


#-------------------Save data-------------------
#save.image(file.path(file, 'LPCAT1.RData'))

