library(tidyverse)
library(plyr)
library(gplots)
library(ggsci)
library(ggpubr)
library(ggvenn)
library(gridExtra)
library(ggrepel)
library(ggtext)
library(igraph)
library(visNetwork)
library(xlsx)
library(gtools)
library(data.table)

file <- dirname(rstudioapi::getSourceEditorContext()$path)
#load(file.path(file,'LPCAT1.RData'))

#-------------------Data upload-------------------

load(file.path(file,'Required_data/Required_function.RData'))
load(file.path(file,'Required_data/Required_data.RData'))

each_FA <- str_extract_all(char$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char <- char %>% mutate(each_FA=each_FA)


#-------------------Analysis for unprocessed data-------------------
#Remove Free fatty acids
lipid_exp <- exp[-c(81:84),]
lipid_char <- char[-c(81:84),]

no_sub_t <- non_processed_data_test(lipid_exp, lipid_char, 't.test', 'adj_p_value', 2:6,  7:11)


#-------------------Decompose lipids into substructure-------------------

sub_species <- lipid_species_substructure_transform(lipid_char, FA_substructure, lipid_substructure)


species_substructure <- sub_species[[2]] %>% map(.f = function(x){names(x) <- 1:length(x);return(x)}) %>% 
  plyr::ldply(rbind) %>% mutate(Lipid=sub_species[[1]]) %>% 
  dplyr::select(Lipid, everything())

species_substructure[is.na(species_substructure)] <- ''

colnames(species_substructure) <- c('Lipid', str_c('Unit', 1:(ncol(species_substructure)-1)))

#-------------------original data provides FA composition for PS----------------
#so that we use this information to correct PS-LPS connections

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

LPS_m <- LPS_modify(lipid_char, lipid_original_name, species_substructure)
LPS_m <- rbindlist(LPS_m)

species_substructure <- species_substructure %>% filter(!str_detect(Lipid,'LPS')) %>% 
  rbind(LPS_m)


#-------------------extract substructure using fold changes-------------------

species_t <- no_sub_t[[2]] %>% filter(type=='species')


species_sub_stop <- lipid_substructure_w_stop(species_substructure, species_t, 'species', 0.1)

species_sub_stop <- list(species_sub_stop[[1]], apply(species_sub_stop[-1], MARGIN = 1, FUN = function(x){x[x!='']}))

#-------------------transform lipid exp data into substructue exp data-------------------

species_sub_exp <- lipid_substructure_matrix(lipid_exp, species_sub_stop, 'species')


#-------------------differential expression analysis for substructures-------------------

species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:5, 6:10, 't.test', 'adj_p_value')

#-------------------substructures improve statistical power-------------------

species_t <- no_sub_t[[2]] %>% filter(type=='species')

stat_improve_fig <- species_sub_exp_t %>% 
  left_join(species_t, by=c('lipid')) %>% 
  .[c('lipid', 'mlog10padj.y', 'mlog10padj.x')] %>% 
  `colnames<-`(c('lipid','Lipid species','Substructure')) %>% 
  filter(!is.na(`Lipid species`)) %>% 
  filter(`Lipid species`>-log10(0.05)| Substructure>-log10(0.05)) %>% 
  ggpaired(cond1 = "Lipid species", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = '-log10(padj)')+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed')+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.24, label.y = 7.3)+
  theme(legend.position = 'none')




#-------------------biosynthetic network data transformation-------------------


species_network <- split(species_substructure[-1], seq(nrow(species_substructure)))
species_network <- species_network %>% map(.f = function(x){x[x!='']})

species_network <- species_network %>% map(.f = function(x){head(rep(x, each=2)[-1],-1)})
species_network <-  matrix(unlist(species_network), ncol=2, byrow = T,dimnames = list(NULL, c('S1','P1'))) %>% 
  as.data.frame() %>%  unique()

species_network <- species_network %>% mutate(S1=str_replace(S1, '_FA\\d',''),
                                              P1=str_replace(P1, '_FA\\d',''))


species_network_data <- draw_network(species_network, species_sub_exp_t, 'no','adj_p_value')

#draw network

#species_network_edge <- species_network_data[[2]]

#species_network_node <- species_network_data[[1]
                                             
#species_network_edge$color <- 'gray'

#visNetwork(species_network_node, species_network_edge)


#-------------------essential pathway analysis for substructures-------------------

species_net <- species_network_data[[2]][1:2]
colnames(species_net) <- c('S1','P1')


set.seed(1)
path_score_species <-  calculate_path_activation_node(species_net, species_sub_exp_t,
                                                      calibrate = T, if_FA = 'no')


rep_lipid <- path_score_species%>% 
  #filter(Significant=='yes') %>% 
  .$path %>% 
  str_extract_all('\\d+:\\d+') %>% 
  map(.f = function(x){k=table(x) %>% sort(decreasing = T);
  if(length(k)==0){''}
  else if(sum(k!=1)==0){sort(names(k), decreasing = T)[1]}
  else if(sum(k!=1)!=0){names(k)[1]}}) %>% 
  unlist


path_score_species <- path_score_species%>% 
  mutate(rep_sub_path=rep_lipid)


top5_path <- path_score_species %>% 
  filter(Type=='Suppressed') %>% 
  arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,]

top5_path_fig <- top5_path
top5_path_fig$path[1:5] <- top5_path_fig$path[1:5] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], tail(x,2)[1],tail(x,2)[2],sep=' --> ')})


top5_path_fig <- top5_path_fig %>% 
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

#write.csv(path_score_species, file.path(file, 'result/path_score_subspecies_LPCAT1.csv'), row.names = F)


#-------------------essential edges (reactions) analysis for substructures-------------------

species_net_w_rev <- add_rev_rection(network_edge, species_net)



#calculate perturbation score for each edge (reaction)
perturbation_score_species <- calculate_perturbation_point(species_net_w_rev,
                                                           species_sub_exp[[3]],
                                                           species_sub_exp_t,
                                                           ctrl=1:5, exp=6:10,
                                                           stat = 'p')



top5_rep_path <- path_score_species %>% filter(Type=='Suppressed') %>% 
  arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path

top5_rep_path <- path_score_species %>% filter(rep_sub_path%in%top5_rep_path, Significant=='yes') %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top5_path <- perturbation_score_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top5_rep_path && x[2] %in% top5_rep_path})
edge_in_top5_path <- perturbation_score_species[edge_in_top5_path,]$edge_name



top5_edge_lands_cycle <- perturbation_score_species%>% 
  filter(edge_name %in%edge_in_top5_path) %>% 
  filter(str_detect(edge_name, 'L.+ --> P.+|P.+ --> L.+')) %>% 
  filter(!str_detect(edge_name, 'LPA')) %>% 
  .[c((nrow(.)-4):nrow(.)),] %>% 
  filter(p_value<0.05) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color, FA_change))


top5_edge_lands_cycle_fig <- top5_edge_lands_cycle %>% ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
                          fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(top5_edge_lands_cycle$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', legend.box = "vertical",
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = pal_lancet()(2))+
  scale_color_manual(values = c('white','gold'))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 suppressed edges (Lands cycle)',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))

#write.csv(perturbation_score_species, file.path(file, 'result/perturbation_score_subspecies_LPCAT1.csv'), row.names = F)


#-------------------construct biosynthetic network-------------------



top5_net_edge <- species_net_w_rev %>% filter(S1 %in% top5_rep_path, P1 %in% top5_rep_path)


top5_net_edge <- top5_net_edge %>% mutate(color='gray', arrows='to',length=100)

colnames(top5_net_edge)[1:2] <- c('from','to')





top5_net_node <- species_network_data[[1]] %>% 
  filter(id %in% unlist(top5_net_edge))


top5_net_path_col <- path_score_species %>% filter(Type=='Suppressed') %>% 
  arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% 
  mutate(rank=c(6:10))

top5_net_reaction_col <- perturbation_score_species%>% 
  filter(edge_name %in%edge_in_top5_path) %>% 
  filter(str_detect(edge_name, 'L.+ --> P.+|P.+ --> L.+')) %>% 
  filter(!str_detect(edge_name, 'LPA')) %>% 
  .[c((nrow(.)-4):nrow(.)),] %>% 
  filter(p_value<0.05)%>% 
  mutate(rank=c(6:10))


top5_net_node$label <- str_replace_all(top5_net_node$label , ';0','')

net_edge <- path_color(top5_net_edge, top5_net_path_col,top5_net_reaction_col)


net_node <- top5_net_node %>% filter(id %in%c(net_edge$from, net_edge$to))



#building network

visNetwork(net_node, net_edge)



species_sub_t_net <- species_sub_exp_t
maxinf <- ceiling(max(species_sub_t_net$log2FC[is.finite(species_sub_t_net$log2FC)]))
mininf <- ceiling(min(species_sub_t_net$log2FC[is.finite(species_sub_t_net$log2FC)]))
species_sub_t_net$log2FC[species_sub_t_net$log2FC>0 & is.infinite(species_sub_t_net$log2FC)] <- maxinf
species_sub_t_net$log2FC[species_sub_t_net$log2FC<0 & is.infinite(species_sub_t_net$log2FC)] <- mininf


# path_color(top5_net_edge, top5_net_path_col,top5_net_reaction_col) %>% 
#    left_join(species_sub_t_net, by=c('from'='lipid')) %>% 
#    left_join(species_sub_t_net, by=c('to'='lipid')) %>% 
#    .[c(1:9,16,21,22,29,34,35)] %>% 
#    write.csv(file.path(file, 'network/subspecies_network_LPCAT1.csv'),
#              row.names = F, na = '')


#BIOPAN's algorithm cannot build the connections involving FA transfer if no free FA data provided
net_edge_biopan <- net_edge %>% filter(from %in% char$feature, to %in% char$feature)

net_edge_biopan <- net_edge_biopan %>% left_join(char[c('feature','class')], by=c('from'='feature')) %>% 
  left_join(char[c('feature','class')], by=c('to'='feature')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.x'='Abbreviation')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.y'='Abbreviation')) %>% 
  filter(FA.x==FA.y)

net_node_biopan <- net_node  %>% filter(id %in% char$feature)


visNetwork(net_node_biopan, net_edge_biopan)

#-------------------free fatty acid analysis-------------------
options(scipen = -1)
FA_exp <- exp[81:84,]
FA_exp_fig <- FA_exp %>% 
  gather(-feature, key='label',value='value') %>% 
  mutate(group=str_extract(label,'[A-Za-z]+')) %>% 
  ggboxplot(x = "feature", y = "value",
            color = "group", palette = 'jco',
            add = "jitter")+ 
  stat_compare_means(aes(group=group),method = "t.test",label.y = 440, 
                     label.x = 1.18)+
  theme_classic() +
  #scale_y_continuous(limits = c(0,25))+
  labs(x='', title='', y='FA abundance')

#-------------------essential pathway analysis for raw data-------------------


raw_species_net <- species_network_data[[2]][1:2]
colnames(raw_species_net) <- c('S1','P1')

raw_species_net <- raw_species_net%>% filter(S1%in%species_t$lipid,
                                             P1%in%species_t$lipid)

set.seed(1)
path_score_raw_species <-  calculate_path_activation_node(raw_species_net, species_t,
                                                          calibrate = T, if_FA = 'no')



rep_lipid <- path_score_raw_species%>% 
  .$path %>% 
  str_extract_all('\\d+:\\d+') %>% 
  map(.f = function(x){k=table(x) %>% sort(decreasing = T);
  if(length(k)==0){''}
  else if(sum(k!=1)==0){sort(names(k), decreasing = T)[1]}
  else if(sum(k!=1)!=0){names(k)[1]}}) %>% 
  unlist


path_score_raw_species <- path_score_raw_species%>% 
  mutate(rep_sub_path=rep_lipid)


top5_path_raw <- path_score_raw_species %>% 
  filter(Type=='Suppressed') %>% 
  arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,]


top5_path_raw_fig <- top5_path_raw %>% 
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


#write.csv(path_score_raw_species, file.path(file, 'result/path_score_rawspecies_LPCAT1.csv'), row.names = F)


#-------------------essential edges (reactions) analysis for raw data-------------------

raw_species_net_w_rev <- add_rev_rection(network_edge, raw_species_net)


#calculate perturbation score for each edge (reaction)
perturbation_score_raw_species <- calculate_perturbation_point(raw_species_net_w_rev,
                                                           column_to_rownames(lipid_exp, var = 'feature'),
                                                           species_t,
                                                           ctrl=1:5, exp=6:10,
                                                           stat = 'p')



top5_rep_path_raw <- path_score_raw_species %>% filter(Type=='Suppressed') %>% 
  arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path

top5_rep_path_raw <- path_score_raw_species %>% filter(rep_sub_path%in%top5_rep_path_raw) %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top5_path_raw <- perturbation_score_raw_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top5_rep_path_raw && x[2] %in% top5_rep_path_raw})
edge_in_top5_path_raw <- perturbation_score_raw_species[edge_in_top5_path_raw,]$edge_name



top5_edge_lands_cycle_raw <- perturbation_score_raw_species%>% 
  filter(edge_name %in%edge_in_top5_path_raw) %>% 
  filter(str_detect(edge_name, 'L.+ --> P.+|P.+ --> L.+')) %>% 
  filter(!str_detect(edge_name, 'LPA')) %>% 
  .[c((nrow(.)-4):nrow(.)),] %>% 
  filter(p_value<0.05) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color, FA_change))


top5_edge_lands_cycle_raw_fig <- top5_edge_lands_cycle_raw %>% ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
                                                                  fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(top5_edge_lands_cycle_raw$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', legend.box = "vertical",
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = pal_lancet()(2))+
  scale_color_manual(values = c('white','gold'))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 suppressed edges (Lands cycle)',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))

#write.csv(perturbation_score_raw_species, file.path(file, 'result/perturbation_score_rawspecies_LPCAT1.csv'), row.names = F)


#-------------------construct network for raw data-------------------


top5_net_edge_raw <- raw_species_net_w_rev %>% filter(S1 %in% top5_rep_path_raw, P1 %in% top5_rep_path_raw)


top5_net_edge_raw <- top5_net_edge_raw %>% mutate(color='gray', arrows='to',length=100)

colnames(top5_net_edge_raw)[1:2] <- c('from','to')


top5_net_node_raw <- species_network_data[[1]] %>% 
  filter(id %in% unlist(top5_net_edge_raw))


top5_net_path_col_raw <- path_score_raw_species %>% filter(Type=='Suppressed') %>% 
  arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% 
  mutate(rank=c(6:10))

top5_net_reaction_col_raw <- perturbation_score_raw_species%>% 
  filter(edge_name %in%edge_in_top5_path_raw) %>% 
  filter(str_detect(edge_name, 'L.+ --> P.+|P.+ --> L.+')) %>% 
  filter(!str_detect(edge_name, 'LPA')) %>% 
  .[c((nrow(.)-4):nrow(.)),] %>% 
  filter(p_value<0.05)%>% 
  mutate(rank=c(6:10))


top5_net_node_raw$label <- str_replace_all(top5_net_node_raw$label , ';0','')

net_edge_raw <- path_color(top5_net_edge_raw, top5_net_path_col_raw,top5_net_reaction_col_raw)


net_node_raw <- top5_net_node_raw %>% filter(id %in%c(net_edge_raw$from, net_edge_raw$to))


#building network

visNetwork(net_node_raw, net_edge_raw)




species_raw_t_net <- species_t
maxinf <- ceiling(max(species_raw_t_net$log2FC[is.finite(species_raw_t_net$log2FC)]))
mininf <- ceiling(min(species_raw_t_net$log2FC[is.finite(species_raw_t_net$log2FC)]))
species_raw_t_net$log2FC[species_raw_t_net$log2FC>0 & is.infinite(species_raw_t_net$log2FC)] <- maxinf
species_raw_t_net$log2FC[species_raw_t_net$log2FC<0 & is.infinite(species_raw_t_net$log2FC)] <- mininf


# path_color(top5_net_edge_raw, top5_net_path_col_raw,top5_net_reaction_col_raw) %>% 
#   left_join(species_raw_t_net, by=c('from'='lipid')) %>% 
#    left_join(species_raw_t_net, by=c('to'='lipid')) %>% 
#   .[c(1:9,17,22,23,31,36,37)] %>% 
#   write.csv(file.path(file, 'network/rawspecies_network_LPCAT1.csv'),
#               row.names = F, na = '')

