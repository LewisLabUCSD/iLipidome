library(tidyverse)
library(plyr)
library(gplots)
library(ggsci)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(ggtext)
library(igraph)
library(visNetwork)
library(xlsx)
library(gtools)
library(data.table)

dir_name <- dirname(rstudioapi::getSourceEditorContext()$path)
#load(file.path(dir_name, 'CVD.RData'))
#-------------------Data upload-------------------

load(file.path(dir_name, 'Required_data/Required_data_CVD.RData'))
load(file.path(dir_name,'Required_data/Required_function.RData'))

each_FA <- str_extract_all(char$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char <- char %>% mutate(each_FA=each_FA)


#-------------------Analysis for unprocessed data-------------------


no_sub_t <- non_processed_data_test(exp, char, 't.test', 'adj_p_value', 2:86,  87:147)


#-------------------Decompose lipids into substructure-------------------
sub_species <- lipid_species_substructure_transform(char, FA_substructure, lipid_substructure)


species_substructure <- sub_species[[2]] %>% map(.f = function(x){names(x) <- 1:length(x);return(x)}) %>% 
  plyr::ldply(rbind) %>% mutate(Lipid=sub_species[[1]]) %>% 
  dplyr::select(Lipid, everything())

species_substructure[is.na(species_substructure)] <- ''

colnames(species_substructure) <- c('Lipid', str_c('Unit', 1:(ncol(species_substructure)-1)))

#-------------------extract substructure using fold changes-------------------

species_t <- no_sub_t[[2]] %>% filter(type=='species')


species_sub_stop <- lipid_substructure_w_stop(species_substructure, species_t, 'species', 0.3)
species_sub_stop <- list(species_sub_stop[[1]], apply(species_sub_stop[-1], MARGIN = 1, FUN = function(x){x[x!='']}))

#-------------------transform lipid profile into substructure profile-------------------



species_sub_exp <- lipid_substructure_matrix(exp, species_sub_stop, 'species')

#-------------------differential expression analysis for substructures-------------------

species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:85, 86:146, 't.test', 'adj_p_value')

#-------------------substructures improve statistical power-------------------


options(scipen = 0)

stat_improve_fig <- species_sub_exp_t %>% 
  left_join(species_t, by=c('lipid')) %>% 
  .[c('lipid', 'mlog10padj.y', 'mlog10padj.x')] %>% 
  `colnames<-`(c('lipid','Lipid species','Substructure')) %>%
  filter(!is.na(`Lipid species`)) %>% 
  filter(`Lipid species`> -log10(0.05) | `Substructure`> -log10(0.05)) %>% 
  ggpaired(cond1 = "Lipid species", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = '-log10(padj)')+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed')+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.24, label.y = 7.3)+
  theme(legend.position = 'none')


#-------------------reconstruct biosynthetic network-------------------


species_network <- split(species_substructure[-1], seq(nrow(species_substructure)))
species_network <- species_network %>% map(.f = function(x){x[x!='']})

species_network <- species_network %>% map(.f = function(x){head(rep(x, each=2)[-1],-1)})
species_network <-  matrix(unlist(species_network), ncol=2, byrow = T,dimnames = list(NULL, c('S1','P1'))) %>% 
  as.data.frame() %>%  unique()

species_network <- species_network %>% mutate(S1=str_replace(S1, '_FA\\d',''),
                                              P1=str_replace(P1, '_FA\\d',''))


species_network_data <- draw_network(species_network, species_sub_exp_t, 'no','adj_p_value')

#-------------------essential pathway analysis for substructures-------------------

species_net <- species_network_data[[2]][1:2]
colnames(species_net) <- c('S1','P1')


set.seed(1)
path_score_species <-  calculate_path_activation_node(species_net, species_sub_exp_t,
                                                      calibrate = T,if_FA = 'no')


#write.csv(path_score_species, file.path(dir_name, 'result/path_score_subspecies_CVD.csv'), row.names = F)

top5_path <- rbind(path_score_species %>% 
                     filter(Type=='Active') %>%
                     arrange(desc(cal_score)) %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,])

top5_path_fig <- top5_path

top5_path_fig$path[c(1:3,5,7,8,9)] <- top5_path_fig$path[c(1:3,5,7,8,9)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], tail(x,2)[1],tail(x,2)[2],sep=' --> ')})

top5_path_fig$path[c(6,10)] <- top5_path_fig$path[c(6,10)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], x[2], tail(x,2)[1],tail(x,2)[2],sep=' --> ')})


top5_path_fig$path <- str_replace_all(top5_path_fig$path, ';0','')
top5_path_fig <- top5_path_fig %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  #scale_y_continuous(limits = c(-5,5.5))+
  scale_fill_manual(values = bluered(100)[c(100,90,80,70,60,1,10,20,30,40)])+
  scale_y_continuous(limits = c(-4,6))+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')



#-------------------essential edges (reactions) analysis for substructures-------------------
#consider all reversible reactions

species_net_w_rev <- add_rev_rection(network_edge, species_net)

#calculate pertubation score for each edge (reaction)


perturbation_score_species <- calculate_perturbation_point(species_net_w_rev,
                                                           species_sub_exp[[3]],
                                                           species_sub_exp_t,
                                                           1:85,86:146)

#write.csv(perturbation_score_species, file.path(dir_name, 'result/perturbation_score_subspecies_CVD.csv'), row.names = F)


top5_rep_path <- c(path_score_species %>% filter(Type=='Active') %>% 
                     arrange(desc(cal_score)) %>% 
                     .[!duplicated(.$rep_lipid),] %>% .[1:5,] %>% .$rep_lipid
                   ,path_score_species %>% filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_lipid),] %>% .[1:5,] %>% .$rep_lipid)

top5_rep_path <- path_score_species %>% filter(rep_lipid%in%top5_rep_path, Significant=='yes') %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top5_path <- perturbation_score_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top5_rep_path && x[2] %in% top5_rep_path})
edge_in_top5_path <- perturbation_score_species[edge_in_top5_path,]$edge_name



top5_edge <- perturbation_score_species%>% 
  filter(edge_name %in%edge_in_top5_path) %>% 
  mutate(edge_name=str_replace_all(edge_name, ';0','')) %>% 
  .[c(1:5,(nrow(.)-4):nrow(.)),] %>% 
  filter(p_value<0.05) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color, FA_change))


top5_edge_fig <- top5_edge %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c('Same as\nreaction',''))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase','Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
                                          fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(top5_edge$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  #geom_vline(xintercept = c(3.92,-3.92), linetype='dashed',color='gray')+
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

grid.arrange(top5_path_fig, top5_edge_fig, 
             nrow=1,layout_matrix=rbind(c(1,2)))


#-------------------biosynthetic network for top 5 pathways and reactions-------------------


top5_net_edge <- species_net_w_rev %>% filter(S1 %in% top5_rep_path, P1 %in% top5_rep_path)

top5_net_edge <- top5_net_edge %>% mutate(color='gray', arrows='to',length=100)

colnames(top5_net_edge)[1:2] <- c('from','to')


top5_net_node <- species_network_data[[1]] %>% 
  filter(id %in% unlist(top5_net_edge))

top5_net_path_col <- rbind(path_score_species %>% filter(Type=='Active') %>% 
                             .[!duplicated(.$rep_lipid),] %>% .[1:5,] %>% 
                             mutate(rank=c(1:5)),
                           path_score_species %>% filter(Type=='Suppressed') %>% 
                             arrange(cal_score) %>% .[!duplicated(.$rep_lipid),] %>% .[1:5,] %>% 
                             mutate(rank=c(6:10)))


top5_net_reaction_col <- perturbation_score_species%>% 
  filter(edge_name %in%edge_in_top5_path) %>% 
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  filter(p_value<0.05)%>% 
  mutate(rank=c(1:10))


#building network

top5_net_node$label <- str_replace_all(top5_net_node$label, ';0','')

net_edge <- path_color(top5_net_edge, top5_net_path_col,top5_net_reaction_col)


net_node <- top5_net_node %>% filter(id %in%c(net_edge$from, net_edge$to))

visNetwork(net_node, net_edge)


#Cytoscape data
#net_edge %>% 
#  left_join(species_sub_exp_t[c('lipid', 'log2FC','mlog10padj','sig')], by=c('from'='lipid')) %>% 
#  left_join(species_sub_exp_t[c('lipid', 'log2FC','mlog10padj','sig')], by=c('to'='lipid')) %>% 
#  mutate(from=str_replace_all(from, ';0','')) %>% 
#  mutate(to=str_replace_all(to, ';0','')) %>% 
#  write.csv(file.path(dir_name, 'network/subspecies_network_CVD.csv'), row.names = F)


#BIOPAN's algorithm cannot build the connections involving FA transfer if no free FA data provided
net_edge_biopan <- net_edge %>% filter(from %in% char$feature, to %in% char$feature)

net_edge_biopan <- net_edge_biopan %>% left_join(char[c('feature','class')], by=c('from'='feature')) %>% 
  left_join(char[c('feature','class')], by=c('to'='feature')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.x'='Abbreviation')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.y'='Abbreviation')) %>% 
  filter(FA.x==FA.y)

net_node_biopan <- net_node  %>% filter(id %in% char$feature)


visNetwork(net_node_biopan, net_edge_biopan)


#-------------------biomarker analysis-------------------
stat_data <- species_sub_exp_t %>% 
  left_join(species_t, by=c('lipid')) %>% 
  .[c('lipid', 'mlog10p.y', 'mlog10p.x','log2FC.y','log2FC.x')] %>% 
  `colnames<-`(c('lipid','Lipid species','Substructure','log2FC_raw','log2FC_sub')) %>% 
  filter(!is.na(`Lipid species`)) 

top5_sig_fig <- rbind(arrange(stat_data, desc(Substructure))[1:5,],
                 arrange(stat_data, desc(`Lipid species`))[1:5,]) %>% unique() %>% 
  mutate(Significance='Decrease') %>% 
  mutate(Improve=ifelse(Substructure>`Lipid species`, 'yes','no')) %>% 
  mutate(lipid=str_replace_all(lipid, ';0','')) %>% 
  ggplot(aes(x=`Lipid species`, y=`Substructure`, fill=Significance, color=Improve))+
  geom_point(pch=21,stroke = 1.5, size=3)+
  geom_abline(slope = 1, linetype='dashed', color='gray')+
  scale_color_manual(values = c('white','orange'))+
  scale_fill_manual(values = c('blue'))+
  theme_bw()+
  scale_x_continuous(limits = c(5,10))+
  scale_y_continuous(limits = c(5,10))+
  geom_text_repel(aes(label=lipid), color='black')+
  labs(x='Lipid species: -log10 (p-value)', 
       y='Substructure: -log10 (p-value)',
       color='Increased\npower',
       title = 'Top 5 significant lipid species and substructures')+
  theme(plot.title = element_text(hjust = 0.5))



PC_18_2 <- which(rownames(species_sub_exp[[1]])=='PC_18:2;0_20:2;0')

PC_18_2_sub <- data.frame(group=names(species_sub_exp[[3]][PC_18_2,]),
                value=species_sub_exp[[3]][PC_18_2,]) %>% 
  mutate(group=str_extract(group, '[A-z]+')) %>% 
  mutate(group=ifelse(group=='Ctrl','Control',group)) %>% 
  ggboxplot(x = "group", y = "value",
            color = "group", palette = c('#999999','#E69F00'),
            add = "jitter")+ 
  stat_compare_means(method = "t.test",label.y = 8, 
                     label.x = 1.23, method.args = list(var.equal=T))+
  theme_classic() +
  #scale_y_continuous(limits = c(0,3.2))+
  labs(x='', title="PC_18:2_20:2' = PC_18:2_20:2 + LPC_18:2", y="PC_18:2_20:2'")+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(size=12))

PC_18_2_raw <- no_sub_t[[1]] %>% filter(feature=='PC_18:2;0_20:2;0') %>% .[-c(1:2)] %>% 
  gather(key='group', value='value') %>% 
  mutate(group=str_extract(group, '[A-z]+')) %>% 
  mutate(group=ifelse(group=='Ctrl','Control',group)) %>% 
  ggboxplot(x = "group", y = "value",
            color = "group", palette = c('#999999','#E69F00'),
            add = "jitter")+ 
  stat_compare_means(method = "t.test",label.y = 6, 
                     label.x = 1.23, method.args = list(var.equal=T))+
  theme_classic() +
  #scale_y_continuous(limits = c(0,3.2))+
  labs(x='', title="PC_18:2_20:2", y="PC_18:2_20:2")+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(size=12))


PE_18_2 <- which(rownames(species_sub_exp[[1]])=='PE_18:2;0_18:2;0')

PE_18_2_sub <- data.frame(group=names(species_sub_exp[[3]][PE_18_2,]),
                value=species_sub_exp[[3]][PE_18_2,]) %>% 
  mutate(group=str_extract(group, '[A-z]+')) %>% 
  mutate(group=ifelse(group=='Ctrl','Control',group)) %>% 
  ggboxplot(x = "group", y = "value",
            color = "group", palette = c('#999999','#E69F00'),
            add = "jitter")+ 
  stat_compare_means(method = "t.test",label.y = 43, 
                     label.x = 1.23, method.args = list(var.equal=T))+
  theme_classic() +
  #scale_y_continuous(limits = c(0,3.2))+
  labs(x='', title="PE_18:2_18:2' = PC_18:2_18:2 +\nPE_18:2_18:2 + LPC_18:2 + LPE_18:2", y="PE_18:2_18:2'")+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(size=12))

PE_18_2_raw <- no_sub_t[[1]] %>% filter(feature=='PE_18:2;0_18:2;0') %>% .[-c(1:2)] %>% 
  gather(key='group', value='value') %>% 
  mutate(group=str_extract(group, '[A-z]+')) %>% 
  mutate(group=ifelse(group=='Ctrl','Control',group)) %>% 
  ggboxplot(x = "group", y = "value",
            color = "group", palette = c('#999999','#E69F00'),
            add = "jitter")+ 
  stat_compare_means(method = "t.test",label.y = 2, 
                     label.x = 1.23, method.args = list(var.equal=T))+
  theme_classic() +
  #scale_y_continuous(limits = c(0,3.2))+
  labs(x='', title="PE_18:2_18:2", y="PE_18:2_18:2")+
  theme(legend.position = 'none',plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(size=12))

options(scipen = 0)
grid.arrange(PC_18_2_raw, PC_18_2_sub, PE_18_2_raw, PE_18_2_sub, 
             nrow=2,layout_matrix=rbind(c(1,2),
                                        c(3,4)))
#-------------------essential pathway analysis for raw data-------------------


raw_species_net <- species_network_data[[2]][1:2]
colnames(raw_species_net) <- c('S1','P1')

raw_species_net <- raw_species_net%>% filter(S1%in%species_t$lipid,
                                             P1%in%species_t$lipid)

set.seed(2)
path_score_raw_species <-  calculate_path_activation_node(raw_species_net, species_t,
                                                          calibrate = T, if_FA = 'no')



top5_path_raw <- rbind(path_score_raw_species %>% 
                         filter(Type=='Suppressed') %>% 
                         arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                       path_score_raw_species %>% 
                         filter(Type=='Active') %>% 
                         .[!duplicated(.$rep_sub_path),]%>% .[1:5,])

top5_path_raw_fig <- top5_path_raw

top5_path_raw_fig$path[c(1,2,4)] <- top5_path_raw_fig$path[c(1,2,4)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1],tail(x,2)[1],tail(x,2)[2],sep=' --> ')})



top5_path_raw_fig <- top5_path_raw_fig %>% 
  mutate(path=str_replace_all(path, ';0','')) %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  scale_fill_manual(values = c(bluered(100)[c(1,10,20,30,40)],bluered(100)[c(100,90,80,70)], 'gray'))+
  scale_y_continuous(limits = c(-4,4))+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')

#write.csv(path_score_raw_species, file.path(dir_name, 'result/path_score_rawspecies_CVD.csv'), row.names = F)


#-------------------essential edges (reactions) analysis for raw data-------------------

raw_species_net_w_rev <- add_rev_rection(network_edge, raw_species_net)


#calculate perturbation score for each edge (reaction)
perturbation_score_raw_species <- calculate_perturbation_point(raw_species_net_w_rev,
                                                               column_to_rownames(exp, var = 'feature'),
                                                               species_t,
                                                               ctrl=1:85, exp=86:146,
                                                               stat = 'p')



top5_rep_path_raw <- c(path_score_raw_species %>% filter(Type=='Suppressed') %>% 
                         arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path,
                       path_score_raw_species %>% filter(Type=='Active')  %>% 
                         .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path)


top5_rep_path_raw <- path_score_raw_species %>% filter(rep_sub_path%in%top5_rep_path_raw, Significant=='yes') %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()


edge_in_top5_path_raw <- perturbation_score_raw_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top5_rep_path_raw && x[2] %in% top5_rep_path_raw})
edge_in_top5_path_raw <- perturbation_score_raw_species[edge_in_top5_path_raw,]$edge_name


top5_edge_raw <- perturbation_score_raw_species%>% 
  filter(edge_name %in%edge_in_top5_path_raw) %>% 
  filter(p_value<0.05) %>%
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color, FA_change))


top5_edge_raw_fig <- top5_edge_raw %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase','Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(top5_edge_raw$edge_color)
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

grid.arrange(top5_path_raw_fig,top5_path_fig,
             layout_matrix=rbind(c(1,1,1,1,1,1,2,2,2,2,2,2)))

grid.arrange(top5_edge_raw_fig,top5_edge_fig,
             layout_matrix=rbind(c(1,1,1,1,1,1,1,2,2,2,2,2,2,2)))
#write.csv(perturbation_score_raw_species, file.path(dir_name, 'result/perturbation_score_rawspecies_CVD.csv'), row.names = F)


#-------------------construct network for raw species data-------------------


top5_rep_path_raw <- c(path_score_raw_species %>% filter(Type=='Suppressed') %>% 
                         arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path,
                       path_score_raw_species %>% filter(Type=='Active')  %>% 
                         .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path)



top5_rep_path_raw <- path_score_raw_species %>% filter(rep_sub_path%in%top5_rep_path_raw) %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top5_path_raw <- perturbation_score_raw_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top5_rep_path_raw && x[2] %in% top5_rep_path_raw})
edge_in_top5_path_raw <- perturbation_score_raw_species[edge_in_top5_path_raw,]$edge_name



top5_net_edge_raw <- raw_species_net_w_rev %>% filter(S1 %in% top5_rep_path_raw, P1 %in% top5_rep_path_raw)

top5_net_edge_raw <- top5_net_edge_raw %>% mutate(color='gray', arrows='to',length=100)

colnames(top5_net_edge_raw)[1:2] <- c('from','to')


top5_net_node_raw <- species_network_data[[1]] %>% 
  filter(id %in% unlist(top5_net_edge_raw))

top5_net_path_col_raw <- rbind(path_score_raw_species %>% filter(Type=='Active') %>% 
                                 .[!duplicated(.$rep_sub_path),] %>% .[1:5,], 
                               path_score_raw_species %>% filter(Type=='Suppressed') %>% 
                                 arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% 
                                 .[1:5,]) %>% mutate(rank=c(1:10))

top5_net_reaction_col_raw <- perturbation_score_raw_species%>% 
  filter(edge_name %in% edge_in_top5_path_raw) %>% 
  filter(p_value<0.05) %>% 
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  filter(p_value<0.05)%>% 
  mutate(rank=c(1:10))


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
#    left_join(species_raw_t_net, by=c('from'='lipid')) %>% 
#    left_join(species_raw_t_net, by=c('to'='lipid')) %>% 
#     .[c(1:9,17,22,23,31,36,37)] %>% 
#     write.csv(file.path(dir_name, 'network/rawspecies_network_CVD.csv'),
#               row.names = F, na = '')



#-------------------Explore reactions related to ferropotsis-----------------------
Ferr_reaction_sub_fig <- perturbation_score_species %>% 
  mutate(`Ferroptosis-related reactions`=ifelse(FA_change==' (20:4)' & str_detect(edge_name, 'L.+ --> P.+'), 'yes','no')) %>% 
  mutate(label=str_replace_all(edge_name,';0','')) %>% 
  mutate(label=str_replace_all(label,' ','')) %>% 
  mutate(label=ifelse(`Ferroptosis-related reactions`=='yes',label,'')) %>% 
  ggplot(aes(x=perturbation_score, y=mlog10p))+
  geom_point(aes(color=`Ferroptosis-related reactions`, alpha=`Ferroptosis-related reactions`))+
  geom_text_repel(aes(label=label), max.overlaps = 10000,
                  size  = 3)+
  scale_color_manual(values = c('gray','red'))+
  scale_alpha_manual(values = c(0.3,1))+
  scale_y_continuous(limits = c(0,12.5))+
  scale_x_continuous(limits = c(-20,25))+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed', )+
  labs(y='-log10 (p-value)', x='Perturbation score')+
  theme_bw()+
  theme(legend.position = 'top')


Ferr_reaction_raw_fig <- perturbation_score_raw_species %>% 
  mutate(`Ferroptosis-related reactions`=ifelse(FA_change==' (20:4)' & str_detect(edge_name, 'L.+ --> P.+'), 'yes','no')) %>% 
  mutate(label=str_replace_all(edge_name,';0','')) %>% 
  mutate(label=str_replace_all(label,' ','')) %>% 
  mutate(label=ifelse(`Ferroptosis-related reactions`=='yes',label,'')) %>% 
  ggplot(aes(x=perturbation_score, y=mlog10p))+
  geom_point(aes(color=`Ferroptosis-related reactions`, alpha=`Ferroptosis-related reactions`))+
  geom_text_repel(aes(label=label), max.overlaps = 10000,
                  size  = 3)+
  scale_y_continuous(limits = c(0,12.5))+
  scale_x_continuous(limits = c(-20,25))+
  scale_color_manual(values = c('gray','red'))+
  scale_alpha_manual(values = c(0.3,1))+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed', )+
  labs(y='-log10 (p-value)', x='Perturbation score')+
  theme_bw()+
  theme(legend.position = 'top')


grid.arrange(Ferr_reaction_raw_fig,Ferr_reaction_sub_fig,
             layout_matrix=rbind(c(1,1,1,1,1,1,1,2,2,2,2,2,2,2)))
