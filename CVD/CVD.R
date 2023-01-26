library(ggsci)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(ggtext)

file <- dirname(rstudioapi::getSourceEditorContext()$path)

#-------------------Data upload-------------------

source(file.path(dirname(file),'Required_function/required_function.R'))
load(file.path(dirname(file),'Required_data/required_data.RData'))


exp <- read.csv(file.path(file,'lipidome_data/exp_CVD_raw.csv'))
char <- read.csv(file.path(file,'lipidome_data/char_CVD.csv'))

each_FA <- str_extract_all(char$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char <- char %>% mutate(each_FA=each_FA)


#-------------------Analysis for unprocessed data-------------------

no_sub_t <- unprocessed_data_test(exp, char, 't.test', 'adj_p_value', 1:85,  86:146)

#-------------------Lipid species substructure analysis-------------------

#Decompose lipids into substructures

species_substructure <- species_sub_transform(char, lipid_substructure, 
                                              network_node)

#Extract species substructures using fold change

species_t <- no_sub_t[[2]] %>% filter(type=='species')

species_sub_stop <- species_sub_extract(species_substructure, species_t, 'species', pct_limit = 0.3)


#Transform lipid exp into substructure exp

species_sub_exp <- lipid_sub_matrix(exp, species_sub_stop, 'species')


#Differential expression analysis for substructures

species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:85, 86:146, 't.test', 'adj_p_value')

#Species biosynthetic network data transformation

species_network <- build_species_net(species_substructure)

#Essential pathway analysis for species substructures


set.seed(1)

path_score_species <-  path_scoring(species_network, species_sub_exp_t,
                                    calibrate = T, data_type = 'species')



SF4c_data <- rbind(path_score_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_species %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
  mutate(path=str_replace_all(path, ';0',''))



SF4c_data$path[c(2:4,6,7,8,10)] <- SF4c_data$path[c(2:4,6,7,8,10)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], tail(x,2)[1],tail(x,2)[2],sep=' --> ')})

SF4c_data$path[c(1,5)] <- SF4c_data$path[c(1,5)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], x[2], tail(x,2)[1],tail(x,2)[2],sep=' --> ')})


#write.csv(SF4c_data, file.path(file, 'source_data/SF4c.csv'), row.names = F)

SF4c <- SF4c_data %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  #scale_y_continuous(limits = c(-5,5.5))+
  scale_fill_manual(values = bluered(100)[c(1,10,20,30,40,100,90,80,70,60)])+
  scale_y_continuous(limits = c(-4,6))+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')


#Essential edges (reactions) analysis for species substructures

species_net_w_rev <- add_rev_rection(network_edge, species_network)


reaction_score_species <- reaction_scoring(species_net_w_rev,
                                           species_sub_exp[[3]],
                                           species_sub_exp_t,
                                           ctrl=1:85, exp=86:146,
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



SF4e_data <- reaction_score_species%>% 
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

#write.csv(SF4e_data[-c(20:22)], file.path(file, 'source_data/SF4e.csv'), row.names = F)

SF4e <- SF4e_data %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase','Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(SF4e_data$edge_color)
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


#Construct species biosynthetic network


species_network_data <- draw_network(species_net_w_rev, species_sub_exp_t,
                                     if_species = T,significant = 'adj_p_value',
                                     path_scoring_result = path_score_species,
                                     reaction_scoring_result = reaction_score_species,
                                     top_n = 5, path_type = 'both')


F6b <- visNetwork(species_network_data[[1]], species_network_data[[2]])

species_sub_t_net <- species_sub_exp_t

maxinf <- ceiling(max(species_sub_t_net$log2FC[is.finite(species_sub_t_net$log2FC)]))
mininf <- floor(min(species_sub_t_net$log2FC[is.finite(species_sub_t_net$log2FC)]))
species_sub_t_net$log2FC[species_sub_t_net$log2FC>0 & is.infinite(species_sub_t_net$log2FC)] <- maxinf
species_sub_t_net$log2FC[species_sub_t_net$log2FC<0 & is.infinite(species_sub_t_net$log2FC)] <- mininf


F6b_data <- species_network_data[[2]]%>% left_join(species_sub_t_net, by=c('from'='lipid')) %>% 
  left_join(species_sub_t_net, by=c('to'='lipid')) %>% 
  mutate(from=str_replace_all(from, ';0','')) %>% 
  mutate(to=str_replace_all(to, ';0','')) %>% 
  .[c(1:9,16,21,22,29,34,35)]

#write.xlsx(F6b_data, file.path(file, 'source_data/F6b.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)



#BIOPAN's algorithm cannot build the connections involving FA transfer without free FA data provided
biopan_network_edge <- species_network_data[[2]] %>% filter(from %in% char$feature, to %in% char$feature)

biopan_network_edge <- biopan_network_edge %>% left_join(char[c('feature','class')], by=c('from'='feature')) %>% 
  left_join(char[c('feature','class')], by=c('to'='feature')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.x'='Abbreviation')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.y'='Abbreviation')) %>% 
  filter(FA.x==FA.y)

biopan_network_node <- species_network_data[[1]]  %>% filter(id %in% char$feature)

SF5c <- visNetwork(biopan_network_node, biopan_network_edge)

SF5c_data1 <- biopan_network_node %>% 
  left_join(species_sub_t_net, by=c('id'='lipid')) %>% 
  mutate(id=str_replace_all(id, ';0','')) %>% 
  mutate(label=str_replace_all(label, ';0','')) %>% 
  .[c(1:7,14,19,20)]

SF5c_data2 <- biopan_network_edge[1:9] %>% left_join(species_sub_t_net, by=c('from'='lipid')) %>% 
  left_join(species_sub_t_net, by=c('to'='lipid')) %>% 
  mutate(from=str_replace_all(from, ';0','')) %>% 
  mutate(to=str_replace_all(to, ';0','')) %>% 
  .[c(1:9,16,21,22,29,34,35)]

#write.xlsx(SF5c_data1, file.path(file, 'source_data/SF5c1.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)
#write.xlsx(SF5c_data2, file.path(file, 'source_data/SF5c2.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)

#-------------------Raw lipid species analysis-------------------

#Essential pathway analysis for raw species data

raw_species_net <- species_network

raw_species_net <- raw_species_net%>% filter(S1%in%species_t$lipid,
                                             P1%in%species_t$lipid)


set.seed(1)
path_score_raw_species <-  path_scoring(raw_species_net, species_t,
                                        calibrate = T, data_type = 'species')





SF4b_data <- rbind(path_score_raw_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_raw_species %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
  mutate(path=str_replace_all(path, ';0',''))

SF4b_data$path[c(1,2,4)] <- SF4b_data$path[c(1,2,4)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1],tail(x,2)[1],tail(x,2)[2],sep=' --> ')})



#write.csv(SF4b_data, file.path(file, 'source_data/SF4b.csv'), row.names = F)

SF4b <- SF4b_data %>% 
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


#Essential edges (reactions) analysis for raw species data


raw_species_net_w_rev <- add_rev_rection(network_edge, raw_species_net)



reaction_score_raw_species <- reaction_scoring(raw_species_net_w_rev,
                                               column_to_rownames(exp, var = 'feature'),
                                               species_t,
                                               ctrl=1:85, exp=86:146,
                                               Species = 'human')


top5_rep_path_raw <- c(path_score_raw_species %>% filter(Type=='Active') %>% 
                         arrange(desc(cal_score)) %>% 
                         .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path
                       ,path_score_raw_species %>% filter(Type=='Suppressed') %>% 
                         arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path)


top5_rep_path_raw <- path_score_raw_species %>% filter(rep_sub_path%in%top5_rep_path_raw, Significant=='yes') %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top5_path_raw <- reaction_score_raw_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top5_rep_path_raw && x[2] %in% top5_rep_path_raw})
edge_in_top5_path_raw <- reaction_score_raw_species[edge_in_top5_path_raw,]$edge_name



SF4d_data <- reaction_score_raw_species%>% 
  filter(edge_name %in%edge_in_top5_path_raw) %>% 
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

#write.csv(SF4d_data[-c(20:22)], file.path(file, 'source_data/SF4d.csv'), row.names = F)

SF4d <- SF4d_data %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase','Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(SF4d_data$edge_color)
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


#-------------------Species substructures improve statistical power-------------------



SF4a_data <- species_sub_exp_t %>% 
  left_join(species_t, by=c('lipid')) %>% 
  .[c('lipid', 'mlog10padj.y', 'mlog10padj.x')] %>% 
  `colnames<-`(c('lipid','Lipid species','Substructure')) %>% 
  filter(`Lipid species`>-log10(0.05)| Substructure>-log10(0.05)) 

#write.csv(SF4a_data, file.path(file, 'source_data/SF4a.csv'), row.names = F, na = '')

options(scipen = 0)

SF4a <- SF4a_data %>% 
  filter(!is.na(`Lipid species`)) %>% 
  ggpaired(cond1 = "Lipid species", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = '-log10(padj)')+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed')+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.24, label.y = 7.3)+
  theme(legend.position = 'none')



#-------------------Biomarker analysis-------------------
stat_data <- species_sub_exp_t %>% 
  left_join(species_t, by=c('lipid')) %>% 
  .[c('lipid', 'mlog10p.y', 'mlog10p.x','log2FC.y','log2FC.x')] %>% 
  `colnames<-`(c('lipid','Lipid species','Substructure','log2FC_raw','log2FC_sub')) %>% 
  filter(!is.na(`Lipid species`))

F6c_data <- rbind(arrange(stat_data, desc(Substructure))[1:5,],
      arrange(stat_data, desc(`Lipid species`))[1:5,]) %>% unique() %>% 
  mutate(Significance='Decrease') %>% 
  mutate(Improve=ifelse(Substructure>`Lipid species`, 'yes','no')) %>% 
  mutate(lipid=str_replace_all(lipid, ';0','')) 

#write.csv(F6c_data, file.path(file, 'source_data/F6c.csv'), row.names = F, na='')


F6c <- F6c_data %>% 
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


F6d1_data <- which(rownames(species_sub_exp[[1]])=='PC_18:2;0_20:2;0')

F6d1_data <- data.frame(group=names(species_sub_exp[[3]][F6d1_data,]),
           value=species_sub_exp[[3]][F6d1_data,]) %>% 
  mutate(group=str_extract(group, '[A-z]+')) %>% 
  mutate(group=ifelse(group=='Ctrl','Control',group)) %>% 
  mutate(Data='Substructure', feature='PC_18:2;0_20:2;0')

#write.csv(F6d1_data, file.path(file, 'source_data/F6d1.csv'), row.names = F, na='')

F6d1 <- F6d1_data %>% 
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


F6d2_data <- no_sub_t[[1]] %>% filter(feature=='PC_18:2;0_20:2;0') %>% .[-c(1:2)] %>% 
  gather(key='group', value='value') %>% 
  mutate(group=str_extract(group, '[A-z]+')) %>% 
  mutate(group=ifelse(group=='Ctrl','Control',group)) %>% 
  mutate(Data='Lipid species', feature='PC_18:2;0_20:2;0')


#write.csv(F6d2_data, file.path(file, 'source_data/F6d2.csv'), row.names = F, na='')

F6d2 <- F6d2_data %>% 
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



F6d3_data <- which(rownames(species_sub_exp[[1]])=='PE_18:2;0_18:2;0')

F6d3_data <- data.frame(group=names(species_sub_exp[[3]][F6d3_data,]),
                        value=species_sub_exp[[3]][F6d3_data,]) %>% 
  mutate(group=str_extract(group, '[A-z]+')) %>% 
  mutate(group=ifelse(group=='Ctrl','Control',group)) %>% 
  mutate(Data='Substructure', feature='PC_18:2;0_18:2;0')

#write.csv(F6d3_data, file.path(file, 'source_data/F6d3.csv'), row.names = F, na='')



F6d3 <- F6d3_data %>% 
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

F6d4_data <- no_sub_t[[1]] %>% filter(feature=='PE_18:2;0_18:2;0') %>% .[-c(1:2)] %>% 
  gather(key='group', value='value') %>% 
  mutate(group=str_extract(group, '[A-z]+')) %>% 
  mutate(group=ifelse(group=='Ctrl','Control',group)) %>% 
  mutate(Data='Lipid species', feature='PE_18:2;0_18:2;0')


#write.csv(F6d4_data, file.path(file, 'source_data/F6d4.csv'), row.names = F, na='')


F6d4 <- F6d4_data %>% 
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
grid.arrange(F6d2, F6d1, F6d4, F6d3, 
             nrow=2,layout_matrix=rbind(c(1,2),
                                        c(3,4)))
#-------------------Explore reactions related to ferropotsis-----------------------
SF4g_data <- reaction_score_species %>% 
  mutate(`Ferroptosis-related reactions`=ifelse(FA_change==' (20:4)' & str_detect(edge_name, 'L.+ --> P.+'), 'yes','no')) %>% 
  mutate(label=str_replace_all(edge_name,';0','')) %>% 
  mutate(label=str_replace_all(label,' ','')) %>% 
  mutate(label=ifelse(`Ferroptosis-related reactions`=='yes',label,''))

#write.csv(filter(SF4g_data, label!=''), file.path(file, 'source_data/SF4g.csv'), row.names = F)

SF4g <- SF4g_data %>% 
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


SF4f_data <- reaction_score_raw_species %>% 
  mutate(`Ferroptosis-related reactions`=ifelse(FA_change==' (20:4)' & str_detect(edge_name, 'L.+ --> P.+'), 'yes','no')) %>% 
  mutate(label=str_replace_all(edge_name,';0','')) %>% 
  mutate(label=str_replace_all(label,' ','')) %>% 
  mutate(label=ifelse(`Ferroptosis-related reactions`=='yes',label,''))

#write.csv(filter(SF4f_data, label!=''), file.path(file, 'source_data/SF4f.csv'), row.names = F)

SF4f <- SF4f_data %>% ggplot(aes(x=perturbation_score, y=mlog10p))+
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
#-------------------Save data-------------------

#save.image(file.path(file, 'CVD.RData'))
