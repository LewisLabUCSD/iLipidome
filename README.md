# iLipidome

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Quick Example](#quick-example)
- [License](#license)

# Overview
Here, we present ``iLipidome``, a novel method based on R (v 4.1.0) platform for analyzing lipidomics data in the context of the lipid biosynthetic network, thus accounting for the interdependence of measured lipids. iLipidome includes a series of functions for lipid substructure transformation, extraction, lipid biosynthetic network reconstruction, and essential pathway and reaction scoring. With three demo datasets, we demonstrated that iLipidome could enhance statistical power, enable reliable clustering and lipid enrichment analysis, and link lipidomic changes to their genetic origins. To summarize, iLipidome facilitates systems-level comparison of lipid profiles with lipid biosynthetic information efficiently and provides a deeper insight into complex lipidomic alterations across samples.


# System Requirements
## Hardware requirements
To run demo datasets with iLipidome, it requires only a standard computer with enough RAM and installed R software version over 4.0.0.

## Software requirements
### OS Requirements
The functions and demo datasets has been tested on the following systems:
+ macOS: Monterey (12.4)


### R Dependencies
The version information about R, the OS and attached or loaded packages for `iLipidome` are listed below.

```
R version 4.1.0 (2021-05-18)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 12.4

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] zh_TW.UTF-8/zh_TW.UTF-8/zh_TW.UTF-8/C/zh_TW.UTF-8/zh_TW.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] cowplot_1.1.1         fpc_2.2-9             ComplexHeatmap_2.11.1
 [4] data.table_1.14.0     gtools_3.9.2          xlsx_0.6.5           
 [7] visNetwork_2.1.0      igraph_1.2.6          ggtext_0.1.1         
[10] ggrepel_0.9.1         gridExtra_2.3         ggvenn_0.1.9         
[13] ggpubr_0.4.0          ggsci_2.9             gplots_3.1.1         
[16] plyr_1.8.6            forcats_0.5.1         stringr_1.4.0        
[19] dplyr_1.0.7           purrr_0.3.4           readr_1.4.0          
[22] tidyr_1.1.3           tibble_3.1.2          ggplot2_3.3.5        
[25] tidyverse_1.3.1      

loaded via a namespace (and not attached):
 [1] colorspace_2.0-2    ggsignif_0.6.3      rjson_0.2.20        ellipsis_0.3.2     
 [5] class_7.3-19        modeltools_0.2-23   mclust_5.4.9        circlize_0.4.13    
 [9] GlobalOptions_0.1.2 fs_1.5.0            gridtext_0.1.4      clue_0.3-60        
[13] rstudioapi_0.13     flexmix_2.3-17      fansi_0.5.0         lubridate_1.7.10   
[17] xml2_1.3.3          codetools_0.2-18    doParallel_1.0.16   robustbase_0.93-9  
[21] jsonlite_1.7.2      rJava_1.0-6         broom_0.7.8         cluster_2.1.2      
[25] kernlab_0.9-29      dbplyr_2.1.1        png_0.1-7           compiler_4.1.0     
[29] httr_1.4.2          backports_1.2.1     assertthat_0.2.1    cli_3.1.0          
[33] htmltools_0.5.1.1   tools_4.1.0         gtable_0.3.0        glue_1.6.0         
[37] Rcpp_1.0.8.3        carData_3.0-4       cellranger_1.1.0    vctrs_0.3.8        
[41] iterators_1.0.13    xlsxjars_0.6.1      rvest_1.0.3         lifecycle_1.0.0    
[45] rstatix_0.7.0       DEoptimR_1.0-10     MASS_7.3-54         scales_1.1.1       
[49] hms_1.1.0           parallel_4.1.0      RColorBrewer_1.1-2  stringi_1.6.2      
[53] S4Vectors_0.32.3    foreach_1.5.1       caTools_1.18.2      BiocGenerics_0.40.0
[57] shape_1.4.6         rlang_1.0.4         pkgconfig_2.0.3     prabclus_2.3-2     
[61] matrixStats_0.61.0  bitops_1.0-7        lattice_0.20-44     htmlwidgets_1.5.4  
[65] tidyselect_1.1.1    magrittr_2.0.1      R6_2.5.0            IRanges_2.28.0     
[69] generics_0.1.0      DBI_1.1.1           pillar_1.6.1        haven_2.4.1        
[73] withr_2.4.2         abind_1.4-5         nnet_7.3-16         modelr_0.1.8       
[77] crayon_1.4.1        car_3.0-12          KernSmooth_2.23-20  utf8_1.2.1         
[81] GetoptLong_1.0.5    readxl_1.3.1        reprex_2.0.0        digest_0.6.27      
[85] diptest_0.76-0      stats4_4.1.0        munsell_0.5.0     
```

# Installation Guide
Users only need to install R software (https://www.r-project.org/) and required packages (This process usually takes within 30 mins).

# Summary for source codes and demo datasets
Three lipidomics datasets with their analysis codes and required data were deposited in the corresponding folders. Here, we only provided the processed lipidomics data. Detailed data preprocessing approach can be found in the manuscript. The code (.R file) for each dataset has a comprehensive analysis pipeline with step-by-step comments explaining each process. All figures and tables in the manuscipt can also be generated in the codes.
- DHA (docosahexaenoic acid) dataset:
  - It contains 13 membrane lipid profiles from 7 controls and 6 DHA treated rat basophilic leukemia (RBL) cells and covers 444 individual lipid species across 16 lipid classes.
- LPCAT1 (lysophosphatidylcholine acyltransferase 1) knockout dataset:
  - The second dataset contains 10 lipid profiles of human glioblastoma cell line (U87EGFRvIII), expressing LPCAT1 shRNA or a non-targeting control. There are 82 lipid species across 6 lipid classes and 4 free FAs (16:0, 18:1, 18:0, and 20:4) in this dataset.
- CVD (cardiovascular disease) dataset:
  - In this dataset, 146 plasma lipid profiles including 85 controls and 61 potential CVD patients were collected. 


# Quick Example
Typically, each example dataset can be finished in 10 minutes. Here is a quick run using the DHA dataset for lipid species substructure analysis. Users are able to reproduce the Fig. 4f in the manuscript.
```diff
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
library(ComplexHeatmap)
library(fpc)
library(cowplot)

# Input file path for DHA dataset

file <- file_path


# Data upload

load(file.path(file,'Required_data/Required_function.RData'))
load(file.path(file,'Required_data/Required_data_DHA.Rdata'))

each_FA <- str_extract_all(char$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char <- char %>% mutate(each_FA=each_FA)


# Analysis for unprocessed data

no_sub_t <- non_processed_data_test(exp, char, 't.test', 'adj_p_value', 2:8,  9:14)


# Decompose lipids into species substructures
# We excluded ether lipids since we doidn't know they are Alkyl (O-) or Alkenyl- (P-) linked so that we cannot map them to the pathways

char_sel <- char[!str_detect(char$feature, 'O-'),]

exp_sel <- exp[!str_detect(exp$feature, 'O-'),]

sub_species <- lipid_species_substructure_transform(char_sel, FA_substructure, lipid_substructure)

species_substructure <- sub_species[[2]] %>% map(.f = function(x){names(x) <- 1:length(x);return(x)}) %>% 
  plyr::ldply(rbind) %>% mutate(Lipid=sub_species[[1]]) %>% 
  dplyr::select(Lipid, everything())

species_substructure[is.na(species_substructure)] <- ''

colnames(species_substructure) <- c('Lipid', str_c('Unit', 1:(ncol(species_substructure)-1)))


# extract species substructures using fold changes accross conditons

species_t <- no_sub_t[[2]] %>% filter(type=='species')

species_sub_stop <- lipid_substructure_w_stop(species_substructure, species_t, 'species', pct_limit = 0.3)
species_sub_stop <- list(species_sub_stop[[1]], apply(species_sub_stop[-1], MARGIN = 1, FUN = function(x){x[x!='']}))


# transform lipid exp data into species substructures exp data

species_sub_exp <- lipid_substructure_matrix(exp_sel, species_sub_stop, 'species')


# differential expression analysis for substructures

species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:7, 8:13, 't.test', 'adj_p_value')


# species biosynthetic network data transformation

species_network <- split(species_substructure[-1], seq(nrow(species_substructure)))
species_network <- species_network %>% map(.f = function(x){x[x!='']})
species_network <- species_network %>% map(.f = function(x){head(rep(x, each=2)[-1],-1)})
species_network <-  matrix(unlist(species_network), ncol=2, byrow = T,dimnames = list(NULL, c('S1','P1'))) %>% 
  as.data.frame() %>%  unique()
species_network <- species_network %>% mutate(S1=str_replace(S1, '_FA\\d',''),
                                              P1=str_replace(P1, '_FA\\d',''))

species_network_data <- draw_network(species_network, species_sub_exp_t, 'no','adj_p_value')


# essential pathway analysis for species substructures

species_net <- species_network_data[[2]][1:2]

colnames(species_net) <- c('S1','P1')

set.seed(1)

path_score_species <-  calculate_path_activation_node(species_net, species_sub_exp_t,
                                                      calibrate = T, if_FA = 'no')

top5_path <- rbind(path_score_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_species %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,])

top5_path_fig <- top5_path
top5_path_fig <- top5_path_fig %>% 
  mutate(path=str_replace_all(path, ';0','')) %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  ###scale_y_continuous(limits = c(-7,7))+
  scale_fill_manual(values = bluered(100)[c(1,10,20,30,40,100,90,80,70,60)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')


# essential edges (reactions) analysis for species substructures

species_net_w_rev <- add_rev_rection(network_edge, species_net)


# calculate perturbation score for each edge (reaction)

perturbation_score_species <- calculate_perturbation_point(species_net_w_rev,
                                                           species_sub_exp[[3]],
                                                           species_sub_exp_t,
                                                           ctrl=1:7, exp=8:13,
                                                           stat = 'p')

top5_rep_path <- c(path_score_species %>% filter(Type=='Active') %>% 
                     arrange(desc(cal_score)) %>% 
                     .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path
                   ,path_score_species %>% filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:5,] %>% .$rep_sub_path)
top5_rep_path <- path_score_species %>% filter(rep_sub_path%in%top5_rep_path, Significant=='yes') %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top5_path <- perturbation_score_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top5_rep_path && x[2] %in% top5_rep_path})
edge_in_top5_path <- perturbation_score_species[edge_in_top5_path,]$edge_name

top5_edge <- perturbation_score_species%>% 
  filter(edge_name %in%edge_in_top5_path) %>% 
  filter(p_value<0.05) %>% 
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:###FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:###0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:###FF0000'>", node2, "</i>"),
                            paste0("<i style='color:###0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color, FA_change))

top5_edge_fig <- top5_edge %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase','Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(top5_edge$edge_color)
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


# construct species biosynthetic network. To simplify, we only drew top3 pathways

top3_rep_path <- c(path_score_species %>% filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:3,] %>% .$rep_sub_path,
                   path_score_species %>% filter(Type=='Active')  %>% 
                     .[!duplicated(.$rep_sub_path),] %>% .[1:3,] %>% .$rep_sub_path)

top3_rep_path <- path_score_species %>% filter(rep_sub_path%in%top3_rep_path, Significant=='yes') %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top3_path <- perturbation_score_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top3_rep_path && x[2] %in% top3_rep_path})
edge_in_top3_path <- perturbation_score_species[edge_in_top3_path,]$edge_name

top3_net_edge <- species_net_w_rev %>% filter(S1 %in% top3_rep_path, P1 %in% top3_rep_path)
top3_net_edge <- top3_net_edge %>% mutate(color='gray', arrows='to',length=100)

colnames(top3_net_edge)[1:2] <- c('from','to')

top3_net_node <- species_network_data[[1]] %>% 
  filter(id %in% unlist(top3_net_edge))

top3_net_path_col <- rbind(path_score_species %>% filter(Type=='Active') %>% 
                             .[!duplicated(.$rep_sub_path),] %>% .[1:3,], 
                           path_score_species %>% filter(Type=='Suppressed') %>% 
                             arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% 
                             .[1:3,]) %>% mutate(rank=c(1,2,3,6,7,8))

top3_net_reaction_col <- perturbation_score_species%>% 
  filter(edge_name %in% edge_in_top3_path) %>% 
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  filter(p_value<0.05)%>% 
  mutate(rank=c(1:10))

top3_net_node$label <- str_replace_all(top3_net_node$label , ';0','')

net_edge <- path_color(top3_net_edge, top3_net_path_col,top3_net_reaction_col)

net_node <- top3_net_node %>% filter(id %in%c(net_edge$from, net_edge$to))


# building biosynthetic network for top 3 activa and suppressed pathways@@

visNetwork(net_node, net_edge)

```

# License
