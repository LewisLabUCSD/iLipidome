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

file <- dirname(rstudioapi::getSourceEditorContext()$path)
#load(file.path(file,'DHA.RData'))
#-------------------Data upload-------------------

load(file.path(file,'Required_data/Required_function.RData'))
load(file.path(file,'Required_data/Required_data_DHA.Rdata'))

each_FA <- str_extract_all(char$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char <- char %>% mutate(each_FA=each_FA)


#-------------------Analysis for unprocessed data-------------------

no_sub_t <- non_processed_data_test(exp, char, 't.test', 'adj_p_value', 2:8,  9:14)

#-------------------Conventional lipidomics analysis-------------------

FA_length <- char$FA_split %>% str_extract_all('\\d+') %>% 
  map(.f = function(x){if(length(x)==6){str_c(x[1],x[4],sep = ',')}
    else if(length(x)==3){x[1]}
    else{NA}})

FA_db <- char$FA_split %>% str_extract_all('\\d+') %>% 
  map(.f = function(x){if(length(x)==6){str_c(x[2],x[5],sep = ',')}
    else if(length(x)==3){x[2]}
    else{NA}})


lipid_char_analysis <- function(exp, char, var){
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
  try <- Species2Char(exp, char, char_var = var)
  
  try <- data_process(try, exclude_var_missing=F, 
                      missing_pct_limit=70,
                      replace_zero=F, zero2what='min', xmin=0.5,
                      replace_NA=T, NA2what='min', ymin=0,
                      pct_transform=T,
                      data_transform=F,trans_type='log',
                      centering=F,
                      scaling=F)
  rownames(try) <- try[[1]]
  
  try <- t_test(try[-1], 1:7, 8:13, 't.test', 'adj_value') %>% 
    mutate(star=stars.pval(.$adj_p_value) %>% str_replace_all(' ','') %>% 
             str_replace_all('\\.',''))
  
  try <- data.frame(mean_all=c(try$mean_all,try$mean_all),
                    mean=c(try$mean_ctrl, try$mean_exp),
                    sd=c(try$sd_ctrl, try$sd_exp),
                    group=c(rep('Ctrl',nrow(try)),rep('DHA',nrow(try)))) %>% 
    cbind(try[c(1,7:15)])
  return(try)
}


#FA length analysis
FA_length_analysis <- lipid_char_analysis(exp, mutate(char, FA_length=unlist(FA_length)), 'FA_length')

FA_length_analysis_fig <- FA_length_analysis %>% 
  mutate(star=ifelse(mean<mean_all,'',star)) %>% 
  mutate(star_loc=ifelse(mean<mean_all,0,mean+sd)) %>% 
  ggplot(aes(x=lipid, y=mean, fill=group)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(x=lipid, y=star_loc+3, label=star))+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))+
  scale_y_sqrt()+
  labs(x='Chain length', title='', y='GPL+DAG (mol%)')+
  theme(legend.position = 'top')


#FA double bond analysis

FA_db_analysis <- lipid_char_analysis(exp, mutate(char, FA_db=unlist(FA_db)), 'FA_db')

FA_db_analysis_fig <- FA_db_analysis %>% 
  mutate(star=ifelse(mean<mean_all,'',star)) %>% 
  mutate(star_loc=ifelse(mean<mean_all,0,mean+sd)) %>% 
  ggplot(aes(x=lipid, y=mean, fill=group)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(x=lipid, y=star_loc+3, label=star))+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))+
  scale_y_sqrt()+
  labs(x='Double bond', title='', y='GPL+DAG (mol%)')+
  theme(legend.position = 'top')

#FA double bond-chain length combined analysis

FA_db_length_analysis_fig <- no_sub_t[[2]] %>% filter(type=='FA') %>%
  mutate(star=stars.pval(.$adj_p_value) %>% str_replace_all(' ','') %>% 
           str_replace_all('\\.','')) %>% 
  mutate(char2=str_extract(lipid, '\\d+')) %>% 
  mutate(char1=str_sub(str_extract(lipid, ':\\d+'),start = 2)) %>% 
  mutate(char1=factor(char1, levels = as.character(sort(unique(as.numeric(.$char1)))))) %>% 
  mutate(char2=factor(char2, levels = as.character(sort(unique(as.numeric(.$char2)))))) %>% 
  ggplot(aes(x=char2, y=char1, fill=log2FC))+
  geom_tile(color='black')+
  geom_text(aes(label = star),vjust = 0.75,size=5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",na.value = 'gray', midpoint = 0)+
  theme_bw()+
  theme(axis.text.x =element_text(size=12),
        axis.text.y =element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.title =element_text(size=14))+
  labs(x='Chain length', y='Double bond')

#Differential expression for lipid species

DE_lipid_fig <- no_sub_t[[2]] %>% filter(type=='species') %>%  
  mutate(log2FC=ifelse(is.infinite(log2FC),10*sign(log2FC),log2FC)) %>% 
  mutate(Significance=ifelse(log2FC>0, 'Increase','Decrease')) %>% 
  mutate(Significance=ifelse(sig=='yes', Significance, 'No change')) %>% 
  ggplot(aes(x=log2FC, y=mlog10padj, col=Significance)) +
  geom_point() + 
  scale_x_continuous(limits = c(-10,10), labels = c('-Inf','-5','0','5', 'Inf'))+
  scale_color_manual(values=c("blue", "red", "gray")) +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='dashed')+
  theme_classic()+
  theme(legend.position = 'top')+
  labs(y='-log10(padj)')


DE_lipid_wo_fdr_fig <- no_sub_t[[2]] %>% filter(type=='species') %>%  
  mutate(sig=ifelse(p_value<0.05, 'yes','no')) %>% 
  mutate(log2FC=ifelse(is.infinite(log2FC),10*sign(log2FC),log2FC)) %>% 
  mutate(Significance=ifelse(log2FC>0, 'Increase','Decrease')) %>% 
  mutate(Significance=ifelse(sig=='yes', Significance, 'No change')) %>% 
  ggplot(aes(x=log2FC, y=mlog10p, col=Significance)) +
  geom_point() + 
  scale_x_continuous(limits = c(-10,10), labels = c('-Inf','-5','0','5', 'Inf'))+
  scale_color_manual(values=c("blue", "red", "gray")) +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='dashed')+
  theme_classic()+
  theme(legend.position = 'top')+
  labs(y='-log10(p-value)')

#sample coverage

FA_coverage_fig <- data.frame(sample_coverage=apply(filter(no_sub_t[[1]], type=='FA')[-c(1,2)], MARGIN = 1, FUN = function(x){sum(x!=0)})) %>% 
  gghistogram(x = "sample_coverage", 
    rug = TRUE, alpha = 0.5,
    fill="#E7B800",binwidth = 1)+
  scale_x_continuous(breaks = 4:13, labels = as.character(4:13))+
  labs(x='Sample coverage (1/13)', title='Fatty acid')+
  theme(plot.title = element_text(hjust = 0.5))

species_coverage_fig <- data.frame(sample_coverage=apply(filter(no_sub_t[[1]], type=='species')[-c(1,2)], MARGIN = 1, FUN = function(x){sum(x!=0)})) %>% 
  gghistogram( x = "sample_coverage", 
    rug = TRUE, alpha = 0.5,
    fill="#E7B800",binwidth = 1)+
  scale_x_continuous(breaks = 4:13, labels = as.character(4:13))+
  labs(x='Sample coverage (1/13)', title='Lipid species')+
  theme(plot.title = element_text(hjust = 0.5))

#FA enrichment analysis

FA_enrich_analysis <- function(all_lipid, stat_data, type, p='padj'){
  
  all <- all_lipid %>% str_extract_all('\\d+:\\d+;\\d+') %>% 
    unlist() %>% table()
  if(p!='padj'){
    stat_data <- stat_data %>% mutate(adj_p_value=p_value)
  }
  
  if(type=='up'){
    sig <- stat_data %>%filter(adj_p_value<0.05, log2FC>0) %>% .$lipid
    
    sig <- str_extract_all(sig,'\\d+:\\d+;\\d+') %>% 
      unlist() %>% table()
    
  }
  else{
    sig <- stat_data %>%filter(adj_p_value<0.05, log2FC<0) %>% .$lipid
    
    sig <- str_extract_all(sig,'\\d+:\\d+;\\d+') %>% 
      unlist() %>% table()
  }
  
  lipid <- character()
  pvalue <- numeric()
  
  num <- 1
  for (a in names(all)){
    lipid[num] <- a
    
    if(length(which(names(sig)==a))==0){
      pvalue[num] <- tryCatch({fisher.test(matrix(c(0,
                                                    sum(sig),
                                                    all[names(all)==a],
                                                    sum(all)-sum(sig)-all[names(all)==a]),nrow = 2),alternative = 'greater')$p.value},
                              warning=function(x){NA},
                              error=function(x){NA})
    }
    else{
      pvalue[num] <- tryCatch({fisher.test(matrix(c(sig[names(sig)==a],
                                                    sum(sig),
                                                    all[names(all)==a]-sig[names(sig)==a],
                                                    sum(all)-sum(sig)),nrow = 2),alternative = 'greater')$p.value},
                              warning=function(x){NA},
                              error=function(x){NA})
    }
    
    num <- num+1
  }
  if(type=='up'){
    return(data.frame(FA=lipid, pvalue_up=pvalue))
  }
  else{
    return(data.frame(FA=lipid, pvalue_down=pvalue))
  }
}


species_t <- no_sub_t[[2]] %>% filter(type=='species')


enrich_FDR <- cbind(FA_enrich_analysis(species_t$lipid[-c(1:11)], 
                                       species_t, 'up'),
                    FA_enrich_analysis(species_t$lipid[-c(1:11)], 
                                       species_t,'down')) %>% .[-3] %>% 
  gather(-FA, key = 'group', value='pvalue') %>% 
  mutate(mlogp=-log10(pvalue)) %>% 
  mutate(mlogp=ifelse(group=='pvalue_down',-mlogp, mlogp)) %>% 
  mutate(group=ifelse(pvalue>0.05, 'no', group)) %>% 
  arrange(desc(mlogp)) %>% 
  mutate(Data='Adjusted p-value\n(FDR)')


enrich_p <- cbind(FA_enrich_analysis(filter(no_sub_t[[2]], type=='species')$lipid[-c(1:11)], 
                         species_t, 'up', 'p'),
      FA_enrich_analysis(filter(no_sub_t[[2]], type=='species')$lipid[-c(1:11)], 
                         species_t,'down','p')) %>% .[-3] %>% 
  gather(-FA, key = 'group', value='pvalue') %>% 
  mutate(mlogp=-log10(pvalue)) %>% 
  mutate(mlogp=ifelse(group=='pvalue_down',-mlogp, mlogp)) %>% 
  mutate(group=ifelse(pvalue>0.05, 'no', group)) %>% 
  arrange(desc(mlogp)) %>% 
  mutate(Data='Original p-value')




enrich_FDR_fig <- rbind(enrich_FDR, enrich_p) %>% 
  mutate(FA=str_replace(FA, ';0','')) %>% 
  filter(FA %in%c('22:6', '16:0', '20:1', 
                  '17:1', '20:5', '20:2', '18:1')) %>% 
  mutate(FA=factor(FA, levels =rev(c('22:6', '20:5', '16:0', 
                                     '17:1', '18:1', '20:1',  '20:2')))) %>% 
  group_by(Data) %>% 
  filter(!(pvalue==1& FA!='20:5')) %>% 
  mutate(sig=ifelse(pvalue<0.05, 'yes', 'no'),
         color=ifelse(mlogp>0, 'Enriched','Depleted')) %>% 
  mutate(color=ifelse(pvalue<0.05, color,'Non-sig')) %>% 
  ggplot(aes(x=FA, y=mlogp, fill=Data, color=color))+
  geom_bar(stat='identity',position = position_dodge(), size=0.7)+
  coord_flip()+
  geom_hline(yintercept = c(-log10(0.05),log10(0.05)), linetype='dashed',color='black')+
  geom_hline(yintercept = 0, color='black')+
  theme_bw()+
  theme(legend.position = 'right',
        plot.title = element_text(hjust = 0.5,size = 10))+
  scale_fill_manual(values=c('#E69F00','#999999'), breaks = c('Original p-value','Adjusted p-value\n(FDR)' ))+
  scale_color_manual(values =c('#EE0000FF','blue','white'), breaks = c('Enriched','Depleted','Non-sig'))+
  scale_y_continuous(breaks = c(-5,0,5,10),labels = c('5','0','5','10'))+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))+
  labs(y='-log10(p-value)', x='Fatty acid enrichment', color='For DHA group')

grid.arrange(DE_lipid_fig, species_coverage_fig, enrich_FDR_fig,
             layout_matrix=rbind(c(1,1,1,2,2,2,3,3,3,3)))

#999999','#E69F00
grid.arrange(FA_length_analysis_fig, FA_db_analysis_fig, FA_db_length_analysis_fig,
             layout_matrix=rbind(c(1,1,1,1,2,2,2,2,3,3,3,3,3,3)))

#-------------------extract FA substructure using fold changes-------------------
FA_t <- no_sub_t[[2]] %>% filter(type=='FA')
#18:2 and 20:4 are majorly omega-6 fatty acids, so we only kept omega-6 pathways for them
FA_substructure_sel <- FA_substructure[-c(19,20,42),]

FA_t <- no_sub_t[[2]] %>% filter(type=='FA')

FA_sub_stop <- FA_substructure_w_stop(FA_substructure_sel, FA_t,exact_FA='no', exo_lipid='w3-22:6;0')


FA_sub_stop <- FA_sub_stop %>% mutate(FA=str_extract(FA, '\\d+:\\d+;\\d+'))

FA_sub_stop <- FA_substructure_transform(char, FA_sub_stop)


#-------------------transform lipid FA exp data into substructue exp data-------------------

FA_sub_exp <- lipid_substructure_matrix(exp, FA_sub_stop, 'FA')

#-------------------differential expression analysis for FA substructures-------------------

FA_sub_exp_t <- t_test(FA_sub_exp[[3]], 1:7, 8:13, 't.test', 'adj_value')

#-------------------FA substructures improve statistical power-------------------
options(scipen = -1)
FA_stat_improve_fig <- FA_sub_exp_t %>% 
  mutate(simple=str_extract(lipid, '\\d+:\\d+;\\d+')) %>% 
  left_join(FA_t, by=c('simple'='lipid')) %>% 
  .[c('lipid', 'mlog10padj.y', 'mlog10padj.x')] %>% 
  `colnames<-`(c('lipid','Fatty acid','Substructure')) %>% 
  filter(!is.na(`Fatty acid`)) %>% 
  ggpaired(cond1 = "Fatty acid", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = '-log10(padj)')+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed')+
  scale_y_continuous(limits = c(0,4))+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.2, label.y = 3.8)+
  theme(legend.position = 'none', axis.text.x = element_text(size=10))

FA_stat_improve_sub_fig <- FA_sub_exp_t %>% 
  mutate(simple=str_extract(lipid, '\\d+:\\d+;\\d+')) %>% 
  left_join(FA_t, by=c('simple'='lipid')) %>% 
  .[c('lipid', 'mlog10padj.y', 'mlog10padj.x')] %>% 
  `colnames<-`(c('lipid','Fatty acid','Substructure')) %>% 
  filter(!is.na(`Fatty acid`)) %>% 
  left_join(gather(FA_network, -pathway, key='key', value='lipid')[-2], by='lipid') %>%
  filter(!is.na(pathway)) %>% unique() %>% 
  ggpaired(cond1 = "Fatty acid", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = '-log10(padj)')+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed')+
  facet_wrap(~pathway)+
  scale_y_continuous(limits = c(0,4))+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.15, label.y = 3.6)+
  theme(legend.position = 'none', axis.text.x = element_text(size=10))


#-------------------FA biosynthetic network data transformation-------------------

FA_network_data <- draw_network(FA_network, FA_sub_exp_t, 'no','padj')

visNetwork(FA_network_data[[1]], FA_network_data[[2]])%>% 
  visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                  physics = F, smooth = TRUE, randomSeed =5) 

#-------------------essential pathway analysis for FA substructures-------------------

set.seed(2)

path_score_FA <- calculate_path_activation_node(FA_network, FA_sub_exp_t, calibrate = T, if_FA = T)

top5_path_FA <- rbind(path_score_FA %>% filter(Type=='Active') %>% 
                        .[!duplicated(.$rep_sub_path),] %>% .[1:5,], 
                      path_score_FA %>% filter(Type=='Suppressed') %>% 
                        arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>%.[1:5,])

top5_path_FA_fig <- top5_path_FA
top5_path_FA_fig$path[c(1, 3,4,8)] <- top5_path_FA_fig$path[c(1, 3,4,8)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], last(x),sep=' --> ')})

top5_path_FA_fig$path[c(6,7)] <- top5_path_FA_fig$path[c(6,7)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], x[2], last(x),sep=' --> ')})

top5_path_FA_fig$path[2] <- top5_path_FA_fig$path[2] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], tail(x,2)[1],tail(x,2)[2],sep=' --> ')})


top5_path_FA_fig <- top5_path_FA_fig%>% 
  mutate(path=str_replace_all(path, ';0','')) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(1.96,-1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(limits = c(-4.5,4.5))+
  scale_fill_manual(values = bluered(100)[c(100,90,80,70,60, 1,10,20)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 1),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Significant representative pathways')

#write.csv(path_score_FA, file.path(file, 'result/path_score_subFA_DHA.csv'), row.names = F)

#-------------------essential edges (reactions) analysis for FA substructures-------------------

#calculate pertubation score for each edge (reaction)
perturbation_score_FA <- calculate_perturbation_point(FA_network, 
                                                      FA_sub_exp[[3]], 
                                                      FA_sub_exp_t, 
                                                      ctrl = 1:7, exp = 8:13, 
                                                      stat = 'p')

top_5_edge_FA <-perturbation_score_FA %>%
  filter(!is.na(perturbation_score)) %>% 
  filter(p_value<0.05) %>% 
  filter(!str_detect(edge_name, '22:6')) %>% 
  .[-4,]

top_5_edge_FA_fig <- top_5_edge_FA %>%
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color))


top_5_edge_FA_fig <- top_5_edge_FA_fig %>%
  mutate(Mode=factor(Mode, levels = c("Increase","Decrease"))) %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(top_5_edge_FA_fig$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', 
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = rev(pal_lancet()(2)))+
  scale_color_manual(values = c('gold','white'))+
  scale_x_continuous(limits = c(-5,3))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 significant edges',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))

options(scipen = 0)
grid.arrange(FA_stat_improve_fig, top5_path_FA_fig, top_5_edge_FA_fig,
             layout_matrix=rbind(c(1,1,1,2,2,2,2,3,3,3,3,3)))


#write.csv(perturbation_score_FA, file.path(file, 'result/perturbation_score_subFA_DHA.csv'), row.names = F)

#-------------------construct FA biosynthetic network-------------------
top5_path_FA <- rbind(path_score_FA %>% filter(Type=='Active') %>% 
                        .[!duplicated(.$rep_sub_path),] %>% .[1:5,], 
                      path_score_FA %>% filter(Type=='Suppressed') %>% 
                        arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>%.[1:5,])

sig_path_col <- top5_path_FA %>%  
  filter(Significant=='yes') %>% 
  mutate(rank=c(1:5,6,7,8))


sig_edge_col <- top_5_edge_FA %>%
  mutate(rank=c(1:3, 6:10))

#building network

visNetwork(FA_network_data[[1]], path_color(FA_network_data[[2]], sig_path_col, sig_edge_col)) %>% 
  visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                  physics = F, smooth = TRUE, randomSeed =5) 


# path_color(FA_network_data[[2]], sig_path_col, sig_edge_col) %>%
#    left_join(FA_sub_exp_t, by=c('from'='lipid')) %>% 
#    left_join(FA_sub_exp_t, by=c('to'='lipid')) %>% 
#    .[c(1:9,16,21,22,29,34,35)] %>% 
#    write.csv(file.path(file, 'network/subFA_network_DHA.csv'),
#              row.names = F, na = '')
  

#-------------------essential pathway analysis for raw FA data-------------------

FA_mapping <- FA_substructure %>% apply(MARGIN = 1, FUN = function(x){c(x[1], last(x[x!='']))}) %>% 
  t() %>% as.data.frame() %>% unique()

FA_t2 <- FA_t %>% left_join(FA_mapping, by=c('lipid'='FA')) %>% 
  mutate(lipid=V2) %>% dplyr::select(-V2)

set.seed(7)

path_score_raw_FA <- calculate_path_activation_node(FA_network, FA_t2, calibrate = T, if_FA = T)


top5_path_raw_FA <- rbind(path_score_raw_FA %>% filter(Type=='Active') %>% 
                        .[!duplicated(.$rep_sub_path),] %>% .[1:5,], 
                        path_score_raw_FA %>% filter(Type=='Suppressed') %>% 
                        arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>%.[1:5,])


top5_path_raw_FA_fig <- top5_path_raw_FA
top5_path_raw_FA_fig$path[c(3,4,5,7,8)] <- top5_path_raw_FA_fig$path[c(3,4,5,7,8)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], last(x),sep=' --> ')})


top5_path_raw_FA_fig$path[1] <- top5_path_raw_FA_fig$path[1] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], tail(x,2)[1],tail(x,2)[2],sep=' --> ')})

top5_path_raw_FA_fig <- top5_path_raw_FA_fig%>% 
  mutate(path=str_replace_all(path, ';0','')) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(1.96,-1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  scale_y_continuous(limits = c(-4.5,4.5))+
  scale_fill_manual(values = bluered(100)[c(100,90,70,80,60, 1,10,20)])+
  theme(legend.position='none',
        plot.title = element_text(hjust = 1),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Significant representative pathways')

#write.csv(path_score_raw_FA, file.path(file, 'result/path_score_rawFA_DHA.csv'), row.names = F)

#-------------------essential edges (reactions) analysis for raw FA data-------------------
#calculate perturbation score for each edge (reaction)

FA_exp<- FA_mapping %>% 
  left_join(filter(no_sub_t[[1]], type=='FA'), by=c('FA'='feature')) %>% 
  filter(!is.na(Ctrl1)) %>% 
  column_to_rownames(var='V2') %>% 
  select(-1,-2)


perturbation_score_raw_FA <- calculate_perturbation_point(FA_network,
                                                          FA_exp, 
                                                          FA_t2, 
                                                          ctrl=1:7,exp=8:13,
                                                          stat = 'p')


top_5_edge_raw_FA <-perturbation_score_raw_FA %>%
  filter(!is.na(perturbation_score)) %>% 
  filter(p_value<0.05) %>% 
  filter(!str_detect(edge_name, '22:6')) %>% 
  .[c(1:3,(nrow(.)-4):(nrow(.))),]

top_5_edge_raw_FA_fig <- top_5_edge_raw_FA %>%
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color))


top_5_edge_raw_FA_fig <- top_5_edge_raw_FA_fig %>%
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase', 'Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(top_5_edge_raw_FA_fig$edge_color)
  ) +
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(legend.position = 'right', 
        axis.text.y = element_markdown(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = rev(pal_lancet()(2)))+
  scale_color_manual(values = c('gold','white'))+
  scale_x_continuous(limits = c(-10,10))+
  labs(y='', fill='Reaction', x='Perturbation score', 
       title='Top 5 significant edges',color='Edge type')+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))

grid.arrange(top5_path_raw_FA_fig, top_5_edge_raw_FA_fig,
             layout_matrix=rbind(c(1,1,1,2,2,2,2)))

#write.csv(perturbation_score_raw_FA, file.path(file, 'result/perturbation_score_rawFA_DHA.csv'), row.names = F)

#-------------------construct network for raw FA-------------------

FA_network_data_raw <- draw_network(FA_network, FA_t2, 'no','padj')
sig_path_raw_col <- top5_path_raw_FA %>%  
  filter(Significant=='yes') %>% 
  mutate(rank=c(1:5,6,7,8))


sig_edge_raw_col <- top_5_edge_raw_FA %>%
  mutate(rank=c(1:3, 6:10))


#building network

visNetwork(FA_network_data_raw[[1]], path_color(FA_network_data_raw[[2]], sig_path_raw_col, sig_edge_raw_col)) %>% 
  visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                  physics = F, smooth = TRUE, randomSeed =5) 


# path_color(FA_network_data_raw[[2]], sig_path_raw_col, sig_edge_raw_col) %>%
#   left_join(FA_t2, by=c('from'='lipid')) %>% 
#   left_join(FA_t2, by=c('to'='lipid')) %>% 
#   .[c(1:9,17,22,23,31,36,37)] %>%
#   write.csv(file.path(file, 'network/rawFA_network_DHA.csv'),
#             row.names = F, na = '')

#-------------------FA substructures improve hierarchical clustering-------------------
cluster_compare <- function(data, group){
  dist.fun=function(x){
    x=t(x)
    cor.mat=cor(x,method='spearman',use = 'complete.obs')
    cor.mat=(1-cor.mat)
    cor.dist=as.dist(cor.mat)
    return(cor.dist)
  }
  d <- dist.fun(data)

  fit <- hclust(d, method="complete") 
  groups <- cutree(fit, k=2)
  res.stat <- cluster.stats(d, group, groups)
  return(list(groups, res.stat))
}

adjRand_test=function(A, B, perm=1000) {
  if (length(A)!=length(B)) { stop("A and B should have the same length") }
  # Make sure that the two groups of partitions have the same length
  
  ARIo=mclust::adjustedRandIndex(A, B)
  # Observed adjusted Rand index
  
  Aperm=lapply(seq_len(perm), function(X) sample(A, length(A), replace=FALSE))
  Bperm=lapply(seq_len(perm), function(X) sample(B, length(B), replace=FALSE))
  # Generate permuted samples
  ARIperm=unlist(lapply(seq_len(perm),
                        function(i) mclust::adjustedRandIndex(Aperm[[i]], Bperm[[i]])))
  # Compute adjusted Rand index for the permuted samples
  m=mean(ARIperm); v=var(ARIperm)
  # compute mean and variance of the permuted samples
  
  NARI=(ARIperm-m)/sqrt(v)
  # compute NARI (normalized ARI)
  
  NARIo=(ARIo-m)/sqrt(v)
  # compute observed NARI
  
  p_value=(length(which(NARI>NARIo))+1)/(perm+1)
  # Compute p value as proportion of permuted NARI larger than the observed
  
  Results=c(ARIo, p_value)
  names(Results)=c("Adjusted_Rand_index", "p_value")
  return(Results)
}

#Clustering for raw FA data
FA_exp_raw <- no_sub_t[[1]] %>% filter(type=='FA') %>% 
  column_to_rownames(var='feature') %>% .[-1]


FA_exp_m <- as.matrix(FA_exp_raw) %>% t() %>% scale()

res.stat <- cluster_compare(FA_exp_m, c(rep(1,7),rep(2,6)))
set.seed(2)
res.stat <- adjRand_test(c(rep(1,7),rep(2,6)) ,res.stat[[1]])

ARI_col <- round(res.stat[1],3)
ARIP_col <- round(res.stat[2],3)

FA_group <- colnames(FA_exp_m) %>% str_extract_all('\\d+') %>% 
  map(.f = function(x){tail(x,2)[1]})%>% as.integer()
FA_group <- ifelse(FA_group>0 & FA_group<3, 1,2)


res.stat <- cluster_compare(t(FA_exp_m), FA_group)
set.seed(1)
res.stat <- adjRand_test(FA_group ,res.stat[[1]])

ARI_row <- round(res.stat[1],3)
ARIP_row <- round(res.stat[2],3)


FA_exp_heatmap <- FA_exp_m %>% t()


row_color <- rownames(FA_exp_heatmap) %>% str_extract_all('\\d+') %>% 
  map(.f = function(x){tail(x,2)[1]})%>% as.integer()
row_color <- cut(row_color, breaks = c(-1,0,2,6), labels = c('black','blue','red')) %>% 
  as.character()


FA_exp_heatmap_fig <- Heatmap(FA_exp_heatmap, clustering_method_rows = 'complete',
               clustering_method_columns = 'complete',
               clustering_distance_rows = 'spearman',
               clustering_distance_columns = 'spearman',
               column_split = 2,
               row_split = 2,
               column_names_gp = gpar(col=c('orange','green')),
               top_annotation = HeatmapAnnotation(foo=anno_block(gp=gpar(fill=c('orange','green')),
                                                                 labels=c('group1','group2'))),
               row_names_gp = gpar(fontsize = 8,
                                   col=row_color),
               heatmap_legend_param = list(title=''),
               show_heatmap_legend = F,
               column_title = str_c('Raw fatty acid\n',
                                    '(ARI = ',ARI_col,', p = ',ARIP_col,')'),
               row_title = str_c('(ARI = ',ARI_row,', p = ',ARIP_row,')'))


#Clustering for FA substructure data

FA_sub_exp_m <- as.matrix(FA_sub_exp[[3]]) %>% 
  t() %>% scale()

res.stat <- cluster_compare(FA_sub_exp_m, c(rep(1,7),rep(2,6)))
set.seed(1)
res.stat <- adjRand_test(c(rep(1,7),rep(2,6)) ,res.stat[[1]])

ARI_col <- round(res.stat[1],3)
ARIP_col <- round(res.stat[2],3)

FA_group <- colnames(FA_sub_exp_m) %>% str_extract_all('\\d+') %>% 
  map(.f = function(x){tail(x,2)[1]})%>% as.integer()
FA_group <- ifelse(FA_group>0 & FA_group<3, 1,2)


res.stat <- cluster_compare(t(FA_sub_exp_m), FA_group)
set.seed(1)
res.stat <- adjRand_test(FA_group ,res.stat[[1]])

ARI_row <- round(res.stat[1],3)
ARIP_row <- round(res.stat[2],3)


FA_sub_exp_heatmap <- FA_sub_exp_m %>% t()


row_color <- rownames(FA_sub_exp_heatmap) %>% str_extract_all('\\d+') %>% 
  map(.f = function(x){tail(x,2)[1]})%>% as.integer()
row_color <- cut(row_color, breaks = c(-1,0,2,6), labels = c('black','blue','red')) %>% 
  as.character()


FA_sub_exp_heatmap_fig <- Heatmap(FA_sub_exp_heatmap, clustering_method_rows = 'complete',
               clustering_method_columns = 'complete',
               clustering_distance_rows = 'spearman',
               clustering_distance_columns = 'spearman',
               column_split = 2,
               row_split = 2,
               column_dend_reorder = F,
               show_row_names = T,
               #column_title_gp = gpar(fill=c('orange','green')),
               column_names_gp = gpar(col=c('green','orange')),
               top_annotation = HeatmapAnnotation(foo=anno_block(gp=gpar(fill=c('green','orange')),
                                                                 labels=c('group2','group1'))),
               row_names_gp = gpar(fontsize = 8,
                                   col=row_color),
               show_heatmap_legend = F,
               column_title = str_c('Fatty acid substructure\n',
                                    '(ARI = ',ARI_col,', p = ',ARIP_col,')'),
               row_title = str_c('(ARI = ',ARI_row,', p = ',ARIP_row,')'))


#Outliar detect
DHA_sub_plot <- FA_sub_exp[[3]] %>% as.data.frame() %>% 
  .[str_detect(rownames(FA_sub_exp[[3]]), '22:6'),] %>% 
  gather(key='label',value='value') %>% 
  mutate(group=str_extract(label,'[A-Za-z]+')) %>% 
  mutate(label=ifelse(label=='DHA4', label,'')) %>% 
  ggboxplot(x = "group", y = "value",
            color = "group", palette = c('#999999','#E69F00'),
            add = "jitter")+ 
  stat_compare_means(method = "t.test",label.y = 22, 
                     label.x = 1.2)+
  theme_classic() +
  scale_y_continuous(limits = c(-1,25))+
  geom_text_repel(aes(label=label))+
  labs(x='', title='Substructure', y='DHA abundance')+
  theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))


DHA_raw_plot <- FA_exp %>% as.data.frame() %>% 
  .[str_detect(rownames(FA_exp), '22:6'),] %>% 
  gather(key='label',value='value') %>% 
  mutate(group=str_extract(label,'[A-Za-z]+')) %>% 
  mutate(label=ifelse(label=='DHA4', label,'')) %>% 
  ggboxplot(x = "group", y = "value",
            color = "group", palette = c('#999999','#E69F00'),
            add = "jitter")+ 
  stat_compare_means(method = "t.test",label.y = 22, 
                     label.x = 1.2)+
  theme_classic() +
  scale_y_continuous(limits = c(-1,25))+
  geom_text_repel(aes(label=label))+
  labs(x='', title='Fatty acid', y='DHA abundance')+
  theme(legend.position = 'right',plot.title = element_text(hjust = 0.5))

#-------------------FA substructures improve sample and network coverage-------------------

#sample coverage
FA_not_in_network <- FA_substructure$FA[77:128]


FA_raw_cov <- data.frame(sample_coverage=apply(no_sub_t[[1]] %>% filter(type=='FA', !feature %in% FA_not_in_network) %>% 
                                                 .[-c(1,2)], MARGIN = 1, 
                                               FUN = function(x){sum(x!=0)}))

FA_sub_cov <- data.frame(sample_coverage=apply(FA_sub_exp[[3]][!rownames(FA_sub_exp[[3]])%in%FA_not_in_network,], 
                                               MARGIN = 1, FUN = function(x){sum(x!=0)}))

FA_raw_cov_avg <- mean(FA_raw_cov$sample_coverage)
FA_sub_cov_avg <- mean(FA_sub_cov$sample_coverage)

FA_cov_improve_fig1 <- gghistogram(
  rbind(FA_raw_cov %>% mutate(Data='Fatty acid'),
        FA_sub_cov %>% mutate(Data='Substructure')), 
  x = "sample_coverage", 
  rug = TRUE, alpha = 0.5,
  fill = "Data", palette = c("#00AFBB", "#E7B800"),binwidth = 1)+
  scale_x_continuous(breaks = 4:13, labels = as.character(4:13))+
  labs(x='Sample coverage (1/13)')+
  annotate('text',x=8,y=25,
           label=str_c('Substrtucture improves ',
                       round((FA_sub_cov_avg-FA_raw_cov_avg)/FA_raw_cov_avg*100,1),'%'), 
           color="red", size=4)

FA_cov_improve_fig2 <- ggdensity(
  rbind(FA_raw_cov %>% mutate(Data='Fatty acid'),
        FA_sub_cov %>% mutate(Data='Substructure')), 
  x = "sample_coverage", 
  color= "Data", palette = c("#00AFBB", "#E7B800"),
  alpha = 0
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right")  +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend")+
  rremove("y.axis")+
  rremove("ylab") +
  rremove("y.text") +
  rremove("y.ticks") 


FA_cov_improve_fig <- align_plots(FA_cov_improve_fig1, FA_cov_improve_fig2, align="hv", axis="tblr")
FA_cov_improve_fig <- ggdraw(FA_cov_improve_fig[[1]]) + draw_plot(FA_cov_improve_fig[[2]])


#network coverage
FA_network_node <- filter(FA_network_data[[1]] , value!=0)$id

FA_raw_net_cov <- data.frame(lipid=FA_network_node) %>% 
  left_join(mutate(FA_exp, lipid=rownames(FA_exp)), by='lipid')

FA_raw_net_cov[is.na(FA_raw_net_cov)] <- 0

FA_raw_net_cov <- data.frame(sample_coverage=apply(FA_raw_net_cov[-1], MARGIN = 2, FUN = function(x){sum(x!=0)/length(x)*100}))


FA_sub_net_cov <- data.frame(lipid=FA_network_node) %>% 
  left_join(mutate(as.data.frame(FA_sub_exp[[3]]), lipid=rownames(FA_sub_exp[[3]])), by='lipid')

FA_sub_net_cov[is.na(FA_sub_net_cov)] <- 0

FA_sub_net_cov <- data.frame(sample_coverage=apply(FA_sub_net_cov[-c(1)], MARGIN = 2, FUN = function(x){sum(x!=0)/length(x)*100}))


improve_net_cov <- round(mean(FA_sub_net_cov$sample_coverage-FA_raw_net_cov$sample_coverage),1)


FA_net_cov_improve_fig <- cbind(FA_raw_net_cov, FA_sub_net_cov) %>% `colnames<-`(c('Fatty acid','Substructure')) %>% 
  ggpaired(cond1 = "Fatty acid", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = 'Network coverage',
           title = '')+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.16, label.y = 102)+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(labels = function(x) paste0(x,'%'))+
  annotate('text',x=1.5,y=108,label=str_c('Substrtucture improves ',improve_net_cov, '% (N=13)'), color="red", size=3.5)



grid.arrange(DHA_raw_plot,DHA_sub_plot,
             FA_cov_improve_fig, FA_net_cov_improve_fig, nrow=2)
#-------------------Decompose lipids into species substructures-------------------

#We excluded ether lipids because we don't know they are Alkyl (O-) or Alkenyl- (P-) linked
#so that we cannot map them to the pathways

char_sel <- char[!str_detect(char$feature, 'O-'),]
exp_sel <- exp[!str_detect(exp$feature, 'O-'),]


sub_species <- lipid_species_substructure_transform(char_sel, FA_substructure, lipid_substructure)


species_substructure <- sub_species[[2]] %>% map(.f = function(x){names(x) <- 1:length(x);return(x)}) %>% 
  plyr::ldply(rbind) %>% mutate(Lipid=sub_species[[1]]) %>% 
  dplyr::select(Lipid, everything())

species_substructure[is.na(species_substructure)] <- ''

colnames(species_substructure) <- c('Lipid', str_c('Unit', 1:(ncol(species_substructure)-1)))


#-------------------extract species substructures using fold changes-------------------

species_t <- no_sub_t[[2]] %>% filter(type=='species')


species_sub_stop <- lipid_substructure_w_stop(species_substructure, species_t, 'species', pct_limit = 0.3)

species_sub_stop <- list(species_sub_stop[[1]], apply(species_sub_stop[-1], MARGIN = 1, FUN = function(x){x[x!='']}))

#-------------------transform lipid exp data into species substructures exp data-------------------

species_sub_exp <- lipid_substructure_matrix(exp_sel, species_sub_stop, 'species')


#-------------------differential expression analysis for substructures-------------------

species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:7, 8:13, 't.test', 'adj_p_value')

#-------------------species substructures improve statistical power-------------------

stat_improve <- species_sub_exp_t %>% 
  left_join(species_t, by=c('lipid')) %>% 
  .[c('lipid', 'mlog10padj.y', 'mlog10padj.x')] %>% 
  `colnames<-`(c('lipid','Lipid species','Substructure')) 


stat_improve_venn <- ggvenn(
  list(`Significant\nLipid species` = stat_improve %>% filter(`Lipid species`> -log10(0.05)) %>% .$lipid, 
       `Significant\nSubstructure`= stat_improve %>% filter(Substructure> -log10(0.05)) %>% .$lipid),
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4,
)

stat_improve_fig <- stat_improve %>% 
  filter(!is.na(`Lipid species`)) %>% 
  filter(`Lipid species`>-log10(0.05)| Substructure>-log10(0.05)) %>% 
  ggpaired(cond1 = "Lipid species", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = '-log10(padj)')+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed')+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.24, label.y = 7.3)+
  theme(legend.position = 'none')


#-------------------species substructures improve enrichment analysis-------------------


#remove cardiolipins because the data does not provide their detailed FA composition
lipid_sub_wo_CL <- stat_improve$lipid[-c(1:11)]


enrich_sub <- cbind(FA_enrich_analysis(lipid_sub_wo_CL, 
                         species_sub_exp_t[-c(1:11),],'up'),
      FA_enrich_analysis(lipid_sub_wo_CL,
                         species_sub_exp_t[-c(1:11),],'down')) %>% .[-3] %>% 
  gather(-FA, key = 'group', value='pvalue') %>% 
  mutate(mlogp=-log10(pvalue)) %>% 
  mutate(mlogp=ifelse(group=='pvalue_down',-mlogp, mlogp)) %>% 
  mutate(group=ifelse(pvalue>0.05, 'no', group)) %>% 
  arrange(desc(mlogp)) %>% 
  mutate(Data='Substructure')

lipid_raw_wo_CL <- stat_improve %>% filter(!is.na(`Lipid species`)) %>% 
  .$lipid %>% .[-c(1:11)]


enrich_raw <- cbind(FA_enrich_analysis(lipid_raw_wo_CL, 
                                       filter(species_t, lipid%in%lipid_raw_wo_CL),
                                       'up'),
                    FA_enrich_analysis(lipid_raw_wo_CL,
                                       filter(species_t, lipid%in%lipid_raw_wo_CL),
                                       'down')) %>% .[-3] %>% 
  gather(-FA, key = 'group', value='pvalue') %>% 
  mutate(mlogp=-log10(pvalue)) %>% 
  mutate(mlogp=ifelse(group=='pvalue_down',-mlogp, mlogp)) %>% 
  mutate(group=ifelse(pvalue>0.05, 'no', group)) %>% 
  arrange(desc(mlogp)) %>% 
  mutate(Data='Lipid species')



FA_enrich_fig <- rbind(enrich_raw, enrich_sub) %>% 
  mutate(FA=str_replace(FA, ';0','')) %>% 
  filter(FA %in%c('22:6', '16:0', '18:0', '20:1', '17:1', '19:1', '20:2', '18:1')) %>% 
  mutate(FA=factor(FA, levels =rev(c('22:6', '16:0', '18:0', '20:1', '17:1', '19:1', '20:2', '18:1')))) %>% 
  group_by(Data) %>% 
  filter(!(pvalue==1 & FA!='20:1')) %>% 
  mutate(sig=ifelse(pvalue<0.05, 'yes', 'no'),
         color=ifelse(mlogp>0, 'Enriched','Depleted')) %>% 
  mutate(color=ifelse(pvalue<0.05, color,'Non-sig')) %>% 
  ggplot(aes(x=FA, y=mlogp, fill=Data, color=color))+
  geom_bar(stat='identity',position = position_dodge(), size=0.7)+
  coord_flip()+
  geom_hline(yintercept = c(-log10(0.05),log10(0.05)), linetype='dashed',color='black')+
  geom_hline(yintercept = 0, color='black')+
  theme_bw()+
  theme(legend.position = 'right',plot.title = element_text(hjust = 0.5,size = 10))+
  scale_fill_manual(values=c('#E69F00','#999999'), breaks = c('Substructure','Lipid species' ))+
  scale_color_manual(values =c('#EE0000FF','blue','white'), breaks = c('Enriched','Depleted','Non-sig'))+
  scale_y_continuous(breaks = c(-5,0,5),labels = c('5','0','5'))+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))+
  labs(y='-log10(p-value)', x='Fatty acid enrichment', color='For DHA group')


#-------------------species biosynthetic network data transformation-------------------


species_network <- split(species_substructure[-1], seq(nrow(species_substructure)))
species_network <- species_network %>% map(.f = function(x){x[x!='']})

species_network <- species_network %>% map(.f = function(x){head(rep(x, each=2)[-1],-1)})
species_network <-  matrix(unlist(species_network), ncol=2, byrow = T,dimnames = list(NULL, c('S1','P1'))) %>% 
  as.data.frame() %>%  unique()

species_network <- species_network %>% mutate(S1=str_replace(S1, '_FA\\d',''),
                                              P1=str_replace(P1, '_FA\\d',''))


species_network_data <- draw_network(species_network, species_sub_exp_t, 'no','adj_p_value')

#-------------------essential pathway analysis for species substructures-------------------

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

top5_path_fig$path[c(1:5,8:10)] <- top5_path_fig$path[c(1:5,8:10)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1],last(x),sep=' --> ')})

top5_path_fig$path[6:7] <- top5_path_fig$path[6:7] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], tail(x,2)[1],tail(x,2)[2],sep=' --> ')})



top5_path_fig <- top5_path_fig %>% 
  mutate(path=str_replace_all(path, ';0','')) %>% 
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

#write.csv(path_score_species, file.path(file, 'result/path_score_subspecies_DHA.csv'), row.names = F)


#-------------------essential edges (reactions) analysis for species substructures-------------------

species_net_w_rev <- add_rev_rection(network_edge, species_net)


#calculate perturbation score for each edge (reaction)
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
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
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

#write.csv(perturbation_score_species, file.path(file, 'result/perturbation_score_subspecies_DHA.csv'), row.names = F)

#-------------------construct species biosynthetic network-------------------

#To simplify, we only drew top3 pathways
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

#building network

visNetwork(net_node, net_edge)



species_sub_t_net <- species_sub_exp_t
maxinf <- ceiling(max(species_sub_t_net$log2FC[is.finite(species_sub_t_net$log2FC)]))
mininf <- ceiling(min(species_sub_t_net$log2FC[is.finite(species_sub_t_net$log2FC)]))
species_sub_t_net$log2FC[species_sub_t_net$log2FC>0 & is.infinite(species_sub_t_net$log2FC)] <- maxinf
species_sub_t_net$log2FC[species_sub_t_net$log2FC<0 & is.infinite(species_sub_t_net$log2FC)] <- mininf


# path_color(top3_net_edge, top3_net_path_col,top3_net_reaction_col) %>% 
#   left_join(species_sub_t_net, by=c('from'='lipid')) %>% 
#   left_join(species_sub_t_net, by=c('to'='lipid')) %>% 
#   .[c(1:9,16,21,22,29,34,35)] %>% 
#   write.csv(file.path(file, 'network/subspecies_network_DHA.csv'),
#             row.names = F, na = '')


#BIOPAN's algorithm cannot build the connections involving FA transfer if no free FA data provided
net_edge_biopan <- net_edge %>% filter(from %in% char$feature, to %in% char$feature)

net_edge_biopan <- net_edge_biopan %>% left_join(char[c('feature','class')], by=c('from'='feature')) %>% 
  left_join(char[c('feature','class')], by=c('to'='feature')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.x'='Abbreviation')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.y'='Abbreviation')) %>% 
  filter(FA.x==FA.y)

net_node_biopan <- net_node  %>% filter(id %in% char$feature)


visNetwork(net_node_biopan, net_edge_biopan)

#-------------------essential pathway analysis for raw species data-------------------

raw_species_net <- species_network_data[[2]][1:2]
colnames(raw_species_net) <- c('S1','P1')
raw_species_net <- raw_species_net%>% filter(S1%in%species_t$lipid,
                                    P1%in%species_t$lipid)


set.seed(1)
path_score_raw_species <-  calculate_path_activation_node(raw_species_net, species_t,
                                                          calibrate = T, if_FA = F)


top5_path_raw <- rbind(path_score_raw_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_raw_species %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,])

top5_path_raw_fig <- top5_path_raw

top5_path_raw_fig$path[c(1:10)] <- top5_path_raw_fig$path[c(1:10)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1],last(x),sep=' --> ')})



top5_path_raw_fig <- top5_path_raw_fig %>% 
  mutate(path=str_replace_all(path, ';0','')) %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  scale_fill_manual(values = c(bluered(100)[c(1,10,20,30)],'gray',bluered(100)[c(100,90,80,70,60)], 'gray'))+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8))+
  labs(x='', y='Path score',
       title='Top 5 representative pathways')

#write.csv(path_score_raw_species, file.path(file, 'result/path_score_rawspecies_DHA.csv'), row.names = F)

#-------------------essential edges (reactions) analysis for raw species data-------------------

raw_species_net_w_rev <- add_rev_rection(network_edge, raw_species_net)


#calculate perturbation score for each edge (reaction)
perturbation_score_raw_species <- calculate_perturbation_point(raw_species_net_w_rev,
                                                               column_to_rownames(exp, var = 'feature'),
                                                               species_t,
                                                               ctrl=1:7, exp=8:13,
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
             top5_edge_raw_fig,top5_edge_fig,
             layout_matrix=rbind(c(1,1,1,1,1,2,2,2,2,2,2),
                                 c(3,3,3,3,3,4,4,4,4,4,4)))

#write.csv(perturbation_score_raw_species, file.path(file, 'result/perturbation_score_rawspecies_DHA.csv'), row.names = F)

#-------------------construct network for raw species data-------------------

#To simplify, we only drew top3 pathways

top3_rep_path_raw <- c(path_score_raw_species %>% filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% .[1:3,] %>% .$rep_sub_path,
                   path_score_raw_species %>% filter(Type=='Active')  %>% 
                     .[!duplicated(.$rep_sub_path),] %>% .[1:3,] %>% .$rep_sub_path)



top3_rep_path_raw <- path_score_raw_species %>% filter(rep_sub_path%in%top3_rep_path_raw) %>% 
  .$path %>% str_split(' --> ') %>% unlist() %>% unique()

edge_in_top3_path_raw <- perturbation_score_raw_species$edge_name %>% str_split(' --> ') %>% 
  map_lgl(.f = function(x){x[1] %in% top3_rep_path_raw && x[2] %in% top3_rep_path_raw})
edge_in_top3_path_raw <- perturbation_score_raw_species[edge_in_top3_path_raw,]$edge_name



top3_net_edge_raw <- raw_species_net_w_rev %>% filter(S1 %in% top3_rep_path_raw, P1 %in% top3_rep_path_raw)

top3_net_edge_raw <- top3_net_edge_raw %>% mutate(color='gray', arrows='to',length=100)

colnames(top3_net_edge_raw)[1:2] <- c('from','to')


top3_net_node_raw <- species_network_data[[1]] %>% 
  filter(id %in% unlist(top3_net_edge_raw))

top3_net_path_col_raw <- rbind(path_score_raw_species %>% filter(Type=='Active') %>% 
                             .[!duplicated(.$rep_sub_path),] %>% .[1:3,], 
                             path_score_raw_species %>% filter(Type=='Suppressed') %>% 
                             arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>% 
                             .[1:3,]) %>% mutate(rank=c(1,2,3,6,7,8))

top3_net_reaction_col_raw <- perturbation_score_raw_species%>% 
  filter(edge_name %in% edge_in_top3_path_raw) %>% 
  .[c(1:5, (nrow(.)-4):nrow(.)),] %>% 
  filter(p_value<0.05)%>% 
  mutate(rank=c(1:3,6:10))


top3_net_node_raw$label <- str_replace_all(top3_net_node_raw$label , ';0','')

net_edge_raw <- path_color(top3_net_edge_raw, top3_net_path_col_raw,top3_net_reaction_col_raw)


net_node_raw <- top3_net_node_raw %>% filter(id %in%c(net_edge_raw$from, net_edge_raw$to))


#building network

visNetwork(net_node_raw, net_edge_raw)


species_raw_t_net <- species_t
maxinf <- ceiling(max(species_raw_t_net$log2FC[is.finite(species_raw_t_net$log2FC)]))
mininf <- ceiling(min(species_raw_t_net$log2FC[is.finite(species_raw_t_net$log2FC)]))
species_raw_t_net$log2FC[species_raw_t_net$log2FC>0 & is.infinite(species_raw_t_net$log2FC)] <- maxinf
species_raw_t_net$log2FC[species_raw_t_net$log2FC<0 & is.infinite(species_raw_t_net$log2FC)] <- mininf


# path_color(top3_net_edge_raw, top3_net_path_col_raw,top3_net_reaction_col_raw) %>% 
#   left_join(species_raw_t_net, by=c('from'='lipid')) %>% 
#   left_join(species_raw_t_net, by=c('to'='lipid')) %>% 
#    .[c(1:9,17,22,23,31,36,37)] %>% 
#    write.csv(file.path(file, 'network/rawspecies_network_DHA.csv'),
#              row.names = F, na = '')


#-------------------species substructures improve hierarchical clustering-------------------

#Clustering for raw FA data

species_exp_raw <- no_sub_t[[1]] %>% filter(type=='species') %>% 
  filter(!str_detect(feature, 'O-')) %>% 
  column_to_rownames(var='feature') %>% .[-1]

species_exp_m <- as.matrix(species_exp_raw) %>% t() %>% scale()

res.stat <- cluster_compare(species_exp_m, c(rep(1,7),rep(2,6)))
set.seed(1)
res.stat <- adjRand_test(c(rep(1,7),rep(2,6)) ,res.stat[[1]])

ARI_col <- round(res.stat[1],3)
ARIP_col <- round(res.stat[2],3)


species_exp_heatmap <- species_exp_m %>% t()

species_exp_heatmap_fig <- Heatmap(species_exp_heatmap, clustering_method_rows = 'complete',
        clustering_method_columns = 'complete',
        clustering_distance_rows = 'spearman',
        clustering_distance_columns = 'spearman',
        show_row_names = F,
        heatmap_legend_param = list(title=''),
        column_names_gp =  gpar(col=c(rep('orange',7),rep('green',6))),
        column_split = 2,
        show_heatmap_legend = F,
        column_title = str_c('Lipid species\n',
                             '(ARI = ',ARI_col,', p = ',ARIP_col,')'))


#Clustering for species substructure data

species_sub_exp_m <- as.matrix(species_sub_exp[[3]]) %>% t() %>% scale()


res.stat <- cluster_compare(species_sub_exp_m, c(rep(1,7),rep(2,6)))
set.seed(1)
res.stat <- adjRand_test(c(rep(1,7),rep(2,6)) ,res.stat[[1]])

ARI_col <- round(res.stat[1],3)
ARIP_col <- round(res.stat[2],3)

species_sub_exp_heatmap <- species_sub_exp_m %>% t()

species_sub_exp_heatmap_fig <- Heatmap(species_sub_exp_heatmap, clustering_method_rows = 'complete',
                                   clustering_method_columns = 'complete',
                                   clustering_distance_rows = 'spearman',
                                   clustering_distance_columns = 'spearman',
                                   show_row_names = F,
                                   heatmap_legend_param = list(title=''),
                                   column_names_gp =  gpar(col=c(rep('orange',7),rep('green',6))),
                                   column_split = 2,
                                   show_heatmap_legend = F,
                                   column_title = str_c('Lipid species\n',
                                                        '(ARI = ',ARI_col,', p = ',ARIP_col,')'))

#-------------------species substructures improve sample and network coverage-------------------

#sample coverage

species_raw_cov <- data.frame(sample_coverage=apply(species_exp_raw, MARGIN = 1, FUN = function(x){sum(x!=0)}))
species_sub_cov <- data.frame(sample_coverage=apply(species_sub_exp[[3]], MARGIN = 1, FUN = function(x){sum(x!=0)}))

species_raw_cov_avg <- mean(species_raw_cov$sample_coverage)
species_sub_cov_avg <- mean(species_sub_cov$sample_coverage)


species_cov_improve_fig1 <- gghistogram(
  rbind(species_raw_cov %>% mutate(Data='Lipid species'),
        species_sub_cov %>% mutate(Data='Substructure')), 
  x = "sample_coverage", 
  rug = TRUE, alpha = 0.5,
  fill = "Data", palette = c("#00AFBB", "#E7B800"),binwidth = 1)+
  scale_x_continuous(breaks = 4:13, labels = as.character(4:13))+
  labs(x='Sample coverage (1/13)')+
  annotate('text',x=8,y=250,
           label=str_c('Substrtucture improves ',
                       round((species_sub_cov_avg-species_raw_cov_avg)/species_raw_cov_avg*100,1),'%'), 
           color="red", size=4)



species_cov_improve_fig2 <- ggdensity(
  rbind(species_raw_cov %>% mutate(Data='Lipid species'),
        species_sub_cov %>% mutate(Data='Substructure')), 
  x = "sample_coverage", 
  color= "Data", palette = c("#00AFBB", "#E7B800"),
  alpha = 0
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right")  +
  theme_half_open(11, rel_small = 1) +
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend")+
  rremove("y.axis")+
  rremove("ylab") +
  rremove("y.text") +
  rremove("y.ticks") 


species_cov_improve_fig <- align_plots(species_cov_improve_fig1, species_cov_improve_fig2, align="hv", axis="tblr")
species_cov_improve_fig <- ggdraw(species_cov_improve_fig[[1]]) + draw_plot(species_cov_improve_fig[[2]])



#network coverage

species_network_node <- species_network_data[[1]]$id

species_raw_net_cov <- data.frame(lipid=species_network_node) %>% 
  left_join(mutate(species_exp_raw, lipid=rownames(species_exp_raw)), by='lipid')

species_raw_net_cov[is.na(species_raw_net_cov)] <- 0

species_raw_net_cov <- data.frame(sample_coverage=apply(species_raw_net_cov[-1], MARGIN = 2, FUN = function(x){sum(x!=0)/length(x)*100}))


species_sub_net_cov <- data.frame(lipid=species_network_node) %>% 
  left_join(mutate(as.data.frame(species_sub_exp[[3]]),
                   lipid=rownames(species_sub_exp[[3]])), by='lipid')

species_sub_net_cov[is.na(species_sub_net_cov)] <- 0

species_sub_net_cov <- data.frame(sample_coverage=apply(species_sub_net_cov[-c(1)], MARGIN = 2, FUN = function(x){sum(x!=0)/length(x)*100}))


improve_net_cov <- round(mean(species_sub_net_cov$sample_coverage-species_raw_net_cov$sample_coverage),1)


species_net_cov_improve_fig <- cbind(species_raw_net_cov, species_sub_net_cov) %>% `colnames<-`(c('Lipid species','Substructure')) %>% 
  ggpaired(cond1 = "Lipid species", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = 'Network coverage',
           title = '')+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.16, label.y = 100)+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(labels = function(x) paste0(x,'%'))+
  annotate('text',x=1.5,y=80,label=str_c('Substrtucture improves ',improve_net_cov, '% (N=13)'), color="red", size=3.5)

