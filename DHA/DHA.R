library(ggsci)
library(ggpubr)
library(ggvenn)
library(gridExtra)
library(ggrepel)
library(ggtext)
library(ComplexHeatmap)
library(fpc)
library(cowplot)

file <- dirname(rstudioapi::getSourceEditorContext()$path)



#-------------------Data upload-------------------
source(file.path(dirname(file),'Required_function/required_function.R'))
load(file.path(dirname(file),'Required_data/required_data.RData'))


raw_data <- read.csv(file.path(file,'lipidome_data/exp_DHA_raw.csv'))
char <- read.csv(file.path(file,'lipidome_data/char_DHA.csv'))

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

exp <- data_process(raw_data, exclude_var_missing = T, missing_pct_limit = 70,replace_zero = F,
                    replace_NA = T,NA2what = 'min',ymin = 0, pct_transform = T, data_transform = F)

exp <- remove_rownames(exp)


each_FA <- str_extract_all(char$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char <- char %>% mutate(each_FA=each_FA)


#-------------------Analysis for unprocessed data-------------------

no_sub_t <- unprocessed_data_test(exp, char, 't.test', 'adj_p_value', 1:7,  8:13)

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
  
  char_exp <- Species2Char(exp, char, char_var = var)
  
  char_exp_pro <- data_process(char_exp, exclude_var_missing=F, 
                               missing_pct_limit=70,
                               replace_zero=F, zero2what='min', xmin=0.5,
                               replace_NA=T, NA2what='min', ymin=0,
                               pct_transform=T,
                               data_transform=F,trans_type='log',
                               centering=F,
                               scaling=F)
  rownames(char_exp_pro) <- char_exp_pro[[1]]
  
  char_exp_t <- t_test(char_exp_pro[-1], 1:7, 8:13, 't.test', 'adj_p_value') %>% 
    mutate(star=stars.pval(.$adj_p_value) %>% str_replace_all(' ','') %>% 
             str_replace_all('\\.',''))
  
  char_exp_t <- data.frame(mean_all=c(char_exp_t$mean_all,char_exp_t$mean_all),
                           mean=c(char_exp_t$mean_ctrl, char_exp_t$mean_exp),
                           sd=c(char_exp_t$sd_ctrl, char_exp_t$sd_exp),
                           group=c(rep('Ctrl',nrow(char_exp_t)),rep('DHA',nrow(char_exp_t)))) %>% 
    cbind(char_exp_t[c(1,7:15)])
  
  char_exp_pro <- char_exp_pro %>% 
    gather(-1, key='sample', value='value') %>% 
    mutate(group=str_extract(sample,'[A-z]+'))
  return(list(char_exp_pro, char_exp_t))
}


#FA length analysis
FA_length_analysis <- lipid_char_analysis(exp, mutate(char, FA_length=unlist(FA_length)), 'FA_length')


F1d <- FA_length_analysis[[2]] %>% 
  mutate(star=ifelse(mean<mean_all,'',star)) %>% 
  mutate(star_loc=ifelse(mean<mean_all,0,mean+sd)) %>% 
  ggplot(aes(x=lipid, y=mean, fill=group)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_jitter(data = FA_length_analysis[[1]],
              aes(x=FA_length, y=value,color=group),
              position = position_dodge(width = 0.9),
              size=0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(x=lipid, y=star_loc+3, label=star),color='red')+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))+
  scale_color_manual(values=c('black','black'))+
  scale_y_sqrt()+
  labs(x='Chain length', title='', y='GPL+DAG (mol%)')+
  theme(legend.position = 'top')

#write.csv(FA_length_analysis[[1]], file.path(file, 'source_data/F1d1.csv'), row.names = F)
#write.csv(FA_length_analysis[[2]], file.path(file, 'source_data/F1d2.csv'), row.names = F)

#FA double bond analysis

FA_db_analysis <- lipid_char_analysis(exp, mutate(char, FA_db=unlist(FA_db)), 'FA_db')

F1e <- FA_db_analysis[[2]] %>% 
  mutate(star=ifelse(mean<mean_all,'',star)) %>% 
  mutate(star_loc=ifelse(mean<mean_all,0,mean+sd)) %>% 
  ggplot(aes(x=lipid, y=mean, fill=group)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_jitter(data = FA_db_analysis[[1]],
              aes(x=FA_db, y=value,color=group),
              position = position_dodge(width = 0.9),
              size=0.5)+
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_text(aes(x=lipid, y=star_loc+3, label=star), color='red')+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))+
  scale_color_manual(values=c('black','black'))+
  scale_y_sqrt()+
  labs(x='Double bond', title='', y='GPL+DAG (mol%)')+
  theme(legend.position = 'top')

#write.csv(FA_db_analysis[[1]], file.path(file, 'source_data/F1e1.csv'), row.names = F)
#write.csv(FA_db_analysis[[2]], file.path(file, 'source_data/F1e2.csv'), row.names = F)


#FA double bond-chain length combined analysis
F1f_data <- no_sub_t[[2]] %>% filter(type=='FA') %>%
  mutate(star=stars.pval(.$adj_p_value) %>% str_replace_all(' ','') %>% 
           str_replace_all('\\.','')) %>% 
  mutate(char2=str_extract(lipid, '\\d+')) %>% 
  mutate(char1=str_sub(str_extract(lipid, ':\\d+'),start = 2)) %>% 
  mutate(char1=factor(char1, levels = as.character(sort(unique(as.numeric(.$char1)))))) %>% 
  mutate(char2=factor(char2, levels = as.character(sort(unique(as.numeric(.$char2))))))

#write.csv(F1f_data, file.path(file, 'source_data/F1f.csv'), row.names = F)

F1f <- F1f_data %>% 
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

DE_lipid_data <- no_sub_t[[2]] %>% filter(type=='species') %>%  
  mutate(log2FC=ifelse(is.infinite(log2FC),10*sign(log2FC),log2FC)) %>% 
  mutate(Significance=ifelse(log2FC>0, 'Increase','Decrease')) %>% 
  mutate(Significance=ifelse(sig=='yes', Significance, 'No change'))

#write.csv(DE_lipid_data, file.path(file, 'source_data/F1a.csv'), row.names = F)


F1a <- DE_lipid_data %>% 
  ggplot(aes(x=log2FC, y=mlog10padj, col=Significance)) +
  geom_point() + 
  scale_x_continuous(limits = c(-10,10), labels = c('-Inf','-5','0','5', 'Inf'))+
  scale_color_manual(values=c("blue", "red", "gray")) +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='dashed')+
  theme_classic()+
  theme(legend.position = 'top')+
  labs(y='-log10(padj)')



#sample coverage

F1b_data <- data.frame(feature=filter(no_sub_t[[1]], type=='species')$feature,
                       sample_coverage=apply(filter(no_sub_t[[1]], type=='species')[-c(1,2)], MARGIN = 1, FUN = function(x){sum(x!=0)}))


#write.csv(F1b_data, file.path(file, 'source_data/F1b.csv'), row.names = F)

F1b <- F1b_data %>% 
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

lipid_exact_FA <- char %>% filter(FA_split!='') %>% .$feature
lipid_exact_FA <- species_t$lipid[species_t$lipid%in% lipid_exact_FA]

enrich_FDR <- cbind(FA_enrich_analysis(lipid_exact_FA,  species_t, 'up'),
                    FA_enrich_analysis(lipid_exact_FA,  species_t,'down')) %>% .[-3] %>% 
  gather(-FA, key = 'group', value='pvalue') %>% 
  mutate(mlogp=-log10(pvalue)) %>% 
  mutate(mlogp=ifelse(group=='pvalue_down',-mlogp, mlogp)) %>% 
  mutate(sig=ifelse(pvalue>0.05, 'no', 'yes')) %>% 
  arrange(desc(mlogp)) %>% 
  mutate(Data='Adjusted p-value\n(FDR)')


enrich_p <- cbind(FA_enrich_analysis(lipid_exact_FA,  species_t, 'up', 'p'),
                  FA_enrich_analysis(lipid_exact_FA, species_t,'down','p')) %>% .[-3] %>% 
  gather(-FA, key = 'group', value='pvalue') %>% 
  mutate(mlogp=-log10(pvalue)) %>% 
  mutate(mlogp=ifelse(group=='pvalue_down',-mlogp, mlogp)) %>% 
  mutate(sig=ifelse(pvalue>0.05, 'no', 'yes')) %>% 
  arrange(desc(mlogp)) %>% 
  mutate(Data='Original p-value')



F1c_data <- rbind(enrich_FDR, enrich_p) %>% 
  mutate(FA=str_replace(FA, ';0','')) %>% 
  mutate(color=ifelse(mlogp>0, 'Enriched','Depleted')) %>% 
  mutate(color=ifelse(pvalue<0.05, color,'Non-sig'))

#write.xlsx(F1c_data, file.path(file, 'source_data/F1c.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)

sig_enriched_feature <- F1c_data %>% filter(sig=='yes') %>% .$FA %>% unique()

dup <- F1c_data %>% arrange(pvalue) %>% .[c('FA','Data')] %>% duplicated()

F1c <- F1c_data %>% arrange(pvalue) %>% .[!dup,] %>%
  filter(FA %in% sig_enriched_feature) %>% 
  ggplot(aes(x=reorder(FA, mlogp), y=mlogp, fill=Data, color=color))+
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

#-------------------FA substructure analysis-------------------

#FA biosynthetic network data transformation

FA_network <- build_FA_net(FA_network, no_sub_t)


#Decompose lipids into FA substructures
#18:2 and 20:4 are majorly omega-6 FAs, so we only kept omega-6 forms of them
unmapped_FA <-  c('w9-18:2;0','w3-20:4;0')
FA_substructure <- FA_sub_transform(FA_network, no_sub_t,
                                    unmapped_FA =unmapped_FA)


#Extract FA substructure using fold changes


FA_t <- no_sub_t[[2]] %>% filter(type=='FA')

FA_sub_stop <- FA_sub_extract(char, FA_substructure, FA_t,
                              exact_FA='no', exo_lipid='w3-22:6;0')


#Transform lipid FA exp data into substructure exp data

FA_sub_exp <- lipid_sub_matrix(exp, FA_sub_stop, 'FA')

#Differential expression analysis for FA substructures

FA_sub_exp_t <- t_test(FA_sub_exp[[3]], 1:7, 8:13, 't.test', 'adj_p_value')

#Essential pathway analysis for FA substructures
set.seed(1)
path_score_FA <- path_scoring(FA_network, FA_sub_exp_t,
                              calibrate = T, data_type ='FA')


F3d_data <- rbind(path_score_FA %>% filter(Type=='Active') %>% 
                    .[!duplicated(.$rep_sub_path),] %>% .[1:5,], 
                  path_score_FA %>% filter(Type=='Suppressed') %>% 
                    arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>%.[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


F3d_data$path[c(1,3,4,5,8)] <- F3d_data$path[c(1,3,4,5,8)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], last(x),sep=' --> ')})

F3d_data$path[c(2,6,7)] <- F3d_data$path[c(2,6,7)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], x[2],last(x),sep=' --> ')})

#write.csv(F3d_data, file.path(file, 'source_data/F3d.csv'), row.names = F)

F3d <- F3d_data %>% 
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

#Essential edges (reactions) analysis for FA substructures

reaction_score_FA <- reaction_scoring(FA_network, 
                                      FA_sub_exp[[3]], FA_sub_exp_t, 
                                      ctrl = 1:7, exp = 8:13, 
                                      Species = 'rat')


F3e_data <-reaction_score_FA %>%
  filter(!is.na(perturbation_score)) %>% 
  filter(p_value<0.05) %>% 
  filter(!str_detect(edge_name, '22:6')) %>% 
  .[-4,] %>% 
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color))

#write.csv(F3e_data[-c(20:22)], file.path(file, 'source_data/F3e.csv'), row.names = F)


F3e <- F3e_data %>%
  mutate(Mode=factor(Mode, levels = c("Increase","Decrease"))) %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(F3e_data$edge_color)
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


#Construct FA biosynthetic network

FA_network_data <- draw_network(FA_network, FA_sub_exp_t,
                                if_species = 'no', significant = 'adj_p_value',
                                path_scoring_result = path_score_FA,
                                reaction_scoring_result = 
                                  filter(reaction_score_FA, !str_detect(edge_name, '22:6')),
                                top_n = 5)

F3b <- visNetwork(FA_network_data[[1]],FA_network_data[[2]]) %>% 
  visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                  physics = F, smooth = TRUE, randomSeed =5) 


F3b_data <- FA_network_data[[2]]%>% left_join(FA_sub_exp_t, by=c('from'='lipid')) %>% 
  left_join(FA_sub_exp_t, by=c('to'='lipid')) %>% 
  mutate(from=str_replace_all(from, ';0','')) %>% 
  mutate(to=str_replace_all(to, ';0','')) %>% 
  .[c(1:9,16,21,22,29,34,35)]

#write.xlsx(F3b_data, file.path(file, 'source_data/F3b.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)


#-------------------Raw FA analysis-------------------

#Essential pathway analysis for raw FA data

FA_mapping <- FA_substructure %>% apply(MARGIN = 1, FUN = function(x){c(x[1], last(x[x!='']))}) %>% 
  t() %>% as.data.frame() %>% unique() %>% 
  rbind(data.frame(FA=str_extract(unmapped_FA, '\\d+:\\d+;\\d+'),
                   V2=unmapped_FA))


FA_t2 <- FA_t %>% left_join(FA_mapping, by=c('lipid'='FA')) %>% 
  mutate(lipid=V2) %>% dplyr::select(-V2)

set.seed(7)

path_score_raw_FA <- path_scoring(FA_network, FA_t2, 
                                  calibrate = T, data_type = 'FA')


SF1a_data <- rbind(path_score_raw_FA %>% filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),] %>% .[1:5,], 
                   path_score_raw_FA %>% filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),] %>%.[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))


SF1a_data$path[c(2:5,7,8)] <- SF1a_data$path[c(2:5,7,8)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], last(x),sep=' --> ')})

SF1a_data$path[c(1,6)] <- SF1a_data$path[c(1,6)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], tail(x,2)[1], tail(x,2)[2],sep=' --> ')})

#write.csv(SF1a_data, file.path(file, 'source_data/SF1a.csv'), row.names = F)

SF1a <- SF1a_data %>% 
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

#Essential edges (reactions) analysis for raw FA data

FA_exp <- FA_mapping %>% 
  left_join(filter(no_sub_t[[1]], type=='FA'), by=c('FA'='feature')) %>% 
  filter(!is.na(Ctrl1)) %>% 
  column_to_rownames(var='V2') %>% 
  dplyr::select(-1,-2)



reaction_score_raw_FA <- reaction_scoring(FA_network, 
                                          FA_exp, FA_t2, 
                                          ctrl = 1:7, exp = 8:13, 
                                          Species = 'rat')


SF1b_data <-reaction_score_raw_FA %>%
  filter(!is.na(perturbation_score)) %>% 
  filter(p_value<0.05) %>% 
  filter(!str_detect(edge_name, '22:6')) %>% 
  .[c(1:3,6:10),] %>% 
  mutate(edge_name=str_replace_all(edge_name,  ';0','')) %>% 
  mutate(Edge_direction=ifelse(edge_type==Mode, 'Same as\nreaction', '')) %>% 
  mutate(node1=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[1]}),
         node2=str_split(.$edge_name, ' --> ') %>% map_chr(.f = function(x){x[2]})) %>% 
  mutate(node1_color=ifelse(node1_log2FC>0, paste0("<i style='color:#FF0000'>", node1, ' --> ',"</i>"),
                            paste0("<i style='color:#0000FF'>", node1, ' --> ',"</i>"))) %>% 
  mutate(node2_color=ifelse(node2_log2FC>0, paste0("<i style='color:#FF0000'>", node2, "</i>"),
                            paste0("<i style='color:#0000FF'>", node2, "</i>"))) %>% 
  mutate(edge_color=paste0(node1_color,node2_color))

#write.csv(SF1b_data[-c(20:22)], file.path(file, 'source_data/SF1b.csv'), row.names = F)


SF1b <- SF1b_data %>%
  mutate(Mode=factor(Mode, levels = c("Increase","Decrease"))) %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(SF1b_data$edge_color)
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

#Construct network for raw FA


FA_network_data_raw <- draw_network(FA_network, FA_t2,
                                    if_species = 'no', significant = 'adj_p_value',
                                    path_scoring_result = path_score_raw_FA,
                                    reaction_scoring_result = 
                                      filter(reaction_score_raw_FA, !str_detect(edge_name, '22:6')),
                                    top_n = 5)

F3a <- visNetwork(FA_network_data_raw[[1]],FA_network_data_raw[[2]]) %>% 
  visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                  physics = F, smooth = TRUE, randomSeed =5) 


F3a_data <- FA_network_data_raw[[2]]%>% left_join(FA_t2, by=c('from'='lipid')) %>% 
  left_join(FA_t2, by=c('to'='lipid')) %>% 
  mutate(from=str_replace_all(from, ';0','')) %>% 
  mutate(to=str_replace_all(to, ';0','')) %>% 
  .[c(1:9,17,22,23,31,36,37)]

#write.xlsx(F3a_data, file.path(file, 'source_data/F3a.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)

#-------------------FA substructures improve following analysis-------------------


#FA substructures improve statistical power
F3c_data <- FA_sub_exp_t %>% 
  mutate(simple=str_extract(lipid, '\\d+:\\d+;\\d+')) %>% 
  left_join(FA_t, by=c('simple'='lipid')) %>% 
  .[c('lipid', 'mlog10padj.y', 'mlog10padj.x')] %>% 
  `colnames<-`(c('lipid','Fatty acid','Substructure')) %>% 
  filter(!is.na(`Fatty acid`)) 

#write.csv(F3c_data, file.path(file, 'source_data/F3c.csv'), row.names = F)

options(scipen = -1)


F3c <- F3c_data %>% 
  ggpaired(cond1 = "Fatty acid", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = '-log10(padj)')+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed')+
  scale_y_continuous(limits = c(0,4))+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.2, label.y = 3.8)+
  theme(legend.position = 'none', axis.text.x = element_text(size=10))

#FA substructures improve hierarchical clustering



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
set.seed(1)
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

F3f_data <- as.data.frame(FA_exp_heatmap)


#write.csv(F3f_data, file.path(file, 'source_data/F3f.csv'))

rownames(FA_exp_heatmap) <- str_replace_all(rownames(FA_exp_heatmap),';0','')

F3f <- Heatmap(FA_exp_heatmap, clustering_method_rows = 'complete',
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

F3g_data <- as.data.frame(FA_sub_exp_heatmap)

#write.csv(F3g_data, file.path(file, 'source_data/F3g.csv'))

rownames(FA_sub_exp_heatmap) <- str_replace_all(rownames(FA_sub_exp_heatmap),';0','')
F3g <- Heatmap(FA_sub_exp_heatmap, clustering_method_rows = 'complete',
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


#DHA outlier
SF1d_data <- FA_sub_exp[[3]] %>% as.data.frame() %>% 
  .[str_detect(rownames(FA_sub_exp[[3]]), '22:6'),] %>% 
  gather(key='label',value='value') 

#write.csv(SF1d_data, file.path(file, 'source_data/SF1d.csv'), row.names = F)


SF1d <- SF1d_data %>% 
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

SF1c_data <-  FA_exp %>% as.data.frame() %>% 
  .[str_detect(rownames(FA_exp), '22:6'),] %>% 
  gather(key='label',value='value')

#write.csv(SF1c_data, file.path(file, 'source_data/SF1c.csv'), row.names = F)


SF1c <- SF1c_data %>% 
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

#FA substructures improve sample and network coverage

#sample coverage
FA_not_in_network <- FA_substructure$FA[!FA_substructure$FA %in% 
                                          str_extract(c(FA_network$S1, FA_network$P1), '\\d+:\\d+;\\d+')]

FA_raw_cov <- data.frame(sample_coverage=apply(no_sub_t[[1]] %>% filter(type=='FA', !feature %in% FA_not_in_network) %>% 
                                                 column_to_rownames(var='feature') %>% .[-1], MARGIN = 1, 
                                               FUN = function(x){sum(x!=0)})) %>% 
  mutate(Data='Fatty acid', feature=rownames(.)) %>% 
  dplyr::select(feature, everything())

FA_sub_cov <- data.frame(sample_coverage=apply(FA_sub_exp[[3]][!rownames(FA_sub_exp[[3]])%in%FA_not_in_network,], 
                                               MARGIN = 1, FUN = function(x){sum(x!=0)})) %>% 
  mutate(Data='Substructure', feature=rownames(.)) %>% 
  dplyr::select(feature, everything())


FA_raw_cov_avg <- mean(FA_raw_cov$sample_coverage)
FA_sub_cov_avg <- mean(FA_sub_cov$sample_coverage)


SF1e_data <- rbind(FA_raw_cov, FA_sub_cov)

#write.csv(SF1e_data, file.path(file, 'source_data/SF1e.csv'), row.names = F)

SF1e1 <- gghistogram(SF1e_data, 
                     x = "sample_coverage", 
                     rug = TRUE, alpha = 0.5,
                     fill = "Data", palette = c("#00AFBB", "#E7B800"),binwidth = 1)+
  scale_x_continuous(breaks = 4:13, labels = as.character(4:13))+
  labs(x='Sample coverage (1/13)')+
  annotate('text',x=8,y=25,
           label=str_c('Substrtucture improves ',
                       round((FA_sub_cov_avg-FA_raw_cov_avg)/FA_raw_cov_avg*100,1),'%'), 
           color="red", size=4)

SF1e2 <- ggdensity(SF1e_data,
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


SF1e <- align_plots(SF1e1, SF1e2, align="hv", axis="tblr")
SF1e <- ggdraw(SF1e[[1]]) + draw_plot(SF1e[[2]])


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

SF1f_data <- cbind(FA_raw_net_cov, FA_sub_net_cov) %>% `colnames<-`(c('Fatty acid','Substructure')) %>% 
  mutate(sample=rownames(.)) %>% 
  dplyr::select(sample, everything())

#write.csv(SF1f_data, file.path(file, 'source_data/SF1f.csv'), row.names = F)


SF1f <- SF1f_data %>% 
  ggpaired(cond1 = "Fatty acid", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = 'Network coverage',
           title = '')+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.16, label.y = 102)+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(labels = function(x) paste0(x,'%'))+
  annotate('text',x=1.5,y=108,label=str_c('Substrtucture improves ',improve_net_cov, '% (N=13)'), color="red", size=3.5)



#-------------------Lipid species substructure analysis-------------------

#Decompose lipids into species substructures
#We excluded ether lipids since we cannot differentiate Alkyl (O-) or Alkenyl- (P-) linked ether lipids

char_sel <- char[!str_detect(char$feature, 'O-'),]
exp_sel <- exp[!str_detect(exp$feature, 'O-'),]

species_substructure <- species_sub_transform(char_sel, lipid_substructure, 
                                              network_node)

#Extract species substructures using fold change

species_t <- no_sub_t[[2]] %>% filter(type=='species')

species_sub_stop <- species_sub_extract(species_substructure, species_t,
                                        'species', pct_limit = 0.3)

#Transform lipid exp into substructure exp

species_sub_exp <- lipid_sub_matrix(exp_sel, species_sub_stop, 'species')

#Differential expression analysis for substructures

species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:7, 8:13,
                            't.test', 'adj_p_value')


#Species biosynthetic network data transformation

species_network <- build_species_net(species_substructure)

#Essential pathway analysis for species substructures

set.seed(1)
path_score_species <-  path_scoring(species_network, species_sub_exp_t,
                                    calibrate = T, data_type = 'species')

SF2i_data <- rbind(path_score_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_species %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
  filter(Significant=='yes') %>% 
  mutate(path=str_replace_all(path, ';0',''))



SF2i_data$path[c(1:5,8:10)] <- SF2i_data$path[c(1:5,8:10)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1],last(x),sep=' --> ')})

SF2i_data$path[6:7] <- SF2i_data$path[6:7] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1], tail(x,2)[1],tail(x,2)[2],sep=' --> ')})


#write.csv(SF2i_data, file.path(file, 'source_data/SF2i.csv'), row.names = F)



SF2i <- SF2i_data %>% 
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


species_net_w_rev <- add_rev_rection(network_edge, species_network)


reaction_score_species <- reaction_scoring(species_net_w_rev,
                                           species_sub_exp[[3]],
                                           species_sub_exp_t,
                                           ctrl=1:7, exp=8:13,
                                           Species = 'rat')


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



SF2k_data <- reaction_score_species%>% 
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

#write.csv(SF2k_data[-c(20:22)], file.path(file, 'source_data/SF2k.csv'), row.names = F)

SF2k <- SF2k_data %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase','Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(SF2k_data$edge_color)
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
                                     top_n = 3)


F4f <- visNetwork(species_network_data[[1]], species_network_data[[2]])

species_sub_t_net <- species_sub_exp_t

maxinf <- ceiling(max(species_sub_t_net$log2FC[is.finite(species_sub_t_net$log2FC)]))
mininf <- floor(min(species_sub_t_net$log2FC[is.finite(species_sub_t_net$log2FC)]))
species_sub_t_net$log2FC[species_sub_t_net$log2FC>0 & is.infinite(species_sub_t_net$log2FC)] <- maxinf
species_sub_t_net$log2FC[species_sub_t_net$log2FC<0 & is.infinite(species_sub_t_net$log2FC)] <- mininf

F4f_data <- species_network_data[[2]]%>% left_join(species_sub_t_net, by=c('from'='lipid')) %>% 
  left_join(species_sub_t_net, by=c('to'='lipid')) %>% 
  mutate(from=str_replace_all(from, ';0','')) %>% 
  mutate(to=str_replace_all(to, ';0','')) %>% 
  .[c(1:9,16,21,22,29,34,35)]

F4f_data$font.size[str_c(F4f_data$from, ' --> ',F4f_data$to) %in% SF2j_data$edge_name] <- 40
F4f_data$label[str_c(F4f_data$from, ' --> ',F4f_data$to) %in% SF2j_data[SF2j_data$perturbation_score>0,]$edge_name] <- 'Increase'
F4f_data$label[str_c(F4f_data$from, ' --> ',F4f_data$to) %in% SF2j_data[SF2j_data$perturbation_score<0,]$edge_name] <- 'Decrease'
F4f_data$font.color[str_c(F4f_data$from, ' --> ',F4f_data$to) %in% SF2j_data[SF2j_data$perturbation_score>0,]$edge_name] <- 'red'
F4f_data$font.color[str_c(F4f_data$from, ' --> ',F4f_data$to) %in% SF2j_data[SF2j_data$perturbation_score<0,]$edge_name] <- 'blue'


#write.xlsx(F4f_data, file.path(file, 'source_data/F4f.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)



#BIOPAN's algorithm cannot build the connections involving FA transfer without free FA data provided
biopan_network_edge <- species_network_data[[2]] %>% filter(from %in% char$feature, to %in% char$feature)

biopan_network_edge <- biopan_network_edge %>% left_join(char[c('feature','class')], by=c('from'='feature')) %>% 
  left_join(char[c('feature','class')], by=c('to'='feature')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.x'='Abbreviation')) %>% 
  left_join(network_node[c('Abbreviation','FA')], by=c('class.y'='Abbreviation')) %>% 
  filter(FA.x==FA.y)

biopan_network_node <- species_network_data[[1]]  %>% filter(id %in% char$feature)

SF5a <- visNetwork(biopan_network_node, biopan_network_edge)

SF5a_data1 <- biopan_network_node %>% 
  left_join(species_sub_t_net, by=c('id'='lipid')) %>% 
  mutate(id=str_replace_all(id, ';0','')) %>% 
  mutate(label=str_replace_all(label, ';0','')) %>% 
  .[c(1:7,14,19,20)]

SF5a_data2 <- biopan_network_edge[1:9] %>% left_join(species_sub_t_net, by=c('from'='lipid')) %>% 
  left_join(species_sub_t_net, by=c('to'='lipid')) %>% 
  mutate(from=str_replace_all(from, ';0','')) %>% 
  mutate(to=str_replace_all(to, ';0','')) %>% 
  .[c(1:9,16,21,22,29,34,35)]

#write.xlsx(SF5a_data1, file.path(file, 'source_data/SF5a1.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)
#write.xlsx(SF5a_data2, file.path(file, 'source_data/SF5a2.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)

#-------------------Raw lipid species analysis-------------------

#Essential pathway analysis for raw species data

raw_species_net <- species_network

raw_species_net <- raw_species_net%>% filter(S1%in%species_t$lipid,
                                             P1%in%species_t$lipid)


set.seed(1)
path_score_raw_species <-  path_scoring(raw_species_net, species_t,
                                        calibrate = T, data_type = 'species')



SF2h_data <- rbind(path_score_raw_species %>% 
                     filter(Type=='Suppressed') %>% 
                     arrange(cal_score) %>% .[!duplicated(.$rep_sub_path),]%>% .[1:5,],
                   path_score_raw_species %>% 
                     filter(Type=='Active') %>% 
                     .[!duplicated(.$rep_sub_path),]%>% .[1:5,]) %>% 
  mutate(path=str_replace_all(path, ';0',''))



SF2h_data$path[c(1:10)] <- SF2h_data$path[c(1:10)] %>% 
  str_split(' --> ') %>% 
  map_chr(.f = function(x){str_c(x[1],last(x),sep=' --> ')})


#write.csv(SF2h_data, file.path(file, 'source_data/SF2h.csv'), row.names = F)

SF2h <- SF2h_data %>% 
  mutate(path=factor(.$path, levels = .$path)) %>% 
  ggplot(aes(x=reorder(path, cal_score), y=cal_score, fill=path))+
  geom_bar(stat='identity')+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-1.96,1.96), linetype='dashed',color='gray')+
  coord_flip()+
  theme_bw()+
  scale_fill_manual(values = c(bluered(100)[c(1,10,20,30)],'gray',bluered(100)[c(100,90,80,70,60)]))+
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
                                               ctrl=1:7, exp=8:13,
                                               Species = 'rat')


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



SF2j_data <- reaction_score_raw_species%>% 
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

#write.csv(SF2j_data[-c(20:22)], file.path(file, 'source_data/SF2j.csv'), row.names = F)

SF2j <- SF2j_data %>% 
  mutate(Edge_direction=factor(Edge_direction, levels = c("Same as\nreaction",""))) %>% 
  mutate(Mode=factor(Mode, levels = c('Increase','Decrease'))) %>% 
  ggplot(aes(x=perturbation_score, y=reorder(edge_name, perturbation_score), 
             fill=Mode, color=Edge_direction))+
  geom_bar(stat='identity', size=0.8)+
  scale_y_discrete(
    labels=rev(SF2j_data$edge_color)
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

#Construct network for raw species data


raw_species_network_data <- draw_network(raw_species_net_w_rev, species_t,
                                         if_species = T,significant = 'adj_p_value',
                                         path_scoring_result = path_score_raw_species,
                                         reaction_scoring_result = reaction_score_raw_species,
                                         top_n = 3)


F4e <- visNetwork(raw_species_network_data[[1]], raw_species_network_data[[2]])

species_raw_t_net <- species_t

maxinf <- ceiling(max(species_raw_t_net$log2FC[is.finite(species_raw_t_net$log2FC)]))
mininf <- floor(min(species_raw_t_net$log2FC[is.finite(species_raw_t_net$log2FC)]))
species_raw_t_net$log2FC[species_raw_t_net$log2FC>0 & is.infinite(species_raw_t_net$log2FC)] <- maxinf
species_raw_t_net$log2FC[species_raw_t_net$log2FC<0 & is.infinite(species_raw_t_net$log2FC)] <- mininf


F4e_data <- raw_species_network_data[[2]]%>% left_join(species_raw_t_net, by=c('from'='lipid')) %>% 
  left_join(species_raw_t_net, by=c('to'='lipid')) %>% 
  mutate(from=str_replace_all(from, ';0','')) %>% 
  mutate(to=str_replace_all(to, ';0','')) %>% 
  .[c(1:9,17,22,23,31,36,37)]

F4e_data$font.size[str_c(F4e_data$from, ' --> ',F4e_data$to) %in% SF2i_data$edge_name] <- 40
F4e_data$label[str_c(F4e_data$from, ' --> ',F4e_data$to) %in% SF2i_data[SF2i_data$perturbation_score>0,]$edge_name] <- 'Increase'
F4e_data$label[str_c(F4e_data$from, ' --> ',F4e_data$to) %in% SF2i_data[SF2i_data$perturbation_score<0,]$edge_name] <- 'Decrease'
F4e_data$font.color[str_c(F4e_data$from, ' --> ',F4e_data$to) %in% SF2i_data[SF2i_data$perturbation_score>0,]$edge_name] <- 'red'
F4e_data$font.color[str_c(F4e_data$from, ' --> ',F4e_data$to) %in% SF2i_data[SF2i_data$perturbation_score<0,]$edge_name] <- 'blue'

#write.xlsx(F4e_data, file.path(file, 'source_data/F4e.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)



#-------------------Species substructures improve following analysis-------------------

#Species substructures improve statistical power

F4a_SF2a_data <- species_sub_exp_t %>% 
  left_join(species_t, by=c('lipid')) %>% 
  .[c('lipid', 'mlog10padj.y', 'mlog10padj.x')] %>% 
  `colnames<-`(c('lipid','Lipid species','Substructure')) %>% 
  filter(`Lipid species`>-log10(0.05)| Substructure>-log10(0.05)) 

#write.csv(F4a_SF2a_data, file.path(file, 'source_data/F4a_SF2a.csv'), row.names = F, na = '')

F4a <- ggvenn(
  list(`Significant\nLipid species` = F4a_SF2a_data %>% filter(`Lipid species`> -log10(0.05)) %>% .$lipid, 
       `Significant\nSubstructure`= F4a_SF2a_data %>% filter(Substructure> -log10(0.05)) %>% .$lipid),
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4,
)

SF2a <- F4a_SF2a_data %>% 
  filter(!is.na(`Lipid species`)) %>% 
  ggpaired(cond1 = "Lipid species", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = '-log10(padj)')+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed')+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.24, label.y = 4.5)+
  theme(legend.position = 'none')


#Species substructures improve enrichment analysis
#remove Cardiolipins since data did not provide their exact FA


lipid_sub_wo_CL <- species_sub_exp_t$lipid[!str_detect(species_sub_exp_t$lipid, 'CL')]

enrich_sub <- cbind(FA_enrich_analysis(lipid_sub_wo_CL, 
                                       filter(species_sub_exp_t, lipid%in%lipid_sub_wo_CL),
                                       'up'),
                    FA_enrich_analysis(lipid_sub_wo_CL,
                                       filter(species_sub_exp_t, lipid%in%lipid_sub_wo_CL),
                                       'down')) %>% 
  .[-3] %>% 
  gather(-FA, key = 'group', value='pvalue') %>% 
  mutate(mlogp=-log10(pvalue)) %>% 
  mutate(mlogp=ifelse(group=='pvalue_down',-mlogp, mlogp)) %>% 
  mutate(sig=ifelse(pvalue>0.05, 'no', 'yes')) %>% 
  arrange(desc(mlogp)) %>% 
  mutate(Data='Substructure')

lipid_raw_wo_CL <- species_t$lipid[!str_detect(species_t$lipid, '(CL)|(O-)')]



enrich_raw <- cbind(FA_enrich_analysis(lipid_raw_wo_CL, 
                                       filter(species_t, lipid%in%lipid_raw_wo_CL),
                                       'up'),
                    FA_enrich_analysis(lipid_raw_wo_CL,
                                       filter(species_t, lipid%in%lipid_raw_wo_CL),
                                       'down')) %>% 
  .[-3] %>% 
  gather(-FA, key = 'group', value='pvalue') %>% 
  mutate(mlogp=-log10(pvalue)) %>% 
  mutate(mlogp=ifelse(group=='pvalue_down',-mlogp, mlogp)) %>% 
  mutate(sig=ifelse(pvalue>0.05, 'no', 'yes')) %>% 
  arrange(desc(mlogp)) %>% 
  mutate(Data='Lipid species')


F4b_data <- rbind(enrich_sub, enrich_raw) %>% 
  mutate(FA=str_replace(FA, ';0','')) %>% 
  mutate(color=ifelse(mlogp>0, 'Enriched','Depleted')) %>% 
  mutate(color=ifelse(pvalue<0.05, color,'Non-sig'))

#write.xlsx(F4b_data, file.path(file, 'source_data/F4b.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)

sig_enriched_feature <- F4b_data %>% filter(sig=='yes') %>% .$FA %>% unique()

dup <- F4b_data %>% arrange(pvalue) %>% .[c('FA','Data')] %>% duplicated()

F4b <- F4b_data %>% arrange(pvalue) %>% .[!dup,] %>%
  filter(FA %in% sig_enriched_feature) %>% 
  ggplot(aes(x=reorder(FA, mlogp), y=mlogp, fill=Data, color=color))+
  geom_bar(stat='identity',position = position_dodge(), size=0.7)+
  coord_flip()+
  geom_hline(yintercept = c(-log10(0.05),log10(0.05)), linetype='dashed',color='black')+
  geom_hline(yintercept = 0, color='black')+
  theme_bw()+
  theme(legend.position = 'right',
        plot.title = element_text(hjust = 0.5,size = 10))+
  scale_fill_manual(values=c('#E69F00','#999999'), breaks = c('Substructure','Lipid species' ))+
  scale_color_manual(values =c('#EE0000FF','blue','white'), breaks = c('Enriched','Depleted','Non-sig'))+
  scale_y_continuous(breaks = c(-5,0,5,10),labels = c('5','0','5','10'))+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))+
  labs(y='-log10(p-value)', x='Fatty acid enrichment', color='For DHA group')


#Species substructures improve hierarchical clustering

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

SF2f_data <- as.data.frame(species_exp_heatmap)
#write.csv(SF2f_data, file.path(file, 'source_data/SF2f.csv'))


SF2f <- Heatmap(species_exp_heatmap, clustering_method_rows = 'complete',
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

SF2g_data <- as.data.frame(species_sub_exp_heatmap)
#write.csv(SF2g_data, file.path(file, 'source_data/SF2g.csv'))


SF2g <- Heatmap(species_sub_exp_heatmap, clustering_method_rows = 'complete',
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

#Species substructures improve sample and network coverage

#sample coverage

species_raw_cov <- data.frame(sample_coverage=apply(species_exp_raw, MARGIN = 1, FUN = function(x){sum(x!=0)})) %>% 
  mutate(Data='Lipid species', feature=rownames(.)) %>% 
  dplyr::select(feature, everything())

species_sub_cov <- data.frame(sample_coverage=apply(species_sub_exp[[3]], MARGIN = 1, FUN = function(x){sum(x!=0)})) %>% 
  mutate(Data='Substructure', feature=rownames(.)) %>% 
  dplyr::select(feature, everything())

species_raw_cov_avg <- mean(species_raw_cov$sample_coverage)
species_sub_cov_avg <- mean(species_sub_cov$sample_coverage)


F4c_data <- rbind(species_raw_cov, species_sub_cov)

#write.csv(F4c_data, file.path(file, 'source_data/F4c.csv'), row.names = F)


F4c1 <- gghistogram(F4c_data, 
                    x = "sample_coverage", 
                    rug = TRUE, alpha = 0.5,
                    fill = "Data", palette = c("#00AFBB", "#E7B800"),binwidth = 1)+
  scale_x_continuous(breaks = 4:13, labels = as.character(4:13))+
  labs(x='Sample coverage (1/13)')+
  annotate('text',x=8,y=250,
           label=str_c('Substrtucture improves ',
                       round((species_sub_cov_avg-species_raw_cov_avg)/species_raw_cov_avg*100,1),'%'), 
           color="red", size=4)



F4c2 <- ggdensity(F4c_data, 
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


F4c <- align_plots(F4c1, F4c2, align="hv", axis="tblr")
F4c <- ggdraw(F4c[[1]]) + draw_plot(F4c[[2]])



#network coverage

species_network_node <- unlist(species_network) %>% unique()

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

F4d_data <- cbind(species_raw_net_cov, species_sub_net_cov) %>% 
  `colnames<-`(c('Lipid species','Substructure')) %>% 
  mutate(sample=rownames(.)) %>% 
  dplyr::select(sample, everything())

#write.csv(F4d_data, file.path(file, 'source_data/F4d.csv'), row.names = F)
options(scipen = -1)
F4d <- F4d_data %>% 
  ggpaired(cond1 = "Lipid species", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = 'Network coverage',
           title = '')+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.16, label.y = 80)+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(labels = function(x) paste0(x,'%'))+
  annotate('text',x=1.5,y=90,label=str_c('Substrtucture improves ',improve_net_cov, '% (N=13)'), color="red", size=3.5)


#-------------------Save data-------------------
#save.image(file.path(file, 'DHA.RData'))
