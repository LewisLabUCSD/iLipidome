library(ggsci)
library(ggpubr)
library(ggvenn)
library(gridExtra)
library(ggrepel)
library(ggtext)
library(fgsea)

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
                    replace_NA = T,NA2what = 'min',ymin = 0.5, pct_transform = T, data_transform = F)

exp <- remove_rownames(exp)


each_FA <- str_extract_all(char$FA_split, '\\d+:\\d+;\\d+') %>% map(.f = function(x){x[x!='0:0;0']})

char <- char %>% mutate(each_FA=each_FA)

#-------------------Analysis for unprocessed data-------------------

no_sub_t <- unprocessed_data_test(exp, char, 'mod.t.test', 'adj_p_value', 1:7,  8:13)

#-------------------Conventional lipidomics analysis-------------------

#Differential expression for lipid species
DE_lipid_data <- no_sub_t[[2]] %>% filter(type=='species') %>%  
  mutate(log2FC=ifelse(is.infinite(log2FC),10*sign(log2FC),log2FC)) %>% 
  mutate(Significance=ifelse(log2FC>0, 'Increase','Decrease')) %>% 
  mutate(Significance=ifelse(sig=='yes', Significance, 'No change'))


DE_lipid_data %>% 
  ggplot(aes(x=log2FC, y=mlog10padj, col=Significance)) +
  geom_point() + 
  scale_x_continuous(limits = c(-7,7))+
  scale_color_manual(values=c("blue", "red", "gray")) +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='dashed')+
  theme_classic()+
  theme(legend.position = 'top')+
  labs(y='-log10(padj)')

#-------------------Lipid species substructure analysis-------------------

#Decompose lipids into species substructures
#We excluded ether lipids since we cannot differentiate Alkyl (O-) or Alkenyl- (P-) linked ether lipids

char_sel <- char[!str_detect(char$feature, 'O-'),]
exp_sel <- exp[!str_detect(exp$feature, 'O-'),]

species_substructure <- species_sub_transform(char_sel, lipid_substructure, 
                                              network_node)

#Extract species substructures using fold changes

species_t <- no_sub_t[[2]] %>% filter(type=='species')

species_sub_stop <- species_sub_extract(species_substructure, species_t, 'species', pct_limit = 0.3)

#Transform lipid exp into substructure exp

species_sub_exp <- lipid_sub_matrix(exp_sel, species_sub_stop, 'species')

#Differential expression analysis for substructures

species_sub_exp_t <- t_test(species_sub_exp[[3]], 1:7, 8:13, 'mod.t.test', 'adj_p_value')

#-------------------Species substructures improve following analysis-------------------

#Species substructures improve statistical power

SF2b_SF2c_data <- species_sub_exp_t %>% 
  left_join(species_t, by=c('lipid')) %>% 
  .[c('lipid', 'mlog10padj.y', 'mlog10padj.x')] %>% 
  `colnames<-`(c('lipid','Lipid species','Substructure')) %>% 
  filter(`Lipid species`>-log10(0.05)| Substructure>-log10(0.05)) 

#write.csv(SF2b_SF2c_data, file.path(file, 'source_data/SF2b_SF2c_data.csv'), row.names = F, na = '')

SF2b <- ggvenn(
  list(`Significant\nLipid species` = SF2b_SF2c_data %>% filter(`Lipid species`> -log10(0.05)) %>% .$lipid, 
       `Significant\nSubstructure`= SF2b_SF2c_data %>% filter(Substructure> -log10(0.05)) %>% .$lipid),
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4,
)
options(scipen = -1)
SF2c <- SF2b_SF2c_data %>% 
  filter(!is.na(`Lipid species`)) %>% 
  ggpaired(cond1 = "Lipid species", cond2 = "Substructure", 
           fill = "condition", line.color = "gray", line.size = 0.4,
           palette = "jco",xlab = '', ylab = '-log10(padj)')+
  geom_hline(yintercept = -log10(0.05), color='red', linetype='dashed')+
  stat_compare_means(paired = T, method = 't.test',label.x = 1.24, label.y = 4.5)+
  theme(legend.position = 'none')


#Species substructures improve enrichment analysis

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


SF2d_data <- rbind(enrich_sub, enrich_raw) %>% 
  mutate(FA=str_replace(FA, ';0','')) %>% 
  mutate(color=ifelse(mlogp>0, 'Enriched','Depleted')) %>% 
  mutate(color=ifelse(pvalue<0.05, color,'Non-sig'))

#write.xlsx(SF2d_data, file.path(file, 'source_data/SF2d.xlsx'), sheetName = "Sheet1", col.names = T, row.names = F, append = F,showNA = F)

sig_enriched_feature <- SF2d_data %>% filter(sig=='yes') %>% .$FA %>% unique()

dup <- SF2d_data %>% arrange(pvalue) %>% .[c('FA','Data')] %>% duplicated()

SF2d <- SF2d_data %>% arrange(pvalue) %>% .[!dup,] %>%
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

#Species substructures improve LSEA


lipidset_sub <- data.frame(feature=rep(lipid_sub_wo_CL, (map_int(str_split(lipid_sub_wo_CL, '_'), ~length(.x))-1)),
                       class=unlist(str_extract_all(lipid_sub_wo_CL, '\\d+:\\d+;\\d+'))) %>% 
  group_by(class) %>% 
  dplyr::summarise(lipid=list(feature))

lipidset_sub <- lipidset_sub$lipid %>% `names<-`(lipidset_sub$class)

lipid.stat_sub <- filter(species_sub_exp_t, lipid%in%lipid_sub_wo_CL) %>% 
  mutate(stat = statistics) %>%
  arrange(desc(stat))

lipidrank_sub <- lipid.stat_sub$stat
names(lipidrank_sub) <- lipid.stat_sub$lipid

lipid.fgseaRes_sub <- fgsea(pathways = lipidset_sub, stats = lipidrank_sub)



lipidset_raw <- data.frame(feature=rep(lipid_raw_wo_CL, (map_int(str_split(lipid_raw_wo_CL, '_'), ~length(.x))-1)),
                       class=unlist(str_extract_all(lipid_raw_wo_CL, '\\d+:\\d+;\\d+'))) %>% 
  group_by(class) %>% 
  dplyr::summarise(lipid=list(feature))
lipidset_raw <- lipidset_raw$lipid %>% `names<-`(lipidset_raw$class)

lipid.stat_raw <- filter(species_t, lipid%in%lipid_raw_wo_CL) %>% 
  mutate(stat = statistics) %>%
  arrange(desc(stat))

lipidrank_raw <- lipid.stat_raw$stat
names(lipidrank_raw) <- lipid.stat_raw$lipid

lipid.fgseaRes_raw <- fgsea(pathways = lipidset_raw, stats = lipidrank_raw)

SF2e_data <- rbind(mutate(lipid.fgseaRes_raw, Data='Lipid species'),
      mutate(lipid.fgseaRes_sub, Data='Substructure')) %>% 
  filter(pathway %in%c(filter(lipid.fgseaRes_sub, padj<0.05)$pathway, filter(lipid.fgseaRes_raw, padj<0.05)$pathway)) %>% 
  mutate(FA=str_replace(pathway, ';0','')) %>% 
  group_by(Data) %>% 
  mutate(color=ifelse(NES>0, 'Enriched','Depleted')) %>% 
  mutate(color=ifelse(padj<0.05, color,'Non-sig')) 

#write.csv(SF2e_data[-8], file.path(file, 'source_data/SF2e.csv'), row.names = F)


SF2e <- SF2e_data%>% 
  ggplot(aes(x=reorder(FA, NES), y=NES, fill=Data, color=color))+
  geom_bar(stat='identity',position = position_dodge(), size=0.7)+
  coord_flip()+
  geom_hline(yintercept = 0, color='black')+
  theme_bw()+
  theme(legend.position = 'right',plot.title = element_text(hjust = 0.5,size = 10))+
  scale_fill_manual(values=c('#E69F00','#999999'), breaks = c('Substructure','Lipid species' ))+
  scale_color_manual(values =c('#EE0000FF','blue','white'), breaks = c('Enriched','Depleted','Non-sig'))+
  guides(fill=guide_legend(order=1),
         color=guide_legend(order=2))+
  labs(y='NES', x='Fatty acid enrichment', color='For DHA group')



plotEnrichment(lipidset_sub[[lipid.fgseaRes_sub$pathway[28]]], lipidrank_sub) + 
  labs(title=paste0(lipid.fgseaRes_sub$pathway[28], '\npvalue = ', 
                    round(lipid.fgseaRes_sub$pval[which(lipid.fgseaRes_sub$pathway == lipid.fgseaRes_sub$pathway[28])], digits = 10), 
                    ', NES = ', round(lipid.fgseaRes_sub$NES[which(lipid.fgseaRes_sub$pathway == lipid.fgseaRes_sub$pathway[28])], digits = 3))) +
  theme(text = element_text(size = 15))

#-------------------Save data-------------------
#save.image(file.path(file, 'DHA_mod_t_lsea.RData'))
