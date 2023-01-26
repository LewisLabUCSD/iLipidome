library(tidyverse)
library(xlsx)
library(biomaRt)
library(KEGGREST)
library(data.table)

file <- dirname(rstudioapi::getSourceEditorContext()$path)

#biomart data
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ENSG_human <- getBM(attributes=c('ensembl_gene_id','entrezgene_id',
                                 'uniprotswissprot','external_gene_name'),
                    mart = ensembl) %>% 
  mutate(entrezgene_id=as.character(entrezgene_id))

ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ENSG_mouse <- getBM(attributes=c('ensembl_gene_id','entrezgene_id',
                                 'uniprotswissprot','external_gene_name'),
                    mart = ensembl) %>% 
  mutate(entrezgene_id=as.character(entrezgene_id))


ensembl = useEnsembl(biomart="ensembl", dataset="rnorvegicus_gene_ensembl")
ENSG_rat <- getBM(attributes=c('ensembl_gene_id','entrezgene_id',
                               'uniprotswissprot','external_gene_name'),
                  mart = ensembl) %>% 
  mutate(entrezgene_id=as.character(entrezgene_id))

#RHEA data

rhea_SwissProt <- read.delim(file.path(file,'rhea2uniprot_sprot.tsv'),header = T)

#Lipid class reaction gene

get_reaction_gene <- function(Reaction_data, rhea_data, ENSG_mapping_table, species='human'){
  if(species=='human'){
    ec_id_prefix <- 'HSA: '
  }
  else if(species=='mouse'){
    ec_id_prefix <- 'MMU: '
  }
  else if(species=='rat'){
    ec_id_prefix <- 'RNO: '
  }
  
  gene_id <- vector(mode = 'list', length = nrow(Reaction_data))
  from <- character(nrow(Reaction_data))
  to <- character(nrow(Reaction_data))
  
  for(num in 1:nrow(Reaction_data)){
    from[num] <- Reaction_data$iLipidome_from[num]
    to[num] <- Reaction_data$iLipidome_to[num]
    
    #mapping ec number
    if(Reaction_data$EC.numbers[num]!=''){
      ec_id <- str_c('ec:', Reaction_data$EC.numbers[num])
      
      ec_id <- tryCatch({keggGet(ec_id)[[1]]$GENES},
                        error=function(e){return(NULL)})
      if(is.null(ec_id)){
        ec_id <- data.frame()
        print(num)
      }
      else{
        ec_id_loc <- which(str_detect(ec_id, ec_id_prefix))
        ec_id_list <- str_extract_all(ec_id[ec_id_loc], '\\d+\\(') %>% 
          unlist() %>% str_sub(end = -2)
        ec_id <- ENSG_mapping_table %>% filter(entrezgene_id %in% ec_id_list) %>% 
          dplyr::select(ensembl_gene_id, external_gene_name)
      }
    }
    else{
      ec_id <- data.frame()
      
    }
    
    #mapping rhea id
    if(Reaction_data$RHEA[num]!=''){
      rhea_id <-  str_extract(Reaction_data$RHEA[num], '\\d+')
      rhea_id <- rhea_data %>% filter(RHEA_ID %in% as.integer(rhea_id)) %>%
        filter(ID %in% ENSG_mapping_table$uniprotswissprot) %>% 
        left_join(ENSG_mapping_table, by=c('ID'='uniprotswissprot')) %>% 
        dplyr::select(ensembl_gene_id, external_gene_name)
    }
    else{
      rhea_id <- data.frame()
    }
    
    #mapping entrez id
    
    if(Reaction_data$Entrez[num]!=''){
      Entrez_id <- str_split(Reaction_data$Entrez[num], '_') %>% unlist()
      Entrez_id <- ENSG_mapping_table %>% filter(entrezgene_id %in% Entrez_id) %>% 
        dplyr::select(ensembl_gene_id, external_gene_name)
    }
    else{
      Entrez_id <- data.frame()
    }
    gene_id[[num]] <- rbind(rhea_id, Entrez_id, ec_id) %>% unique()
    print(num)
  }
  reaction_gene <- data.frame(From=rep(from, map_int(gene_id, ~nrow(.x))),
                              To=rep(to, map_int(gene_id, ~nrow(.x)))) %>% 
    cbind(rbindlist(gene_id)) %>% unique()
  return(reaction_gene)
}

Lipid_net <- read.xlsx2(file.path(file, 'lipid_reaction_gene.xlsx'), sheetIndex = 3)

reaction_gene_human <- get_reaction_gene(Reaction_data = Lipid_net, rhea_data = rhea_SwissProt,
                                         ENSG_mapping_table = ENSG_human,  species='human')

Lipid_net <- read.xlsx2(file.path(file, 'lipid_reaction_gene.xlsx'), sheetIndex = 4)

reaction_gene_mouse <- get_reaction_gene(Reaction_data = Lipid_net, rhea_data = rhea_SwissProt,
                                         ENSG_mapping_table = ENSG_mouse,  species='mouse')

Lipid_net <- read.xlsx2(file.path(file, 'lipid_reaction_gene.xlsx'), sheetIndex = 5)

reaction_gene_rat <- get_reaction_gene(Reaction_data = Lipid_net, rhea_data = rhea_SwissProt,
                                       ENSG_mapping_table = ENSG_rat,  species='rat')


reaction_gene <- rbind(reaction_gene_human %>% mutate(species='human'),
      reaction_gene_mouse %>% mutate(species='mouse'),
      reaction_gene_rat %>% mutate(species='rat'))

#write.csv(reaction_gene, file.path(file, 'reaction_gene.csv'), row.names = F)

#FA reaction gene

FA_network <- read.csv(file.path(file, 'FA_network.csv'), header = T)
FA_net <- read.xlsx2(file.path(file, 'lipid_reaction_gene.xlsx'), sheetIndex = 6)

FA_reaction_gene_human <- FA_net %>% map(~filter(ENSG_human, entrezgene_id %in% as.character(.x)) %>% 
                                           dplyr::select(ensembl_gene_id, external_gene_name) %>% unique)

FA_reaction_gene_human <- data.frame(From=rep(FA_network$S1, map_int(FA_reaction_gene_human, ~nrow(.))),
                                     To=rep(FA_network$P1, map_int(FA_reaction_gene_human, ~nrow(.)))) %>% 
  cbind(rbindlist(FA_reaction_gene_human))

FA_net <- read.xlsx2(file.path(file, 'lipid_reaction_gene.xlsx'), sheetIndex = 7)

FA_reaction_gene_mouse <- FA_net %>% map(~filter(ENSG_mouse, entrezgene_id %in% as.character(.x)) %>% 
                                           dplyr::select(ensembl_gene_id, external_gene_name) %>% unique)

FA_reaction_gene_mouse <- data.frame(From=rep(FA_network$S1, map_int(FA_reaction_gene_mouse, ~nrow(.))),
                                     To=rep(FA_network$P1, map_int(FA_reaction_gene_mouse, ~nrow(.)))) %>% 
  cbind(rbindlist(FA_reaction_gene_mouse))

FA_net <- read.xlsx2(file.path(file, 'lipid_reaction_gene.xlsx'), sheetIndex = 8)

FA_reaction_gene_rat <- FA_net %>% map(~filter(ENSG_rat, entrezgene_id %in% as.character(.x)) %>% 
                                         dplyr::select(ensembl_gene_id, external_gene_name) %>% unique)

FA_reaction_gene_rat <- data.frame(From=rep(FA_network$S1, map_int(FA_reaction_gene_rat, ~nrow(.))),
                                   To=rep(FA_network$P1, map_int(FA_reaction_gene_rat, ~nrow(.)))) %>% 
  cbind(rbindlist(FA_reaction_gene_rat))

FA_reaction_gene <- rbind(FA_reaction_gene_human %>% mutate(species='human'),
      FA_reaction_gene_mouse %>% mutate(species='mouse'),
      FA_reaction_gene_rat %>% mutate(species='rat'))

#write.csv(FA_reaction_gene, file.path(file, 'FA_reaction_gene.csv'), row.names = F)

reaction_gene_mapping <- 
  rbind(FA_reaction_gene, reaction_gene) %>%
  mutate(Reaction=str_c(From, '_',To)) %>% 
  dplyr::select(Reaction, everything()) %>% 
  group_by(Reaction, species) %>% 
  summarise(gene=str_c(unique(sort(external_gene_name)), collapse = ', '))

#write.csv(reaction_gene_mapping, file.path(file, 'reaction_gene_mapping.csv'), row.names = F)

#-------------------Save data-------------------

#save.image(file = file.path(file, 'Lipid_reaction_gene.RData'))
