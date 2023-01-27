#-------------------Source function and required data-------------------

file <- dirname(rstudioapi::getSourceEditorContext()$path)

# Source function
suppressWarnings(suppressPackageStartupMessages(source(file.path(file,'required_function.R'))))

# Load required data
load(file.path(file,'required_data.RData'))

#-------------------Data upload and process-------------------


exp <- read.csv(file.path(file, 'exp.csv'))

exp_sel <- build_char_table(exp, network_node = network_node)[[1]]

char_sel <- build_char_table(exp, network_node = network_node)[[2]]

#-------------------Analysis for unprocessed data-------------------

no_sub_t <- unprocessed_data_test(exp_data = exp_sel,
                                  char_table = char_sel,
                                  method = 't.test',
                                  significant='adj_p_value',
                                  ctrl_group = 1:7, exp_group = 8:13)

#-------------------FA substructure analysis-------------------

#FA biosynthetic network data transformation

FA_network_new <- build_FA_net(FA_network = FA_network,
                           unprocessed_data_result = no_sub_t)
FA_network_new <- build_FA_net(FA_network = FA_network,
                               unprocessed_data_result = no_sub_t)
#Decompose lipids into FA substructures
#18:2 and 20:4 are majorly omega-6 FAs, so we only kept omega-6 forms of them

FA_substructure <- FA_sub_transform(FA_network = FA_network_new,
                                    unprocessed_data_result = no_sub_t,
                                    unmapped_FA = c('w9-18:2;0','w3-20:4;0'))

#Extract FA substructures using fold changes

FA_sub_stop <- FA_sub_extract(char_table = char_sel,
                              FA_substructure = FA_substructure,
                              unprocessed_data_result = no_sub_t,
                              exact_FA='no', exo_lipid='w3-22:6;0')

#Transform FA exp into substructure exp

FA_sub_exp <- lipid_sub_matrix(exp_data = exp_sel, sub_data = FA_sub_stop,
                               sub_type = 'FA')


#Differential expression analysis for FA substructures

FA_sub_exp_t <- t_test(data = FA_sub_exp[[3]], ctrl = 1:7, exp = 8:13,
                       method = 't.test', significant = 'adj_p_value')


#Essential pathway analysis for FA substructures

set.seed(1)

path_score_FA <- path_scoring(network = FA_network_new, sub_t = FA_sub_exp_t, 
                              calibrate = T, data_type = 'FA')

#Essential edges (reactions) analysis for FA substructures

reaction_score_FA <- reaction_scoring(network = FA_network_new, 
                                      sub_exp = FA_sub_exp[[3]],
                                      sub_t = FA_sub_exp_t, 
                                      ctrl = 1:7, exp = 8:13, 
                                      Species = 'rat')


#FA biosynthetic network construction
FA_network_data <- draw_network(network_data = FA_network_new,
                                DE_data = FA_sub_exp_t,
                                if_species = F, significant = 'adj_p_value',
                                path_scoring_result = path_score_FA,
                                reaction_scoring_result = reaction_score_FA,
                                top_n = 5, path_type = 'both')


visNetwork(FA_network_data[[1]],FA_network_data[[2]]) %>% 
  visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                  physics = F, smooth = TRUE, randomSeed =5) 


#-------------------Lipid species substructure analysis-------------------

#We excluded ether lipids since we cannot differentiate Alkyl (O-) or Alkenyl- (P-) linked ether lipids

char_wo_EL <- char_sel[!str_detect(char_sel$feature, 'O-'),]
exp_wo_EL <- exp_sel[!str_detect(exp_sel$feature, 'O-'),]

#Decompose lipids into species substructures

species_substructure <- species_sub_transform(char = char_wo_EL,
                                              lipid_substructure = lipid_substructure,
                                              network_node = network_node)

#Extract species substructures using fold changes

species_sub_stop <- species_sub_extract(lipid_substructure = species_substructure,
                                        unprocessed_data_result =  no_sub_t,
                                        type = 'species', pct_limit = 0.3,
                                        exo_lipid=NULL)

#Transform lipid exp into substructure exp

species_sub_exp <- lipid_sub_matrix(exp_data = exp_wo_EL, 
                                    sub_data = species_sub_stop,
                                    sub_type = 'species')


#Differential expression analysis for substructures

species_sub_exp_t <- t_test(data = species_sub_exp[[3]], ctrl = 1:7, exp = 8:13,
                            method = 't.test', significant = 'adj_p_value')


#Species biosynthetic network data transformation

species_network <- build_species_net(species_substructure = species_substructure)


#Essential pathway analysis for species substructures

set.seed(1)
path_score_species <-  path_scoring(network = species_network,
                                    sub_t = species_sub_exp_t,
                                    calibrate = T, data_type = 'Species')


#Essential edges (reactions) analysis for species substructures

species_net_w_rev <- add_rev_rection(network_edge = network_edge,
                                     species_net = species_network)

reaction_score_species <- reaction_scoring(network = species_net_w_rev,
                                           sub_exp = species_sub_exp[[3]],
                                           sub_t = species_sub_exp_t,
                                           ctrl=1:7, exp=8:13,
                                           Species = 'rat')

#Lipid species biosynthetic network construction

species_network_data <- draw_network(network_data = species_net_w_rev,
                                     DE_data = species_sub_exp_t,
                                     if_species = T,significant = 'adj_p_value',
                                     path_scoring_result = path_score_species,
                                     reaction_scoring_result = reaction_score_species,
                                     top_n = 3, path_type = 'both')


visNetwork(species_network_data[[1]], species_network_data[[2]])


#-------------------Lipid class substructure analysis-------------------

#Extract class substructures using fold changes

class_sub_stop <- species_sub_extract(lipid_substructure =lipid_substructure,
                                      unprocessed_data_result = no_sub_t,
                                      type = 'class', pct_limit = 0.01,
                                      exo_lipid=NULL)

#Transform lipid exp into substructures exp

class_exp <- no_sub_t[[1]] %>% filter(type=='class') %>% 
  dplyr::select(-type)


class_sub_exp <- lipid_sub_matrix(exp_data = class_exp, 
                                  sub_data = class_sub_stop,
                                  sub_type = 'Class')



#Differential expression analysis for substructures

class_sub_exp_t <- t_test(data = class_sub_exp[[3]], ctrl = 1:7, exp = 8:13,
                          method = 't.test', significant = 'adj_p_value')


#Class biosynthetic network data transformation

class_network <- network_edge[c('S1','P1')] %>% 
  filter(S1 %in% class_sub_exp_t$lipid, P1 %in% class_sub_exp_t$lipid)

#Essential pathway analysis for class substructures

set.seed(1)
path_score_class <-  path_scoring(network = class_network,
                                  sub_t = class_sub_exp_t,
                                  calibrate = T, data_type = 'Class')

#Essential edges (reactions) analysis for class substructures

reaction_score_class <- reaction_scoring(network = class_network,
                                         sub_exp = class_sub_exp[[3]],
                                         sub_t = class_sub_exp_t,
                                         ctrl=1:7, exp=8:13,
                                         Species = 'rat')

#Lipid class biosynthetic network construction

class_network_data <- draw_network(network_data = class_network,
                                   DE_data = class_sub_exp_t,
                                   if_species = F,significant = 'adj_p_value',
                                   path_scoring_result = path_score_class,
                                   reaction_scoring_result = reaction_score_class,
                                   top_n = 3, path_type = 'both')

visNetwork(class_network_data[[1]], class_network_data[[2]])

#-------------------Save data-------------------
#save.image(file.path(file, 'quick_example.RData'))

