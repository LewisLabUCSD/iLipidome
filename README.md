# iLipidome

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Quick Example](#quick-example)
- [License](#license)

# Overview
Here, we present ``iLipidome``, a method for analyzing lipidomics data in the context of the lipid biosynthetic network, thus accounting for the interdependence of measured lipids. Currently, iLipidome only supports “two-group comparison”, enabling users to identify essential altered lipid pathways and link lipidomic changes to their genetic origins. The tutorial describes a series of iLipidome functions to facilitate systems-level comparison of fatty acid, lipid species, and lipid class profiles using a novel substructure-based approach. We hope it can provide researchers a deeper insight into complex lipidomic alterations across samples.


# System Requirements
## Hardware requirements
To run example datasets with iLipidome, it requires only a standard computer with enough RAM and installed R software version over 4.0.0.

## Software requirements
### OS Requirements
The functions and example datasets has been tested on the following systems:
+ macOS: Ventura (13.0)


### R Dependencies
The version information about R, the OS and attached or loaded packages for `iLipidome` are listed below.

```
R version 4.1.0 (2021-05-18)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 13.0

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] zh_TW.UTF-8/zh_TW.UTF-8/zh_TW.UTF-8/C/zh_TW.UTF-8/zh_TW.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gtools_3.9.2      gplots_3.1.1      MKmisc_1.8        data.table_1.14.0
 [5] xlsx_0.6.5        visNetwork_2.1.0  igraph_1.2.6      forcats_0.5.1    
 [9] stringr_1.4.0     dplyr_1.0.7       purrr_0.3.4       readr_1.4.0      
[13] tidyr_1.1.3       tibble_3.1.2      ggplot2_3.3.5     tidyverse_1.3.1  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.8.3       lubridate_1.7.10   xlsxjars_0.6.1     assertthat_0.2.1  
 [5] digest_0.6.27      utf8_1.2.1         R6_2.5.0           cellranger_1.1.0  
 [9] backports_1.2.1    reprex_2.0.0       httr_1.4.2         pillar_1.6.1      
[13] rlang_1.0.4        readxl_1.3.1       rstudioapi_0.13    htmlwidgets_1.5.4 
[17] munsell_0.5.0      broom_0.7.8        compiler_4.1.0     modelr_0.1.8      
[21] pkgconfig_2.0.3    htmltools_0.5.1.1  tidyselect_1.1.1   fansi_0.5.0       
[25] crayon_1.4.1       dbplyr_2.1.1       withr_2.4.2        bitops_1.0-7      
[29] grid_4.1.0         jsonlite_1.7.2     gtable_0.3.0       lifecycle_1.0.0   
[33] DBI_1.1.1          magrittr_2.0.1     scales_1.1.1       KernSmooth_2.23-20
[37] cli_3.1.0          stringi_1.6.2      fs_1.5.0           limma_3.48.3      
[41] robustbase_0.93-9  xml2_1.3.3         ellipsis_0.3.2     generics_0.1.0    
[45] vctrs_0.3.8        RColorBrewer_1.1-2 tools_4.1.0        glue_1.6.0        
[49] DEoptimR_1.0-10    hms_1.1.0          colorspace_2.0-2   caTools_1.18.2    
[53] rvest_1.0.3        rJava_1.0-6        haven_2.4.1
```

# Installation Guide
Users only need to install R software (https://www.r-project.org/) and required packages (This process usually takes within 30 mins).


# Quick Example
Here is a quick run for fatty acid, lipid species, and lipid class substructure analysis using iLipidome. Before analysis, please download the files in "Documentation" folder. We also provide detailed inforamtion for each function in the documentation.


## Source function and load required data
Required function and data files can be found in the  "Required_function" and "Required_data" folder, respectively.
 
```{r Source function and load required data}

file <- dirname(rstudioapi::getSourceEditorContext()$path)

#Source function
source(file.path(file,'required_function.R'))

#Load required data
load(file.path(file,'required_data.RData'))

```

## Upload lipidomics data

iLipidome only requires users to upload one processed lipid expression table (data.frame) where lipids are rows and samples are columns for analysis. Lipid names should be in the first column named as "feature", and sample names are in the first row (see example below).  At least two samples in each group are required to calculate statistics. Also, data processing or normalization methods, such as missing value imputation or log transformation, may be required based on data source to achieve better results before analysis.

iLipidome only requires users to upload one processed lipid expression table (data.frame) where lipids are rows and samples are columns for analysis. Lipid names should be in the first column named as “feature”, and sample names are in the first row (see example below). At least two samples in each group are required to calculate statistics. Also, data processing or normalization methods, such as missing value imputation or log transformation, may be required based on data source to achieve better results before analysis.

Lipid names can be represented as:
1. [LipidClass]_[sum of FA chain length] : [sum of FA double bonds] ; [sum of FA oxygens]
e.g., PC_34:1;0 or TAG_52:1;0 when the exact identity of FAs is unknown.
2. [LipidClass]_[FA1 chain length] : [FA1 double bonds] ; [FA1 oxygens]_[FA2 chain length] : [FA2 double bonds] ; [FA2 oxygens]…
e.g., PC_16:0;0_18:1;0 or TAG_16:0;0_18:0;0_18:1;0 when the exact identity of FAs is known.

Supported lipid classes, abbreviations, and corresponding FA numbers can be found in the “supported_lipid_class.csv” file. Note that lipid classes with same FA numbers (e.g., PC, PE) in same pathways (e.g., Glycerophospholipid) should have consistent lipid naming format (e.g., PC_36:0;0 and PE_34:0;0 or PC_18:0;0_18:0;0 and PE_16:0;0_18:0;0). Further, dihydrosphingolipids (dh-) specify the sphingolipids with sphingoid bases of 18:0:2 instead of 18:1:2.

```{r Upload lipidomics data and process format}

#Expression table of example lipidomics dataset
exp <- read.csv(file.path(file, 'exp.csv'))

head(exp)

```

## Process data for iLipidome inputs

"build_char_table" transforms lipid expression table ("exp") into two iLipidome inputs: selected lipid expression table  ("exp_sel") and selected lipid characteristics table ("char_sel"). Note that it only considers the lipid classes recorded in the "network_node" table.

```{r Upload lipidomics data and process format2}

exp_sel <- build_char_table(raw_data=exp, network_node = network_node)[[1]]

#selected lipid expression table
head(exp_sel)

char_sel <- build_char_table(exp, network_node = network_node)[[2]]

#selected lipid characteristics table
head(char_sel)

```

## Analysis for unprocessed data
"unprocessed_data_test" uses the output of "build_char_table" to perform differential expression for three types of data: (1) lipid species, (2) fatty acids, and (3) lipid classes.

```{r Analysis for unprocessed data}

no_sub_t <- unprocessed_data_test(exp_data = exp_sel,
                                  char_table = char_sel,
                                  method = 't.test',
                                  significant='adj_p_value',
                                  ctrl_group = 1:7, exp_group = 8:13)

#Expression tables for lipid species, fatty acids, and lipid classes
no_sub_t[[1]] %>% head()

#Statistical result table for lipid species, fatty acids, and lipid classes
no_sub_t[[2]] %>% head()

```

## 1. FA substructure analysis
Here, we provide a step-by-step process to perform FA substructure analysis using the data above and a series of functions.

### 1-1. FA biosynthetic network transformation
Firstly, the reference FA biosynthetic network is trimmed by users' data.
 
```{r FA substructure analysis 1}

FA_network_new <- build_FA_net(FA_network = FA_network,
                           unprocessed_data_result = no_sub_t)

#Trimmed FA biosynthetic network
FA_network_new %>% head()

```

### 1-2. Decompose FAs into FA substructures
"FA_sub_transform" decomposes FAs into FA substructures based on the FA biosynthetic network.
 
```{r FA substructure analysis 2}

#18:2 and 20:4 are majorly omega-6 FAs, so we only kept omega-6 forms of them.

FA_substructure <- FA_sub_transform(FA_network = FA_network_new,
                                    unprocessed_data_result = no_sub_t,
                                    unmapped_FA = c('w9-18:2;0','w3-20:4;0'))

#FA substructure table
FA_substructure %>% head()

```

### 1-3. Extract FA substructures using fold changes
"FA_sub_extract" maps FA substructures in each pathway with fold changes from the "unprocessed_data_test" result and extracts them through a backpropagated process. Specifically, the checking process starts from the last substructure (target FA) and would not stop until it meets a substructure with an opposite fold change along the biosynthetic route. One exception is the endogenous biosynthesis pathway for FAs in the upstream of palmitate (e.g., 14:0 or 12:0). Since they are synthesized as a group (2:0 to 16:0), we do not check their fold change and keep all substructures.

```{r FA substructure analysis 3}

FA_sub_stop <- FA_sub_extract(char_table = char_sel,
                              FA_substructure = FA_substructure,
                              unprocessed_data_result = no_sub_t,
                              exact_FA='no', exo_lipid='w3-22:6;0')

#lipid species
FA_sub_stop[[1]] %>% head()

#Extracted FA substructures for lipid species
FA_sub_stop[[2]] %>% head()

```

### 1-4. Transform FA exp into substructure exp

The function converts expression of FAs to expression of FA substructures.
 
```{r FA substructure analysis 4}

FA_sub_exp <- lipid_sub_matrix(exp_data = exp_sel, sub_data = FA_sub_stop,
                               sub_type = 'FA')

#FA substructure matrix encoding the frequency of each substructure 
FA_sub_exp[[1]][1:5, 1:5]

#Lipid profile
FA_sub_exp[[2]]%>% head()

#FA substructure profile
FA_sub_exp[[3]] %>% head()

```

### 1-5. Differential expression analysis for FA substructures

```{r FA substructure analysis 5}

FA_sub_exp_t <- t_test(data = FA_sub_exp[[3]], ctrl = 1:7, exp = 8:13,
                       method = 't.test', significant = 'adj_p_value')

#Statistical result table for FA substructures
FA_sub_exp_t %>% head()

```

### 1-6. Essential pathway analysis for FA substructures

"path_scoring" use FA substructures to score pathways in FA biosynthetic network.

```{r FA substructure analysis 6}

set.seed(1)
path_score_FA <- path_scoring(network = FA_network_new, sub_t = FA_sub_exp_t, 
                              calibrate = T, data_type = 'FA')

#Pathway scoring result table
path_score_FA %>% head()

```

### 1-7. Essential edges (reactions) analysis for FA substructures

"reaction_scoring" evaluates each reaction in FA biosynthetic network using FA substructures.

```{r FA substructure analysis 7}

reaction_score_FA <- reaction_scoring(network = FA_network_new, 
                                      sub_exp = FA_sub_exp[[3]],
                                      sub_t = FA_sub_exp_t, 
                                      ctrl = 1:7, exp = 8:13, 
                                      Species = 'rat')

#Reaction scoring result table
reaction_score_FA %>% head()

```

### 1-8. FA biosynthetic network construction

Build the FA biosynthetic network using FA substructures, pathway and reaction scoring results.

```{r FA substructure analysis 8}

FA_network_data <- draw_network(network_data = FA_network_new,
                                DE_data = FA_sub_exp_t,
                                if_species = F, significant = 'adj_p_value',
                                path_scoring_result = path_score_FA,
                                reaction_scoring_result = reaction_score_FA,
                                top_n = 5, path_type = 'both')

#FA biosynthetic network node
FA_network_data[[1]] %>% head()

#FA biosynthetic network edge
FA_network_data[[2]] %>% head()

#FA biosynthetic network
visNetwork(FA_network_data[[1]],FA_network_data[[2]]) %>% 
  visIgraphLayout(layout = "layout_with_sugiyama", type='square',
                  physics = F, smooth = TRUE, randomSeed =5) 

```

## 2. Lipid species substructure analysis

A similar approach can be used to analyze lipid species substructures.

### 2-1. Decompose lipids into species substructures

"species_sub_transform" decomposes lipids into species substructures based on the lipid biosynthetic network.
 
```{r Lipid species substructure analysis 1}

#We excluded ether lipids since we cannot differentiate Alkyl (O-) or Alkenyl- (P-) linked ether lipids

char_wo_EL <- char_sel[!str_detect(char_sel$feature, 'O-'),]
exp_wo_EL <- exp_sel[!str_detect(exp_sel$feature, 'O-'),]

species_substructure <- species_sub_transform(char = char_wo_EL,
                                              lipid_substructure = lipid_substructure,
                                              network_node = network_node)


#Lipid species substructure table
species_substructure %>% head()

```

### 2-2. Extract species substructures using fold changes
"species_sub_extract" maps species substructures in each pathway with fold changes from the "unprocessed_data_test" result and extracts them through a backpropagated process. Specifically, the checking process starts from the last substructure (target species) and would not stop until it meets a substructure with an opposite fold change along the biosynthetic route.

```{r Lipid species substructure analysis 2}

species_sub_stop <- species_sub_extract(lipid_substructure = species_substructure,
                                        unprocessed_data_result =  no_sub_t,
                                        type = 'species', pct_limit = 0.3,
                                        exo_lipid=NULL)

#Lipid species
species_sub_stop[[1]] %>% head()

#Extracted species substructures for lipid species
species_sub_stop[[2]] %>% head()

```

### 2-3. Transform lipid exp into substructure exp

The function converts expression of lipid species to expression of species substructures.
 
```{r Lipid species substructure analysis 3}

species_sub_exp <- lipid_sub_matrix(exp_data = exp_wo_EL, 
                                    sub_data = species_sub_stop,
                                    sub_type = 'Species')


#Species substructure matrix encoding the frequency of each substructure 
species_sub_exp[[1]][1:5, 1:5]

#Lipid profile
species_sub_exp[[2]] %>% head()

#Species substructure profile
species_sub_exp[[3]] %>% head()

```

### 2-4. Differential expression analysis for species substructures

```{r Lipid species substructure analysis 4}

species_sub_exp_t <- t_test(data = species_sub_exp[[3]], ctrl = 1:7, exp = 8:13,
                            method = 't.test', significant = 'adj_p_value')


#Statistical result table for species substructures
species_sub_exp_t %>% head()

```

### 2-5. Lipid species biosynthetic network transformation

"build_species_net" uses species substructures to contruct lipid biosynthetic network.

```{r Lipid species substructure analysis 5}

#species_substructure: Output of "species_sub_transform".

species_network <- build_species_net(species_substructure = species_substructure)

#Lipid species biosynthetic network
species_network %>% head()

```

### 2-6. Essential pathway analysis for species substructures

"path_scoring" use species substructures to score pathways in lipid species biosynthetic network.

```{r Lipid species substructure analysis 6}


set.seed(1)
path_score_species <-  path_scoring(network = species_network,
                                    sub_t = species_sub_exp_t,
                                    calibrate = T, data_type = 'Species')


#Pathway scoring result table
path_score_species %>% head()

```

### 2-7. Essential edges (reactions) analysis for species substructures

"add_rev_rection" completes all reversible reactions in lipid species biosynthetic network, where "reaction_scoring" evaluates each reaction using species substructures.

```{r Lipid species substructure analysis 7}


species_net_w_rev <- add_rev_rection(network_edge = network_edge,
                                     species_net = species_network)

#Lipid species biosynthetic network with complete reversible reactions
species_net_w_rev %>% head()

reaction_score_species <- reaction_scoring(network = species_net_w_rev,
                                           sub_exp = species_sub_exp[[3]],
                                           sub_t = species_sub_exp_t,
                                           ctrl=1:7, exp=8:13,
                                           Species = 'rat')

#Reaction scoring result table
reaction_score_species %>% head()

```

### 2-8. Lipid species biosynthetic network construction

Build the lipid species biosynthetic network using species substructures, pathway and reaction scoring results.

```{r Lipid species substructure analysis 8}


species_network_data <- draw_network(network_data = species_net_w_rev,
                                     DE_data = species_sub_exp_t,
                                     if_species = T,significant = 'adj_p_value',
                                     path_scoring_result = path_score_species,
                                     reaction_scoring_result = reaction_score_species,
                                     top_n = 3, path_type = 'both')



#Lipid species biosynthetic network node
species_network_data[[1]] %>% head()

#Lipid species biosynthetic network edge
species_network_data[[2]] %>% head()

#Lipid species biosynthetic network
visNetwork(species_network_data[[1]], species_network_data[[2]])

```

## 3. Lipid class substructure analysis

<font size="3"> A similar approach can be used to analyze lipid class substructures.

### 3-1. Extract class substructures using fold changes

```{r Lipid class substructure analysis 1}

class_sub_stop <- species_sub_extract(lipid_substructure =lipid_substructure,
                                      unprocessed_data_result = no_sub_t,
                                      type = 'class', pct_limit = 0.01,
                                      exo_lipid=NULL)

#Lipid classes
class_sub_stop[[1]] %>% head()

#Extracted class substructures for lipid classes
class_sub_stop[[2]] %>% head()

```

### 3-2. Transform lipid class exp into substructure exp

The function converts expression of lipid classes to expression of class substructures.
 
```{r Lipid class substructure analysis 2}

#Lipid class expression table.
class_exp <- no_sub_t[[1]] %>% filter(type=='class') %>% 
  dplyr::select(-type)

class_sub_exp <- lipid_sub_matrix(exp_data = class_exp, 
                                    sub_data = class_sub_stop,
                                    sub_type = 'Class')


#Class substructure matrix encoding the frequency of each substructure 
class_sub_exp[[1]][1:5, 1:5]

#Lipid class profile
class_sub_exp[[2]] %>% head()

#Class substructure profile
class_sub_exp[[1]] %>% head()

```

### 3-3. Differential expression analysis for lipid class substructures

```{r Lipid class substructure analysis 3}

class_sub_exp_t <- t_test(data = class_sub_exp[[3]], ctrl = 1:7, exp = 8:13,
                          method = 't.test', significant = 'adj_p_value')

#Statistical result table for class substructures
class_sub_exp_t %>% head()

```

### 3-4. Lipid class biosynthetic network transformation

The reference lipid biosynthetic network in iLipdiome is trimmed by class substructures to build lipid class network.

```{r Lipid class substructure analysis 4}

class_network <- network_edge[c('S1','P1')] %>% 
  filter(S1 %in% class_sub_exp_t$lipid, P1 %in% class_sub_exp_t$lipid)

#Lipid class biosynthetic network
class_network %>% head()

```

### 3-5. Essential pathway analysis for species substructures

"path_scoring" use class substructures to score pathways in lipid class biosynthetic network.

```{r Lipid class substructure analysis 5}

set.seed(1)
path_score_class <-  path_scoring(network = class_network,
                                    sub_t = class_sub_exp_t,
                                    calibrate = T, data_type = 'Class')


#Pathway scoring result table
path_score_class %>% head()

```

### 3-6. Essential edges (reactions) analysis for species substructures

“reaction_scoring” evaluates each reaction in lipid class biosynthetic network using class substructures.

```{r Lipid class substructure analysis 6}

reaction_score_class <- reaction_scoring(network = class_network,
                                         sub_exp = class_sub_exp[[3]],
                                         sub_t = class_sub_exp_t,
                                         ctrl=1:7, exp=8:13,
                                         Species = 'rat')



#Reaction scoring result table
reaction_score_class %>% head()

```

### 3-7. Lipid class biosynthetic network construction

Build the lipid class biosynthetic network using class substructures, pathway and reaction scoring results.

```{r Lipid class substructure analysis 7}

class_network_data <- draw_network(network_data = class_network,
                                     DE_data = class_sub_exp_t,
                                     if_species = F,significant = 'adj_p_value',
                                     path_scoring_result = path_score_class,
                                     reaction_scoring_result = reaction_score_class,
                                     top_n = 3, path_type = 'both')


#Lipid class biosynthetic network node
class_network_data[[1]] %>% head()

#Lipid class biosynthetic network edge
class_network_data[[2]] %>% head()

#Lipid class biosynthetic network
visNetwork(class_network_data[[1]], class_network_data[[2]])

```

# License

