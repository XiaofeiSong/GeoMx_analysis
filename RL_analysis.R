


library(NanoStringNCTools) 
library(GeomxTools) 
library(GeoMxWorkflows) 
library(gridExtra)
library(scales) 
library(GeoDiff)
library(tidyverse)

load('target_dcc_set_filtered_normed_2023-10-26.Rdata')

##prepare inputs

##expression matrix
expr_mat <- t(assayData(target_dcc_set)[["log_q"]])

##clincial data
smpl_anno <- pData(target_dcc_set) %>%
  select(Sample_ID2, region_type)


# RL database -------------------------------------------------------------


#test correlation in each type of region ------------------------
source("S_geomx_RL_analysis.R")

## run the analysis
Tu_IDs <- smpl_anno %>%
  filter(region_type == 'Tu') %>%
  pull(Sample_ID2)
Tu_expr_mat <- expr_mat[Tu_IDs, ]
correl_RL_avg_Tu <- geomx_RL_test(Tu_expr_mat)

St_IDs <- smpl_anno %>%
  filter(region_type == 'St') %>%
  pull(Sample_ID2)
St_expr_mat <- expr_mat[St_IDs, ]
correl_RL_avg_St <- geomx_RL_test(St_expr_mat)

Separative_IDs <- smpl_anno %>%
  filter(region_type == 'Separative') %>%
  pull(Sample_ID2)
Separative_expr_mat <- expr_mat[Separative_IDs, ]
correl_RL_avg_Separative <- geomx_RL_test(Separative_expr_mat)

Mixture_IDs <- smpl_anno %>%
  filter(region_type == 'Mixture') %>%
  pull(Sample_ID2)
Mixture_expr_mat <- expr_mat[Mixture_IDs, ]
correl_RL_avg_Mixture <- geomx_RL_test(Mixture_expr_mat)


## compare the results
Tu <- correl_RL_avg_Tu %>%
  mutate(Tu_cor = Correl, Tu_sig = if_else(Tu_cor > 0.75, 1, 0)) %>%
  select(ligand_target, source, database, Tu_cor, Tu_sig)

St <- correl_RL_avg_St %>%
  mutate(St_cor = Correl, St_sig = if_else(St_cor > 0.75, 1, 0)) %>%
  select(ligand_target, St_cor, St_sig)

Separative <- correl_RL_avg_Separative %>%
  mutate(Separative_cor = Correl, Separative_sig = if_else(Separative_cor > 0.75, 1, 0)) %>%
  select(ligand_target, Separative_cor, Separative_sig)

Mixture <- correl_RL_avg_Mixture %>%
  mutate(Mixture_cor = Correl, Mixture_sig = if_else(Mixture_cor > 0.75, 1, 0)) %>%
  select(ligand_target, Mixture_cor, Mixture_sig)

RL_results_sum <- Tu %>%
  full_join(St) %>%
  full_join(Separative) %>%
  full_join(Mixture) %>%
  filter(!database %in% c('ppi_prediction', 'ppi_prediction_go')) 
  
write.csv(RL_results_sum, file = 'RL_results_sum.csv', row.names = F)

# visualization ------------------------------------------------------------


# prepare input for cc plot----------------------------------------------------

ccplot_Tu <- correl_RL_avg_Tu %>%
  filter(!database %in% c('ppi_prediction', 'ppi_prediction_go')) %>%
  filter(Correl > 0.75) %>%
  mutate(source = 'Tu_sender', 
         target = 'Tu_receiver',
         ligand = from,
         receptor = to,
         score = Correl) %>%
  arrange(score)

nrow(ccplot_Tu)

ccplot_St <- correl_RL_avg_St %>%
  filter(!database %in% c('ppi_prediction', 'ppi_prediction_go')) %>%
  filter(Correl > 0.75) %>%
  mutate(source = 'St_sender', 
         target = 'St_receiver',
         ligand = from,
         receptor = to,
         score = Correl) %>%
  arrange(score)

nrow(ccplot_St)
ccplot_Separative <- correl_RL_avg_Separative %>%
  filter(!database %in% c('ppi_prediction', 'ppi_prediction_go')) %>%
  filter(Correl > 0.75) %>%
  mutate(source = 'Separative_sender', 
         target = 'Separative_receiver',
         ligand = from,
         receptor = to,
         score = Correl) %>%
  arrange(score)
nrow(ccplot_Separative)

ccplot_Mixture <- correl_RL_avg_Mixture %>%
  filter(!database %in% c('ppi_prediction', 'ppi_prediction_go')) %>%
  filter(Correl > 0.75) %>%
  mutate(source = 'Mixture_sender', 
         target = 'Mixture_receiver',
         ligand = from,
         receptor = to,
         score = Correl) %>%
  arrange(score)
nrow(ccplot_Mixture)


# create a function to retrieve express data -------------------------------

##expression matrix
expr_mat <- t(assayData(target_dcc_set)[["log_q"]])

##clincial data
smpl_anno <- pData(target_dcc_set) %>%
  select(Sample_ID2, region_type)

get_expr_mean <- function(type, top_n, ccplot_type){
  type_IDs <- smpl_anno %>%
    filter(region_type == type) %>%
    pull(Sample_ID2)
  
  type_mean_expr <- colMeans(expr_mat[type_IDs, ]) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'gene')
  
  colnames(type_mean_expr)[2] <- 'mean_exp'

  gene_pair <- ccplot_type %>%
    slice(1:top_n)
  
  temp1 <- gene_pair %>%
    select(source, from)
  
  temp2 <- gene_pair %>%
    select(target, to)
  
  colnames(temp1) <- c('cell_type', 'gene')
  colnames(temp2) <- c('cell_type', 'gene')
  
  
  type_mean_expr_top <- rbind(temp1, temp2) %>%
    left_join(type_mean_expr) %>%
    unique()
  
  return(type_mean_expr_top)
}



# plot --------------------------------------------------------------------
library(CCPlotR)
library(circlize)
library(ComplexHeatmap)
source('S_ccplotr_circos.R')

## region_color 

region_color_ccplot <- c('#c96587', '#adc6e8',  '#a57fe6',  '#fcc563','#802d49', '#4d82cd','#7f47dc', '#faad16')
names(region_color_ccplot) <- c(paste0(c('Mixture', 'Separative', 'St', 'Tu'), '_sender'),
                                paste0(c('Mixture', 'Separative', 'St', 'Tu'), '_receiver'))



##Mixture

top_n = nrow(ccplot_Mixture)
Mixture_mean_expr_top <- get_expr_mean('Mixture', top_n, ccplot_Mixture)

pdf('cc_circos_Mixture.pdf', width = 10, height = 10)

cc_circos_2(ccplot_Mixture, option = "C", n_top_ints = top_n, exp_df = Mixture_mean_expr_top,
            cell_cols = region_color_ccplot, palette = "PuRd")

dev.off()





##TU

top_n = nrow(ccplot_Tu)
Tu_mean_expr_top <- get_expr_mean('Tu', top_n, ccplot_Tu)

pdf('cc_circos_Tu.pdf', width = 10, height = 10)

cc_circos_2(ccplot_Tu, option = "C", n_top_ints = top_n, exp_df = Tu_mean_expr_top,
          cell_cols = region_color_ccplot, palette = "PuRd")

dev.off()



##St
top_n = nrow(ccplot_St)
St_mean_expr_top <- get_expr_mean('St', top_n, ccplot_St)

pdf('cc_circos_St.pdf', width = 10, height = 10)

cc_circos_2(ccplot_St, option = "C", n_top_ints = top_n, exp_df = St_mean_expr_top,
          cell_cols = region_color_ccplot, palette = "PuRd")

dev.off()

##Separative
top_n = nrow(ccplot_Separative)
Separative_mean_expr_top <- get_expr_mean('Separative', top_n, ccplot_Separative)
pdf('cc_circos_Separative.pdf', width = 10, height = 10)

cc_circos_2(ccplot_Separative, option = "C", n_top_ints = top_n, exp_df = Separative_mean_expr_top, 
            cell_cols = region_color_ccplot, palette = "PuRd")
dev.off()



pdf('cc_circos_Mixture_Separative.pdf', width = 10, height = 10)
top_n = nrow(ccplot_Mixture)
cc_circos_2(ccplot_Mixture, option = "C", n_top_ints = top_n, exp_df = Mixture_mean_expr_top,
            cell_cols = region_color_ccplot, palette = "PuRd")

top_n = nrow(ccplot_Separative)
cc_circos_2(ccplot_Separative, option = "C", n_top_ints = top_n, exp_df = Separative_mean_expr_top, 
            cell_cols = region_color_ccplot, palette = "PuRd")

dev.off()


