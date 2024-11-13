


library(NanoStringNCTools) 
library(GeomxTools) 
library(GeoMxWorkflows) 
library(gridExtra)
library(scales) 
library(GeoDiff)
library(tidyverse)
source('S_geomx_de_analysis.R')


outputdir <- file.path("DE_analysis/")
mydate <- Sys.Date()

load('target_dcc_set_filtered_normed_2023-10-26.Rdata')


# compute DE --------------------------------------------------------------


#### Mixture vs Separative
my_compare <- 'Mixture_vs_Separative'

pData(target_dcc_set)$region_type <- factor(pData(target_dcc_set)$region_type, levels = c("St","Tu","Mixture", "Separative"))

ind <- pData(target_dcc_set)$region_type %in% c("Mixture", "Separative")

mixedOutmc <- mixedModelDE(target_dcc_set[, ind],
                           elt = "log_q",
                           modelFormula = ~ region_type + (1 + region_type | slide),
                           groupVar = "region_type",
                           nCores = parallel::detectCores(),
                           multiCore = FALSE)

outputdir_sub <- file.path(outputdir, my_compare)
dir.create(outputdir_sub)
compute_de(mixedOutmc = mixedOutmc, fc_cutoff = 1.5, p_cutoff = 0.05, outputdir = outputdir_sub)
save(mixedOutmc, file = file.path(outputdir_sub, 'mixedOutmc.Rdata'))


#### Up vs Lw Tumor
my_compare <- 'UpTu_vs_LwTu'

pData(target_dcc_set)$region_type <- factor(pData(target_dcc_set)$region_type, levels = c("St","Tu","Mixture", "Separative"))
pData(target_dcc_set)$ID_site <- factor(pData(target_dcc_set)$ID_site, levels = c("Up", "Lw"))

ind <- pData(target_dcc_set)$region_type == "Tu"

mixedOutmc <- mixedModelDE(target_dcc_set[, ind],
                           elt = "log_q",
                           modelFormula = ~ ID_site + (1 + ID_site | slide),
                           groupVar = "ID_site",
                           nCores = parallel::detectCores(),
                           multiCore = FALSE)

outputdir_sub <- file.path(outputdir, my_compare)
dir.create(outputdir_sub)
compute_de(mixedOutmc = mixedOutmc, fc_cutoff = 1.5, p_cutoff = 0.05, outputdir = outputdir_sub)
save(mixedOutmc, file = file.path(outputdir_sub, 'mixedOutmc.Rdata'))

