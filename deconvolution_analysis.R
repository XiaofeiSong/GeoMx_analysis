
# set up evrn -------------------------------------------------------------

library(tidyverse)
library(NanoStringNCTools) ##1.4.0
library(GeomxTools) ##3.0.1
library(GeoMxWorkflows) ##1.2.0
library(gridExtra)
library(scales) 
library(patchwork)
library(SpatialDecon)
library(ggpubr)
library(ComplexHeatmap)


##output folder
outputdir <- file.path("deconvolution/")

mydate <- Sys.Date()

##load predefined colors
load('color_book.Rdata')

# Run SpatialDecon -------------------------------------------------------------

##input the tcc after normalizaiaton - target_dcc_set
load('target_dcc_set_filtered_normed_2023-10-26.Rdata')

##add meta data column to indicate whether it's a tumor or not

target_dcc_set$istumor = target_dcc_set$region == "Tu"

##Run spatialDecon
restils = runspatialdecon(object = target_dcc_set,
                          norm_elt = "q_norm",                
                          raw_elt = "exprs",                      
                          X = safeTME,                           
                          cellmerges = safeTME.matches,           
                          is_pure_tumor = target_dcc_set$istumor,   
                          n_tumor_clusters = 1)                  


##extract deconvolution result matrix

### include tumors
res_mat <- t(restils$beta.granular)
res_pct_mat <- sweep(res_mat, 2, colSums(res_mat), "/")
rownames(res_pct_mat) <- str_replace(rownames(res_pct_mat), 'tumor.1', 'tumor')
res_pct_mat <- apply(res_pct_mat,2,rev)

### not include tumors
res_TME <- t(restils$prop_of_nontumor)
res_TME[is.na(res_TME)] <- 0

# Visualization -------------------------------------------------------------


##plot complex heatmap and stacked barplot - clustering columns


mymeta <- pData(restils)
all(rownames(mymeta) == colnames(res_pct_mat))

mymeta$region <- factor(mymeta$region, levels = c('Tu', 'TuSt', 'St'))

stacked_barplot = anno_barplot(t(res_pct_mat), 
                               gp=gpar(fill = cellcols[rownames(res_pct_mat)]),
                               height = unit(6, "cm"))

ha = HeatmapAnnotation(ID_site = mymeta$ID_site,
                       immu_status = mymeta$immu_status,
                       region = mymeta$region,
                       composition = stacked_barplot,
                       col = list(region = region_color,
                                  immu_status = immu_status_color,
                                  ID_site = ID_site_color))

colnames(res_TME) <- str_remove(colnames(res_TME), 'DSP-1001660008738-A-|\\.dcc')

H1 <- Heatmap(res_pct_mat,
              top_annotation = ha,
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = F,
              height = unit(6, "cm")
)


H2 <- Heatmap(res_TME,
              show_column_dend = FALSE,
              row_names_gp = gpar(fontsize = 8),
              show_row_dend = FALSE,
              height = unit(6, "cm")
)

pdf(file = paste0(outputdir, 'spatialDecon_result_231108.pdf'), width = 12, height = 12)
H1 %v% H2

dev.off()


pdf(file = paste0(outputdir, 'spatialDecon_result_celltype_legend_231108.pdf'), width = 4, height = 8)

lgd = Legend(labels = names(cellcols[rownames(res_pct_mat)]), 
             title = "celltype", 
             legend_gp = gpar(fill = cellcols[rownames(res_pct_mat)]))
draw(lgd)

dev.off()

