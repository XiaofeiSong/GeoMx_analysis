library(tidyverse)
library(cowplot)
library(patchwork)
library(NanoStringNCTools) 
library(GeomxTools) 
library(GeoMxWorkflows) 


outputdir <- "main_figures/"

load('color_book.Rdata')

mydate <- Sys.Date()

metadata <- pData(target_dcc_set)

# umaps -------------------------------------------------------------------
metadata$region_type <- factor(metadata$region_type, levels = c('Tu', 'Separative', 'Mixture', 'St'))
metadata$ID_site <- factor(metadata$ID_site, levels = c('Up', 'Lw'))

p1 <- ggplot(metadata,
             aes(x = UMAP1, y = UMAP2, color = ID_site)) +
  geom_point(size = 3) +
  scale_color_manual(values = ID_site_color) +
  theme_bw() +
  theme(legend.position = "bottom")  

p2 <- ggplot(metadata,
             aes(x = UMAP1, y = UMAP2, color = region_type)) +
  scale_color_manual(values = region_type_color) +
  geom_point(size = 3) +
  theme_bw() +
  theme(legend.position = "bottom")  

p3 <- ggplot(metadata,
             aes(x = UMAP1, y = UMAP2, color = immu_status)) +
  scale_color_manual(values =immu_status_color) +
  geom_point(size = 3) +
  theme_bw() +
  theme(legend.position = "bottom")  


pdf(paste0(outputdir, "UMAP_", mydate, '.pdf'), width=8, height=3)

(p1|p2|p3) 

dev.off()

pdf(paste0(outputdir, "UMAP_legend_", mydate, '.pdf'), width=15, height=3)

(p1|p2|p3) 

dev.off()


pdf(paste0(outputdir, 'UMAP_region_type.pdf'), width=6, height=4)

ggplot(metadata,
       aes(x = UMAP1, y = UMAP2, color = region_type)) +
  scale_color_manual(values = region_type_color) +
  geom_point(size = 3) +
  theme_bw()

dev.off()





# barplot for deconvolution results ---------------------------------------
library(ComplexHeatmap)

metadata <- metadata %>%
  arrange(region_type, ID_site)

res_pct_mat <- res_pct_mat[,metadata$Sample_ID2]

all(rownames(metadata) == colnames(res_pct_mat))

stacked_barplot = anno_barplot(t(res_pct_mat), 
                               gp=gpar(fill = cellcols[rownames(res_pct_mat)]),
                               height = unit(6, "cm"))

all(colnames(res_pct_mat) == metadata$Sample_ID2)

ha = HeatmapAnnotation(Tissue_site = metadata$ID_site,
                       Immu_status = metadata$immu_status,
                       ROI_region = metadata$region_type,
                       Composition = stacked_barplot,
                       col = list(ROI_region = region_type_color,
                                  Immu_status = immu_status_color,
                                  Tissue_site = ID_site_color),
                       annotation_name_side = "left")


pdf(paste0(outputdir, "Deconvolution_cellular_composition_", mydate, '.pdf'), width=12, height=5)

H1 <- Heatmap(
  matrix(nr = 0, nc = ncol(res_pct_mat)),
  top_annotation = ha,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  row_names_gp = gpar(fontsize = 8),
  cluster_rows = F,
  height = unit(6, "cm")
  
)

draw(H1,
     annotation_legend_side = "bottom")

dev.off()


pdf(paste0(outputdir, "Deconvolution_cellular_composition_legend_", mydate, '.pdf'), width=15, height=5)

lgd = Legend(labels = names(cellcols[rownames(res_pct_mat)]), 
             title = "celltype", 
             legend_gp = gpar(fill = cellcols[rownames(res_pct_mat)]),
             ncol = 1)


draw(lgd)

dev.off()


# compare cc pct between two methods --------------------------------------
library(ggpubr)

##spatialdecon results

decon_pct <- t(res_pct_mat) %>%
  as_tibble(rownames = 'Sample_ID') %>%
  mutate(Sample_ID = str_remove(Sample_ID, '\\.dcc'),
         cc_pct_decon = tumor) %>%
  select(Sample_ID, cc_pct_decon)


##image results

cc_pct_sum <- metadata %>%
  left_join(decon_pct)

cc_pct_sum$region_type <- factor(cc_pct_sum$region_type, levels = c('Tu', 'Separative', 'Mixture', 'St'))

ggplot(cc_pct_sum, aes(x = cc_pct_decon, y = cc_pct_image)) +            
  geom_smooth(method = "lm", formula = y ~ x) +
  geom_point(data = cc_pct_sum, aes(color=region_type)) +
  scale_color_manual(values = region_type_color) +
  stat_cor(label.y = 1) +
  xlab("cancer cell% - Deconvolution") +
  ylab("cancer cell% - Image") +
  theme_bw() 


ggsave(paste0(outputdir, 'cc_pct_correlation.pdf'), width = 5, height = 3)


cc_pct_sum$region <- factor(cc_pct_sum$region, levels = c('Tu', 'TuSt', 'St'))


quartz(pointsize = 6) # define point size
par(mar=c(3,3,1,1), family = "Arial")  


ggplot(cc_pct_sum, aes(x =  cc_pct_image, y =cc_pct_decon)) +            
  geom_smooth(method = "lm", formula = y ~ x) +
  geom_point(data = cc_pct_sum, aes(color=region)) +
  scale_color_manual(values = region_color) +
  stat_cor(label.y = 1) +
  xlab("cancer cell% - ROICellTrack") +
  ylab("cancer cell% - Deconvolution") +
  theme_classic() +
  theme(legend.position="bottom")


ggsave(paste0(outputdir, 'cc_pct_correlation_region.pdf'), width = 2.73, height = 4.095)

# barplot for deconvolution results -subset ---------------------------------------


library(ComplexHeatmap)


##select samples

cc_pct_sum_selected <- cc_pct_compare %>%
  filter(region != 'TuSt') %>%
  filter((region == 'Tu' & abs(image_decon_diff) > 0.2) | (region == 'St' & abs(image_decon_diff) > 0.25)) %>%
  arrange(desc(region), desc(abs(image_decon_diff)))

res_pct_mat <- res_pct_mat[,paste0(cc_pct_sum_selected$Sample_ID,'.dcc')]

metadata <- metadata[paste0(cc_pct_sum_selected$Sample_ID,'.dcc'),] %>%
  filter(Sample_ID2 %in% paste0(cc_pct_sum_selected$Sample_ID,'.dcc'))

all(rownames(metadata) == colnames(res_pct_mat))


stacked_barplot = anno_barplot(t(res_pct_mat), 
                               gp=gpar(fill = cellcols[rownames(res_pct_mat)]),
                               height = unit(6, "cm"))

all(colnames(res_pct_mat) == metadata$Sample_ID2)

ha = HeatmapAnnotation(
                       Composition = stacked_barplot,
                    
                       ROI_region = metadata$region_type,
                       col = list(ROI_region = region_color),
                       annotation_name_side = "left")



H1 <- Heatmap(
  matrix(nr = 0, nc = ncol(res_pct_mat)),
  bottom_annotation = ha,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  row_names_gp = gpar(fontsize = 8),
  cluster_rows = F,
  height = unit(6, "cm")
  
)
H1 

pdf(paste0(outputdir, 'Deconvolution_cellular_composition_selected_ROI_.pdf'), width=5, height=25)

draw(H1,
     annotation_legend_side = "bottom")


dev.off()


# plot spatial clustering -------------------------------------------------

##plot spatial clustering

plot_list <- list() 

metadata <- metadata %>%
  arrange(factor(region, levels = c('Tu', 'TuSt', 'St')), auc_value)

for(i in 1:nrow(metadata)){
  print(i)
  ## two inputs - input csv file and cutoff
  cell_stat_one_file <- metadata$cell_stat_file[i]
  
  G_Int_cutoff <- ifelse(metadata$PTID[i] == "Patient1", 15, 20)
  
  ## Process
  
  cell_stat <- read.csv(cell_stat_one_file) %>%
    mutate(tu_status = if_else(G_Int > G_Int_cutoff, 'Tu', 'NotTu'))
  
  
  plot_list[[i]] <- ggplot(cell_stat, aes(x = X_coordinate, y = Y_coordinate, color = tu_status)) + 
    geom_point() +
    scale_color_manual(values = c('Tu' = 'red', 'NotTu'= 'darkgrey')) + 
    theme_void() +
    scale_y_reverse() +
    #ggtitle(round(metadata$auc_value[i], 3)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none")
  
}


names(plot_list) <- metadata$region

p1 <- plot_grid(plotlist=plot_list[names(plot_list)=='Tu'], ncol=6)
p2 <- plot_grid(plotlist=plot_list[names(plot_list)=='TuSt'], ncol=6)
p3 <- plot_grid(plotlist=plot_list[names(plot_list)=='St'], ncol=6)

#plot_grid(p1, p, p3, ncol = 1, rel_heights = c(2,4,3))
#ggsave(file = file.path(outputdir, 'spatial_pattern_sum_split_by_region_type.pdf'), width = 14, height = 18)

plot_grid(p1, NULL, p3, NULL, p2, ncol = 1, rel_heights = c(2,1, 3,1, 5))
ggsave(file = file.path(outputdir, 'spatial_pattern_sum_split_by_region_type_v2.pdf'), width = 9.4, height = 18.8)

ggsave(file = file.path(outputdir, 'spatial_pattern_sum_split_by_region_type_v2.png'), width = 9.4, height = 18.8)



# pathway analysis result -------------------------------------------------


GSEA.result.df <- GSEA.result.df %>%
  slice(unique(c(1:5, n() - 4:0)) ) %>%
  mutate(Description = str_replace_all(str_replace(Description, 'HALLMARK_', ''),
                                                   '_', ' ')) %>%
  mutate(region = if_else(NES > 0, 'Mixture', 'Separative')) %>%
  mutate(NES.abs = abs(NES)) %>%
  arrange(region, desc(NES.abs))

## modify abb locally


ggplot(GSEA.result.df, aes(x=NES.abs, y=reorder(Pathway, NES.abs), fill = region)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = region_type_color) +
  labs(y = "",
       x = "abs(NES)") +
  theme(axis.text.y = element_text(face = 'bold')) +
  theme_classic() + 
  theme(legend.position="none")


ggsave(file = file.path(outputdir, 'selected_mix_vs_sep_pathways_240715.pdf'), width = 3, height = 3)


























