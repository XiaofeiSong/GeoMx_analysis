

# set up evrn -------------------------------------------------------------

library(tidyverse)
library(NanoStringNCTools) 
library(GeomxTools) 
library(GeoMxWorkflows) 
library(gridExtra)
library(scales) 
library(patchwork)
library(umap)


##output folder
outputdir <- file.path("unsupervised")

mydate <- Sys.Date()

##input the tcc after normalizaiaton - target_dcc_set
load('target_dcc_set_filtered_normed_2023-10-26.Rdata')


# Dimension reduction -----------------------------------------------------
norm_method = "q_norm"


# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

# run UMAP
umap_out <- umap(t(log2(assayDataElement(target_dcc_set , elt = norm_method))),  
                 config = custom_umap)

pData(target_dcc_set)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]

# plot

## define colors 
col25_cols = cols25(n = 10)
slide_color <- col25_cols[1:3]
names(slide_color) <- unique(pData(target_dcc_set)$`slide name`)

region_color <- col25_cols[4:6]
names(region_color) <- unique(pData(target_dcc_set)$region)

immu_status_color <- col25_cols[7:8]
names(immu_status_color) <- unique(pData(target_dcc_set)$immu_status)

ID_site_color <- col25_cols[9:10]
names(ID_site_color) <- unique(pData(target_dcc_set)$ID_site)


pdf(paste0(outputdir, "unsupervised_umap", mydate, '.pdf'), width=12, height=10)

p1 <- ggplot(pData(target_dcc_set),
             aes(x = UMAP1, y = UMAP2, color = `slide name`)) +
  scale_color_manual(values = slide_color) +
  geom_point(size = 3) +
  theme_bw()

p2 <- ggplot(pData(target_dcc_set),
       aes(x = UMAP1, y = UMAP2, color = region)) +
  scale_color_manual(values = region_color) +
  geom_point(size = 3) +
  theme_bw()

p3 <- ggplot(pData(target_dcc_set),
             aes(x = UMAP1, y = UMAP2, color = immu_status)) +
  scale_color_manual(values =immu_status_color) +
  geom_point(size = 3) +
  theme_bw()  

p4 <- ggplot(pData(target_dcc_set),
             aes(x = UMAP1, y = UMAP2, color = ID_site)) +
  geom_point(size = 3) +
  scale_color_manual(values = ID_site_color) +
  theme_bw()

((p1/p3) + plot_layout(guides = 'collect')) |
  ((p2/p4) + plot_layout(guides = 'collect'))

dev.off()


