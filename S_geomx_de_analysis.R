
##load in curated gene information 
load('msigdb_hallmark_cosmic_hallmark_230710.Rdata')
load('maker_annotation_sum_230518.Rdata')

library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)



compute_de <- function(mixedOutmc, fc_cutoff = 1.5, p_cutoff = 0.05, outputdir){
  
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  r_test$Gene <- unlist(lapply(colnames(mixedOutmc), rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Contrast", "Estimate", "Pr(>|t|)", "FDR")]
  
  # add some gene info
  r_test <- r_test %>%
    left_join(gene_anno, by = c('Gene' ='gene.name')) %>%
    left_join(marker_sum, by = c('Gene' = 'marker_genes')) 
  
  write.csv(r_test, file = file.path(outputdir, 'DEG_results.csv'), row.names = F)
  
  # filter DEG
  r_test_deg0.05 <- r_test %>%
    filter(abs(Estimate) > log2(fc_cutoff) & `Pr(>|t|)` < p_cutoff & FDR < 0.05)
  
  write.csv(r_test_deg0.05, file = file.path(outputdir, paste0('DEG_results_FC',fc_cutoff, '_p', p_cutoff, '_fdr', 0.05, '.csv')), row.names = F)
  
  r_test_deg0.25 <- r_test %>%
    filter(abs(Estimate) > log2(fc_cutoff) & `Pr(>|t|)` < p_cutoff & FDR < 0.25)
  
  write.csv(r_test_deg0.25, file = file.path(outputdir, paste0('DEG_results_FC',fc_cutoff, '_p', p_cutoff, '_fdr', 0.25, '.csv')), row.names = F)
  
  r_test_deg1 <- r_test %>%
    filter(abs(Estimate) > log2(fc_cutoff) & `Pr(>|t|)` < p_cutoff & FDR < 1)
  
  write.csv(r_test_deg1, file = file.path(outputdir, paste0('DEG_results_FC',fc_cutoff, '_p', p_cutoff, '_fdr', 1, '.csv')), row.names = F)
  
  # plot
  library(ggrepel) 
  results <- r_test
  #genes of interset
  
  if(nrow(r_test_deg0.25) > 0){
    selected_genes <- r_test_deg0.25 %>%
      filter(!is.na(hallmark_pathways)) %>%
      pull(Gene)
  }else{
    selected_genes <- r_test_deg1 %>%
      filter(!is.na(hallmark_pathways)) %>%
      pull(Gene)
  }
 
  
  # Categorize Results based on P-value & FDR for plotting
  results$Color <- "NS or FC < 0.58"
  results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
  results$Color[results$FDR < 0.25] <- "FDR < 0.25"
  results$Color[results$FDR < 0.05] <- "FDR < 0.05"
  results$Color[abs(results$Estimate) < log2(1.5)] <- "NS or FC < 0.58"
  results$Color <- factor(results$Color,
                          levels = c("NS or FC < 0.58", "P < 0.05",
                                     "FDR < 0.25", "FDR < 0.05"))
  
  # pick top genes for either side of volcano to label
  # order genes for convenience:
  results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
  results <- results[, -1*ncol(results)] # remove invert_P from matrix
  
  # Graph results
  
  if(length(selected_genes) > 0){
    ggplot(results, aes(x = Estimate, y = -log10(`Pr(>|t|)`), color = Color, label = Gene)) +
      geom_vline(xintercept = c(0.58, -0.58), lty = "dashed") +
      geom_hline(yintercept = -log10(0.05), lty = "dashed") +
      geom_point() +
      labs(x = results$Contrast[1],
           y = "Significance, -log10(P)",
           color = "Significance") +
      scale_color_manual(values = c(`FDR < 0.05` = "dodgerblue",
                                    `FDR < 0.25` = "lightblue",
                                    `P < 0.05` = "orange2",
                                    `NS or FC < 0.58` = "gray"),
                         guide = guide_legend(override.aes = list(size = 4))) +
      scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
      geom_text_repel(data = subset(results, Gene %in% selected_genes),
                      size = 4, point.padding = 0.15, color = "black",
                      min.segment.length = .1, box.padding = .2, lwd = 2,
                      max.overlaps = 50) +
      theme_bw(base_size = 16) +
      theme(legend.position = "bottom")
  }else{
    ggplot(results, aes(x = Estimate, y = -log10(`Pr(>|t|)`), color = Color, label = Gene)) +
      geom_vline(xintercept = c(0.58, -0.58), lty = "dashed") +
      geom_hline(yintercept = -log10(0.05), lty = "dashed") +
      geom_point() +
      labs(x = results$Contrast[1],
           y = "Significance, -log10(P)",
           color = "Significance") +
      scale_color_manual(values = c(`FDR < 0.05` = "dodgerblue",
                                    `FDR < 0.25` = "lightblue",
                                    `P < 0.05` = "orange2",
                                    `NS or FC < 0.58` = "gray"),
                         guide = guide_legend(override.aes = list(size = 4))) +
      scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
      theme_bw(base_size = 16) +
      theme(legend.position = "bottom")
  }
  
  
  ggsave(file.path(outputdir, paste0('DEG_results_FC',fc_cutoff, '_p', p_cutoff, '_volcano.pdf')), width = 8, height = 8)
  
  # GSEA analysis
  
  gene.converter = bitr(r_test$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  DEG_df2 <- r_test %>%
    left_join(gene.converter, by=c("Gene"="SYMBOL")) %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::arrange(desc(Estimate))
  
  geneList.gsea <- as.numeric(DEG_df2[,"Estimate"])
  names(geneList.gsea) <- as.character(DEG_df2[,"ENTREZID"])
  
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)
  
  
  mygsea <- GSEA(geneList.gsea, TERM2GENE = m_t2g, pvalueCutoff = 0.25, pAdjustMethod = "BH")
  
  GSEA.result.df <- mygsea@result %>%
    mutate(Sig = -1*log10(qvalue)) %>%
    arrange(NES) 
  
  save(GSEA.result.df, file = file.path(outputdir, 'GSEA_bar_chart.Rdata'))
  
  ggplot(GSEA.result.df, aes(x=NES, y=reorder(Description, NES), fill = Sig)) + 
    geom_bar(stat = "identity") + 
    scale_fill_viridis_c(option = "plasma", begin =0.4, end = 1) +
    labs(y = "Hallmark pathway",
         x = "normalized enrichment score",
         fill = "-log10(qvalues)") +
    theme(axis.text.y = element_text(face = 'bold'))
  
  ggsave(file.path(outputdir, 'GSEA_bar_chart.pdf'), width = 7, height = 6)
  
}