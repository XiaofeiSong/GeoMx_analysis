
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds")) %>%
  group_by(from, to) %>%
  summarise(source = paste(source, collapse = ';'),
            database = paste(database, collapse = ';')) %>%
  unique()

##expr_mat row-sample names, col-gene

geomx_RL_test <- function(expr_mat, r_thred = 0.2){
  
  sender <- expr_mat
  receiver <- expr_mat
  
  ### Correlations of all genes between Sender and receiver 
  LR_correl <- cor(sender, receiver)
  # convert to data frame
  LR_correl <- as.data.frame(LR_correl)
  ## Convert rownames to column
  LR_correl <- rownames_to_column(LR_correl, var = "from")
  ## stack correlations
  correl_stack <- LR_correl %>%
    gather("to", "Correl", -1) %>%
    filter(from != to) 
  ## select RL genes only
  correl_RL <- correl_stack %>%
    inner_join(lr_network)
  
  ## Average Correlation of ligands with target genes
  correl_RL_avg <- correl_RL %>% 
    group_by(from) %>% 
    summarize(Mean.Correl =mean(Correl))
  
  ## annotate the RL pairs with the average information
  
  correl_RL <- correl_RL %>%
    left_join(correl_RL_avg) %>%
    arrange(-Mean.Correl, from, -Correl) %>%
    mutate(pass_thred = if_else(Correl>r_thred, 1, 0)) %>%
    mutate(ligand_target =  paste0(from, '_', to)) %>%
    select(ligand_target, from, to, Mean.Correl, Correl, pass_thred, source, database)
  
  return(correl_RL)
  
}




