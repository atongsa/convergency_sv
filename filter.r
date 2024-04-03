#!/usr/bin/env Rscript

df_val <- read_tsv(val_f, col_names=T) %>% replace(is.na(.), 0)
df_sp <- read_tsv(sp_f, col_names=T)

df_val_m <- df_val 
df_val_t <- as.data.frame(t(df_val_m)) %>% filter(V1 < 0.05)
df_val_s <- df_val[,c('Sample', rownames(df_val_t))]

df_val_l <- w2l(df_val_s) 
df_sp_l <- w2l(df_sp)

i_j <- inner_join(df_val_l, df_sp_l, by=c('Sample', 'wind'))
i_j_r <- i_j %>% select(wind) %>% unique()

