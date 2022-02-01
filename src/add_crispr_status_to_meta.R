library(tidyverse)

bcG15 <- scan('https://raw.githubusercontent.com/davemcg/chick_analysis_collab/main/data/M007_stopCodon.bc', what = 'character') %>% gsub("CB:Z:","",.) %>% unique() %>% paste0(., "_2")
bcG16 <- scan('https://raw.githubusercontent.com/davemcg/chick_analysis_collab/main/data/M008_stopCodon.bc', what = 'character') %>% gsub("CB:Z:","",.) %>% unique() %>% paste0(., "_3")

MD <- combined@meta.data
MD <- MD %>% as_tibble(rownames = 'BC') %>% mutate(CRISPR = case_when(BC %in% c(bcG15, bcG16) ~ 'Yes', TRUE ~ 'No'))
combined@meta.data$CRISPR <- MD$CRISPR

DimPlot(combined, split.by = 'CRISPR')

MD$SALL1pos <- MD %>% mutate(SALL1pos = case_when(BC %in% c(bcG15, bcG16) ~ 'Yes', TRUE ~ 'No'))
