library(tidyverse)

t2g <- read_tsv('/data/mcgaugheyd/projects/outside/chick_MGT/ref/t2g.txt', col_names = FALSE)
t2g_u <- t2g %>% select(X1, X3, X2, X4, X5, X6, X7, X8) %>%  mutate(X3 = case_when(is.na(X3) ~ X2, TRUE ~ X3))
write_tsv(t2g_u, file = '/data/mcgaugheyd/projects/outside/chick_MGT/ref/t2g_named.txt', col_names = FALSE)
