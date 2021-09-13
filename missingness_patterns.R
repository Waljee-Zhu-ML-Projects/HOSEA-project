# exploring missingness patterns within the data
# to check whether it would be viable to do reduced models

# some dataset for development: dplyr::starwars
# df = starwars
missingness_patterns = function(df){
  colnms <- colnames(df)
  counts = df %>%
    filter_at(vars(all_of(colnms)), any_vars(is.na(.))) %>%
    is.na %>% as_tibble %>%
    group_by_all() %>% summarise(count = n())
  return(counts)
}

# load data
load('R_data/subsample/sub_complete_data_impute.RData')

# check patterns
counts = missingness_patterns(train_data_impute$clean)
counts = counts %>% arrange()

df = starwars 

library(tidyverse)
# build tree
counts = missingness_patterns(df)
pattern_tree = function(df){
  # observed patterns
  patterns = df %>%
    filter_at(vars(all_of(colnms)), any_vars(is.na(.))) %>%
    is.na %>% (function(x) !x) %>% as_tibble
  # diana clustering
  dmat = stats::dist(patterns, method="binary")
  hc = cluster::diana(dmat, diss=T)
  plot(hc)
  rect.hclust(hc, h=2, border = "red")
  rect.hclust(hc, h=1.5, border = "blue")
  rect.hclust(hc, h=1.2, border = "green")
  rect.hclust(hc, h=1, border = "orange")
  
  # counts per patter
  counts = patterns %>%
    group_by_all() %>% summarise(count = n())
  x = counts %>% select(!"count")
  y = counts %>% select("count")
  root = sapply(x, prod) > 0
}

split = root
grow_tree = function(x, y, split, candidates){
  sapply(candidates, function(col){
    proposal = split
    proposal[col] = TRUE
    left = apply(x, 1, function(x) superset(x, proposal))
    right = !left
    total_left = y %>% filter(left)
    total = sum(y)
  })
}

superset = function(x, y){ # is x a superset of y?
  return(all(x>y))
}
