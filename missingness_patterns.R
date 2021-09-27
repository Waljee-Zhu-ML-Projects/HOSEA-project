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

df = train_data_impute$clean
df = starwars 
counts = missingness_patterns(df)
n_groups = 10
min_size = 10

library(tidyverse)
library(purrrlyr)
# build tree
pattern_tree = function(df, n_groups=1, min_size=10){
  # observed patterns
  colnms = colnames(df)
  patterns = df %>%
    # filter_at(vars(all_of(colnms)), any_vars(is.na(.))) %>%
    is.na %>% (function(x) !x) %>% as_tibble
  
  # counts per pattern
  counts = patterns %>%
    group_by_all() %>% summarise(count = n(), .groups="drop")
  common = counts %>% summarise_all(all)
  
  # keep largest
  counts[nrow(counts)+1, ] = common
  x = counts %>% select(!count)
  y = counts %>% select(count)
  
  # size of each
  counts$total = 0
  for(i in seq(nrow(counts))){
    pattern = x[i, ]
    counts[i, "total"] = sum(apply(patterns, 1, function(row) all(row >= pattern)))
  }
  
  # keep largests and min size
  large_patterns = counts %>% 
    slice_max(order_by=total, n=n_groups) %>%
    filter(total>min_size)
  
  # assignments
  assignments = matrix(NA, nrow(patterns), nrow(large_patterns))
  for(i in seq(nrow(large_patterns))){
    pattern = large_patterns[i, 1:ncol(df)]
    assignments[, i] = apply(patterns, 1, function(row) all(row >= pattern))
  }
  
  hist(apply(assignments, 1, sum))
  hist(apply(assignments, 2, sum))
  
  return(large_patterns)
}

large_patterns = pattern_tree(df, 20)
