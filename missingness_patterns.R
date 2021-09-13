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
n_groups = 5

library(tidyverse)
library(purrrlyr)
# build tree
pattern_tree = function(df, n_groups){
  # observed patterns
  colnms = colnames(df)
  patterns = df %>%
    filter_at(vars(all_of(colnms)), any_vars(is.na(.))) %>%
    is.na %>% (function(x) !x) %>% as_tibble
  
  # counts per pattern
  counts = patterns %>%
    group_by_all() %>% summarise(count = n(), .groups="drop")
  common = counts %>% summarise_all(all)
  
  # keep largest
  large_patterns = counts %>% slice_max(order_by=count, n=n_groups)
  common$count = sum(counts$count) - sum(large_patterns$count)
  large_patterns[nrow(large_patterns)+1, ] = common
  x = large_patterns %>% select(!count)
  y = large_patterns %>% select(count)
  
  #
  large_patterns$total = 0
  for(i in range(nrow(large_patterns))){
    pattern = x[i, ]
    large_patterns[i, "total"] = sum(apply(patterns, 1, function(row) all(row >= pattern)))
  }
  
  return(large_patterns)
}

large_patterns = pattern_tree(df, 20)
