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
