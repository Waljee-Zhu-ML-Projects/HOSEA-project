mmean <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    mean(x,na.rm=TRUE)
  }
}

mmax <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    max(x,na.rm=TRUE)
  }
}

mmin <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    min(x,na.rm=TRUE)
  }
}

ssum <- function(x){
  if(all(is.na(x))){
    NA
  }
  else{
    sum(x,na.rm=TRUE)
  }
}

# apply median imputation
fill_by_median <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}
# apply zero imputation
fill_by_zero <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}

expit <- function(x){
  exp(x) / (1+exp(x))
}

