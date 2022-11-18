
loadpkg <- list("caret", "data.table")
lapply(loadpkg, require, character.only = T)




data_partition <- function(dt, target_name, trainFrac, seed,  list=FALSE) {
  
  
  
  setDT(dt)
  set.seed(seed)
  trainIndex <- caret::createDataPartition(y=dt[, get(target_name)], times=1, p=trainFrac, list=list)
  
  
  #  class label distribution 
  message(paste0(target_name, ":"))
  message("class label distribution")
  
  
  message("In training set: ")
  print(prop.table(table(dt[,get(target_name)][trainIndex])))
  
  
  message("In testing set: ")
  print(prop.table(table(dt[,get(target_name)][-trainIndex])))
  
  
  # rm NaN label
  train <- dt[ trainIndex][!is.na(get(target_name))]
  test  <- dt[-trainIndex][!is.na(get(target_name))]

  
  # check
  if (!identical(train[,get(target_name)] , dt[,get(target_name)][trainIndex] ) ){
    warning("NaN label removed in training")
  }
  if(!identical( test[,get(target_name)] , dt[,get(target_name)][-trainIndex] ) ){
    warning("NaN label removed in testing")
  }
  
  
  output.list <- list(train=train, test=test)
  message("Returned a list of 2 data.table elements: train, test")
  return(output.list) 
  
}


