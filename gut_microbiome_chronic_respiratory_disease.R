#
#  gut_microbiome_respiratory_disease
#
# 
#--------------------------------
# 
# taxon-specific Cox
# 
#--------------------------------

require("survival")
require("data.table")


# ds <- 'ASTHMA'
# ds <- 'COPD'


disease.name <- paste0("INCIDENT_",ds);disease.name
time2event <- paste0(gsub(x = disease.name, pattern = "INCIDENT_", replacement = ""), "_AGEDIFF");time2event


coxres <- 
  lapply(c('P','C','O','F','G','S'), function(rk){
    
    
    x <- readRDS(paste0('file_',rk,'.rds'))
    factorcol <- c('MEN')
    x$metadata[, (factorcol) := lapply(.SD, as.factor), .SDcols = factorcol]
    formulas <- sapply(x$taxanames, function(taxon) as.formula(paste('Surv(',time2event, ',', disease.name,')~', taxon,'+BL_AGE+MEN')))
    mds <- lapply( formulas, function(fml){coxph(fml, data = x$metadata)})
    phtest <- lapply( mds, function(md){cox.zph(md)})
    
    
    res <- lapply(mds,
                  function(x){ 
                    x <- summary(x)
                    P_value <- x$coef[1,"Pr(>|z|)"]
                    beta <- x$coef[1,"coef"]  
                    expbeta <- x$coef[1,"exp(coef)"]
                    lower_95 <- x$conf.int[1,"lower .95"]  
                    upper_95 <- x$conf.int[1,"upper .95"]     
                    HR <- paste0(round(expbeta, digits=2), " (", 
                                 round(lower_95,2), "-", 
                                 round(upper_95,2), ")")
                    output <- data.table(rownames( x$coef)[1], P_value, beta, HR, expbeta,lower_95,upper_95)
                    names(output)<-c("Taxa","P_value", "Beta", "HR", "expbeta","lower_95","upper_95")
                    return(output)
                  })
    
    
    res <- rbindlist(res)
    res <- res[order(-abs(Beta))]
    res[,phtest:= unlist(lapply(phtest, function(test){test$table['GLOBAL','p']}))]
    res[,FDR:= p.adjust(p = P_value, method = "fdr")]
    return(res)
  })

names(coxres) <- c('P','C','O','F','G','S')
bindres <- rbindlist(coxres,idcol = 'TaxoLevel')
setcolorder(x = bindres, neworder = c("TaxoLevel", "Taxa","Beta", "HR", "P_value", "FDR"))
bindres[,.(TaxoLevel,Taxa,Beta,HR,P_value,FDR)]






#--------------------------------
# 
# data partitioning
# 
#--------------------------------


loadpkg <- list("caret", "data.table")
lapply(loadpkg, require, character.only = T)



data_partition <- function(dt, target_name, trainFrac=0.7, seed=12345, nfold=5) {
  
  
  setDT(dt)
  FoldIndex <- caret::createDataPartition(y=dt[, get(target_name)], times=1, p=trainFrac, list=FALSE)
  
  
  train <- dt[ FoldIndex]
  test  <- dt[-FoldIndex]
  
  
  message(paste0(target_name, ":"))
  message("class label distribution")
  
  
  message("In training set: ")
  print(prop.table(table(train[,get(target_name)])))
  
  
  message("In testing set: ")
  print(prop.table(table(test[,get(target_name)]) ))
  
  
  train.ind <- train[, ind]
  test.ind <- test[, ind]
  
  
  # fold indices
  set.seed(seed)
  FoldIndex <- createFolds(train[, get(target_name)], k = nfold, list = F)
  set.seed(seed)
  FoldIndexList <- createFolds(train[, get(target_name)], k = nfold, list = T)
  
  output <- list(train.ind=train.ind, 
                 test.ind=test.ind,
                 train.FoldIndex=FoldIndex,
                 train.FoldIndexList=FoldIndexList)
  return(output) 
  
}





#--------------------------------
# 
# feature - univariate analysis
# 
#--------------------------------



unlink(".RData")
packages <- c("data.table","survival")
invisible(lapply(packages, require, character.only = T))




factorcol <- c('MEN','INCIDENT_ASTHMA','INCIDENT_COPD')
train[, (factorcol) := lapply(.SD, as.factor), .SDcols = factorcol]
lapply(train[,.(BL_AGE,MEN,INCIDENT_ASTHMA,INCIDENT_COPD)], summary)



# Spearman 
spearman <- sapply(taxanames, function(nm){
  cor.test(as.numeric(as.character(train[,get(disease.name)])), 
           train[,get(nm)], method = "spearman", exact=FALSE)
}) 

spearman.pval <-  apply(spearman, 2, (function(x){ x$p.value }) )  
spearman.rho <-  apply(spearman, 2, (function(x){ x$estimate }) )  
spearman.pval <- as.data.table(spearman.pval, keep.rownames="taxa")
spearman.rho <- as.data.table(spearman.rho, keep.rownames="taxa")
spearmanCorr <- merge(spearman.pval, spearman.rho, sort = F)




# logistic regression
var.adjust <- c("MEN","BL_AGE")

formulas_adjusted <- sapply(taxanames, function(nm){
  as.formula(paste0(disease.name, "~",  paste(c(nm, var.adjust), collapse = "+")))
})


models_adjusted <- lapply(formulas_adjusted, function(form){
  md <- glm(form, data = train, family = "binomial")
})


beta <- sapply(models_adjusted, function(model){
  coef(summary(model))[2,1]
})
beta<- as.data.table(beta, keep.rownames="taxa")


pval <- sapply(models_adjusted, function(model){
  coef(summary(model))[2,4]
})
pval<- as.data.table(pval, keep.rownames="taxa")


OR <- sapply(models_adjusted, function(model){
  exp(coef(summary(model))[2,1])
})
OR<- as.data.table(OR, keep.rownames="taxa")


mergefunc = function(x,y) merge(x, y, by="taxa", sort = F) 
logistic_adjusted <- Reduce(mergefunc,list(pval, beta, OR))
setcolorder(logistic_adjusted, c("taxa","pval"))




#  cox regression 
var.adjust <- c("MEN","BL_AGE")

formulas_adjusted <- sapply(taxanames, function(nm){
  fml <- as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", paste(c(nm, var.adjust), collapse = "+")))
})


models_adjusted <- lapply(formulas_adjusted, function(form){
  sv <- coxph(form, data = train)
})

phtest <- lapply(models_adjusted, function(md){
  cox.zph(md)
})

beta <- sapply(models_adjusted, function(model){
  summary(model)$coef[1,1]
})
beta<- as.data.table(beta, keep.rownames="taxa")


pval <- sapply(models_adjusted, function(model){
  summary(model)$coefficients[1,5]
})
pval<- as.data.table(pval, keep.rownames="taxa")


waldpval <- sapply(models_adjusted, function(model){
  summary(model)$wald["pvalue"][[1]]
})
waldpval<- as.data.table(waldpval, keep.rownames="taxa")


HR <- sapply(models_adjusted, function(model){
  summary(model)$coef[1,2] 
})
HR<- as.data.table(HR, keep.rownames="taxa")

# HR.lower <- sapply(models_adjusted, function(model){
#   summary(model)$conf.int[1,"lower .95"]
# })
# HR.lower<- as.data.table(HR.lower, keep.rownames="taxa")
# 
# HR.upper <- sapply(models_adjusted, function(model){
#   summary(model)$conf.int[1,"upper .95"]
# })
# HR.upper<- as.data.table(HR.upper, keep.rownames="taxa")


mergefunc = function(x,y) merge(x, y, by="taxa", sort = F) 
cox_adjusted <- Reduce(mergefunc,list(beta, pval, HR))
setcolorder(cox_adjusted, c("taxa","pval"))



#--------------------------------
# 
# xgboost
# 
#--------------------------------



pkgnames <- c( "data.table","mltools","xgboost","smoof","mlrMBO","rgenoud")
for(pkg in pkgnames){
  if(!require(pkg,character.only = TRUE)) {
    install_decision <- readline(prompt=paste0("Package ",pkg, " hasn't been found. Install package? Y or N: "))
    if(install_decision == "Y"){
      install.packages(pkg)
    }
  }
  invisible(require(pkg, character.only = TRUE))
}



xgb_bayes     <-           function(train, 
                                    target_name, 
                                    feature_names, 
                                    tree_method = "hist",  
                                    grow_policy = "depthwise", # depthwise, lossguide
                                    eval_metric = "auc", 
                                    ParamSet = makeParamSet(
                                      makeNumericParam("eta",              lower = 0.001, upper = 0.5),
                                      makeIntegerParam("max_depth",        lower= 3,      upper = 10),
                                      makeIntegerParam("min_child_weight", lower= 5,      upper = 100),
                                      makeNumericParam("subsample",        lower = 0.8,   upper = 0.95),
                                      makeNumericParam("colsample_bytree", lower = 0.6,   upper = 0.85),
                                      makeNumericParam("lambda",           lower = 1e-4,     upper = 1),  
                                      makeNumericParam("gamma",            lower = 0,     upper = 5),  
                                      makeNumericParam("max_delta_step",   lower = 0,     upper = 8)),
                                    nround=200, 
                                    early_stopping_rounds=10, 
                                    nFolds = 5,      
                                    FoldList = NULL, 
                                    savepath = NULL, 
                                    savename = NULL,
                                    genPnts = 40,  
                                    bayesIter = 40, 
                                    seed = NULL
){
  
  
  
  if(!is.null(savepath)){
    dir.create(savepath, showWarnings = TRUE, recursive = TRUE)
  }
  
  
  setDT(train)
  train <- train[, c(feature_names,target_name), with = F]
  
  
  factor_col <- colnames(train)[which(as.vector(train[, lapply(.SD, is.factor), .SDcols = feature_names])==TRUE)]
  if (length(factor_col)>0){
    train <- mltools::one_hot(train, cols = factor_col, 
                              naCols = TRUE, 
                              dropCols = TRUE, 
                              dropUnusedLevels = FALSE)   
    
    feature_names <- setdiff(colnames(train), target_name)
  }
  
  
  dtrain = xgb.DMatrix(data = as.matrix(train[,..feature_names]), label=as.matrix(train[ ,get(target_name)]))
  
  
  obj  <- smoof::makeSingleObjectiveFunction(
    name = "xgbcv",
    fn =   function(x){
      set.seed(seed)
      cv <- xgb.cv(
        params = list(
          booster          = "gbtree",
          eta              = x["eta"],
          max_depth        = x["max_depth"],            
          min_child_weight = x["min_child_weight"], 
          max_delta_step   = x["max_delta_step"],    
          subsample        = x["subsample"],            
          colsample_bytree = x["colsample_bytree"],     
          lambda           = x["lambda"],   
          gamma            = x["gamma"],               
          # alpha            = x["alpha"],
          # lambda_bias      = x["lambda_bias"],        
          # colsample_bylevel= x["colsample_bylevel"],  
          # colsample_bynode = x["colsample_bynode"],  
          # scale_pos_weight = x["scale_pos_weight"],     
          # max_leaves       = x["max_leaves"],
          objective        = 'binary:logistic', 
          eval_metric      = eval_metric),
        data             = dtrain,
        tree_method      = tree_method,
        grow_policy      = grow_policy,
        nfold            = nFolds, 
        folds            = FoldList,  
        nround           = nround,
        early_stopping_rounds = early_stopping_rounds,     
        prediction       = TRUE,
        verbose          = 1, 
        print_every_n    = 5)
      cv$evaluation_log[, max(test_auc_mean)]
    },
    description = "gbtree booster", 
    par.set = ParamSet,
    minimize = FALSE
  )
  
  
  des = generateDesign(n=genPnts,
                       par.set = getParamSet(obj), 
                       fun = lhs::randomLHS)  
  
  
  control = makeMBOControl()
  control = setMBOControlTermination(control, iters = bayesIter)
  
  
  runBayes = mbo(fun = obj, 
                 design = des,  
                 control = control, 
                 show.info = TRUE)
  
  
  message("AUC:")
  message(runBayes$y)
  
  
  
  if(!is.null(savepath) & !is.null(savename)){
    saveRDS(runBayes, file = file.path(savepath, savename))
  }
  
  
  return(runBayes)
  
}




#--------------------------------
# 
# glmnet
# 
#--------------------------------


loadpkg <- list("data.table", "glmnet", "doParallel", "foreach")
invisible(lapply(loadpkg, require, character.only = T))



run_glmnet <- function(train, 
                       target_name, 
                       predictor_names, 
                       alpha = 0, 
                       lambda = exp(seq(log(1e-4), log(200), length.out=200)), 
                       params=list(type.measure="auc", family = "binomial", maxit =  10^4),
                       nfolds=5,
                       foldid = foldid,
                       pe=TRUE, 
                       seed=NULL
){
  
  
  setDT(train)
  
  
  use_parallel <- FALSE
  if(pe){  
    registerDoParallel(cores = nfolds)  
    use_parallel <- TRUE
  }
  
  
  options(na.action="na.pass")
  set.seed(seed)
  
  
  
  cvglmnet <- do.call(cv.glmnet,
                      c(list(x = model.matrix( ~ .-1, train[,predictor_names, with=F]),
                             y = train[,get(target_name)], 
                             parallel = use_parallel,
                             alpha = alpha,  
                             nfolds = nfolds, 
                             foldid = foldid,
                             lambda = lambda), params))
  
  
  
  coef.1se <- coef(cvglmnet, s = "lambda.1se")  
  coef.1se <- data.frame(coef.1se=coef.1se[which(coef.1se!=0),])
  setDT(coef.1se, keep.rownames="feature")
  coef.1se <- coef.1se[order(-abs(coef.1se))]
  
  
  coef.min <- coef(cvglmnet, s = "lambda.min")  
  coef.min <- data.frame(coef.min=coef.min[which(coef.min!=0),])
  setDT(coef.min, keep.rownames="feature")
  coef.min <- coef.min[order(-abs(coef.min))]
  
  
  tmp <- list(cvglmnet=cvglmnet,
              lambda.min = cvglmnet$lambda.min,
              lambda.1se = cvglmnet$lambda.1se,
              coef.1se = coef.1se,
              coef.min = coef.min) 
  
  return(tmp)
}





#--------------------------------
# 
#  Cox analysis in validation dataset
# 
#--------------------------------


require(data.table)
require(survival)
require(Publish)




factorcol <- c('MEN','CURR_SMOKE','EX_SMOKE','smoking','age_60plus','bmi_25plus','INCIDENT_ASTHMA','INCIDENT_COPD')
metadata[, (factorcol) := lapply(.SD, as.factor), .SDcols = factorcol]
disease.name <- paste0("INCIDENT_",ds);disease.name
time2event <- paste0(gsub(x = disease.name, pattern = "INCIDENT_", replacement = ""), "_AGEDIFF");time2event
vars <- c(disease.name,time2event,"MEN","BL_AGE","BMI","CURR_SMOKE","EX_SMOKE",'smoking','age_60plus','bmi_25plus')
lapply(metadata[,..vars], summary)
metadata <- na.omit(metadata, cols=vars, invert=FALSE)



# set ref level
setRefVar <- c("MEN", "age_60plus", "bmi_25plus","CURR_SMOKE","EX_SMOKE")
metadata[, (setRefVar) := lapply(.SD, relevel, ref="0"), .SDcols = setRefVar]
metadata[, smoking:=relevel(x = metadata$smoking, ref = "Never smoker")]



# univariate
vars <-  c("scale(score)","MEN", "BL_AGE", "BMI","CURR_SMOKE")
univ_formulas <- sapply(vars,function(x) as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", x )))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = metadata)})
univ_zph <- lapply( univ_models, function(x){cox.zph(x)})
univ_stats <- lapply(univ_models, function(sv){
  res.tmp <- as.data.table(summary(sv)$conf.int, keep.rownames=T)
  res.tmp[, pvalue:=summary(sv)$coefficients[,"Pr(>|z|)"]]
  res.tmp[, concordance:=summary(sv)$concordance[['C']]]
  return(res.tmp)
})
univ_stats <- rbindlist(univ_stats)
univ_stats[pvalue<=0.001, P:="<=0.001"]
univ_stats[pvalue>0.001, P:=as.character(round(pvalue,3))]
univ_stats[,HR:=paste0( round(get('exp(coef)'),2), " (", round(get('lower .95'),2), "-", round(get('upper .95'),2), ")")]



# multiple
multiv_formulas <- list(
  md1 = as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", 
                          paste0(c("MEN","BL_AGE","BMI"), collapse = "+") )),
  md2 = as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", 
                          paste0(c("MEN","BL_AGE","BMI","CURR_SMOKE"), collapse = "+") )),
  md3 = as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", 
                          paste0(c("MEN","BL_AGE","BMI","CURR_SMOKE","scale(score)"), collapse = "+") ))
)

multiv_models <- lapply( multiv_formulas, function(x){coxph(x, data = metadata)})
multiv_zph <- lapply( multiv_models, function(x){cox.zph(x)})
multiv_stats <- lapply(multiv_models, function(sv){
  res.tmp <- as.data.table(summary(sv)$conf.int, keep.rownames=T)
  res.tmp[,pvalue:=summary(sv)$coefficients[,"Pr(>|z|)"]]
  res.tmp[, concordance:=summary(sv)$concordance[['C']]]
  res.tmp[pvalue<=0.001, P:="<=0.001"]
  res.tmp[pvalue>0.001, P:=as.character(round(pvalue,3))]
  res.tmp[,HR:=paste0( round(get('exp(coef)'),2), " (", round(get('lower .95'),2), "-", round(get('upper .95'),2), ")")]
  return(res.tmp)
})
multiv_stats



# age_60plus
fml <- list(as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", paste0('scale(score) * age_60plus')  )),
            as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", paste0('scale(score) * age_60plus + MEN + CURR_SMOKE + BMI')  ))
)
md <- lapply(fml, function(x){coxph(x, data = metadata)})
phtest <- lapply(md, function(x){cox.zph(x)})
res <- lapply(md, function(x){publish(x)})


# SEX
fml <- list(as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", paste0('scale(score) * MEN')  )),
            as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", paste0('scale(score) * MEN + BL_AGE + CURR_SMOKE + BMI')  ))
)
md <- lapply(fml, function(x){coxph(x, data = metadata)})
phtest <- lapply(md, function(x){cox.zph(x)})
res <- lapply(md, function(x){publish(x)})


# bmi_25plus
fml <- list(as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", paste0('scale(score) * bmi_25plus')  )),
            as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", paste0('scale(score) * bmi_25plus + MEN + CURR_SMOKE + BMI')  ))
)
md <- lapply(fml, function(x){coxph(x, data = metadata)})
phtest <- lapply(md, function(x){cox.zph(x)})
res <- lapply(md, function(x){publish(x)})


# CURR_SMOKE
fml <- list(as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", paste0('scale(score) * CURR_SMOKE')  )),
            as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", paste0('scale(score) * CURR_SMOKE + BL_AGE + BMI + MEN')  ))
)
md <- lapply(fml, function(x){coxph(x, data = metadata)})
phtest <- lapply(md, function(x){cox.zph(x)})
res <- lapply(md, function(x){publish(x)})


# non-SMOKE
nonsmoke <- metadata[CURR_SMOKE==0]
setRefVar <- c("MEN","EX_SMOKE")
nonsmoke[, (setRefVar) := lapply(.SD, relevel, ref="0"), .SDcols = setRefVar]


fml <- as.formula(paste0("Surv(", time2event, ",as.numeric( as.character(", disease.name, "))) ~ ", 
                         paste0(c("scale(score)*EX_SMOKE","MEN","BL_AGE","BMI"), collapse = "+") ))
md <- coxph(fml, data = nonsmoke)
phtest <- cox.zph(md)
res <- as.data.table(summary(md)$conf.int, keep.rownames=T)
res[,pvalue:=summary(md)$coefficients[,"Pr(>|z|)"]]
res[, concordance:=summary(md)$concordance[['C']]]
res[pvalue<=0.001, P:="<=0.001"]
res[pvalue>0.001, P:=as.character(round(pvalue,3))]
res[,HR:=paste0( round(get('exp(coef)'),2), 
                 " (", round(get('lower .95'),2), "-", round(get('upper .95'),2), 
                 ")")]
resPub <- publish(md)









