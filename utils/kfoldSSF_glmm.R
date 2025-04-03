kfoldSSF <- function(mod, strata_name = "unique_step", k = 5, nrepet = 100, jitter = FALSE,
                     reproducible = FALSE, details = FALSE)
{
  require(glmmTMB)
  require(tidyverse)
  
  # Get model frame
  dt <- model.frame(mod)
  
  # Rename response variable
  names(dt)[1] <- "case"
  
  # Function for updating model fit
  modUpdate <- function(mod, data){
    mod_up <- update(mod, case ~ ., data = data)
    mod_up$parameters$theta[1] <- log(1e3)
    ntheta = length(mod_up$parameters$theta)
    mod_up$mapArg <- list(theta=factor(c(NA, rep(1,ntheta-1))))
    mod_up_fit <- glmmTMB:::fitTMB(mod_up)
    return(mod_up_fit)
  }
  
  kfold <- rd <- warn <- numeric(length = nrepet)
  if (details)
    dbg <- list()
  
  # Loop through folds and refit model
  print(paste("start:", Sys.time()))
  for (i in 1:nrepet) {
    print(paste("loop position:", i, "of", nrepet))
    
    flag = 1
    while(flag == 1){
      dt$sets <- "train"
      
      # Take a random subset of strata (1/k of total) and set as test set. There rest are train set
      dt$sets[dt[, strata_name] %in% sample(unique(dt[, strata_name]),
                                            length(unique(dt[, strata_name]))/k)] <- "test"
      
      # Refit model with training set only
      reg = NULL
      tryCatch({
        reg <- modUpdate(mod, data = filter(dt, sets == 'train'))
      }, error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")}
      )
      
      if(!is.null(reg)) flag = flag +1
      
      print(paste("while loop flag = ", flag))
    }
    
    # Predict on test set
    dtest <- droplevels(subset(dt, sets == "test"))
    coefs <- fixef(reg)$cond
    if("dev:slope" %in% names(coefs)) dtest$`dev:slope` <- dtest$dev * dtest$slope
    covs <- dtest %>% select(names(coefs)) %>% as.matrix
    dtest$predall <- exp(covs %*% coefs)
    
    if (jitter) {
      if (reproducible)
        set.seed(i)
      dtest$predall <- jitter(dtest$predall)
    }
    samplepred <- function(df) {
      nrand <- sum(df$case == 0)
      obs <- rank(df$predall)[df$case == 1]
      if (reproducible)
        set.seed(i)
      rand <- sample(rank(df$predall[df$case == 0]),
                     1)
      return(data.frame(obs = obs, rand = rand, nrand = nrand))
    }
    
    colTest <- dtest %>% group_by(!!sym(strata_name)) %>% summarise(n1 = length(which(case==1)), n0 = length(which(case==0)))
    if(length(which(colTest$n1 < 1)) > 0 | length(which(colTest$n0 < 1)) > 0) next
    
    ranks <- do.call(rbind, by(dtest, dtest[, strata_name], samplepred))
    nrand <- unique(ranks$nrand)
    if (length(nrand) != 1) {
      nrand <- max(nrand)
      warn[i] <- 1
    }
    kfold[i] <- cor(1:(nrand+1), table(factor(ranks$obs,
                                              levels = 1:(nrand+1))), method = "spearman")
    rd[i] <- cor(1:(nrand), table(factor(ranks$rand, levels = 1:(nrand))),
                 method = "spearman")
    if (details)
      dbg[[i]] <- ranks
  }
  res <- data.frame(kfold = c(kfold, rd), type = rep(c("obs",
                                                       "rand"), each = nrepet))
  if (details)
    attr(res, "details") <- dbg
  if (sum(warn) > 0)
    warning(paste("The number of controls was not equal among stratas for",
                  sum(warn), "repetitions. Correlations might be biased.",
                  "Use 'details = TRUE' to get more details."))
  print(paste("end:", Sys.time()))
  return(res)
}
