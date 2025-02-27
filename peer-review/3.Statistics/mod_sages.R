mod_sages <- function(model, dat){
  
  
  if (model == 'lin' | model == 'quad'){
    ages_all <- na.omit(dat$age)
    sages <- (c(ages_all, 50) - 40) / sd(ages_all)
    sage_50 <- sages[length(sages)]
  } else if (model == 'log'){
    sage_50 <- log(50+1)
  } else if (model == 'asym'){
    sage_50 <- 50
  }
  
  return(sage_50)
  
}