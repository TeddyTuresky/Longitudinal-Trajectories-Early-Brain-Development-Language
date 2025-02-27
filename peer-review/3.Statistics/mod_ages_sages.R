mod_ages <- function(model, dat, sd_age, inv){
  
  
  if (inv == 1){
    if (model == 'lin' | model == 'quad'){
      out <- (dat$age * sd_age) + 40
    } else if (model == 'log'){
      out <- exp((dat$age)) - 1
    } else if (model == 'asym'){
      out <- dat$age
    }
  } else if (inv == 0) {
    if (model == 'lin' | model == 'quad'){
      out <- (dat$age - 40) / sd_age
    } else if (model == 'log'){
      out <- log(dat$age + 1)
    } else if (model == 'asym'){
      out <- dat$age
    }
  }
    
  
  return(out)
}



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