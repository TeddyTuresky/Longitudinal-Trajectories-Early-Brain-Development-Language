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
