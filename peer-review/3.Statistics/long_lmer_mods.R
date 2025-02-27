long_lmer_mods <- function(dir_in, model, meas, trk, dat, cov_str, cov_var, brain_name){
  
  # Runs linear mixed effects models and outputs results
  
  if (model == 'lin'){
  
    # Random intercepts and random slopes model ---- 
    FormulaM <- paste0('brain1 ~ 1 + age + ', cov_str, ' + (1 + age|id)')
    M1_lin <- lmer(formula=FormulaM, data=dat, REML=FALSE, 
                   control=lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=5e5)))
    M1 <- M1_lin
    if (meas != 'dti_fa' & meas != 'dti_md') {
      save(M1_lin, file = paste0(dir_in, '/M1_lin_', brain_name, '.RData'))
    } else {
      save(M1_lin, file = paste0(dir_in, '/M1_lin_', trk, '_', brain_name, '.RData'))
    }
    
    # for testing models without random slopes
    FormulaM <- paste0('brain1 ~ 1 + age + ', cov_str, ' + (1|id)')
    M1_lin_no_slope <- lmer(formula=FormulaM, data=dat, REML=FALSE, 
                            control=lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=5e5)))
    if (meas != 'dti_fa' & meas != 'dti_md') {
      save(M1_lin_no_slope, file = paste0(dir_in, '/M1_lin_no_slope_', brain_name, '.RData'))
    } else {
      save(M1_lin_no_slope, file = paste0(dir_in, '/M1_lin_no_slope_', trk, '_', brain_name, '.RData'))
    }
    
    # for testing models without HLE or FHD covariates (includes random slopes)
    if (grepl(cov_var, cov_str) == 1){
      cov_str <- str_remove(cov_str, fixed(paste0(' + ', cov_var)))
      FormulaM <- paste0('brain1 ~ 1 + age + ', cov_str, ' + (1 + age|id)')
      M1_lin_no_cov_var <- lmer(formula=FormulaM, data=dat, REML=FALSE, 
                     control=lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=5e5)))
      if (meas != 'dti_fa' & meas != 'dti_md') {
        save(M1_lin_no_cov_var, file = paste0(dir_in, '/M1_lin_no_', cov_var, '_', brain_name, '.RData'))
      } else {
        save(M1_lin_no_cov_var, file = paste0(dir_in, '/M1_lin_no_', cov_var, '_', trk, '_', brain_name, '.RData'))
      }
    }
    
  } else if (model == 'log'){
      
      # Random intercepts and random slopes model ---- 
      FormulaM <- paste0('brain1 ~ 1 + age + ', cov_str, ' + (1 + age|id)')
      M1_log <- lmer(formula=FormulaM, data=dat, REML=FALSE, 
                     control=lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=5e5)))
      M1 <- M1_log      
      if (meas != 'dti_fa' & meas != 'dti_md') {
        save(M1_log, file = paste0(dir_in, '/M1_log_', brain_name, '.RData'))
      } else {
        save(M1_log, file = paste0(dir_in, '/M1_log_', trk, '_', brain_name, '.RData'))
      }
      
      # for testing models without random slopes
      FormulaM <- paste0('brain1 ~ 1 + age + ', cov_str, ' + (1|id)')
      M1_log_no_slope <- lmer(formula=FormulaM, data=dat, REML=FALSE, 
                              control=lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=5e5)))
      
      if (meas != 'dti_fa' & meas != 'dti_md') {
        save(M1_log_no_slope, file = paste0(dir_in, '/M1_log_no_slope_', brain_name, '.RData'))
      } else {
        save(M1_log_no_slope, file = paste0(dir_in, '/M1_log_no_slope_', trk, '_', brain_name, '.RData'))
      }
    
      # for testing models without HLE or FHD covariates (includes random slopes)
      if (grepl(cov_var, cov_str) == 1){
        cov_str <- str_remove(cov_str, fixed(paste0(' + ', cov_var)))
        FormulaM <- paste0('brain1 ~ 1 + age + ', cov_str, ' + (1 + age|id)')
        M1_log_no_cov_var <- lmer(formula=FormulaM, data=dat, REML=FALSE, 
                              control=lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=5e5)))
        if (meas != 'dti_fa' & meas != 'dti_md') {
          save(M1_log_no_cov_var, file = paste0(dir_in, '/M1_log_no_', cov_var, '_', brain_name, '.RData'))
        } else {
          save(M1_log_no_cov_var, file = paste0(dir_in, '/M1_log_no_', cov_var, '_', trk, '_', brain_name, '.RData'))
        }
      }
    
  } else if (model == 'quad'){
  
      # Random intercepts and random slopes model ---- 
      FormulaM <- paste0('brain1 ~ 1 + age + age_sq + ', cov_str, ' + (1 + age|id)')
      M1_sq <- lmer(formula=FormulaM, data=dat, REML=FALSE, 
                    control=lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=5e5)))
      M1 <- M1_sq
      if (meas != 'dti_fa' & meas != 'dti_md') {
        save(M1_sq, file = paste0(dir_in, '/M1_sq_', brain_name, '.RData'))
      } else {
        save(M1_sq, file = paste0(dir_in, '/M1_sq_', trk, '_', brain_name, '.RData'))
      }
      
      # for testing models without random slopes
      FormulaM <- paste0('brain1 ~ 1 + age + age_sq + ', cov_str, ' + (1|id)')
      M1_sq_no_slope <- lmer(formula=FormulaM, data=dat, REML=FALSE, 
                             control=lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=5e5)))
      if (meas != 'dti_fa' & meas != 'dti_md') {
        save(M1_sq_no_slope, file = paste0(dir_in, '/M1_sq_no_slope_', brain_name, '.RData'))
      } else {
        save(M1_sq_no_slope, file = paste0(dir_in, '/M1_sq_no_slope_', trk, '_', brain_name, '.RData'))
      }
      
      # for testing models without HLE or FHD covariates (includes random slopes)
      if (grepl(cov_var, cov_str) == 1){
        cov_str <- str_remove(cov_str, fixed(paste0(' + ', cov_var)))
        FormulaM <- paste0('brain1 ~ 1 + age + age_sq + ', cov_str, ' + (1 + age|id)')
        M1_lin_no_cov_var <- lmer(formula=FormulaM, data=dat, REML=FALSE, 
                              control=lmerControl(optimizer='bobyqa', optCtrl=list(maxfun=5e5)))
        if (meas != 'dti_fa' & meas != 'dti_md') {
          save(M1_lin_no_cov_var, file = paste0(dir_in, '/M1_lin_no_', cov_var, '_', brain_name, '.RData'))
        } else {
          save(M1_lin_no_cov_var, file = paste0(dir_in, '/M1_lin_no_', cov_var, '_', trk, '_', brain_name, '.RData'))
        }
      }
      
  }
  
  return(M1)
  
}