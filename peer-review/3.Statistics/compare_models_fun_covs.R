compare_models_fun_covs <- function(model, meas, hem, trk, trk_sz, dir_in, cov_var, regions){

  # compares models generated with lme4
  
quarters <- c(1:4)
nodes <- c(1:100)

if (meas == 'dti_fa' | meas == 'dti_md') {
  if (trk_sz == 'nodes'){
    brain_names <- paste0('n', nodes, '_', meas)
  } else if (trk_sz == 'quarters'){
    brain_names <- paste0('q', quarters, '_', meas)
  } else if (trk_sz == 'average') {
    brain_names <- meas
  }
} else if (meas == 'wm'){
  brain_names <- c('eTIV', paste0(meas, '.', hem, '.', regions),
                   paste0(hem, 'CerebralWhiteMatterVol'), 'CerebralWhiteMatterVol')
} else if (meas == 'area') {
  brain_names <- c('eTIV', paste0(hem, '_', regions, '_', meas), paste0(hem, '_WhiteSurfArea_area'))
} else if (meas == 'thickness') {
  brain_names <- c('eTIV', paste0(hem, '_', regions, '_', meas), paste0(hem, '_MeanThickness_thickness'))
}  else {
  brain_names <- c('eTIV', paste0(hem, '_', regions, '_', meas))
}

model_names <- c("Covariates", "No Covariates")

bic_comp <- matrix(ncol=length(model_names), nrow=length(brain_names))
p_comp <- matrix(ncol=length(model_names), nrow=length(brain_names))

for (i in (1:length(brain_names))){
  
  print(brain_names[i])
  if (meas == 'dti_fa' | meas == 'dti_md'){
    load(paste0(dir_in, '/M1_log_', trk, '_', brain_names[i], '.RData'))
    load(paste0(dir_in, '/M1_log_no_', cov_var, '_', trk, '_', brain_names[i], '.RData'))
  } else {
    load(paste0(dir_in, '/M1_', model, '_', brain_names[i], '.RData'))
    load(paste0(dir_in, '/M1_', model, '_no_', cov_var, '_', brain_names[i], '.RData'))
  }
  #mods_comp <- anova(M1_log, M1_log_no_cov_var)
  mods_comp <- anova(eval(as.symbol(paste0('M1_', model))), eval(as.symbol(paste0('M1_', model, '_no_cov_var'))))
  print(mods_comp)
  bic_comp[i,] <- mods_comp$BIC
  p_comp[i,] <- mods_comp$`Pr(>Chisq)`
}

if ( meas == 'dti_fa' | meas == 'dti_md' ){
  sz_in <- 1
  sz_out <- length(brain_names)
} else {
  sz_in <- 2
  sz_out <- length(regions) + 1
}

all_comp <- cbind(bic_comp[sz_in:sz_out,],p_comp[sz_in:sz_out,])
if (meas != 'dti_fa' & meas != 'dti_md') {
  write.csv(all_comp, paste0(dir_in, '/', hem, '_', meas, '_mod_comp_no_', cov_var, '.csv'))
} else {
  write.csv(all_comp, paste0(dir_in, '/', hem, '_', meas, '_', trk, '_mod_comp_no_', cov_var, '.csv'))
}

}