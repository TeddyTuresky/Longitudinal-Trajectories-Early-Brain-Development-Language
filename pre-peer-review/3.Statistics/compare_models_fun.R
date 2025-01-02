compare_models_fun <- function(meas, hem, trk, trk_sz, dir_in, regions){

  # compares models generated with lme4
  
quarters <- c(1:4)
nodes <- c(1:100)

if (meas == 'dti_fa' | meas == 'dti_md') {
  if (trk_sz == 'quarters'){
    brain_names <- paste0('q', quarters, '_', meas)
  } else if (trk_sz == 'nodes'){
    brain_names <- paste0('n', nodes, '_', meas)
  } else if (trk_sz == 'average') {
    brain_names <- paste0(trk, '_', meas)
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

model_names <- c("Logarithmic: Random Intercepts + Slopes", "Quadratic: Random Intercepts + Slopes", 
              "Linear: Random Intercepts + Slopes", "Logarithmic: Random Intercepts",
              "Quadratic: Random Intercepts", "Linear: Random Intercepts")

aic_comp <- matrix(ncol=length(model_names), nrow=length(brain_names))
p_comp <- matrix(ncol=length(model_names), nrow=length(brain_names))

for (i in c(1:length(brain_names))){
  
  print(brain_names[i])
  load(paste0(dir_in, '/M1_log_', brain_names[i], '.RData'))
  load(paste0(dir_in, '/M1_sq_', brain_names[i], '.RData'))
  load(paste0(dir_in, '/M1_lin_', brain_names[i], '.RData'))
  load(paste0(dir_in, '/M1_log_no_slope_', brain_names[i], '.RData'))
  load(paste0(dir_in, '/M1_sq_no_slope_', brain_names[i], '.RData'))
  load(paste0(dir_in, '/M1_lin_no_slope_', brain_names[i], '.RData'))
  mods_comp <- anova(M1_log, M1_sq, M1_lin, M1_log_no_slope, M1_sq_no_slope, M1_lin_no_slope)
  print(mods_comp)
  aic_comp[i,] <- mods_comp$AIC
  p_comp[i,] <- mods_comp$`Pr(>Chisq)`
}

if ( meas == 'dti_fa' | meas == 'dti_md' ){
  sz_in <- 1
  sz_out <- length(brain_names)
} else {
  sz_in <- 2
  sz_out <- length(regions) + 1
}

all_comp <- cbind(aic_comp[sz_in:sz_out,],p_comp[sz_in:sz_out,])
if (meas != 'dti_fa' & meas != 'dti_md') {
  write.csv(all_comp, paste0(dir_in, '/', hem, '_', meas, '_mod_comp.csv'))
} else {
  write.csv(all_comp, paste0(dir_in, '/', hem, '_', meas, '_', trk, '_mod_comp.csv'))
}

}