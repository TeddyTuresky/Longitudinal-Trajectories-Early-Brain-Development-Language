long_model_report_stats <- function(meas, hem, trk_sz, dir_in, regions,
                             beh_names, covariates, stats_sum){
  
  # reports stats for dev_struct_dwi analysis
  

  # specify brain names as output by FS or py(Baby)AFQ
  if (meas == 'dti_fa' | meas == 'dti_md') {
    if (trk_sz == 'nodes'){
      nodes <- c(1:100)
      brain_names <- regions <- paste0('n', nodes, '_', meas)
    } else if (trk_sz == 'quarters'){
      quarters <- c(1:4)
      brain_names <- regions <- paste0('q', quarters, '_', meas)
    } else if (trk_sz == 'average'){
      brain_names <- regions <- meas
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
  
  # print stats
  if ( meas == 'dti_fa' | meas == 'dti_md' ){
    sz_in <- 1
    sz_out <- length(brain_names)
  } else {
    sz_in <- 2
    sz_out <- length(regions) + 1
  }
  
  # reduce to regions/tracts of interest
  x_est_covs <- matrix(stats_sum$est_covs[sz_in:sz_out,], nrow = length(regions))
  x_est_beh <- matrix(stats_sum$est_beh[sz_in:sz_out,], nrow = length(regions))
  x_est_ageXbeh <- matrix(stats_sum$est_ageXbeh[sz_in:sz_out,], nrow = length(regions))
  x_est_brain_beh <- matrix(stats_sum$est_brain_beh[sz_in:sz_out,], nrow = length(regions))
  x_se_covs <- matrix(stats_sum$se_covs[sz_in:sz_out,], nrow = length(regions))
  x_p_covs <- matrix(stats_sum$p_covs[sz_in:sz_out,], nrow = length(regions))
  x_p_beh <- matrix(stats_sum$p_beh[sz_in:sz_out,], nrow = length(regions))
  x_p_ageXbeh <- matrix(stats_sum$p_ageXbeh[sz_in:sz_out,], nrow = length(regions))
  x_p_brain_beh <- matrix(stats_sum$p_brain_beh[sz_in:sz_out,], nrow = length(regions))

  print('Uncorrected stats...')
  print(x_est_covs)
  print(x_est_beh)
  print(x_est_ageXbeh)
  print(x_est_brain_beh)
  print(x_se_covs)
  print(x_p_covs)
  print(x_p_beh)
  print(x_p_ageXbeh)
  print(x_p_brain_beh)
  
  if (trk_sz != 'nodes'){
    
  print('pFDR: covariates')
  p_fdr_cov <- p.adjust(as.numeric(c(x_p_covs)), 'fdr') 
  mat_cov <- matrix(p_fdr_cov, nrow = length(regions), ncol = length(covariates))
  print(mat_cov)
  
  print('pFDR: int-behavior')
  p_fdr_beh <- p.adjust(as.numeric(c(x_p_beh)), 'fdr') 
  mat_beh <- matrix(p_fdr_beh, nrow = length(regions), ncol = length(beh_names))
  print(mat_beh)
  
  print('pFDR: slope-behavior')
  p_fdr_ageXbeh <- p.adjust(as.numeric(c(x_p_ageXbeh)), 'fdr') 
  mat_ageXbeh <- matrix(p_fdr_ageXbeh, nrow = length(regions), ncol = length(beh_names))
  print(mat_ageXbeh)
  
  print('pFDR: 50mo-behavior')
  p_fdr_brain_beh <- p.adjust(as.numeric(c(x_p_brain_beh)), 'fdr') 
  mat_brain_beh <- matrix(p_fdr_brain_beh, nrow = length(regions), ncol = length(beh_names))
  print(mat_brain_beh)
  
  print('int-slope correlations')
  print(stats_sum$re_cor)
  
  # Index FDR-corrected regions and behavioral measures by brain measure and hemisphere
  # and their corresponding uncorrected p-values
  indx_cov <- which(mat_cov < 0.05, arr.ind = T)
  indx_beh <- which(mat_beh < 0.05, arr.ind = T)
  indx_ageXbeh <- which(mat_ageXbeh < 0.05, arr.ind = T)
  indx_brain_beh <- which(mat_brain_beh < 0.05, arr.ind = T)
  
  indx_cov_txt <- matrix(nrow = nrow(indx_cov), ncol = 4)
  indx_beh_txt <- matrix(nrow = nrow(indx_beh), ncol = 4)
  indx_ageXbeh_txt <- matrix(nrow = nrow(indx_ageXbeh), ncol = 4)
  indx_brain_beh_txt <- matrix(nrow = nrow(indx_brain_beh), ncol = 4)
  
  if (nrow(indx_cov_txt) > 0){
    indx_cov_txt[,1] <- regions[indx_cov[,1]]
    indx_cov_txt[,2] <- covariates[indx_cov[,2]]
    indx_cov_txt[,3] <- x_p_covs[indx_cov]
    indx_cov_txt[,4] <- 'covs'
  }
  
  if (nrow(indx_beh_txt) > 0){
    indx_beh_txt[,1] <- regions[indx_beh[,1]]
    indx_beh_txt[,2] <- beh_names[indx_beh[,2]]
    indx_beh_txt[,3] <- x_p_beh[indx_beh]
    indx_beh_txt[,4] <- 'int'
  }
  
  if (nrow(indx_ageXbeh_txt) > 0){
    indx_ageXbeh_txt[,1] <- regions[indx_ageXbeh[,1]]
    indx_ageXbeh_txt[,2] <- beh_names[indx_ageXbeh[,2]]
    indx_ageXbeh_txt[,3] <- x_p_ageXbeh[indx_ageXbeh]
    indx_ageXbeh_txt[,4] <- 'slope'
  }
  
  if (nrow(indx_brain_beh_txt) > 0){
    indx_brain_beh_txt[,1] <- regions[indx_brain_beh[,1]]
    indx_brain_beh_txt[,2] <- beh_names[indx_brain_beh[,2]]
    indx_brain_beh_txt[,3] <- x_p_ageXbeh[indx_brain_beh]
    indx_brain_beh_txt[,4] <- '50mo'
  }
  
  # consolidate fdr stats for ggseg visualizations
  indx_fdr <<- rbind(indx_cov_txt, indx_beh_txt, indx_ageXbeh_txt, indx_brain_beh_txt)
  
  print('Numbers of participants and observations')
  print(stats_sum$n_sub_obs)
  
  print('Numbers of participants for behavioral analyses')
  print(stats_sum$n_beh)
  
  }
  
  print('Minimum number of participants and observations by measure')
  print(apply(stats_sum$n_sub_obs,2,min))
  
  print('Minimum number of participants and observations by measure for behavioral analyses')
  print(apply(stats_sum$n_beh,2,min))
  
}



long_model_report_stats_asym <- function(meas, hem, trk_sz, dir_in, regions,
                                             beh_names, covariates, stats_sum){
  
  # reports stats for dev_struct_dwi analysis
  
  
  # specify brain names as output by FS or py(Baby)AFQ
  if (meas == 'dti_fa' | meas == 'dti_md') {
    if (trk_sz == 'nodes'){
      nodes <- c(1:100)
      brain_names <- paste0('n', nodes, '_', meas)
    } else if (trk_sz == 'quarters'){
      quarters <- c(1:4)
      brain_names <- paste0('q', quarters, '_', meas)
    } else if (trk_sz == 'average'){
      brain_names <- regions <- meas
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
  
  # print stats
  if ( meas == 'dti_fa' | meas == 'dti_md' ){
    sz_in <- 1
    sz_out <- length(brain_names)
  } else {
    sz_in <- 2
    sz_out <- length(regions) + 1
  }
  
  # reduce to regions/tracts of interest
  x_est_cov_asym <- matrix(stats_sum$est_cov_asym[sz_in:sz_out,], nrow = length(regions))
  x_est_cov_r0 <- matrix(stats_sum$est_cov_r0[sz_in:sz_out,], nrow = length(regions))
  x_est_asym <- matrix(stats_sum$est_asym[sz_in:sz_out,], nrow = length(regions))
  x_est_r0 <- matrix(stats_sum$est_r0[sz_in:sz_out,], nrow = length(regions))
  x_p_cov_asym <- matrix(stats_sum$p_cov_asym[sz_in:sz_out,], nrow = length(regions))
  x_p_cov_r0 <- matrix(stats_sum$p_cov_r0[sz_in:sz_out,], nrow = length(regions))
  x_p_asym <- matrix(stats_sum$p_asym[sz_in:sz_out,], nrow = length(regions))
  x_p_r0 <- matrix(stats_sum$p_r0[sz_in:sz_out,], nrow = length(regions))
  
  if (trk_sz != 'nodes'){
    
  print('Uncorrected stats...')
  print(x_est_cov_asym)
  print(x_est_cov_r0)
  print(x_est_asym)
  print(x_est_r0)
  print(x_p_cov_asym)
  print(x_p_cov_r0)
  print(x_p_asym)
  print(x_p_r0)
  
  print('pFDR: covariates')
  p_fdr_cov_asym <- p.adjust(as.numeric(c(x_p_cov_asym)), 'fdr') 
  mat_cov_asym <- matrix(p_fdr_cov_asym, nrow = length(regions), ncol = length(covariates))
  print(mat_cov_asym)
  
  p_fdr_cov_r0 <- p.adjust(as.numeric(c(x_p_cov_r0)), 'fdr') 
  mat_cov_r0 <- matrix(p_fdr_cov_r0, nrow = length(regions), ncol = length(covariates))
  print(mat_cov_r0)
  
  print('pFDR: asym-behavior')
  p_fdr_asym <- p.adjust(as.numeric(c(x_p_asym)), 'fdr') 
  mat_asym <- matrix(p_fdr_asym, nrow = length(regions), ncol = length(beh_names))
  print(mat_asym)
  
  print('pFDR: r0-behavior')
  p_fdr_r0 <- p.adjust(as.numeric(c(x_p_r0)), 'fdr') 
  mat_r0 <- matrix(p_fdr_r0, nrow = length(regions), ncol = length(beh_names))
  print(mat_r0)
  
  print('asym-r0 correlations')
  
  print(stats_sum$re_cor)
  
  
  # Index FDR-corrected regions and behavioral measures by brain measure and hemisphere
  # and their corresponding uncorrected p-values
  indx_cov_asym <- which(mat_cov_asym < 0.05, arr.ind = T)
  indx_cov_r0 <- which(mat_cov_r0 < 0.05, arr.ind = T)
  indx_asym <- which(mat_asym < 0.05, arr.ind = T)
  indx_r0 <- which(mat_r0 < 0.05, arr.ind = T)
  
  indx_cov_asym_txt <- matrix(nrow = nrow(indx_cov_asym), ncol = 4)
  indx_cov_r0_txt <- matrix(nrow = nrow(indx_cov_r0), ncol = 4)
  indx_asym_txt <- matrix(nrow = nrow(indx_asym), ncol = 4)
  indx_r0_txt <- matrix(nrow = nrow(indx_r0), ncol = 4)
  
  if (nrow(indx_cov_asym_txt) > 0){
    indx_cov_asym_txt[,1] <- regions[indx_cov_asym[,1]]
    indx_cov_asym_txt[,2] <- covariates[indx_cov_asym[,2]]
    indx_cov_asym_txt[,3] <- x_p_cov_asym[indx_cov_asym]
    indx_cov_asym_txt[,4] <- 'cov_asym'
  }
  
  if (nrow(indx_cov_r0_txt) > 0){
    indx_cov_r0_txt[,1] <- regions[indx_cov_r0[,1]]
    indx_cov_r0_txt[,2] <- covariates[indx_cov_r0[,2]]
    indx_cov_r0_txt[,3] <- x_p_cov_r0[indx_cov_r0]
    indx_cov_r0_txt[,4] <- 'cov_r0'
  }
  
  if (nrow(indx_asym_txt) > 0){
    indx_asym_txt[,1] <- regions[indx_asym[,1]]
    indx_asym_txt[,2] <- beh_names[indx_asym[,2]]
    indx_asym_txt[,3] <- x_p_asym[indx_asym]
    indx_asym_txt[,4] <- 'asym'
  }
  
  if (nrow(indx_r0_txt) > 0){
    indx_r0_txt[,1] <- regions[indx_r0[,1]]
    indx_r0_txt[,2] <- beh_names[indx_r0[,2]]
    indx_r0_txt[,3] <- x_p_r0[indx_r0]
    indx_r0_txt[,4] <- 'r0'
  }
  
  # consolidate fdr stats for ggseg visualizations
  indx_fdr <<- rbind(indx_cov_asym_txt, indx_cov_r0_txt, indx_asym_txt, indx_r0_txt)
  
  print('Numbers of participants and observations')
  print(stats_sum$n_sub_obs)
  
  print('Numbers of participants for behavioral analyses')
  print(stats_sum$n_beh)
  
  }
  
  print('Minimum number of participants and observations by measure')
  print(apply(stats_sum$n_sub_obs,2,min))
  
  print('Minimum number of participants and observations by measure for behavioral analyses')
  print(apply(stats_sum$n_beh,2,min))
  
}