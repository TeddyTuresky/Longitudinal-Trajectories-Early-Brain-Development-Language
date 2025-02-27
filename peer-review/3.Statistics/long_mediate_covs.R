long_mediate <- function(med_elem, stats_sum, indx_fdr_bin, meas, hem, trk, trk_sz, dir_in, 
                         other_dat, covariates_mediation, regions, beh_names, med_names, ctr_str){

# performs mediation on significant (uncorrected) longitudinal associations
  
  #cov_str_med <- paste(covariates_mediation, collapse = ' + ')
  
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
  if (trk_sz == 'nodes'){
    sz_in <- 6
    sz_out <- 95
    regions <- regions[sz_in:sz_out]
  } else if (trk_sz == 'quarters'){
    sz_in <- 1
    sz_out <- 4
  } else if (trk_sz == 'average'){
    sz_in <- 1
    sz_out <- 1
  }
} else {
  sz_in <- 2
  sz_out <- length(regions) + 1
}


dat4_int_mat1 <- med_elem$dat4_int_mat1
dat4_slope_mat1 <- med_elem$dat4_slope_mat1
refs_eTIV <- med_elem$refs_eTIV

  if ((meas == 'dti_fa' | meas == 'dti_md') & trk_sz == 'nodes') {
    x_p_int <- sig_nodes_med$sig_nodes_int_mat
    x_p_slope <- sig_nodes_med$sig_nodes_slope_mat
  } else {
    x_p_int <- indx_fdr_bin$indx_fdr_int
    x_p_slope <- indx_fdr_bin$indx_fdr_slope
  }

# incorporate dependent variables and covariates
med_behav_covs <- other_dat[, c('id', med_names, covariates_mediation)]
med_behav_covs <- distinct(med_behav_covs)


#covs <- other_dat[, c('id', covariates)]
#covs <- distinct(covs)
#med_behav_covs <- left_join(med_behav, covs, by = 'id')

k <- 0
indx_med_all <- list()
for (m in (1:length(med_names))){
  brain_meds <- matrix(ncol=length(beh_names), nrow=length(regions)) 
  est_meds <- matrix(ncol=length(beh_names), nrow=length(regions)) 
  ci1_meds <- matrix(ncol=length(beh_names), nrow=length(regions))
  ci2_meds <- matrix(ncol=length(beh_names), nrow=length(regions))
  p_meds <- matrix(ncol=length(beh_names), nrow=length(regions))
  prop_med <- matrix(ncol=length(beh_names), nrow=length(regions))
  
  for (feat in c('int', 'slope')){
  
    dat4_mat1 <- eval(as.symbol(paste0('dat4_', feat, '_mat1')))
    x_p <- eval(as.symbol(paste0('x_p_', feat)))
    dat4_mat2 <- list()
    
    if (ctr_str == 'eTIV'){
      eTIV_as_cov <- refs_eTIV[, c('id', paste0(feat, '_eTIV'))]
      med_behav_covs_fin <- left_join(med_behav_covs, eTIV_as_cov, by = 'id', multiple = 'first')
      covariates_mediation_fin <- c(covariates_mediation, paste0(feat, '_eTIV'))
    } else {
      med_behav_covs_fin <- med_behav_covs
      covariates_mediation_fin <- covariates_mediation
    }
    cov_str_med <<- paste(covariates_mediation_fin, collapse = ' + ')
    
    # reorganize data for mediation
    for (j in (1:length(beh_names))){
      dat4_mat2_reg <- list()
      
      for (i in (sz_in:sz_out)){ 
        dat4_mat2_reg[[i-(sz_in-1)]] <- dat4_mat1[[i]][[j]]
      }
      dat4_mat2[[j]] <- dat4_mat2_reg
    }
    
    
    for (j in (1:length(beh_names))){
      for (i in (1:(sz_out - sz_in + 1))){

        if (x_p[i,j] == TRUE){
          
          dat4_mat3 <- left_join(dat4_mat2[[j]][[i]], med_behav_covs_fin, by = "id")
          dat4_mat4 <- data.frame(brain3 = dat4_mat3$brain3, beh = dat4_mat3$beh, med = dat4_mat3[, med_names[[m]]])
          dat4_mat4 <- cbind(dat4_mat4, dat4_mat3[, c(covariates_mediation_fin)])
          dat4_mat4 <- na.omit(dat4_mat4)
          dat4_mat4_covs <- dat4_mat4[, c(covariates_mediation_fin)]

          # generate linear regression models
          print(paste0('mediation covariates: ', cov_str_med))
          model.X <- lm(paste("med", "~", "brain3", "+", cov_str_med, sep =' '), dat4_mat4)
          model.M <- lm(paste("beh", "~", "brain3", "+", cov_str_med, sep =' '), dat4_mat4)
          model.Y <- lm(paste("med", "~", "beh", "+", "brain3", "+", cov_str_med, sep =' '), dat4_mat4)
          model.B <- lm(paste("med", "~", "beh", "+", cov_str_med, sep =' '), dat4_mat4)
          
          # also estimate standardized beta coefficients
          model.X.beta <- lm.beta(model.X)
          model.M.beta <- lm.beta(model.M)
          model.Y.beta <- lm.beta(model.Y)
          model.B.beta <- lm.beta(model.B)
          
        
          #print(paste0(beh_names[[j]], ' - ', regions[[i]], ' - ', feat))
          #print(model.M.beta)
          
          #print(paste0(beh_names[[j]], ' - ', regions[[i]], ' - ', feat, ' - ', med_names[[m]]))
          #print(model.Y.beta)
          
          #print(paste0(beh_names[[j]], ' - ', med_names[[m]]))
          #print(model.B.beta)
          
          results <- mediate(model.M, model.Y, treat='brain3', mediator='beh', boot=TRUE, sims=100)
          
          #print(summary(results))
          s_mod <- summary(model.X)
          #print(s_mod$coefficients[2,1])
          #print(s_mod$coefficients[2,4])
          brain_meds[i,j] <- s_mod$coefficients[2,4]
          est_meds[i,j] <- results$d0 # d = ACME, z = ADE, n = proportion mediated
          ci1_meds[i,j] <- results$d0.ci[[1]]
          ci2_meds[i,j] <- results$d0.ci[[2]]
          p_meds[i,j] <- results$d0.p
          prop_med[i,j] <- results$n0
        } else {
          est_meds[i,j] <- NA 
          ci1_meds[i,j] <- NA
          ci2_meds[i,j] <- NA
          p_meds[i,j] <- NA
          prop_med[i,j] <- NA
        }
      }
      
      # consolidate stats
      indx_med <- which(p_meds[,j] < 0.05, arr.ind = T)
      indx_med_txt <- matrix(nrow = length(indx_med), ncol = 4)
      if (length(indx_med) > 0){
        k = k + 1
        indx_med_txt[,1] <- regions[indx_med]
        indx_med_txt[,2] <- beh_names[[j]]
        indx_med_txt[,3] <- med_names[[m]]
        indx_med_txt[,4] <- feat
        #print(indx_med_txt)
        indx_med_all[[k]] <- indx_med_txt
      }
      
    }
    
    rownames(est_meds) <- regions
    rownames(ci1_meds) <- regions
    rownames(ci2_meds) <- regions
    rownames(p_meds) <- regions
    rownames(prop_med) <- regions
    
    colnames(est_meds) <- beh_names
    colnames(ci1_meds) <- beh_names
    colnames(ci2_meds) <- beh_names
    colnames(p_meds) <- beh_names
    colnames(prop_med) <- beh_names
    
    print(paste0('Mediation results for ', med_names[[m]], ' and ', feat))
    print('estimates')
    print(est_meds)
    print('ci lower')
    print(ci1_meds)
    print('ci upper')
    print(ci2_meds)
    print('p-values')
    print(p_meds)
    #print('proportion mediated')
    #print(prop_med)
  }
  
}

indx_med_all <<- indx_med_all

}