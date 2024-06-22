long_mediate <- function(med_elem, stats_sum, indx_fdr_bin, meas, hem, trk, trk_sz, dir_in, 
                         other_dat, regions, beh_names, med_names){

# performs mediation on significant (uncorrected) longitudinal associations

  
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


#x_p_int <- matrix(stats_sum$p_beh[sz_in:sz_out,], nrow = length(regions))
#x_p_slope <- matrix(stats_sum$p_ageXbeh[sz_in:sz_out,], nrow = length(regions))
dat4_int_mat1 <- med_elem$dat4_int_mat1
dat4_slope_mat1 <- med_elem$dat4_slope_mat1

if ((meas == 'dti_fa' | meas == 'dti_md') & trk_sz == 'nodes') {
    x_p_int <- sig_nodes_med$sig_nodes_int_mat
    x_p_slope <- sig_nodes_med$sig_nodes_slope_mat
  } else {
    x_p_int <- indx_fdr_bin$indx_fdr_int
    x_p_slope <- indx_fdr_bin$indx_fdr_slope
  }



  
med_behav <- other_dat[, c('id', med_names)]
med_behav <- distinct(med_behav)


for (m in (1:length(med_names))){
  
  est_meds <- matrix(ncol=length(beh_names), nrow=length(regions)) 
  ci1_meds <- matrix(ncol=length(beh_names), nrow=length(regions))
  ci2_meds <- matrix(ncol=length(beh_names), nrow=length(regions))
  p_meds <- matrix(ncol=length(beh_names), nrow=length(regions))

  for (feat in c('int', 'slope')){
  
    dat4_mat1 <- eval(as.symbol(paste0('dat4_', feat, '_mat1')))
    x_p <- eval(as.symbol(paste0('x_p_', feat)))
    dat4_mat2 <- list()
    
    
    # reorganize data for mediation
    for (j in (1:length(beh_names))){
      dat4_mat2_reg <- list()
      
      for (i in (sz_in:sz_out)){ 
        dat4_mat2_reg[[i-(sz_in-1)]] <- dat4_mat1[[i]][[j]]
      }
      dat4_mat2[[j]] <- dat4_mat2_reg
    }
    
    
    for (i in (1:(sz_out - sz_in + 1))){
      for (j in (1:length(beh_names))){
        
        dat4_mat3 <- left_join(dat4_mat2[[j]][[i]], med_behav, by = "id")
        dat4_mat4 <- data.frame(brain3 = dat4_mat3$brain3, beh = dat4_mat3$beh, med = dat4_mat3[, med_names[[m]]])

        if (x_p[i,j] == TRUE){
          dat4_mat4 <- na.omit(dat4_mat4)
          model.M <- lm(beh ~ brain3, dat4_mat4)
          model.Y <- lm(med ~ beh + brain3, dat4_mat4)
          
          results <- mediate(model.M, model.Y, treat='brain3', mediator='beh', boot=TRUE, sims=5000)
          
          est_meds[i,j] <- results$n0
          ci1_meds[i,j] <- results$n0.ci[[1]]
          ci2_meds[i,j] <- results$n0.ci[[2]]
          p_meds[i,j] <- results$n0.p
        } else {
          est_meds[i,j] <- NA
          ci1_meds[i,j] <- NA
          ci2_meds[i,j] <- NA
          p_meds[i,j] <- NA
        }
    }
    }
    
    rownames(est_meds) <- regions
    rownames(ci1_meds) <- regions
    rownames(ci2_meds) <- regions
    rownames(p_meds) <- regions
    
    colnames(est_meds) <- beh_names
    colnames(ci1_meds) <- beh_names
    colnames(ci2_meds) <- beh_names
    colnames(p_meds) <- beh_names
    
    print(paste0('Mediation results for ', med_names[[m]], ' and ', feat))
    print(est_meds)
    print(ci1_meds)
    print(ci2_meds)
    print(p_meds)
  }
}
}