long_mod_gen <- function(model, meas, hem, trk, trk_sz, dir_in, other_dat, 
                         covariates, covariates_interaction, cov_var, regions,
                         beh_names, rating, dV, ctr_str, use_inf, use3, 
                         vqc_enigma, rating_vqc){
  
  
  
  # clarify covariate terms
  if (is.null(covariates_interaction) == 1){
    cov_str <- paste(covariates, collapse = ' + ')
    n_covs <- length(covariates)
  } else {
    cov_str <- paste(c(covariates, paste0('age*', covariates_interaction)), collapse = ' + ')
    n_covs <- length(covariates) + length(covariates_interaction)
  }
  
  print(paste0('covariate string: ', cov_str))
  
  # load mri data
  if (meas == 'wm'){
    mri <- read.csv( paste0( dir_in, '/reorg/volume_wm_output.csv' ))
    colnames(mri)[colnames(mri) == 'EstimatedTotalIntraCranialVol'] = 'eTIV'
    colnames(other_dat)[colnames(other_dat) == 'qc_s'] = 'rating'
    colnames(other_dat)[colnames(other_dat) == paste0('qc_parc_', hem)] = 'qc_parc'
  } else if (meas == 'dti_fa' | meas == 'dti_md') {
    if (trk_sz == 'quarters'){
      mri <- read.csv( paste0( dir_in, '/reorg/dwi_trk_quarts_', trk, '.csv' ) )
    } else if (trk_sz == 'nodes'){
      mri <- read.csv( paste0( dir_in, '/reorg/dwi_trk_nodes_', trk, '.csv' ) )
    } else if (trk_sz == 'average'){
      mri <- read.csv( paste0( dir_in, '/reorg/dwi_trk_avg_', trk, '.csv' ) )
    }
    colnames(other_dat)[colnames(other_dat) == paste0('qc_', trk)] = 'rating'
  } else {
    mri <- read.csv( paste0( dir_in, '/reorg/', hem, '_', meas, '_output.csv' ))
    colnames(other_dat)[colnames(other_dat) == 'qc_s'] = 'rating'
    colnames(other_dat)[colnames(other_dat) == paste0('qc_parc_', hem)] = 'qc_parc'
  }
  
  dat <- inner_join(mri, other_dat, by = c('id', 'timepoint', 'timepoint_n'))
  dat <- dat[order(dat$id, dat$age), ]
  dat <- dat[ dat$rating > rating, ]
  
  # save combined dataset to generate descriptive stats and longitudinal plots
  if (meas == 'dti_fa' | meas == 'dti_md'){
    save(dat, file = paste0(dir_in, '/', hem, '_', meas, '_', trk, '_full_dataset_qc.RData'))
  } else {
    save(dat, file = paste0(dir_in, '/', hem, '_', meas, '_full_dataset_qc.RData'))
  }
  
  
  # separate behavior, sex, and qc_parc into their own dataframes
  behav <- dat[, c('id', beh_names)]
  sexes <- dat[, c('id', 'Sex')]
  behav <- distinct(behav)
  sexes <- distinct(sexes)
  if (meas != 'dti_fa' & meas != 'dti_md'){
    qc_parc <- dat[, c('id', 'timepoint', 'timepoint_n', 'qc_parc')]
    qc_parc <- na.omit(qc_parc)
  }
  
  # some parameters for diffusion
  quarters <- c(1:4)
  nodes <- c(1:100)
  dat4_int_mat1 <- list()
  dat4_slope_mat1 <- list()
  dat4_50mo_mat1 <- list()
  dat1.typ_tract <- list()
  dat2.typ.low_tract <- list()
  dat2.typ.high_tract <- list()
  dat2.typ.mid_tract <- list()
  dat1_copy <- list()
  
  
  # Specify regions to examine
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
  
  
  # Graph labels and colors
  labels_colors <- graph_labels_colors(meas)
  y_label <- labels_colors$y_label
  y_ax <- labels_colors$y_ax
  y_scal <- labels_colors$y_scal
  y_lims <- labels_colors$y_lims
  meas_color <- labels_colors$meas_color
  meas_color_low <- labels_colors$meas_color_low
  meas_color_high <- labels_colors$meas_color_high
  
  y_ax_int <- bquote(Intercept ~ .(y_ax))
  y_ax_slope <- bquote(Slope ~ .(y_ax))
  y_ax_pred <- bquote(Predicted ~ .(y_ax))
  
  
  # generate matrices for collecting outputs
  re_cor <- n_obs_removed <- est_int <- est_age <- est_age_sq <- se_int <- se_age <- p_int <- p_age <- p_age_sq <- matrix(ncol=1, nrow=length(brain_names))
  n_sub_obs <- matrix(ncol=2, nrow=length(brain_names))
  est_covs <- se_covs <- p_covs <- matrix(ncol=n_covs, nrow=length(brain_names))
  n_beh <- est_beh <- est_ageXbeh <- est_brain_beh <- p_beh <- p_ageXbeh <- p_brain_beh <- matrix(ncol=length(beh_names), nrow=length(brain_names))
  
  
  # scaled 50 mo age measure
  sage_50 <- mod_sages(model, dat)
  
  
  for (i in (1:length(brain_names))){
    
    if (meas == 'wm'){
      brain_name1 <- gsub(paste0(meas, '.', hem, '.'), '', brain_names[[i]])
    } else if ( meas == 'dti_fa' | meas == 'dti_md'){
      brain_name1 <- gsub(paste0('_', meas), '', brain_names[[i]])
    } else {
      brain_name1 <- gsub(paste0('_', meas), '', gsub(paste0(hem, '_'),'', brain_names[[i]]) )
    }
    
    print(brain_name1)
    
    # reduce table to necessary variables
    if (meas == 'dti_fa' | meas == 'dti_md'){
      dat1 <- dat[, c("id", "timepoint", "timepoint_n", "age", covariates)]
    } else {
      dat1 <- dat[, c("id", "timepoint", "timepoint_n", "age", "eTIV", "euler", covariates)]
    }
    dat1$brain1 <- dat[, brain_names[[i]]]
    
    # incorporate visual eliminations
    if (meas != 'dti_fa' & meas != 'dti_md') {
      qc_parc1 <- qc_parc[ qc_parc$brain_area == brain_name1, ]
      dat1 <- dat1 %>% anti_join( qc_parc1, by=c('id', 'timepoint', 'timepoint_n'))
    }
    
    if (vqc_enigma == 1){
      if (brain_name1 == 'bankssts' | brain_name1 == 'superiortemporal'){
        qc_parcel <- read.csv( paste0(dir_in, '/', 'qc_', hem, '_', brain_name1, '.csv') )
        qc_parcel <- qc_parcel[ qc_parcel$rating > rating_vqc, ]
        dat1 <- inner_join(dat1, qc_parcel, by = c("id", "timepoint", "timepoint_n"))
      }
    }
    
    # record initial number of observations to compare to final number of observations
    n_obs_init <- nrow(dat1[dat1$id %in% dat1$id[duplicated(dat1$id)], ])
    print(n_obs_init)
    
    # average by timepoint
    dat1_tp <- group_by(dat1, timepoint)
    dat1_tp_sum <- summarize(dat1_tp, brain1_tp_mean = mean(brain1), brain1_tp_sd = sd(brain1))
    
    # remove observations with > dV% timepoint-2-timepoint change after 25 mo
    dat1_late_diff_qc <- data.frame( vec = 1 )
    nloops_late <- 0
    while ( nrow(na.omit(dat1_late_diff_qc > 0)) ){
      dat1_late <- dat1[ dat1$age > 25, ]
      dat1_late <- dat1_late[dat1_late$id %in% dat1_late$id[duplicated(dat1_late$id)], ]
      
      dat1_late_diff <- dat1_late %>% 
        group_by(id) %>% 
        mutate(diff = brain1 - lag(brain1, default = first(brain1), order_by = age)) %>%
        mutate(diff_age = age - lag(age, default = first(brain1), order_by = age)) %>%
        mutate(diff_tp = lag(timepoint, order_by = age)) %>%
        mutate(diff_tp_n = lag(timepoint_n, order_by = age))
      
      dat1_late_diff$dec <- abs(dat1_late_diff$diff/(dat1_late_diff$brain1 - dat1_late_diff$diff)) / (dat1_late_diff$diff_age / 12)
      dat1_late_diff <- left_join( dat1_late_diff, dat1_tp_sum, by = 'timepoint')
      dat1_late_diff$brain1_tp_z <- abs( ( dat1_late_diff$brain1 - dat1_late_diff$brain1_tp_mean ) / dat1_late_diff$brain1_tp_sd )
      
      dat1_late_diff <- dat1_late_diff %>% 
        group_by(id) %>% 
        mutate(diff_brain1_tp_z = brain1_tp_z - lag(brain1_tp_z, default = first(brain1_tp_z), order_by = age))
      dat1_late_diff_qc <- dat1_late_diff[ dat1_late_diff$dec > dV, ]
      
      dat1_late_diff_qc_earlier <- dat1_late_diff_qc[ dat1_late_diff_qc$diff_brain1_tp_z < 0, ] # earlier timepoint is farther from mean
      dat1_late_diff_qc_later <- dat1_late_diff_qc[ dat1_late_diff_qc$diff_brain1_tp_z > 0, ] # later timepoint is farther from mean 
      
      # consolidate ids and timepoints
      dat1_late_diff_qc_earlier2 <- dat1_late_diff_qc_earlier[, c('id', 'diff_tp', 'diff_tp_n')]
      colnames(dat1_late_diff_qc_earlier2)[2] <- 'timepoint'
      colnames(dat1_late_diff_qc_earlier2)[3] <- 'timepoint_n'
      dat1_late_diff_qc_later2 <- dat1_late_diff_qc_later[, c('id', 'timepoint', 'timepoint_n')]
      dat1 <- dat1 %>% anti_join( dat1_late_diff_qc_earlier2, by=c('id','timepoint', 'timepoint_n'))
      dat1 <- dat1 %>% anti_join( dat1_late_diff_qc_later2, by=c('id','timepoint', 'timepoint_n'))
      nloops_late <- nloops_late + 1
    }
    print(nloops_late)
    print(nrow(dat1[dat1$id %in% dat1$id[duplicated(dat1$id)], ]))
    
    #--
    
    # remove participants with decrease (increase for md) from INF or TOD to PRE, BEG, or REA, except for thickness
    if ( meas != 'thickness' & meas != 'meancurv') {
      for (tp1 in c('INF', 'TOD')){
        dat1_diff_qc <- data.frame( vec = 1 )
        nloops_early <- 0
        while ( nrow(na.omit(dat1_diff_qc > 0)) ){
          dat1_diff <- dat1 %>% 
            group_by(id) %>% 
            mutate(diff_brain = brain1 - lag(brain1, order_by = age)) %>%
            mutate(diff_tp = lag(timepoint, order_by = age)) 
          
          dat1_diff <- left_join( dat1_diff, dat1_tp_sum, by = 'timepoint')
          dat1_diff$brain1_tp_z <- abs( ( dat1_diff$brain1 - dat1_diff$brain1_tp_mean ) / dat1_diff$brain1_tp_sd )
          
          dat1_diff <- dat1_diff %>% 
            group_by(id) %>% 
            mutate(diff_brain1_tp_z = brain1_tp_z - lag(brain1_tp_z, default = first(brain1_tp_z), order_by = age))
          if (meas == 'dti_md'){
            dat1_diff <- dat1_diff[ dat1_diff$diff_brain > 0, ]
          } else {
            dat1_diff <- dat1_diff[ dat1_diff$diff_brain < 0, ]
          }
          dat1_diff_qc <- dat1_diff[dat1_diff$diff_tp == tp1, ]
          dat1_diff_qc_earlier <- dat1_diff_qc[ dat1_diff_qc$diff_brain1_tp_z < 0, ] # earlier timepoint is farther from mean
          dat1_diff_qc_later <- dat1_diff_qc[ dat1_diff_qc$diff_brain1_tp_z > 0, ] # later timepoint is farther from mean 
          
          # consolidate ids and timepoints
          dat1_diff_qc_earlier2 <- dat1_diff_qc_earlier[, c('id', 'diff_tp')]
          colnames(dat1_diff_qc_earlier2)[2] <- 'timepoint'
          dat1_diff_qc_later2 <- dat1_diff_qc_later[, c('id', 'timepoint')]
          dat1 <- dat1 %>% anti_join( dat1_diff_qc_earlier2, by=c('id','timepoint'))
          dat1 <- dat1 %>% anti_join( dat1_diff_qc_later2, by=c('id','timepoint'))
          nloops_early <- nloops_early + 1
        }
        print(nloops_early)
      }
    }
    print(nrow(dat1[dat1$id %in% dat1$id[duplicated(dat1$id)], ]))
    
    #--
    
    
    dat1 <- na.omit(dat1)
    dat1 <- mutate( dat1, id = as.factor( id ))
    
    # remove participants with 1 observation
    dat1 <- dat1[dat1$id %in% dat1$id[duplicated(dat1$id)], ] 
    
    # ensure at least 1 INF or TOD time point
    if (use_inf == 1){
      dat1_inf_tod <- dat1[dat1$timepoint == 'INF' | dat1$timepoint == 'TOD', ]
      u_inf_tod <- unique(dat1_inf_tod$id)
      du_inf_tod <- data.frame(id = u_inf_tod)
      dat1 <- inner_join(dat1, du_inf_tod, by = 'id')
    }
    
    
    # remove participants with fewer than 3 observations
    if (use3 == 1){
      u <- unique(dat1$id)
      nu <- matrix(ncol=1, nrow=length(u))
      for (k in 1:length(u)){
        nu[k,] <- sum(dat1$id == u[[k]])
      }
      du <- data.frame(id = u, nu)
      du2 <- du[du$nu > 2, ]
      du2 %>% select(-nu)
      dat1 <- inner_join(dat1, du2, by = "id")
    }
    
    # separately, average euler numbers for each individual
    if (meas != 'dti_fa' & meas != 'dti_md'){
      dat1_id <- group_by(dat1, id)
      dat1_id_euler <- summarize(dat1_id, euler = mean(euler)) 
      dat1 <- subset(dat1, select = -c(euler))
    }
    
    #--
    
    
    # Scale age (with centering for linear and quadratic)
    dat1$age_act <- dat1$age
    sd_age <- sd(dat1$age)
    dat1$age <- mod_ages(model, dat1, sd_age, 0)
    dat1$age_sq <- dat1$age^2
   
    
    
    # Run lmer and generate individual growth curves 
    M1 <- long_lmer_mods(dir_in, model, meas, trk, dat1, cov_str, cov_var, brain_names[i])
    
    
    # collect output
    sum_M1 <- summary(M1)
    
    re_cor[i,] <- attr(sum_M1$varcor$id, which = "correlation")[2]
    n_sub_obs[i,1] <- sum_M1$ngrps[[1]]   
    n_sub_obs[i,2] <- nobs(M1)
    n_obs_removed[i,1] <- n_obs_init - n_sub_obs[i,2]
    est_int[i,] <- sum_M1$coefficients[1,1]
    est_age[i,] <- sum_M1$coefficients[2,1]
    se_int[i,] <- sum_M1$coefficients[1,2]
    se_age[i,] <- sum_M1$coefficients[2,2]
    p_int[i,] <- sum_M1$coefficients[1,5]
    p_age[i,] <- sum_M1$coefficients[2,5]
    
    if (model == 'quad'){
      est_age_sq[i,] <- sum_M1$coefficients[3,1]
      est_covs[i,] <- sum_M1$coefficients[4:(3+n_covs),1]
      se_covs[i,] <- sum_M1$coefficients[4:(3+n_covs),2]
      p_age_sq[i,] <- sum_M1$coefficients[3,5]
      p_covs[i,] <- sum_M1$coefficients[4:(3+n_covs),5]
    } else {
      est_covs[i,] <- sum_M1$coefficients[3:(2+n_covs),1]
      se_covs[i,] <- sum_M1$coefficients[3:(2+n_covs),2]
      p_covs[i,] <- sum_M1$coefficients[3:(2+n_covs),5]
      
    }
    
    # Individual growth curves
    dat1$brain_pred <- predict( M1 )
    
    if (meas == 'dti_fa' | meas == 'dti_md'){
      dat_other <- subset( dat1, select = -c(age, age_sq, age_act, timepoint, timepoint_n, brain1, brain_pred))
    } else {
      dat_other <- subset( dat1, select = -c(age, age_sq, age_act, timepoint, timepoint_n, brain1, brain_pred, eTIV))
    }
    
    if (model == 'lin' | model == 'quad'){
      ceeq <- seq(-1, 2.5, 0.1)
    } else if (model == 'log'){
      ceeq <- log(seq(0, 120, 1) + 1)
    }
    synth.dat <- expand.grid( id = unique( dat1$id ), age = ceeq)
    synth.dat$age_act <- mod_ages(model, synth.dat, sd_age, 1)
    synth.dat$age_sq <- (synth.dat$age)^2
    synth.dat1 <- left_join(synth.dat, dat_other, by = 'id', multiple = 'first')
    synth.dat1$brain_pred2 <- predict( M1, newdata = synth.dat1, allow.new.levels = TRUE)
    
    # Generate curve for average child
    dat1.typ <- synth.dat1 %>% 
      group_by(age) %>% 
      summarize_at(vars('brain_pred2'), mean)
    dat1.typ$id <- -1
    dat1.typ$age_act <- mod_ages(model, dat1.typ, sd_age, 1)
    
    
    # resolve title issue for diffusion measures
    if ( meas == 'dti_fa' | meas == 'dti_md' ){
      if (trk_sz == 'average'){
        brain_name_t <- trk
      } else {
        brain_name_t <- paste0(trk, '_', brain_name1)
      }
    } else {
      brain_name_t <- brain_name1
    }
    
    pdf(file= paste0(dir_in, '/', hem, '_', meas, '_', brain_name_t, '_sep.pdf'))
    print(ggplot( data = synth.dat1, aes( x = age_act, y = brain_pred2 ) ) +
            facet_wrap( ~ id ) + 
            geom_line( col = paste0(meas_color), lwd=1, lty=1 ) +
            geom_point( data = dat1, aes( x = age_act, y = brain1 ), size = 2.5, stroke = NA ) +
            scale_x_continuous(name = "Age (Months)", breaks = seq(0, 125, 25)) +
            scale_y_continuous(name = y_ax) + 
            theme(legend.position = "none", axis.text = element_text(size = 6, color = "black"), 
                  axis.title = element_text(size = 10),
                  strip.text = element_text(size = 8, face = "bold", margin = margin(1, 0, 1, 0)),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black")))
    dev.off()
    

    # Consolidate average child data for average (across regions) curves by measure
    if (meas == 'dti_fa' | meas == 'dti_md'){
      modality <- 'dwi'
      tract = trk
      strt <- 1
    } else {
      modality <- 'struct'
      tract <- 'none'
      strt <- 2
    }
    
    # For generating all measures change graph
    if ((meas == 'dti_fa' | meas == 'dti_md') & trk_sz == 'average'){
      dat1.typ.reg <- matrix(nrow = length(dat1.typ$id), ncol = 1)
      dat1.typ.reg[,1] <- dat1.typ$brain_pred2
      dat1.typ.all <- data.frame(id = dat1.typ$id, group = meas, modality = modality, tract = tract, age_act = dat1.typ$age_act, brain_pred_all = rowMeans(dat1.typ.reg))
      save(dat1.typ.all, file = paste0(dir_in, '/dat1.typ.all_', meas, '_', hem, '_', tract, '.RData'))
    } else {
      if (i == strt){
        dat1.typ.reg <- matrix(nrow = length(dat1.typ$id), ncol = length(regions))
        dat1.typ.reg[,1] <- dat1.typ$brain_pred2
      } else if (i == (length(regions)+(strt-1))){
        dat1.typ.reg[,i-1] <- dat1.typ$brain_pred2
        dat1.typ.all <- data.frame(id = dat1.typ$id, group = meas, modality = modality, tract = tract, age_act = dat1.typ$age_act, brain_pred_all = rowMeans(dat1.typ.reg))
        save(dat1.typ.all, file = paste0(dir_in, '/dat1.typ.all_', meas, '_', hem, '_', tract, '.RData'))
      } else if (i != 1 & i < (length(regions)+(strt-1))){
        dat1.typ.reg[,i-1] <- dat1.typ$brain_pred2
      }
    }
    
    pdf(file= paste0(dir_in, '/', hem, '_', meas, '_', brain_name_t, '_full.pdf'), width = 6, height = 6)
    print(ggplot(dat1, aes(x = age_act, y = brain1, group = id)) +
            geom_point(size = 2.5, alpha = 0.1, stroke = NA) +
            geom_path(lwd = 0.5, alpha = 0.1) +
            ggtitle(brain_name_t) +
            scale_x_continuous(name = "Age (Months)", breaks = seq(0, 125, 25)) +
            scale_y_continuous(name = y_ax, breaks = y_scal, limits = y_lims) + 
            theme(legend.position = "none", axis.text = element_text(size = 14, color = "black"), 
                  axis.title = element_text(size = 25),
                  plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                  strip.text = element_text(size = 9, face = "bold"),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  axis.line = element_line(colour = "black")) +
            geom_line(data = dat1.typ, aes(x = age_act, y = brain_pred2), col=paste0(meas_color), lwd=2, lty=1))
    dev.off() 
    
    # Generate graphs if any covariates are significant
    for ( k in (1:length(covariates)) ){
      cov = covariates[k]
      cov.mean <- mean(synth.dat1[, c(cov)])
      
      synth.dat1.cov_low <- synth.dat1[ synth.dat1[, c(cov)] < cov.mean, ]
      synth.dat1.cov_high <- synth.dat1[ synth.dat1[, c(cov)] >= cov.mean, ]
      
      dat1.typ.cov_low <- synth.dat1.cov_low %>% 
        group_by(age) %>% 
        summarize(brain_mean = mean(brain_pred2), brain_sd = sd(brain_pred2), brain_n = n())
      dat1.typ.cov_low$id <- -2
      dat1.typ.cov_low$age_act <- mod_ages(model, dat1.typ.cov_low, sd_age, 1)
      margin_low <- qt(0.975, df = dat1.typ.cov_low$brain_n-1) * dat1.typ.cov_low$brain_sd / sqrt(dat1.typ.cov_low$brain_n)
      dat1.typ.cov_low$lowerinterval <- dat1.typ.cov_low$brain_mean - margin_low
      dat1.typ.cov_low$upperinterval <- dat1.typ.cov_low$brain_mean + margin_low
      
      dat1.typ.cov_high <- synth.dat1.cov_high %>% 
        group_by(age) %>% 
        summarize(brain_mean = mean(brain_pred2), brain_sd = sd(brain_pred2), brain_n = n())
      dat1.typ.cov_high$id <- -3
      dat1.typ.cov_high$age_act <- mod_ages(model, dat1.typ.cov_high, sd_age, 1)
      margin_high <- qt(0.975, df = dat1.typ.cov_high$brain_n-1) * dat1.typ.cov_high$brain_sd / sqrt(dat1.typ.cov_high$brain_n)
      dat1.typ.cov_high$lowerinterval <- dat1.typ.cov_high$brain_mean - margin_high
      dat1.typ.cov_high$upperinterval <- dat1.typ.cov_high$brain_mean + margin_high
      
      if ( brain_name_t != 'eTIV'){
        pdf(file= paste0(dir_in, '/', hem, '_', meas, '_', brain_name_t, '_', cov, '.pdf'), width = 2, height = 2)
        print(ggplot(dat1, aes(x = age_act, y = brain1, group = id)) +
                geom_point(size = 0.8, alpha = 0, stroke = NA) +
                ggtitle(brain_name_t) +
                scale_x_continuous(name = element_blank(), breaks = seq(0, 125, 25)) +  
                scale_y_continuous(name = element_blank()) +
                theme(legend.position = "none", axis.text = element_text(size = 8, color = "black"), 
                      axis.title = element_text(size = 8),
                      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
                      strip.text = element_text(size = 2, face = "bold"),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(colour = "black")) +
                geom_line(data = dat1.typ.cov_low, aes(x = age_act, y = brain_mean), col="firebrick1", lwd=0.75, lty=1) +
                geom_line(data = dat1.typ.cov_high, aes(x = age_act, y = brain_mean), col="gold", lwd=0.75, lty=1) + # lightslateblue
                geom_line(data = dat1.typ.cov_low, aes(x = age_act, y = lowerinterval), col="firebrick1", lwd=0.5, linetype = 'dotted') +
                geom_line(data = dat1.typ.cov_low, aes(x = age_act, y = upperinterval), col="firebrick1", lwd=0.5, linetype = 'dotted') +
                geom_line(data = dat1.typ.cov_high, aes(x = age_act, y = lowerinterval), col="gold", lwd=0.5, linetype = 'dotted') + 
                geom_line(data = dat1.typ.cov_high, aes(x = age_act, y = upperinterval), col="gold", lwd=0.5, linetype = 'dotted'))
        
        dev.off() 
      }
    }
    
    
    # random effects parameters
    refs <- ranef( M1 )$id  # random
    
    refs <- tibble::rownames_to_column(refs, "id")
    colnames(refs)[2] <- "int"
    colnames(refs)[3] <- "slope"
    
    refs_beh <- left_join(refs, behav, by = "id")
    refs_beh <- left_join(refs_beh, sexes, by = "id")
    if (meas != 'dti_fa' & meas != 'dti_md'){
      refs_beh <- left_join(refs_beh, dat1_id_euler, by = 'id')
    }
    
    if (brain_names[[i]] == 'eTIV'){
      refs_eTIV <- refs[, c('id', 'int', 'slope')]
      colnames(refs_eTIV)[colnames(refs_eTIV) == 'int'] = 'int_eTIV'
      colnames(refs_eTIV)[colnames(refs_eTIV) == 'slope'] = 'slope_eTIV'
    } else if (meas != 'dti_fa' & meas != 'dti_md'){
      refs_beh <- left_join(refs_beh, refs_eTIV, by = 'id')
    }
    
    
    # For brain-behavior associations at 50 mo
    dat2 <- dat1[ c('id', covariates)]
    u2 <- unique(dat2$id)
    du3 <- data.frame(id = u2)
    dat3 <- left_join(du3, dat2, by = 'id', multiple = 'first')
    dat3$age = sage_50
    dat3$age_sq <- dat3$age^2
    dat3$brain2 <- predict(M1, dat3)
    dat3_beh <- left_join(dat3, behav, by = 'id')
    if (meas != 'dti_fa' & meas != 'dti_md'){
      dat3_beh <- left_join(dat3_beh, dat1_id_euler, by = 'id')
    }
    
    write.csv(dat3, paste0(dir_in, '/', hem, '_', meas, '_predicted_', brain_name_t, '_50mo.csv'))
    
    if (brain_names[[i]] == 'eTIV'){
      dat3_eTIV <- dat3[, c('id', 'brain2')]
      colnames(dat3_eTIV)[colnames(dat3_eTIV) == 'brain2'] = 'eTIV'
    } else if (meas != 'dti_fa' & meas != 'dti_md'){
      dat3_beh <- left_join(dat3_beh, dat3_eTIV, by = 'id')
    }
    
    dat4_int_mat1_beh <- list()
    dat4_slope_mat1_beh <- list()
    dat4_50mo_mat1_beh <- list()
    dat2.typ.beh_low_tract <- list()
    dat2.typ.beh_high_tract <- list()
    dat2.typ.beh_mid_tract <- list()
    
    for (j in (1:length(beh_names)))
    {
      
      # Associations with intercepts and slopes
      if (ctr_str == 'eTIV' & brain_names[[i]] != 'eTIV'){
        refs_beh1 <- refs_beh[, c("id", "int", "slope", "int_eTIV", "slope_eTIV")]
      }
      else if (ctr_str == 'euler'){
        refs_beh1 <- refs_beh[, c("id", "int", "slope", "euler")]
      }
      else{
        refs_beh1 <- refs_beh[, c("id", "int", "slope")]
      }
      
      refs_beh1$beh <- refs_beh[, beh_names[[j]]]
      refs_beh1 <- na.omit(refs_beh1)
      
      if (ctr_str == 'eTIV' & brain_names[[i]] != 'eTIV'){
        reg_int <- lm(int ~ int_eTIV, refs_beh1) 
        reg_slope <- lm(slope ~ slope_eTIV, refs_beh1)
        dat4_int <- data.frame(brain_res = reg_int$residuals, beh = refs_beh1$beh)
        dat4_slope <- data.frame(brain_res = reg_slope$residuals, beh = refs_beh1$beh)
      }
      else if (ctr_str == 'euler'){
        reg_int <- lm(int ~ euler, refs_beh1)
        reg_slope <- lm(slope ~ euler, refs_beh1)
        dat4_int <- data.frame(brain_res = reg_int$residuals, beh = refs_beh1$beh)
        dat4_slope <- data.frame(brain_res = reg_slope$residuals, beh = refs_beh1$beh)
      }
      else{
        dat4_int <- data.frame(brain_res = refs_beh1$int, beh = refs_beh1$beh)
        dat4_slope <- data.frame(brain_res = refs_beh1$slope, beh = refs_beh1$beh)
        
      }
      
      reg_int_cor <- cor.test(dat4_int$beh, dat4_int$brain_res)
      reg_slope_cor <- cor.test(dat4_slope$beh, dat4_slope$brain_res)
      
      n_beh[i,j] <- nrow(refs_beh1)
      est_beh[i,j] <- reg_int_cor$estimate
      est_ageXbeh[i,j] <- reg_slope_cor$estimate
      p_beh[i,j] <- reg_int_cor$p.value
      p_ageXbeh[i,j] <- reg_slope_cor$p.value
      
      
      # Associations with predicted brain structure at 50 mo
      if (ctr_str != 'none' & brain_names[[i]] != 'eTIV'){
        dat3_beh1 <- dat3_beh[, c("id", "brain2", "Sex", ctr_str)]
      }
      else{
        dat3_beh1 <- dat3_beh[, c("id", "brain2", "Sex")] 
      }
      
      dat3_beh1$beh <- dat3_beh[, beh_names[[j]]]
      dat3_beh1 <- na.omit(dat3_beh1)
      
      dat3_sex_reg <- lm(brain2 ~ Sex, dat3_beh1)
      
      dat4 <- data.frame(brain_res = dat3_sex_reg$residuals, beh = dat3_beh1$beh)
      dat4_cor <- cor.test(dat4$beh, dat4$brain_res)
      
      est_brain_beh[i,j] <- dat4_cor$estimate 
      p_brain_beh[i,j] <- dat4_cor$p.value
      
      
      # Graph significant associations
      if (p_beh[i,j] < 0.05 & is.na(p_beh[i,j]) != 1){
        pdf(file= paste0(dir_in, '/', hem, '_', meas, '_', brain_name_t, '_', beh_names[[j]], '_int.pdf'), width = 6, height = 6)
        print(ggplot(dat4_int, aes(x = beh, y = brain_res)) +
                geom_point(size = 2.5, stroke = NA) +
                geom_smooth(formula = 'y ~ x', method = 'lm', color = paste0(meas_color)) +
                ggtitle(brain_name_t) +
                theme(legend.position = "none") +
                scale_x_continuous(name = paste(beh_names[[j]], sep=' ')) +
                scale_y_continuous(name = y_ax_int) +
                theme(legend.position = "none", axis.text = element_text(size = 11, color = "black"), 
                      axis.title = element_text(size = 20),
                      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                      strip.text = element_text(size = 9, face = "bold"),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(colour = "black")))
        dev.off()
      }
      
      if (p_ageXbeh[i,j] < 0.05 & is.na(p_ageXbeh[i,j]) != 1){
        pdf(file= paste0(dir_in, '/', hem, '_', meas, '_', brain_name_t, '_', beh_names[[j]], '_slope.pdf'), width = 6, height = 6)
        print(ggplot(dat4_slope, aes(x = beh, y = brain_res)) +
                geom_point(size = 2.5, stroke = NA) +
                geom_smooth(formula = 'y ~ x', method = 'lm', color = paste0(meas_color)) +
                ggtitle(brain_name_t) +
                theme(legend.position = "none") +
                scale_x_continuous(name = paste(beh_names[[j]], sep=' ')) +
                scale_y_continuous(name = y_ax_slope) +
                theme(legend.position = "none", axis.text = element_text(size = 11, color = "black"), 
                      axis.title = element_text(size = 20),
                      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                      strip.text = element_text(size = 9, face = "bold"),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(colour = "black")))
        dev.off()
      }
      
      behav1 <- behav[, c('id', beh_names[[j]])]
      synth.dat2 <- left_join(synth.dat1, behav1, by = 'id')
      synth.dat2 <- na.omit(synth.dat2)
      
      synth.dat2.beh_low <- synth.dat2[ synth.dat2[, beh_names[[j]]] <= 85, ]
      synth.dat2.beh_high <- synth.dat2[ synth.dat2[, beh_names[[j]]] >= 115, ] 
      synth.dat2.beh_mid <- synth.dat2[ synth.dat2[, beh_names[[j]]] > 85 & synth.dat2[, beh_names[[j]]] < 115, ] 
      
      dat2.typ.beh_low <- synth.dat2.beh_low %>% 
        group_by(age) %>% 
        summarize_at(vars('brain_pred2'), mean)
      dat2.typ.beh_low$id <- -2
      dat2.typ.beh_low$age_act <- mod_ages(model, dat2.typ.beh_low, sd_age, 1)
      
      dat2.typ.beh_high <- synth.dat2.beh_high %>% 
        group_by(age) %>% 
        summarize_at(vars('brain_pred2'), mean)
      dat2.typ.beh_high$id <- -3
      dat2.typ.beh_high$age_act <- mod_ages(model, dat2.typ.beh_high, sd_age, 1)
      
      dat2.typ.beh_mid <- synth.dat2.beh_mid %>% 
        group_by(age) %>% 
        summarize_at(vars('brain_pred2'), mean)
      dat2.typ.beh_mid$id <- -4
      dat2.typ.beh_mid$age_act <- mod_ages(model, dat2.typ.beh_mid, sd_age, 1)
      
      if ((p_beh[i,j] < 0.05 & is.na(p_beh[i,j]) != 1) | (p_ageXbeh[i,j] < 0.05 & is.na(p_beh[i,j]) != 1)){
        pdf(file= paste0(dir_in, '/', hem, '_', meas, '_', brain_name_t, '_', beh_names[[j]], '_split_traj.pdf'), width = 6, height = 6)
        print(ggplot(dat1, aes(x = age_act, y = brain1, group = id)) +
                geom_point(size = 2.5, alpha = 0.1, stroke = NA) +
                geom_path(lwd = 0.5, alpha = 0.1) +
                ggtitle(brain_name_t) +
                scale_x_continuous(name = "Age (Months)", breaks = seq(0, 125, 25), limits = c(0, 125)) +
                scale_y_continuous(name = y_ax) +
                theme(legend.position = "none", axis.text = element_text(size = 14, color = "black"), 
                      axis.title = element_text(size = 25),
                      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                      strip.text = element_text(size = 9, face = "bold"),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(colour = "black")) +
                geom_line(data = dat2.typ.beh_low, aes(x = age_act, y = brain_pred2), color = paste0(meas_color_low), lwd=2) +
                geom_line(data = dat2.typ.beh_high, aes(x = age_act, y = brain_pred2), color = paste0(meas_color_high), lwd=2) +
                geom_line(data = dat2.typ.beh_mid, aes(x = age_act, y = brain_pred2), color = paste0(meas_color), lwd=2, lty=1))
        dev.off() 
      }
      
      
      if (p_brain_beh[i,j] < 0.05 & is.na(p_brain_beh[i,j]) != 1){
        pdf(file= paste0(dir_in, '/', hem, '_', meas, '_', brain_name_t, '_', beh_names[[j]], '_50mo.pdf'), width = 6, height = 6)
        print(ggplot(dat4, aes(x = beh, y = brain_res)) +
                geom_point(size = 2.5, stroke = NA) +
                geom_smooth(formula = 'y ~ x', method = 'lm', color = paste0(meas_color)) +
                ggtitle(brain_name_t) +
                theme(legend.position = "none") +
                scale_x_continuous(name = paste(beh_names[[j]], sep=' ')) +
                scale_y_continuous(name = y_ax_pred) + 
                theme(legend.position = "none", axis.text = element_text(size = 11, color = "black"), 
                      axis.title = element_text(size = 20),
                      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                      strip.text = element_text(size = 9, face = "bold"),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(colour = "black")))
        dev.off()
      }
      
      # Consolidate data by behavior for node-based diffusion cluster analysis
      dat4_int_mat1_beh[[j]] <- data.frame( id = refs_beh1$id, brain3 = refs_beh1$int, beh3 = refs_beh1$beh)
      dat4_slope_mat1_beh[[j]] <- data.frame( id = refs_beh1$id, brain3 = refs_beh1$slope, beh3 = refs_beh1$beh)
      dat4_50mo_mat1_beh[[j]] <- data.frame( id = dat3_beh1$id, brain3 = dat3_sex_reg$residuals, beh3 = dat3_beh1$beh) # used brain residuals, adjusted for sex
      
      # Save average child data to generate growth curve across all nodes or quarters
      dat2.typ.beh_low_tract[[j]] <- dat2.typ.beh_low
      dat2.typ.beh_high_tract[[j]] <- dat2.typ.beh_high
      dat2.typ.beh_mid_tract[[j]] <- dat2.typ.beh_mid
      
    }
    
    # consolidate data by node for diffusion cluster analysis and whole-tract plots
    dat4_int_mat1[[i]] <- dat4_int_mat1_beh
    dat4_slope_mat1[[i]] <- dat4_slope_mat1_beh
    dat4_50mo_mat1[[i]] <- dat4_50mo_mat1_beh
    dat1.typ_tract[[i]] <- dat1.typ
    dat2.typ.low_tract[[i]] <- dat2.typ.beh_low_tract
    dat2.typ.high_tract[[i]] <- dat2.typ.beh_high_tract
    dat2.typ.mid_tract[[i]] <- dat2.typ.beh_mid_tract
    
    # Save FA and MD estimates for cluster analysis
    if ((meas == 'dti_fa' | meas == 'dti_md') & trk_sz == 'nodes'){
      dat1_copy[[i]] <- dat1
    }
    
  }
  
  # export for node-based diffusion cluster analysis
  diff_clus_elem <<- list('dat4_int_mat1' = dat4_int_mat1, 
                          'dat4_slope_mat1' = dat4_slope_mat1,
                          'dat4_50mo_mat1' = dat4_50mo_mat1,
                          'dat1_copy' = dat1_copy,
                          'dat1.typ_tract' = dat1.typ_tract,
                          'dat2.typ.low_tract' = dat2.typ.low_tract,
                          'dat2.typ.high_tract' = dat2.typ.high_tract,
                          'dat2.typ.mid_tract' = dat2.typ.mid_tract)
  
  # export for mediation analysis
  if (meas == 'dti_fa' | meas == 'dti_md'){
    med_elem <<- list('dat4_int_mat1' = dat4_int_mat1, 
                      'dat4_slope_mat1' = dat4_slope_mat1)
  } else {
    med_elem <<- list('dat4_int_mat1' = dat4_int_mat1, 
                    'dat4_slope_mat1' = dat4_slope_mat1,
                    'refs_eTIV' = refs_eTIV)
  }
  
  # export stats
  stats_sum <<- list('est_int' = est_int, 'est_age' = est_age, 'est_age_sq' = est_age_sq,
                     'est_covs' = est_covs, 'est_beh' = est_beh, 
                     'est_ageXbeh' = est_ageXbeh, 'est_brain_beh' = est_brain_beh, 
                     'se_int' = se_int, 'se_age' = se_age, 'se_covs' = se_covs,
                     'p_int' = p_int, 'p_age' = p_age, 'p_age_sq' = p_age_sq,
                     'p_covs' = p_covs, 'p_beh' = p_beh, 'p_ageXbeh' = p_ageXbeh,
                     'p_brain_beh' = p_brain_beh,
                     're_cor' = re_cor, 'n_beh' = n_beh, 'n_sub_obs' = n_sub_obs, 
                     'n_obs_removed' = n_obs_removed)
  
  

}