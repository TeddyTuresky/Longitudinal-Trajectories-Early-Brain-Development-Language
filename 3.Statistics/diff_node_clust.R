diff_node_clust <- function(diff_clus_elem, stats_sum, meas, hem, trk, dir_in, beh_names){
  
  
  
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
  y_ax_50mo <- bquote(Predicted ~ .(y_ax))
  
  
  x_est_beh <- stats_sum$est_beh[6:95,]
  x_est_ageXbeh <- stats_sum$est_ageXbeh[6:95,]
  x_est_brain_beh <- stats_sum$est_brain_beh[6:95,]
  x_p_beh <- stats_sum$p_beh[6:95,]
  x_p_ageXbeh <- stats_sum$p_ageXbeh[6:95,]
  x_p_brain_beh <- stats_sum$p_brain_beh[6:95,]
  dat4_int_mat1 <- diff_clus_elem$dat4_int_mat1
  dat4_slope_mat1 <- diff_clus_elem$dat4_slope_mat1
  dat4_50mo_mat1 <- diff_clus_elem$dat4_50mo_mat1
  dat1.typ <- diff_clus_elem$dat1.typ_tract[[50]] # for graph area, any node is probably fine, assuming there's data
  dat2.typ.low_tract <- diff_clus_elem$dat2.typ.low_tract
  dat2.typ.high_tract <- diff_clus_elem$dat2.typ.high_tract
  dat2.typ.mid_tract <- diff_clus_elem$dat2.typ.mid_tract
  dat1_copy <- diff_clus_elem$dat1_copy
  dat2.typ.low <- list()
  dat2.typ.high <- list()
  dat2.typ.mid <- list()
  sig_nodes_int_mat <- matrix(0, nrow = 90, ncol = length(beh_names))
  sig_nodes_slope_mat <- matrix(0, nrow = 90, ncol = length(beh_names))
  
  # reorganize low, high, mid
  for (j in (1:length(beh_names))){
    dat2.typ.low_node <- list()
    dat2.typ.high_node <- list()
    dat2.typ.mid_node <- list()
    for (i in (6:95)){ # remove first five and last five nodes
      dat2.typ.low_node[[i-5]] <- dat2.typ.low_tract[[i]][[j]]
      dat2.typ.high_node[[i-5]] <- dat2.typ.high_tract[[i]][[j]]
      dat2.typ.mid_node[[i-5]] <- dat2.typ.mid_tract[[i]][[j]]
    }
    dat2.typ.low[[j]] <- dat2.typ.low_node
    dat2.typ.high[[j]] <- dat2.typ.high_node
    dat2.typ.mid[[j]] <- dat2.typ.mid_node
  }
  
  
  for (feat in c('int', 'slope', '50mo')){
    dat4_mat2 <- list()
    dat4_mat1 <- eval(as.symbol(paste0('dat4_', feat, '_mat1')))
    y_ax <- eval(as.symbol(paste0('y_ax_', feat)))
    
    if (length(dat4_mat1) < 100){
      print('some node values missing')
      dat4_mat1[c((length(dat4_mat1) + 1):100)] <- list(NULL)
    }
    
    # reorganize data for cluster analysis
    for (j in (1:length(beh_names))){
      dat4_mat2_node <- list()
      for (i in (6:95)){ # remove first five and last five nodes
        dat4_mat2_node[[i-5]] <- dat4_mat1[[i]][[j]]
      }
      dat4_mat2[[j]] <- dat4_mat2_node
    }
    
    for (j in (1:length(beh_names))){
      dat4_mat2_filt <- Filter(Negate(is.null), dat4_mat2[[j]])
      for (i in (1:length(dat4_mat2_filt))){
        if (i == 1){
          big_mat <- dat4_mat2_filt[[i]]
        } else {
          big_mat <- inner_join(big_mat, dat4_mat2_filt[[i]], by = 'id')
        }
      }
      dat_id_beh <- data.frame(id = big_mat$id, beh = big_mat$beh3.x)
      big_mat_brain <<- big_mat[c(FALSE, TRUE)]
      
      # run permuco
      clus_cm <<- clusterlm(big_mat_brain ~ beh, dat_id_beh, multcomp = 'clustermass')
      clus_tfce <- clusterlm(big_mat_brain ~ beh, dat_id_beh, P = clus_cm$P, multcomp = 'tfce')
      
      try(print(summary(clus_cm, multcomp = "clustermass")$beh))
      try(print(summary(clus_tfce, multcomp = "tfce")$beh))
      
      pdf(file = paste0(dir_in, '/', hem, '_', meas, '_', trk, '_', beh_names[j], '_clus_tfce_', feat, '.pdf'))
      try(print(plot(clus_tfce, multcomp = "tfce", enhanced_stat = TRUE)))
      dev.off()
      
      # identify significant nodes from tfce
      s_clus_tfce <- summary(clus_tfce)
      
      s <- 1
      sig_nodes <- 0
      for (i in (1:length(s_clus_tfce$beh$start))){
        if (s_clus_tfce$beh$`P(>)`[i] == 'sign'){
          if (s == 1){
            s <- s + 1
            sig_nodes <- c(s_clus_tfce$beh$start[i]:s_clus_tfce$beh$end[i])
          } else {
            sig_nodes <- c(sig_nodes, c(s_clus_tfce$beh$start[i]:s_clus_tfce$beh$end[i]))
          }
        }
      }
      
      if (s > 1){
        
        # Write out significant nodes, factor in 5 node offset
        write.table(sig_nodes+5, paste0(dir_in, '/', hem, '_', meas, '_', trk, '_', feat, '_sig_nodes.txt'), 
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        
        # average significant nodes as designated by permuco and generate plots
        avg_big <- data.frame(id = dat_id_beh$id, beh = dat_id_beh$beh, brain = rowMeans(big_mat_brain[, c(sig_nodes)]))
        
        avg_stats <- matrix(nrow = 2, ncol = 1)
        if (feat == 'int'){
          avg_stats[1,1] <- mean(x_est_beh[sig_nodes])
          avg_stats[2,j] <- mean(x_p_beh[sig_nodes])
        } else if (feat == 'slope'){
          avg_stats[1,1] <- mean(x_est_ageXbeh[sig_nodes])
          avg_stats[2,1] <- mean(x_p_ageXbeh[sig_nodes])
        } else if (feat == '50mo'){
          avg_stats[1,1] <- mean(x_est_brain_beh[sig_nodes])
          avg_stats[2,1] <- mean(x_p_brain_beh[sig_nodes])      
        }
        
        print(paste0('Stats for nodes: ', paste(c(sig_nodes+5), collapse = ' ')))
        print(avg_stats)
        
        for (k in sig_nodes){
          if (k == sig_nodes[1]){
            big_dat1 <- data.frame(id = dat1_copy[[k]]$id, age = dat1_copy[[k]]$age,
                                   age_act = dat1_copy[[k]]$age_act, brain1 = dat1_copy[[k]]$brain1)
            big_dat2.typ.low <- data.frame(id = dat2.typ.low[[j]][[k]]$id, age = dat2.typ.low[[j]][[k]]$age, 
                                           age_act = dat2.typ.low[[j]][[k]]$age_act, brain_pred2 = dat2.typ.low[[j]][[k]]$brain_pred2)
            big_dat2.typ.high <- data.frame(id = dat2.typ.high[[j]][[k]]$id, age = dat2.typ.high[[j]][[k]]$age, 
                                            age_act = dat2.typ.high[[j]][[k]]$age_act, brain_pred2 = dat2.typ.high[[j]][[k]]$brain_pred2)
            big_dat2.typ.mid <- data.frame(id = dat2.typ.mid[[j]][[k]]$id, age = dat2.typ.mid[[j]][[k]]$age, 
                                           age_act = dat2.typ.mid[[j]][[k]]$age_act, brain_pred2 = dat2.typ.mid[[j]][[k]]$brain_pred2)
          } else {
            df <- data.frame(id = dat1_copy[[k]]$id, age = dat1_copy[[k]]$age,
                             age_act = dat1_copy[[k]]$age_act, brain1 = dat1_copy[[k]]$brain1)
            big_dat1 <- left_join(big_dat1, df, by = c('id', 'age', 'age_act'))
            big_dat2.typ.low <- cbind(big_dat2.typ.low, brain_pred2 = dat2.typ.low[[j]][[k]]$brain_pred2)
            big_dat2.typ.high <- cbind(big_dat2.typ.high, brain_pred2 = dat2.typ.high[[j]][[k]]$brain_pred2)
            big_dat2.typ.mid <- cbind(big_dat2.typ.mid, brain_pred2 = dat2.typ.mid[[j]][[k]]$brain_pred2)
          }
        }
        avg_big_dat1 <- data.frame(id = big_dat1$id, age = big_dat1$age, 
                                   age_act = big_dat1$age_act, brain1 = rowMeans(big_dat1[, c(4:ncol(big_dat1))]))
        avg_big_dat2.typ.low <- data.frame(id = big_dat2.typ.low$id, age = big_dat2.typ.low$age, 
                                           age_act = big_dat2.typ.low$age_act, brain = rowMeans(big_dat2.typ.low[, c(4:ncol(big_dat2.typ.low))]))
        avg_big_dat2.typ.high <- data.frame(id = big_dat2.typ.high$id, age = big_dat2.typ.high$age, 
                                            age_act = big_dat2.typ.high$age_act, brain = rowMeans(big_dat2.typ.high[, c(4:ncol(big_dat2.typ.high))]))
        avg_big_dat2.typ.mid <- data.frame(id = big_dat2.typ.mid$id, age = big_dat2.typ.mid$age, 
                                           age_act = big_dat2.typ.mid$age_act, brain = rowMeans(big_dat2.typ.mid[, c(4:ncol(big_dat2.typ.mid))]))
        
        
        
        pdf(file= paste0(dir_in, '/', hem, '_', meas, '_', trk, '_', beh_names[[j]], '_', feat, '_clus_sig_nodes.pdf'), width = 6, height = 6)
        print(ggplot(avg_big, aes(x = beh, y = brain)) +
                geom_point(size = 2.5, stroke = NA) +
                geom_smooth(formula = 'y ~ x', method = 'lm', color = paste0(meas_color)) +
                ggtitle(trk) +
                theme(legend.position = "none") +
                scale_x_continuous(name = paste(beh_names[[j]], sep=' ')) +
                scale_y_continuous(name = y_ax) + 
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
        
        pdf(file= paste0(dir_in, '/', hem, '_', meas, '_', trk, '_', beh_names[[j]], '_', feat, '_clus_sig_nodes_split_traj.pdf'), width = 6, height = 6)
        print(ggplot(avg_big_dat1, aes(x = age_act, y = brain1, group = id)) +
                geom_point(size = 2.5, alpha = 0.1, stroke = NA) +
                geom_path(lwd = 0.5, alpha = 0.1) +
                ggtitle(trk) +
                scale_x_continuous(name = "Age (Months)", breaks = seq(0, 125, 25), limits = c(0, 125)) +
                scale_y_continuous(name = y_label) +
                theme(legend.position = "none", axis.text = element_text(size = 14, color = "black"), 
                      axis.title = element_text(size = 25),
                      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                      strip.text = element_text(size = 9, face = "bold"),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      axis.line = element_line(colour = "black")) +
                geom_line(data = avg_big_dat2.typ.low, aes(x = age_act, y = brain), color = paste0(meas_color_low), lwd=2) +
                geom_line(data = avg_big_dat2.typ.high, aes(x = age_act, y = brain), color = paste0(meas_color_high), lwd=2) +
                geom_line(data = avg_big_dat2.typ.mid, aes(x = age_act, y = brain), color = paste0(meas_color), lwd=2, lty=1))
        dev.off() 
        
      }
      
      # consolidate significant nodes for mediation analysis
      if (feat == 'int'){
        sig_nodes_int_mat[sig_nodes, j] <- 1
      } else if (feat == 'slope'){
        sig_nodes_slope_mat[sig_nodes, j] <- 1
      }
      
    }
  }
  
  # export for mediation analysis
  sig_nodes_med <<- list('sig_nodes_int_mat' = sig_nodes_int_mat,
                         'sig_nodes_slope_mat' = sig_nodes_slope_mat)
  
  
}



diff_node_clust_asym <- function(diff_clus_elem, meas, hem, trk, dir_in, beh_names){


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


dat4_asym_mat2 <- list()
dat4_r0_mat2 <- list()
dat4_asym_mat1 <- diff_clus_elem$dat4_asym_mat1
dat4_r0_mat1 <- diff_clus_elem$dat4_r0_mat1
dat1.typ <- diff_clus_elem$dat1.typ_tract[[50]] # for graph area, any node is probably fine, assuming there's data
dat2.typ.low_tract <- diff_clus_elem$dat2.typ.low_tract
dat2.typ.high_tract <- diff_clus_elem$dat2.typ.high_tract
dat2.typ.mid_tract <- diff_clus_elem$dat2.typ.mid_tract
dat1_copy <- diff_clus_elem$dat1_copy
dat2.typ.low <- list()
dat2.typ.high <- list()
dat2.typ.mid <- list()
  
# reorganize low, high, mid
for (j in (1:length(beh_names))){
  dat2.typ.low_node <- list()
  dat2.typ.high_node <- list()
  dat2.typ.mid_node <- list()
  for (i in (6:95)){ # remove first five and last five nodes
    dat2.typ.low_node[[i-5]] <- dat2.typ.low_tract[[i]][[j]]
    dat2.typ.high_node[[i-5]] <- dat2.typ.high_tract[[i]][[j]]
    dat2.typ.mid_node[[i-5]] <- dat2.typ.mid_tract[[i]][[j]]
  }
  dat2.typ.low[[j]] <- dat2.typ.low_node
  dat2.typ.high[[j]] <- dat2.typ.high_node
  dat2.typ.mid[[j]] <- dat2.typ.mid_node
}


for (feat in c('asym', 'r0')){
  dat4_mat2 <- list()
  dat4_mat1 <- eval(as.symbol(paste0('dat4_', feat, '_mat1')))
  
  if (length(dat4_mat1) < 100){
    print('some node values missing')
    dat4_mat1[c((length(dat4_mat1) + 1):100)] <- list(NULL)
  }
  
  # reorganize data for cluster analysis
  for (j in (1:length(beh_names))){
    dat4_mat2_node <- list()
    for (i in (6:95)){ # remove first five and last five nodes
      dat4_mat2_node[[i-5]] <- dat4_mat1[[i]][[j]]
    }
    dat4_mat2[[j]] <- dat4_mat2_node
  }
  
  for (j in (1:length(beh_names))){
    dat4_mat2_filt <- Filter(Negate(is.null), dat4_mat2[[j]])
    for (i in (1:length(dat4_mat2_filt))){
      if (i == 1){
        big_mat <- dat4_mat2_filt[[i]]
      } else {
        big_mat <- inner_join(big_mat, dat4_mat2_filt[[i]], by = 'id')
      }
    }
    dat_id_beh <- data.frame(id = big_mat$id, beh = big_mat$beh3.x)
    big_mat_brain <<- big_mat[c(FALSE, TRUE)]
    
    # run permuco
    clus_cm <<- clusterlm(big_mat_brain ~ beh, dat_id_beh, multcomp = 'clustermass')
    clus_tfce <- clusterlm(big_mat_brain ~ beh, dat_id_beh, P = clus_cm$P, multcomp = 'tfce')
    
    try(print(summary(clus_cm, multcomp = "clustermass")$beh))
    try(print(summary(clus_tfce, multcomp = "tfce")$beh))
    
    pdf(file = paste0(dir_in, '/', hem, '_', meas, '_', trk, '_', beh_names[j], '_clus_tfce_', feat, '.pdf'))
    try(print(plot(clus_tfce, multcomp = "tfce", enhanced_stat = TRUE)))
    dev.off()
    
    # identify significant nodes from tfce
    s_clus_tfce <- summary(clus_tfce)
    
    s <- 1
    for (i in (1:length(s_clus_tfce$beh$start))){
      if (s_clus_tfce$beh$`P(>)`[i] == 'sign'){
        if (s == 1){
          s <- s + 1
          sig_nodes <- c(s_clus_tfce$beh$start[i]:s_clus_tfce$beh$end[i])
        } else {
          sig_nodes <- c(sig_nodes, c(s_clus_tfce$beh$start[i]:s_clus_tfce$beh$end[i]))
        }
      }
    }
    
    if (s > 1){
      
      # Write out significant nodes, factor in 5 node offset
      write.table(sig_nodes+5, paste0(dir_in, '/', hem, '_', meas, '_', trk, '_', feat, '_sig_nodes.txt'), 
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      
      # average significant nodes as designated by permuco and generate plots
      avg_big <- data.frame(id = dat_id_beh$id, beh = dat_id_beh$beh, brain = rowMeans(big_mat_brain[, c(sig_nodes)]))
      
      avg_stats <- matrix(nrow = 2, ncol = 1)
      if (feat == 'int'){
        avg_stats[1,1] <- mean(x_est_beh[sig_nodes])
        avg_stats[2,j] <- mean(x_p_beh[sig_nodes])
      } else if (feat == 'slope'){
        avg_stats[1,1] <- mean(x_est_ageXbeh[sig_nodes])
        avg_stats[2,1] <- mean(x_p_ageXbeh[sig_nodes])
      } else if (feat == '50mo'){
        avg_stats[1,1] <- mean(x_est_brain_beh[sig_nodes])
        avg_stats[2,1] <- mean(x_p_brain_beh[sig_nodes])      
      }
      
      print(paste0('Stats for nodes: ', paste(c(sig_nodes+5), collapse = ' ')))
      print(avg_stats)
      
      for (k in sig_nodes){
        if (k == sig_nodes[1]){
          big_dat1 <- data.frame(id = dat1_copy[[k]]$id, age = dat1_copy[[k]]$age,
                                brain1 = dat1_copy[[k]]$brain1)
          big_dat2.typ.low <- data.frame(id = dat2.typ.low[[j]][[k]]$id, age = dat2.typ.low[[j]][[k]]$age, 
                                brain_pred2 = dat2.typ.low[[j]][[k]]$brain_pred2)
          big_dat2.typ.high <- data.frame(id = dat2.typ.high[[j]][[k]]$id, age = dat2.typ.high[[j]][[k]]$age, 
                                brain_pred2 = dat2.typ.high[[j]][[k]]$brain_pred2)
          big_dat2.typ.mid <- data.frame(id = dat2.typ.mid[[j]][[k]]$id, age = dat2.typ.mid[[j]][[k]]$age, 
                                brain_pred2 = dat2.typ.mid[[j]][[k]]$brain_pred2)
        } else {
          df <- data.frame(id = dat1_copy[[k]]$id, age = dat1_copy[[k]]$age,
                           brain1 = dat1_copy[[k]]$brain1)
          big_dat1 <- left_join(big_dat1, df, by = c('id', 'age'))
          big_dat2.typ.low <- cbind(big_dat2.typ.low, brain_pred2 = dat2.typ.low[[j]][[k]]$brain_pred2)
          big_dat2.typ.high <- cbind(big_dat2.typ.high, brain_pred2 = dat2.typ.high[[j]][[k]]$brain_pred2)
          big_dat2.typ.mid <- cbind(big_dat2.typ.mid, brain_pred2 = dat2.typ.mid[[j]][[k]]$brain_pred2)
        }
      }
      avg_big_dat1 <- data.frame(id = big_dat1$id, age = big_dat1$age, 
                                 brain1 = rowMeans(big_dat1[, c(4:ncol(big_dat1))]))
      avg_big_dat2.typ.low <- data.frame(id = big_dat2.typ.low$id, age = big_dat2.typ.low$age, 
                                 brain = rowMeans(big_dat2.typ.low[, c(4:ncol(big_dat2.typ.low))]))
      avg_big_dat2.typ.high <- data.frame(id = big_dat2.typ.high$id, age = big_dat2.typ.high$age, 
                                brain = rowMeans(big_dat2.typ.high[, c(4:ncol(big_dat2.typ.high))]))
      avg_big_dat2.typ.mid <- data.frame(id = big_dat2.typ.mid$id, age = big_dat2.typ.mid$age, 
                                brain = rowMeans(big_dat2.typ.mid[, c(4:ncol(big_dat2.typ.mid))]))
      
      
      
      pdf(file= paste0(dir_in, '/', hem, '_', meas, '_', trk, '_', beh_names[[j]], '_', feat, '_clus_sig_nodes.pdf'), width = 6, height = 6)
      print(ggplot(avg_big, aes(x = beh, y = brain)) +
              geom_point(size = 2.5, stroke = NA) +
              geom_smooth(formula = 'y ~ x', method = 'lm', color = paste0(meas_color)) +
              ggtitle(trk) +
              theme(legend.position = "none") +
              scale_x_continuous(name = paste(beh_names[[j]], sep=' ')) +
              scale_y_continuous(name = y_ax) + 
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
      
      pdf(file= paste0(dir_in, '/', hem, '_', meas, '_', trk, '_', beh_names[[j]], '_', feat, '_clus_sig_nodes_split_traj.pdf'), width = 6, height = 6)
      print(ggplot(avg_big_dat1, aes(x = age, y = brain1, group = id)) +
              geom_point(size = 2.5, alpha = 0.1, stroke = NA) +
              geom_path(lwd = 0.5, alpha = 0.1) +
              ggtitle(trk) +
              scale_x_continuous(name = "Age (Months)", breaks = seq(0, 125, 25), limits = c(0, 125)) +
              scale_y_continuous(name = y_label) +
              theme(legend.position = "none", axis.text = element_text(size = 14, color = "black"), 
                    axis.title = element_text(size = 25),
                    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                    strip.text = element_text(size = 9, face = "bold"),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(colour = "black")) +
              geom_line(data = avg_big_dat2.typ.low, aes(x = age, y = brain), color = paste0(meas_color_low), lwd=2) +
              geom_line(data = avg_big_dat2.typ.high, aes(x = age, y = brain), color = paste0(meas_color_high), lwd=2) +
              geom_line(data = avg_big_dat2.typ.mid, aes(x = age, y = brain), color = paste0(meas_color), lwd=2, lty=1))
      dev.off()
        
      }
    }
  }
  
}

