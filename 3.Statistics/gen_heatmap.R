gen_heatmap <- function(stats_sum, meas, dir_in, regions, covariates){

  # stats_sum is a list of statistics from lme models

  
  # Graph labels and colors
  labels_colors <- graph_labels_colors(meas)
  y_label <- labels_colors$y_label
  y_ax <- labels_colors$y_ax
  y_scal <- labels_colors$y_scal
  y_lims <- labels_colors$y_lims
  meas_color <- labels_colors$meas_color
  meas_color_low <- labels_colors$meas_color_low
  meas_color_high <- labels_colors$meas_color_high
  
  # overwrite region names for diffusion
  if (meas == 'dti_fa' | meas == 'dti_md'){
    if (trk_sz == 'nodes'){
      nodes <- c(1:100)
      regions <- paste0('n', nodes, '_', meas)
    } else if (trk_sz == 'quarters'){
      quarters <- c(1:4)
      regions <- paste0('q', quarters, '_', meas)
    } else if (trk_sz == 'average'){
      regions <- meas
    }
  }
  
  # organize covariates
  if ( meas == 'dti_fa' | meas == 'dti_md' ){
    sz_in <- 1
    sz_out <- length(regions)
  } else {
    sz_in <- 2
    sz_out <- length(regions) + 1
  }
  p_covs <- stats_sum$p_covs[c(sz_in:sz_out),]
  p_covs <- matrix(p_covs, nrow = length(regions), ncol = length(covariates))
  rownames(p_covs) <- regions
  colnames(p_covs) <- covariates

  

  p_covs_melt <- melt(p_covs)
  colnames(p_covs_melt) <- c('Region', 'Covariate', 'value')

  pdf(file= paste0(dir_in, '/heatmap_', meas, '.pdf'))
  print(ggplot(p_covs_melt, aes(Region, Covariate)) +
    geom_tile(aes(fill = value), color = 'white', lwd = 1.5) +
    scale_fill_gradientn(colors=c('black', paste0(meas_color)),
                       values=c(1, 0), guide="colorbar") +
    coord_fixed() +
    guides(fill = guide_colorbar(barwidth = 0.5, title = "p-value")) +
    theme(legend.position = "none", axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.background = element_blank(), axis.ticks = element_blank()))
  dev.off()
  
}
  
  