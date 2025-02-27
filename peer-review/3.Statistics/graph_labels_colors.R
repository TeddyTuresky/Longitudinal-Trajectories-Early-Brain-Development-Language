graph_labels_colors <- function(meas){
  
  labels_colors <- list()
  
  if (meas == 'wm'){
    y_label <- 'White Matter Volume'
    y_ax <- bquote(.(y_label) ~ (mm^3))
    y_scal <- waiver()
    y_lims <- NULL
    meas_color <- '#0033FF'
    meas_color_low <- 'deepskyblue'
    meas_color_high <- 'darkblue'
  } else if (meas == 'dti_fa'){
    y_label <- 'Fractional Anisotropy'
    y_ax <- y_label
    y_scal <- seq(0.25, 0.55, 0.05)
    y_lims <- c(0.25, 0.55)
    meas_color <- '#00CC66'
    meas_color_low <- '#00FF99'
    meas_color_high <- 'darkgreen'
  } else if (meas == 'dti_md'){
    y_label <- 'Mean Diffusivity'
    y_ax <- y_label
    y_scal <- seq(0.0007, 0.0012, 0.0001)
    y_lims <- c(0.0007, 0.0012)
    meas_color <- '#00CC66'
    meas_color_low <- '#00FF99'
    meas_color_high <- 'darkgreen'
  } else if (meas == 'volume'){
    y_label <- 'Gray Matter Volume'
    y_ax <- bquote(.(y_label) ~ (mm^3))
    y_scal <- waiver()
    y_lims <- NULL
    meas_color <- '#0033FF'
    meas_color_low <- 'deepskyblue'
    meas_color_high <- 'darkblue'
  } else if ( meas == 'area'){
    y_label <- 'Surface Area'
    y_ax <- bquote(.(y_label) ~ (mm^2))
    y_scal <- waiver()
    y_lims <- NULL
    meas_color <- '#0033FF'
    meas_color_low <- 'deepskyblue'
    meas_color_high <- 'darkblue'
  } else if ( meas == 'thickness'){
    y_label <- 'Cortical Thickness'
    y_ax <- bquote(.(y_label) ~ (mm))
    y_scal <- waiver()
    y_lims <- NULL
    meas_color <- '#0033FF'
    meas_color_low <- 'deepskyblue'
    meas_color_high <- 'darkblue'
  } else if ( meas == 'meancurv'){
    y_label <- 'Mean Curvature'
    y_ax <- bquote(.(y_label) ~ (1/mm))
    y_scal <- waiver()
    y_lims <- NULL
    meas_color <- '#0033FF'
    meas_color_low <- 'deepskyblue'
    meas_color_high <- 'darkblue'
  } else if ( meas == 'gausscurv'){
    y_label <- 'Gaussian Curvature'
    y_ax <- bquote(.(y_label) ~ (1/mm^2))
    y_scal <- waiver()
    y_lims <- NULL
    meas_color <- '#0033FF'
    meas_color_low <- 'deepskyblue'
    meas_color_high <- 'darkblue'
  } else {
    print('Unknown measure')
  }
  
  labels_colors <- list('y_label' = y_label, 'y_ax' = y_ax, 'y_scal' = y_scal, 'y_lims' = y_lims,
                        'meas_color' = meas_color, 'meas_color_low' = meas_color_low, 
                        'meas_color_high' = meas_color_high)
  return(labels_colors)
  
}