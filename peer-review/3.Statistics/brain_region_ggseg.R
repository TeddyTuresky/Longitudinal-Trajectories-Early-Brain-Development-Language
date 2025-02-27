brain_region_ggseg <- function(meas, hem, dir_in, indx, fdr_cor){
  
  # Plots FreeSurfer-labeled brain regions using the ggseg package
  # For questions: theodore_turesky@gse.harvard.edu

if (hem == 'lh'){
  hemi = 'left'
} else {
  hemi = 'right'
}

# Specify color(s) to use
labels_colors <- graph_labels_colors(meas)
meas_color <- labels_colors$meas_color

# extract important elements
regions <- indx[,1]
names_covORbeh <- indx[,2]
p_vals <- indx[,3]
type_stat <- indx[,4]

# rename regions for ggseg
regions_all <- c('bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal',
                 'cuneus', 'entorhinal', 'fusiform', 'inferiorparietal',
                 'inferiortemporal', 'isthmuscingulate', 'lateraloccipital', 
                 'lateralorbitofrontal', 'lingual', 'medialorbitofrontal',
                 'middletemporal', 'parahippocampal', 'paracentral', 'parsopercularis',
                 'parsorbitalis', 'parstriangularis', 'pericalcarine', 'postcentral',
                 'posteriorcingulate', 'precentral', 'precuneus', 'rostralanteriorcingulate',
                 'rostralmiddlefrontal', 'superiorfrontal', 'superiorparietal',
                 'superiortemporal', 'supramarginal', 'frontalpole', 'temporalpole',
                 'transversetemporal', 'insula')

regions_ggseg <- c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal',
                 'cuneus', 'entorhinal', 'fusiform', 'inferior parietal',
                 'inferior temporal', 'isthmus cingulate', 'lateral occipital', 
                 'lateral orbitofrontal', 'lingual', 'medial orbitofrontal',
                 'middle temporal', 'parahippocampal', 'paracentral', 'pars opercularis',
                 'pars orbitalis', 'pars triangularis', 'pericalcarine', 'postcentral',
                 'posterior cingulate', 'precentral', 'precuneus', 'rostral anterior cingulate',
                 'rostral middle frontal', 'superior frontal', 'superior parietal',
                 'superior temporal', 'supramarginal', 'frontal pole', 'temporal pole',
                 'transverse temporal', 'insula')


# Map brain regions individually
if (fdr_cor == 1){
  for (i in c(1:length(regions))){
    reg1 <- regions_ggseg[grepl(paste0('\\b', regions[i], '\\b'), regions_all)]
    someData = data.frame(region = reg1, p = 1, stringsAsFactors = FALSE); 
    someData$hemi = hemi
    ggseg(.data = someData, mapping = aes(fill = factor(p)), hemisphere = hemi, view = 'lateral', color = "black", size = 0.2) +
      theme_void() +
      scale_fill_manual(
        breaks = '1', 
        values = meas_color,
        na.value = "grey") +
      theme(legend.position = "none")
  
    ggsave(paste0(dir_in, '/ggseg_', meas, '_', hem, '_', regions[i], '_', names_covORbeh[i], '_', type_stat[i], '.pdf'), plot = last_plot(), device=cairo_pdf)
  }
}


# Map multiple brain areas and assign estimates or p-values
u_names_covORbeh <- unique(names_covORbeh)
u_type_stat <- unique(type_stat)
  for (s in u_type_stat){
    for (t in u_names_covORbeh){
    type_stat_indx <- indx[indx[,4] == s & indx[,2] == t, ]
    if (length(type_stat_indx) == 0){
      next
    } else if (length(type_stat_indx) == 4){
        type_stat_indx <- matrix(type_stat_indx, nrow = 1, ncol = 4)
    }
    
    regions <- type_stat_indx[,1]
    names_covORbeh <- type_stat_indx[,2] 
    p_vals <- type_stat_indx[,3]
    
    # rename for ggseg
    for (i in c(1:length(regions))){
      reg2a <- regions_ggseg[grepl(paste0('\\b', regions[i], '\\b'), regions_all)]
      if (i == 1){
        reg2b <- reg2a
      } else {
        reg2b <- c(reg2b, reg2a)
      }
    }
  

  someData = data.frame(region = reg2b, p = p_vals, stringsAsFactors = FALSE) 
  #someData$hemi = hemi
  #someData$view = 'lateral'
  
  someData <- someData %>%
    mutate_at(vars(p), as.numeric)

  if (fdr_cor == 0){
    someData$r <- someData$p # p-values are really r-values
    
    ggseg(.data = someData, mapping = aes(fill = r), hemisphere = hemi, view = 'lateral', color = "black", size = 0.2) +
      theme_void() +
      scale_fill_gradientn(colors = c('cyan', 'white', 'deeppink'), # heat.colors(10) purple4 - deeppink
                           limits = c(-1, 1), 
                           #values = c(-1, 0, 1), 
                           na.value = "grey", 
                           guide="colorbar")
    ggsave(paste0(dir_in, '/ggseg_', meas, '_', hem, '_multi-region_', names_covORbeh[i], '_', s, '_stats_est_lateral.pdf'), plot = last_plot(), device=cairo_pdf)
    
    ggseg(.data = someData, mapping = aes(fill = p), hemisphere = hemi, view = 'medial', color = "black", size = 0.2) +
      theme_void() +
      scale_fill_gradientn(colors = c('cyan', 'white', 'deeppink'), # heat.colors(10) purple4 - deeppink
                           limits = c(-1, 1), 
                           #values = c(-1, 0, 1), 
                           na.value = "grey", 
                           guide="colorbar")
    ggsave(paste0(dir_in, '/ggseg_', meas, '_', hem, '_multi-region_', names_covORbeh[i], '_', s, '_stats_est_medial.pdf'), plot = last_plot(), device=cairo_pdf)
  
    } else if (fdr_cor == 1){
    
    ggseg(.data = someData, mapping = aes(fill = p), hemisphere = hemi, view = 'lateral', color = "black", size = 0.2) +
      theme_void() +
      scale_fill_gradientn(colors = c('white', 
                                             paste0(meas_color)), 
                           limits = c(0, 0.05), 
                           values = c(1, 0), 
                           na.value = "grey", 
                           guide="colorbar")
    ggsave(paste0(dir_in, '/ggseg_', meas, '_', hem, '_multi-region_', names_covORbeh[i], '_', s, '_stats_fdr.pdf'), plot = last_plot(), device=cairo_pdf)
  
    }
    } 
}
}


brain_region_ggseg_med <- function(meas, hem, dir_in, indx_med_all){
  
  # Plots FreeSurfer-labeled brain regions using the ggseg package
  # For questions: theodore_turesky@gse.harvard.edu
  
  if (hem == 'lh'){
    hemi = 'left'
  } else {
    hemi = 'right'
  }
  
  # Specify color(s) to use
  labels_colors <- graph_labels_colors(meas)
  meas_color <- labels_colors$meas_color
  
  # rename regions for ggseg
  regions_all <- c('bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal',
                   'cuneus', 'entorhinal', 'fusiform', 'inferiorparietal',
                   'inferiortemporal', 'isthmuscingulate', 'lateraloccipital', 
                   'lateralorbitofrontal', 'lingual', 'medialorbitofrontal',
                   'middletemporal', 'parahippocampal', 'paracentral', 'parsopercularis',
                   'parsorbitalis', 'parstriangularis', 'pericalcarine', 'postcentral',
                   'posteriorcingulate', 'precentral', 'precuneus', 'rostralanteriorcingulate',
                   'rostralmiddlefrontal', 'superiorfrontal', 'superiorparietal',
                   'superiortemporal', 'supramarginal', 'frontalpole', 'temporalpole',
                   'transversetemporal', 'insula')
  
  regions_ggseg <- c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal',
                     'cuneus', 'entorhinal', 'fusiform', 'inferior parietal',
                     'inferior temporal', 'isthmus cingulate', 'lateral occipital', 
                     'lateral orbitofrontal', 'lingual', 'medial orbitofrontal',
                     'middle temporal', 'parahippocampal', 'paracentral', 'pars opercularis',
                     'pars orbitalis', 'pars triangularis', 'pericalcarine', 'postcentral',
                     'posterior cingulate', 'precentral', 'precuneus', 'rostral anterior cingulate',
                     'rostral. middle frontal', 'superior frontal', 'superior parietal',
                     'superior temporal', 'supramarginal', 'frontal pole', 'temporal pole',
                     'transverse temporal', 'insula')
  
  
  
  # Map multiple brain areas and assign p-values
  
  for (k in 1:length(indx_med_all)){
    
    # extract important elements
    indx_med <- indx_med_all[[k]]
    
    names_regions <- indx_med[,1]
    names_beh <- indx_med[,2] 
    names_med <- indx_med[,3]
    names_feat <- indx_med[,4]
    
    # rename for ggseg
    for (i in c(1:length(names_regions))){
      reg2a <- regions_ggseg[grepl(paste0('\\b', names_regions[i], '\\b'), regions_all)]
      if (i == 1){
        reg2b <- reg2a
      } else {
        reg2b <- c(reg2b, reg2a)
      }
    }
    
    
    someData = data.frame(region = reg2b, p = c(1:length(names_regions)), stringsAsFactors = FALSE) 
    someData$hemi = hemi
    someData$view = 'lateral'
    
    someData <- someData %>%
      mutate_at(vars(p), as.numeric)
    
    ggseg(.data = someData, mapping = aes(fill = factor(p)), hemisphere = hemi, view = 'lateral', color = "black", size = 0.2) +
      theme_void() +
      scale_fill_manual(
        breaks = paste0(c(1:length(names_regions))), 
        values = rep(meas_color, length(names_regions)),
        na.value = "grey") +
      theme(legend.position = "none")
    
    ggsave(paste0(dir_in, '/ggseg_', meas, '_', hem, '_multi-region_', names_beh[1], '_', names_med[1], '_', names_feat[1], '_med.pdf'), plot = last_plot(), device=cairo_pdf)
    
  }
}