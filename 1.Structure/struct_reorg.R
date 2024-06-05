rm(list = ls())
library(dplyr)


dir_in <- '<insert_input_directory>' # directory containing FS output from consol_stats.sh script



for (hem in c('lh', 'rh')){
  
for (meas in c('volume', 'area', 'thickness', 'meancurv', 'curvind', 'foldind', 'gauscurv', 'thicknessstd')){

  mri <- read.table(paste0(dir_in, hem, '_', meas, '_aparc_stats.txt'), header = TRUE)
  
  sub <- mri[, c(1)]
  sub_tp <- strsplit(sub, '_')
  sub_tp_mat <- matrix(unlist(sub_tp),ncol=2,byrow=T)
  id <- sub_tp_mat[, c(2)]
  timepoint <- sub_tp_mat[, c(1)]
  dat <- data.frame(id, timepoint, mri[, c(2:length(mri))])
  dat <- dat[order(dat$id, decreasing = FALSE), ]
  
  write.csv( dat, file = paste0(dir_in, '/reorg/', hem, '_', meas, '_output.csv'), row.names = FALSE)

}

}


mri <- read.table(paste0(dir_in, 'volume_wm_stats.txt'), header = TRUE)
sub <- mri[, c(1)]
sub_tp <- strsplit(sub, '_')
sub_tp_mat <- matrix(unlist(sub_tp),ncol=2,byrow=T)
id <- sub_tp_mat[, c(2)]
timepoint <- sub_tp_mat[, c(1)]
dat <- data.frame(id, timepoint, mri[, c(2:length(mri))])
dat <- dat[order(dat$id, decreasing = FALSE), ]

write.csv( dat, file = paste0(dir_in, '/reorg/volume_wm_output.csv'), row.names = FALSE)