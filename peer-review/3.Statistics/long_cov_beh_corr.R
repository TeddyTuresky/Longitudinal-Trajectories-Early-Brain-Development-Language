long_cov_beh_corr <- function(dir_in, hem, covariates, beh_names, med_names){
  
  # performs correlations on behavioral measures after participant list was reduced from QC


ids <- read.table(paste0(dir_in, '/', hem, '_final_ids.txt'))
colnames(ids) <- 'id'
var_names <- c(covariates, beh_names, med_names)

# reduce to necessary columns
dat <- left_join(ids, other_dat, by = 'id', multiple = 'first')
dat_beh_med <- dat[, c(var_names)]

# run correlations
print(cor(dat_beh_med, use = "pairwise.complete.obs"))

# for p-values, load library(Hmisc) and use rcorr(as.matrix(dat_beh_med)), but watch out for dependency issues

}