long_descriptive_stats <- function(dir_in, measures, hem, trks, covariates, beh_names){
  
  
# load full dataframes for each measure and tract
for (meas in measures){
  i <- 0
  j <- 1
  while (i == 0){
    i = i + 1
    if (meas == 'dti_fa' | meas == 'dti_md'){
      if (j < length(trks)){ # needed to loop through tracks
        i <- i - 1 # reset i to maintain while loop
      }
      trk <- trks[j] # changes if dti_fa or dti_md is specified
      j <- j + 1
      load(paste0(dir_in, '/', hem, '_', meas, '_', trk, '_full_dataset_qc.RData'))
    } else {
      load(paste0(dir_in, '/', hem, '_', meas, '_full_dataset_qc.RData'))
    }
    dat <- dat[, c('id', 'timepoint', 'age', 'pre_beh_age', covariates, beh_names)]
    if (meas == measures[1]){
      dat_all_meas <- dat 
    } else {
      dat_all_meas <- full_join(dat_all_meas, dat, by = c('id', 'timepoint', 'age', 'pre_beh_age', covariates, beh_names))
    }
  }
}  

# retain necessary columns
dat1 <- dat_all_meas[, c('id', 'timepoint', 'age', 'pre_beh_age', covariates, beh_names)]

dat1$age <- as.numeric(dat1$age)
dat1 <- dat1[order(dat1$id, dat1$age), ]

dat1_order <- dat1 %>% 
  group_by(id) %>% 
  count()

dat1 <- left_join(dat1, dat1_order, by = 'id')
dat2 <- dat1[dat1$n != 1, ]

# remove rows in which covariates are missing
dat2_cov <- dat2[, c('id', 'timepoint', covariates)]
dat2_cov <- na.omit(dat2_cov)
dat2 <- inner_join(dat2, dat2_cov, by = c('id', 'timepoint', covariates))

# ensure at least one infant/toddler timepoint
dat2_inf_tod <- dat2[dat2$timepoint == 'INF' | dat2$timepoint == 'TOD', ]
u_inf_tod <- unique(dat2_inf_tod$id)
du_inf_tod <- data.frame(id = u_inf_tod)
dat2 <- inner_join(dat2, du_inf_tod, by = 'id')

uid_dat2 <- unique(dat2$id)
duid_dat2 <- data.frame(id = uid_dat2)
u_dat2 <- left_join(duid_dat2, dat2, by = 'id', multiple = 'first')

write.table(u_dat2$id, paste0(dir_in, '/', hem, '_final_ids.txt'), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

u_dat2 <- u_dat2[, -c(1, 2)]


# for other covariates
mat_covs_beh <- matrix(nrow = 6, ncol = length(u_dat2))
for (i in 1:length(u_dat2)){
  mat_vec <- data.frame(vec = u_dat2[, c(i)], n = u_dat2$n)
  mat_vec <- na.omit(mat_vec)
  vec <- mat_vec$vec
  mat_covs_beh[1,i] <- min(vec)
  mat_covs_beh[2,i] <- max(vec)
  mat_covs_beh[5,i] <- length(vec)
  mat_covs_beh[6,i] <- sum(mat_vec$n)
  
  # count binary and average/sd ordinal/continuous variables
  if (length(unique(vec)) == 2){
    mat_covs_beh[3,i] <- length(vec[which(vec %in% mat_covs_beh[1,i])])
    mat_covs_beh[4,i] <- length(vec[which(vec %in% mat_covs_beh[2,i])])
  } else {
    mat_covs_beh[3,i] <- mean(vec)
    mat_covs_beh[4,i] <- sd(vec)
    
  }
  
}

rownames(mat_covs_beh) <- c('minimum', 'maximum', 'n. min. | average', 'n. max. | sd', 'n. participants', 'n. observations')
colnames(mat_covs_beh) <- c('age', 'pre_beh_age', covariates, beh_names, 'n')
print(mat_covs_beh)
write.csv(mat_covs_beh, paste0(dir_in, '/', hem, '_desc_stats.csv'), row.names = TRUE)

}