long_graph <- function(dir_in, measures, hem, trks, covariates){

  

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
    dat <- dat[, c('id', 'timepoint', 'age', covariates)]
    if (meas == measures[1]){
      dat_all_meas <- dat 
    } else {
      dat_all_meas <- full_join(dat_all_meas, dat, by = c('id', 'timepoint', 'age', covariates))
    }
  }
}
  


# retain necessary columns
dat1 <- dat_all_meas[, c('id', 'timepoint', 'age', covariates)]

dat1$age <- as.numeric(dat1$age)
dat1$Sex <- as.factor(dat1$Sex)

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

# add numbered ids
u <- unique(dat2$id)
u_num <- c(1:length(u))
dat_u_num <- data.frame(id = u, id_num = u_num)
dat2 <- left_join(dat2, dat_u_num, by = 'id', multiple = 'all')

print(paste0('Total number of participants graphed: ', length(u)))

pdf(file = paste0(dir_in, '/', hem, '_long_dataset_1.pdf'), width = 3, height = 4)
print(ggplot(dat2, aes(x = age, y = as.numeric(reorder(id_num, age)), color = Sex)) + # to order by number of time points, change 'age' to 'n'
  scale_color_manual(values=c("deepskyblue3", "#ff8c00")) +
  geom_line(aes(group = id_num), linewidth = 0.5) +
  geom_point(aes(group = id_num), size = 0.7) +
  scale_x_continuous(name = "Age (Months)") +
  scale_y_continuous(name = 'Participant ID', limits = c(-1, length(unique(dat2$id_num)) + 2), expand=c(0,0)) + 
  theme(legend.position = "none", axis.text = element_text(size = 6, color = "black"),
        axis.line.x = element_blank(), #element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.line.y = element_blank(), # element_line(linewidth = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"),
        #axis.text.y = element_blank(), # element_text(size = 2), # 
        #axis.ticks.y = element_blank(),
        strip.text = element_text(size = 8, face = "bold", margin = margin(1, 0, 1, 0)),
        panel.background = element_blank(), # element_rect(color = "deepskyblue3", linewidth = 3)
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_line(color = "black")))

dev.off()

}