graph_all_meas_regs <- function(dir_in, measures, hem, trk){
  
  k <- 0
  for (meas in measures){
    
    if (meas == 'dti_fa' | meas == 'dti_md'){
      modality <- 'dwi'
      tract = trk
    } else {
      modality <- 'struct'
      tract <- 'none'
    }
    
    k <- k + 1
    load(paste0(dir_in, '/dat1.typ.all_', meas, '_', hem, '_', tract, '.RData'))
    if (k == 1){
      dat1.typ.all.meas <- dat1.typ.all
    } else {
      dat1.typ.all.meas <- rbind(dat1.typ.all.meas, dat1.typ.all)
    }
  
  }
  
  # average dwi tracts and recompile dataframe
  dat1.typ.struct <- dat1.typ.all.meas[dat1.typ.all.meas$modality == 'struct', ]
  dat1.typ.dwi <- dat1.typ.all.meas[dat1.typ.all.meas$modality == 'dwi', ]
  
  dat1.typ.dwi.fa <- dat1.typ.dwi[dat1.typ.dwi$group == 'dti_fa', ]
  dat1.typ.dwi.md <- dat1.typ.dwi[dat1.typ.dwi$group == 'dti_md', ]
  dat1.typ.dwi.md$brain_pred_all <- 1 / dat1.typ.dwi.md$brain_pred_all
  
  dat1.typ.dwi.all1 <- rbind(dat1.typ.dwi.fa, dat1.typ.dwi.md)
  
  dat1.typ.dwi.all2 <- dat1.typ.dwi.all1 %>%
    group_by(age_act) %>% 
    summarize_at(vars('brain_pred_all'), mean)
  
  dat1.typ.dwi.all3 <- data.frame(id = dat1.typ.dwi.all1$id, group = dat1.typ.dwi.all1$group, modality = dat1.typ.dwi.all1$modality, 
    tract = 'all', age_act = dat1.typ.dwi.all1$age_act, brain_pred_all = dat1.typ.dwi.all1$brain_pred_all)
  
  dat1.typ.all.mod <- rbind(dat1.typ.struct, dat1.typ.dwi.all3)
  
  dat1.typ.all.change <- dat1.typ.all.mod %>%
    group_by(group) %>%
    # mutate("{meas}" := brain_pred_all/first(brain_pred_all) * 100)
    mutate(Percent_Change = (brain_pred_all/first(brain_pred_all) * 100) - 100)
  

  pdf(file= paste0(dir_in, '/', hem, '_all_meas_graph.pdf'), width = 6, height = 6)
  print(ggplot(dat1.typ.all.change, aes(x = age_act, y = Percent_Change, group = group, color = modality)) + # 'modality'
          geom_path(lwd = 1, lty=1) + # col = '#0033FF', '#00CC66' aes(col = modality), 
          scale_x_continuous(name = "Age (Months)", breaks = seq(0, 125, 25)) +
          scale_y_continuous(name = 'Change (%)', breaks = seq(0, 400, 50), limits = c(0, 400)) +  # limits = y_lims
          #scale_color_manual(values = c('gold', 'slateblue2', 'lightblue4', 'cornflowerblue', 'firebrick1', 'deeppink', 'magenta4')) + 
          scale_color_manual(values = c('#00CC66', '#0033FF')) + 
          theme(axis.text = element_text(size = 14, color = "black"), # legend.position = "none",
                axis.title = element_text(size = 25),
                plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                strip.text = element_text(size = 9, face = "bold"),
                panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(colour = "black")))
  dev.off() 
}