gen_histograms_spec_beh <- function(dir_in, hem, other_dat){

# generate histograms for key behavioral variables

ids <- read.table(paste0(dir_in, '/', hem, '_final_ids.txt'))
colnames(ids)[colnames(ids) == 'V1'] = 'id'

other_dat <- inner_join(ids, other_dat, by = 'id', multiple = 'first')

tps <- c('inf', 'tod', 'pre', 'beg', 'rea')
hle_vars <- c('num_books_adult_ord', 'ebooks_adult_ord', 'books_kids_ord',
              'ebooks_kids_ord', 'reading_age_child_ord', 'child_reading_time_ord',
              'child_count', 'child_alphabet', 'fam_readership', 'fam_write_notes',
              'fam_write_creative', 'fam_joke')
raw_vars <- c('pre_word_access_raw',	'pre_word_fluency_raw',	'pre_substitution_raw',  
              'nepsy_phonproc_raw', 'pre_wjiv_oralcomp_raw', 'pre_wjiv_pvocab_raw', 
              'beg_wrmt3_wrd_id_raw','beg_wrmt3_wrdatt_raw', 'kbit2_raw')
stand_vars <- c('pre_phonproc_standard', 'pre_wjiv_orallang_standard', 
                'beg_wrmt3_wrd_id_stnd',	'beg_wrmt3_wrdatt_stnd', 'kbit2_standard')


hle_vars_names <- c('Number of adult books', 'Number of adult ebooks', 
                   'Number of childrens books', 'Number of childrens ebooks',
                   'Age first read to', 'Time read to per week', 
                   'Frequency with which family members teach the child to count',
                   'Frequency with which family members teach the child the alphabet',
                   'Frequency with which family members read books, newspapers, or magazines',
                   'Frequency with which family members write notes',
                   'Frequency with which family members write creatively',
                   'Frequency with which family members share rhymes or jokes with the child')
raw_vars_names <- c('WJ-IV Word access raw score', 'WJ-IV Word fluency raw score', 
                    'WJ-IV Substitution raw score', 'NEPSY-II Phonological processing raw score',
                    'WJ-IV Oral comprehension raw score', 'WJ-IV Picture vocabulary raw score',
                    'WRMT Word id raw score', 'WRMT Word attack raw score', 'KBIT-2 raw score')
stand_vars_names <- c('Phonological processing standard score', 'WJ-IV Oral language standard score',
                      'WRMT Word id standard score', 'WRMT Word attack standard score', 'KBIT-2 standard score')


print('Normality test for HLE')
print(shapiro.test(other_dat[, 'HLE']))
pdf(file = paste0(dir_in, '/hist_', hem, '_HLE_stand.pdf'))
  hist(other_dat[, 'HLE'], xlab = 'Home literacy environment (a.u.)', ylab = NULL, main = NULL, col = 'firebrick1',
       xlim = c(-1.5, 1.5), ylim = c(0, 20), cex.axis = 1.5, cex.lab = 1.5)
dev.off()


for (i in (1:length(hle_vars))){
  for (t in (1:length(tps))){
    print(paste0('Normality test for ', hle_vars[i]))
    print(shapiro.test(as.numeric(other_dat[, paste0(tps[t], '_', hle_vars[i])])))
    #mn1 <- min(na.omit(other_dat[, paste0(tps[t], '_', hle_vars[i])]))
    #mx1 <- max(na.omit(other_dat[, paste0(tps[t], '_', hle_vars[i])]))
    pdf(file = paste0(dir_in, '/hist_', hem, '_', tps[t], '_', hle_vars[i], '.pdf'))
    p1 <- hist(as.numeric(other_dat[, paste0(tps[t], '_', hle_vars[i])]), breaks = c(0, 1, 2, 3, 4, 5, 6, 7), xlim = c(0, 6), ylim = c(0, 40), # ylim = c(0, 300),
               xlab = NULL, ylab = NULL, main = NULL, plot = FALSE) # hle_vars_names[i]
    if (i == 1){
      barplot(p1$counts, col = 'gold', space = 0, xlim = c(0, 6), ylim = c(0, 50), ylab = NULL, width = rep(0.5, 6), cex.axis = 3, cex.lab = 3)
    } else {
      barplot(p1$counts, col = 'gold', space = 0, xlim = c(0, 6), ylim = c(0, 50), ylab = NULL, yaxt = 'n', width = rep(0.5, 6), cex.axis = 3, cex.lab = 3)
    }
    dev.off()
  }
}
  
for (i in (1:length(raw_vars))){
  print(paste0('Normality test for ', raw_vars[i]))
  print(shapiro.test(as.numeric(other_dat[, raw_vars[i]])))
  pdf(file = paste0(dir_in, '/hist_', hem, '_', raw_vars[i], '.pdf'))
  hist(as.numeric(other_dat[, raw_vars[i]]), breaks = seq(0, 40, 5), xlim = c(0, 40), ylim = c(0, 50), # ylim = c(0, 200), 
      col = 'deeppink', xlab = raw_vars_names[i], main = '', cex.axis = 1.5, cex.lab = 1.5)
  dev.off()
}

for (i in (1:length(stand_vars))){
  print(paste0('Normality test for ', stand_vars[i]))
  print(shapiro.test(as.numeric(other_dat[, stand_vars[i]])))
  pdf(file = paste0(dir_in, '/hist_', hem, '_', stand_vars[i], '.pdf'))
  hist(as.numeric(other_dat[, stand_vars[i]]), breaks = seq(40, 160, 5), xlim = c(40, 160), ylim = c(0, 20), # ylim = c(0, 100), 
       col = '#9999FF', xlab = stand_vars_names[i], main = '', cex.axis = 1.5, cex.lab = 1.5)
  dev.off()
} 

}

gen_histograms_spec_brain_refs <- function(dir_in, meas, hem, trk, med_elem, regions, beh_names){
  
  # generate histograms for key variables

# specify brain names as output by FS or py(Baby)AFQ
if (meas == 'dti_fa' | meas == 'dti_md') {
  if (trk_sz == 'nodes'){
    nodes <- c(1:100)
    brain_names <- regions <- paste0('n', nodes, '_', meas)
  } else if (trk_sz == 'quarters'){
    quarters <- c(1:4)
    brain_names <- regions <- paste0('q', quarters, '_', meas)
  } else if (trk_sz == 'average'){
    brain_names <- regions <- meas
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

if ( meas == 'dti_fa' | meas == 'dti_md' ){
  sz_in <- 1
  sz_out <- length(brain_names)
} else {
  sz_in <- 2
  sz_out <- length(regions) + 1
}

# Graph labels and colors
labels_colors <- graph_labels_colors(meas)
meas_color <- labels_colors$meas_color

dat4_int_mat1 <- med_elem$dat4_int_mat1
dat4_slope_mat1 <- med_elem$dat4_slope_mat1
#refs_eTIV <- med_elem$refs_eTIV

# To collect stats from Shapiro-Wilk normality test
W_sw <- matrix(ncol=length(beh_names), nrow=length(regions))
p_sw <- matrix(ncol=length(beh_names), nrow=length(regions))

for (feat in c('int', 'slope')){
  dat4_mat1 <- eval(as.symbol(paste0('dat4_', feat, '_mat1')))
  dat4_mat2 <- list()
  
  # reorganize data for mediation
  for (j in (1:length(beh_names))){
    dat4_mat2_reg <- list()
    
    for (i in (sz_in:sz_out)){ 
      dat4_mat2_reg[[i-(sz_in-1)]] <- dat4_mat1[[i]][[j]]
    }
    dat4_mat2[[j]] <- dat4_mat2_reg
  }

  
  for (j in (1:length(beh_names))){
    for (i in (1:(sz_out - sz_in + 1))){
        sd <- apply(dat4_mat2[[j]][[i]], 2, sd, na.rm = TRUE) #  scaling by the standard deviations without centering (https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/scale)
        sw <- shapiro.test(dat4_mat2[[j]][[i]]$brain3) # agostino.test
        W_sw[i,j] <- sw$statistic # sw$statistic[[2]]
        p_sw[i,j] <- sw$p.value
        if (meas != 'dti_fa' & meas != 'dti_md'){
          pdf(file = paste0(dir_in, '/hist_', hem, '_', meas, '_', regions[i], '_', beh_names[j], '_', feat, '.pdf'))
        } else {
          pdf(file = paste0(dir_in, '/hist_', hem, '_', meas, '_', trk, '_', regions[i], '_', beh_names[j], '_', feat, '.pdf'))
        }
        if (i == 1){
          hist(scale(dat4_mat2[[j]][[i]]$brain3, center = FALSE, scale = sd[2][[1]]),
              col = meas_color, breaks = seq(-6, 6, 1), 
              xlab = regions[i],  xaxt = 'n', xlim = c(-6, 6), # xlab = regions[i]
              ylab = NULL, yaxt = 'n', ylim = c(0, 50), 
              main = NULL, cex.axis = 2, cex.lab = 3)
          axis(side=1, at = 0, cex.axis = 2, pos = 0.5) # at = seq(-5, 5, 5), labels=seq(-5, 5, 5)
          axis(side=2, las = 2, at=seq(0, 50, 50), labels=seq(0, 50, 50), cex.axis = 2)
        } else {
          hist(scale(dat4_mat2[[j]][[i]]$brain3, center = FALSE, scale = sd[2][[1]]),
               col = meas_color, breaks = seq(-6, 6, 1), 
               xlab = regions[i], xaxt = 'n', xlim = c(-6, 6), 
               ylab = NULL, yaxt = 'n', ylim = c(0, 50), 
               main = NULL, cex.axis = 2, cex.lab = 3)
          axis(side=1, at = 0, cex.axis = 2, pos = 0.5) # at=seq(-5, 5, 5), labels=seq(-5, 5, 5)
        }
        dev.off()
    }
  }
  
  # Output normality stats
  print('Tests for normality statistics...')
  print(W_sw)
  print(p_sw)
}





}


