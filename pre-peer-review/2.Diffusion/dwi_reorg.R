rm(list = ls())
library(dplyr)


dir_in <- '<insert_input_directory>'
trk <- 'SLF_L' # specify tract 



subs <- list.files(paste0(dir_in, '/nodes'))
for (i in c(1:length(subs))){
  tryCatch({
  sub <- subs[i]
  id <- strsplit(sub, '_')[[1]][2]
  timepoint <- strsplit(sub, '_')[[1]][1]
  dwi <- read.csv( paste0( dir_in, '/nodes/', sub ) )
  dwi$nodeID <- dwi$nodeID + 1
  dwi_red <- dwi[, c('tractID', 'nodeID', 'dti_fa', 'dti_md')]

  dwi_trk <- dwi_red[dwi_red$tractID == trk, ]
  dwi_trk$id <- id
  dwi_trk$timepoint <- timepoint
  
  # quarter tracts
  dwi_trk$quarter <- 'q1'
  dwi_trk$quarter[dwi_trk$nodeID > 25 & dwi_trk$nodeID <= 50 ] = 'q2'
  dwi_trk$quarter[dwi_trk$nodeID > 50 & dwi_trk$nodeID <= 75 ] = 'q3'
  dwi_trk$quarter[dwi_trk$nodeID > 75 & dwi_trk$nodeID <= 100 ] = 'q4'
  
  col_order <- c('id', 'timepoint', 'nodeID', 'quarter',
                 'dti_fa', 'dti_md')
  dwi_trk <- dwi_trk[, col_order]
  
  if (i == 1){
    dwi_trk_all <- dwi_trk
  }
  else {
    dwi_trk_all <- rbind(dwi_trk_all, dwi_trk)
  }
  
  }, error=function(e){cat("ERROR :", sub, " does not have track: ", trk, "\n", conditionMessage(e), "\n")})
  
}



# generate averages overall and by section
dwi_trk_avg <- dwi_trk_all %>% group_by(id, timepoint) %>% summarize_at(vars('dti_fa', 'dti_md'), mean)
dwi_trk_quart <- dwi_trk_all %>% group_by(id, timepoint, quarter) %>% summarize_at(vars('dti_fa', 'dti_md'), mean)
dwi_trk_all <- dwi_trk_all[, c('id', 'timepoint', 'nodeID', 'dti_fa', 'dti_md')]

# restructure quartered dataset to "wide" format
dwi_trk_quart_q <- list()
for (i in c(1:4)){
  quart <- paste0('q', i)
  dwi_trk_quart_q[[i]] <- dwi_trk_quart[dwi_trk_quart$quarter == quart, ]
  colnames(dwi_trk_quart_q[[i]])[4] <- paste0(quart, '_dti_fa')
  colnames(dwi_trk_quart_q[[i]])[5] <- paste0(quart, '_dti_md')
  dwi_trk_quart_q[[i]] <- subset( dwi_trk_quart_q[[i]], select = -c(quarter))
  if (i == 1){
    dwi_trk_quarts <- dwi_trk_quart_q[[i]]
  } else {
    dwi_trk_quarts <- inner_join( dwi_trk_quarts, dwi_trk_quart_q[[i]], by = c('id', 'timepoint'))
  }
}

# restructure nodes dataset to "wide" format
dwi_trk_node_n <- list()
for (i in c(1:100)){
  node <- paste0('n', i)
  dwi_trk_node_n[[i]] <- dwi_trk_all[dwi_trk_all$nodeID == i, ]
  colnames(dwi_trk_node_n[[i]])[4] <- paste0(node, '_dti_fa')
  colnames(dwi_trk_node_n[[i]])[5] <- paste0(node, '_dti_md')
  dwi_trk_node_n[[i]] <- subset( dwi_trk_node_n[[i]], select = -c(nodeID))
  if (i == 1){
    dwi_trk_nodes <- dwi_trk_node_n[[i]]
  } else {
    dwi_trk_nodes <- inner_join( dwi_trk_nodes, dwi_trk_node_n[[i]], by = c('id', 'timepoint'))
  }
}


write.csv( dwi_trk_avg, paste0( dir_in, '/reorg/dwi_trk_avg_', trk, '.csv'), row.names = FALSE)
write.csv( dwi_trk_quarts, paste0( dir_in, '/reorg/dwi_trk_quarts_', trk, '.csv'), row.names = FALSE)
write.csv( dwi_trk_nodes, paste0( dir_in, '/reorg/dwi_trk_nodes_', trk, '.csv'), row.names = FALSE)


