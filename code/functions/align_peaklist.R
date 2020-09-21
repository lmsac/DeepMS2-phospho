align.peaklists = function(peaklists, tolerance) {
  all.peaks.mz = unlist(sapply(peaklists, function(x) x[, 1]))
  all.peaks.int = unlist(sapply(peaklists, function(x) x[, 2]))
  all.peaks.from = unlist(sapply(1:length(peaklists), function(i) rep(i, nrow(peaklists[[i]]))))
  all.peaks.used = rep(FALSE, length(all.peaks.from))
  
  result = do.call(rbind, lapply(order(all.peaks.int, decreasing = T), function(i) {
    if (!all.peaks.used[i]) {
      all.peaks.used[i] <<- TRUE
      mz = all.peaks.mz[i]
      from = all.peaks.from[i]
      
      result.mz = rep(NA, length(peaklists))
      result.int = rep(0, length(peaklists))
      result.mz[from] = mz
      result.int[from] = all.peaks.int[i]
      
      candidate.indexs = which(
        !all.peaks.used & 
          all.peaks.mz >= mz * (1 - tolerance * 1e-6) & 
          all.peaks.mz <= mz * (1 + tolerance * 1e-6)
      )
      if (length(candidate.indexs) > 0) {
        lapply((1:length(peaklists))[-from], function(j) {
          index.j = candidate.indexs[which(all.peaks.from[candidate.indexs] == j)]
          if (length(index.j) > 1) {
            index.j = index.j[which.min(abs(all.peaks.mz[index.j] - mz))]
          }
          if (length(index.j) > 0) {
            result.mz[j] <<- all.peaks.mz[index.j]
            result.int[j] <<- all.peaks.int[index.j]
            all.peaks.used[index.j] <<- TRUE
          }
        })
      }
      c(result.mz, result.int)
    }
    else {
      NULL
    }
  }))
  
  if (length(peaklists) > 1 && nrow(result) > 1) {
    result = result[order(rowMeans(result[, 1:length(peaklists)], na.rm = TRUE)), ]
  }
  
  colnames(result) = c(
    paste0('mz', 1:length(peaklists)), 
    paste0('int', 1:length(peaklists))
  )
  
  list(
    mz = subset(result, select = 1:length(peaklists)),
    intensity = subset(result, select = (1 + length(peaklists)):(length(peaklists) * 2))
  )
}
