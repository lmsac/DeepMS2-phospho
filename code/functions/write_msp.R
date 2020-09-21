write.msp = function(spectra, file) {
  writeLines(sapply(spectra, function(x) {
    s = paste0('Name: ', x$peptide, '/', x$charge)
    if (!is.null(x$mw)) {
      s = c(s, paste0('MW: ', x$mw))
    }
    
    if (!is.null(x$comment)) {
      comment = x$comment
    }
    else if (!is.null(x$fields)) {
      comment = x$fields
    }
    else {
      comment = c()
    }
    if (!is.null(x$pepMass)) {
      comment['Parent'] = x$pepMass
    }
    if (!is.null(x$mods) && length(x$mods$pos) > 0) {
      comment['Mods'] = paste0(
        c(
          length(x$mods$pos),
          sapply(1:length(x$mods$pos), function(i) {
            pos = x$mods$pos[i]
            name = x$mods$name[i]
            paste0(
              pos - 1, ',', 
              substring(x$peptide, pos, pos), ',',
              name
            )
          })
        ),
        collapse = '/'
      )
    }
    else {
      comment['Mods'] = 0
    }
    
    if (length(comment) > 0) {
      s = c(s, paste0(
        'Comment: ',
        paste0(sapply(1:length(comment), function(i) {
          paste0(names(comment)[i], '=', comment[i])
        }), collapse = ' ')
      ))
    }
    
    s = c(s, paste0('Num peaks: ', nrow(x$ions)))
    if (nrow(x$ions) > 0) {
      if (max(x$ions[, 2]) > 0) {
        x$ions[, 2] = x$ions[, 2] / max(x$ions[, 2]) * 10000
      }
      if (is.null(x$ion.labels)) {
        x$ion.labels = rep('?', nrow(x$ions))
      }
      s = c(s, paste(x$ions[, 1], x$ions[, 2], x$ion.labels, sep = '\t'))
    }
    
    paste0(paste0(s, collapse = '\n'), '\n')
  }), file)
}
