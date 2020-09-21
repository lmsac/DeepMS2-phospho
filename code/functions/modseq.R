get.modified.sequence = function(sequence, modification, symbol = c('Phospho' = '*', 'Oxidation' = '@')) {
  if (is.na(modification) || is.null(modification) || nchar(modification) == 0) {
    return(sequence)
  }
  
  s = strsplit(modification, ';')[[1]]
  name = gsub('[()]', '', stringr::str_extract(s, '\\(.+\\)$'))
  pos = as.integer(gsub('[A-Z(]', '', stringr::str_extract(s, '^[A-Z]+[0-9]+\\(')))
  
  if (length(pos) == 0) {
    return(sequence)
  }
  else {
    aa = strsplit(sequence, '')[[1]]
    sapply(1:length(pos), function(i) {
      pos = pos[i]
      name = name[i]
      sym = symbol[name]
      if (is.na(sym)) {
        stop(paste0('unknown modification "', name, '"'))
      }
      aa[pos] <<- paste0(aa[pos], symbol[name])
    })
    paste0(aa, collapse = '')
  }
}

get.sequence.modification = function(modifiedSequence, symbol = c('*' = 'Phospho', '@' = 'Oxidation')) {
  aa =  stringr::str_extract_all(modifiedSequence, '[A-Z][^A-Z]?')[[1]]
  
  modification = c()
  lapply(1:length(aa), function(i) {
    sym = sub('^[A-Z]', '', aa[i])
    if (nchar(sym) == 0) {
      return()
    }
    name = symbol[sym]
    if (is.na(name)) {
      stop(paste0('unknown modification "', sym, '"'))
    }
    aa[i] <<- substring(aa[i], 1, 1)
    modification <<- append(modification, paste0(aa[i], i, '(', name, ')'))
  })
  
  list(
    strippedSequence = paste0(aa, collapse = ''),
    modification = paste0(modification, collapse = ';')
  )
}