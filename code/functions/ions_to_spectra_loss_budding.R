ions.to.spectra.loss.budding = function(ions, carbamidomethyl = TRUE, charge = NA, lossRatio = 0.9, budding = 5e-3, budding.type = c('b', 'y', 'bp', 'yp'), budding.charge = c(1, 2), ...) {
  lapply(ions, function(x) {
    if (is.null(x$charge)) {
      x$charge = charge
    }
    
    mz = fragment.ions.mz.1(
      x$peptide,
      modification = x$modification,
      type = c('b', 'y'),
      charge = c(1, 2), 
      loss = c('NH3', 'H2O', 'H3PO4'), 
      carbamidomethyl = carbamidomethyl, 
      ...
    )
    
    intensity = sapply(names(mz), function(label) {
      type = sub('\\*', 'n', sub('[0-9]+\\^', '', label))
      if (!is.null(x$ions[[type]])) {
        number = as.integer(sub('[a-z*]+([0-9]+)\\^[0-9]+', '\\1', label))
        x$ions[[type]][number]
      }
      else {
        0
      }
    })
    intensity[is.na(intensity)] = 0
    names(intensity) = names(mz)
    
    sapply(names(intensity), function(label) {
      if (grepl('^[by]p[0-9]+\\^[0-9]+$', label)) {
        noloss = sub('^([by])p', '\\1', label)
        if (noloss %in% names(intensity)) {
          intensitySum = intensity[label] + intensity[noloss]
          intensity[label] <<- intensitySum * lossRatio
          intensity[noloss] <<- intensitySum * (1 - lossRatio)
        }
      }
    })
    
    buddingFragments = grepl(
      paste0(
        '^(', paste0(budding.type, collapse = '|'), ')[0-9]+\\^(',
        paste0(budding.charge, collapse = '|'), ')$'
      ), 
      names(mz)
    )
    
    intensity = ifelse(
      intensity > 0 | buddingFragments, 
      intensity * (1 - budding)  + max(intensity) * budding,
      0
    )
    
    intensity = intensity[order(mz)]
    mz = mz[order(mz)]
    ions = cbind(
      mz = mz[intensity > 0],
      intensity = intensity[intensity > 0]
    )
    
    ion.labels = sapply(rownames(ions), function(label) {
      type = substring(label, 1, 1)
      number = as.integer(sub('[a-z*]+([0-9]+)\\^[0-9]+', '\\1', label))
      loss = local({
        s = substring(label, 2, 2)
        if (s == '*')
          '-17'
        else if (s == 'o')
          '-18'
        else if (s == 'p')
          '-98'
        else
          ''
      })
      charge = as.integer(sub('[a-z*]+[0-9]+\\^([0-9]+)', '\\1', label))
      paste0(type, number, loss, '^', charge)
    })
    rownames(ions) = NULL
    
    pepMass = peptide.mass.1(
      x$peptide, charge = x$charge, 
      modification = x$modification, 
      carbamidomethyl = carbamidomethyl,
      ...
    )
    
    if (!is.null(x$modification) && !is.na(x$modification) && nchar(x$modification) > 0) {
      s = strsplit(x$modification, ';')[[1]]
      name = gsub('[()]', '', stringr::str_extract(s, '\\(.+\\)$'))
      pos = as.integer(gsub('[A-Z(]', '', stringr::str_extract(s, '^[A-Z]+[0-9]+\\(')))
      mods = list(pos = pos, name = name)
    }
    else {
      mods = NULL
    }
    if (carbamidomethyl) {
      cysteine = stringr::str_locate_all(x$peptide, '[Cc]')[[1]][, 1]
      cysteine = setdiff(cysteine, mods$pos[mods$name == 'Carbamidomethyl'])
      mods$pos = c(mods$pos, cysteine)
      mods$name = c(mods$name, rep('Carbamidomethyl', length(cysteine)))
    }
    
    list(
      peptide = x$peptide,
      charge = x$charge,
      mods = mods,
      pepMass = as.numeric(pepMass),
      ions = ions,
      ion.labels = ion.labels
    )
  })
}
