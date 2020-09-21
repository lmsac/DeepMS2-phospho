lapply(list.files(pattern = '\\.pdeep2\\.mgf$'), function(mgf.file) {
  spectra = read.mgf(mgf.file)
  
  fragments = list(
    list('b+1', 'b1', nTerm = TRUE),
    list('b+2', 'b2', nTerm = TRUE),
    list('y+1', 'y1', nTerm = FALSE),
    list('y+2', 'y2', nTerm = FALSE),
    list('b+1-ModLoss', 'bp1', nTerm = TRUE),
    list('b+2-ModLoss', 'bp2', nTerm = TRUE),
    list('y+1-ModLoss', 'yp1', nTerm = FALSE),
    list('y+2-ModLoss', 'yp2', nTerm = FALSE),
    list('b+1-NH3', 'bn1', nTerm = TRUE),
    list('b+2-NH3', 'bn2', nTerm = TRUE),
    list('y+1-NH3', 'yn1', nTerm = FALSE),
    list('y+2-NH3', 'yn2', nTerm = FALSE),
    list('b+1-H2O', 'bo1', nTerm = TRUE),
    list('b+2-H2O', 'bo2', nTerm = TRUE),
    list('y+1-H2O', 'yo1', nTerm = FALSE),
    list('y+2-H2O', 'yo2', nTerm = FALSE)
  )
  
  ions = lapply(spectra, function(x) {
    peptide = sub('|', '', stringr::str_extract(x$title, '^[A-Z]+\\|'), fixed = TRUE)
    charge = as.integer(sub('|', '', stringr::str_extract(x$title, '\\|[0-9]+$'), fixed = TRUE))
    modStr = trimws(gsub('^[A-Z]+\\||\\|[0-9]+$', '', x$title))
    if (nchar(modStr) == 0) {
      modification = NULL
    }
    else {
      modification = paste0(sapply(strsplit(modStr, ';')[[1]], function(m) {
        pos = as.integer(sub(',', '', stringr::str_extract(m, '^[0-9]+,')))
        name = gsub('^[0-9]+,|\\[[A-Z]\\]$', '', m)
        aa = gsub('\\[|\\]', '', stringr::str_extract(m, '\\[[A-Z]\\]$'))
        paste0(aa, pos, '(', name, ')')
      }), collapse = ';')
    }
    
    ions = list()
    lapply(fragments, function(frag) {
      intStr = x$fields[frag[[1]]]
      if (!is.null(intStr) & !is.na(intStr)) {
        intensity = as.numeric(strsplit(intStr, ',')[[1]])
        if (!frag$nTerm) {
          intensity = rev(intensity)
        }
        ions[[frag[[2]]]] <<- intensity
      }
    })
    
    list(
      peptide = peptide,
      modification = modification,
      charge = charge,
      ions = ions
    )
  })
  
  writeLines(
    rjson::toJSON(ions),
    paste0(sub('\\.mgf$', '.ions.json', mgf.file))
  )
})
