lossRatio = 50
budding = 100

lapply(list.files(pattern = '_charge[0-9]+.+\\.ions\\.json$'), function(file) {
  ions = rjson::fromJSON(file = file)
  if (grepl('[ST]phospho1', file)) {
    spectra = ions.to.spectra.loss.budding(
      ions, 
      charge = as.integer(sub('charge', '', stringr::str_extract(file, 'charge[0-9]+'))), 
      lossRatio = lossRatio * 1e-2,
      budding = budding * 1e-3, 
      budding.type = c('b', 'y', 'bp', 'yp'), 
      budding.charge = c(1, 2)
    )
    write.msp(spectra, sub('\\.ions\\.json$', paste0(
      '_loss', lossRatio, '_bud', budding, 'byp2.msp'
    ), file))
  }
  else {
    spectra = ions.to.spectra.budding(
      ions, 
      charge = as.integer(sub('charge', '', stringr::str_extract(file, 'charge[0-9]+'))),       
      budding = budding * 1e-3, 
      budding.type = c('b', 'y'), 
      budding.charge = c(1, 2)
    )
    write.msp(spectra, sub('\\.ions\\.json$', paste0(
      '_bud', budding, 'by2.msp'
    ), file))
  }
})
