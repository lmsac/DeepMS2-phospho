lossRatio = 50

lapply(list.files(pattern = '_charge[0-9]+.+\\.ions\\.json$'), function(file) {
  ions = rjson::fromJSON(file = file)
  if (grepl('[ST]phospho1', file)) {
    spectra = ions.to.spectra.loss(
      ions, 
      charge = as.integer(sub('charge', '', stringr::str_extract(file, 'charge[0-9]+'))), 
      lossRatio = lossRatio * 1e-2
    )
    write.msp(spectra, sub('\\.ions\\.json$', paste0('_loss', lossRatio, '.msp'), file))
  }
  else {
    spectra = ions.to.spectra(
      ions, 
      charge = as.integer(sub('charge', '', stringr::str_extract(file, 'charge[0-9]+')))
    )
    write.msp(spectra, sub('\\.ions\\.json$', '.msp', file))
  }
})
