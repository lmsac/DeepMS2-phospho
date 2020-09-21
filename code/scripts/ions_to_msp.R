lapply(list.files(pattern = '_charge[0-9]+.+\\.ions\\.json$'), function(file) {
  ions = rjson::fromJSON(file = file)
  spectra = ions.to.spectra(ions, charge = as.integer(sub('charge', '', stringr::str_extract(file, 'charge[0-9]+'))))
  write.msp(spectra, sub('\\.ions\\.json$', '.msp', file))
})
