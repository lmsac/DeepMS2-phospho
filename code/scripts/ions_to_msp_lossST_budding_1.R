lossRatio = 50
budding = 100

lapply(list.files(pattern = '_charge[0-9]+.+\\.ions\\.json$'), function(file) {
  ions = rjson::fromJSON(file = file)
  
  spectra = ions.to.spectra.loss.budding.1(
    ions, 
    charge = as.integer(sub('charge', '', stringr::str_extract(file, 'charge[0-9]+'))), 
    loss.ratio = lossRatio * 1e-2,
	loss.site = c('S', 'T'),
    budding = budding * 1e-3, 
    budding.type = c('b', 'y', 'bp', 'yp'), 
    budding.charge = c(1, 2)
  )
  write.msp(spectra, sub('\\.ions\\.json$', paste0(
    '_lossST', lossRatio, '_bud', budding, 'byp2.msp'
  ), file))  
})
