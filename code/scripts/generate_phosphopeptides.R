library(readr)
psm.file = list.files(pattern = '\\.psm\\.csv$')[1]
psm = read_csv(psm.file)
sequences = unique(gsub('[^A-Z]', '', psm$peptide))

sites = c('S', 'T', 'Y')

lapply(sites, function(site) {
  peptides = do.call(rbind, lapply(sequences, function(sequence) {
    modifications = get.candidate.sites(sequence, 'Phospho', site, max.count = 1)
    if (length(modifications) == 0) {
      data.frame(
        sequence = character(),
        modification = character(),
        stringsAsFactors = FALSE
      )
    }
    else {
      data.frame(
        sequence = sequence,
        modification = unlist(modifications),
        stringsAsFactors = FALSE
      )
    }
  }))
  
  peptides = do.call(rbind, c(list(peptides), lapply(1:nrow(peptides), function(i) {
    modifications = get.candidate.sites(peptides$sequence[i], 'Oxidation', 'M', max.count = 5)
    if (length(modifications) == 0) {
      data.frame(
        sequence = character(),
        modification = character(),
        stringsAsFactors = FALSE
      )
    }
    else {
      data.frame(
        sequence = peptides$sequence[i],
        modification = paste0(peptides$modification[i], ';', unlist(modifications)),
        stringsAsFactors = FALSE
      )
    }
  })))
  
  write.csv(
    peptides, 
    paste0(sub('\\.psm\\.csv$', '', psm.file), '_', site, 'phospho1.peptide.csv'), 
    row.names = FALSE
  )
})
