library(readr)

ce = 34
charges = 2:3

lapply(list.files(pattern = '\\.peptide\\.csv$'), function(file) {
  peptides = read_csv(file)
  
  lapply(charges, function(charge) {
    modified_sequence = peptides$sequence
    ox = grep('M[0-9]+\\(Oxidation\\)', peptides$modification)
    if (length(ox) > 0) {
      modified_sequence[ox] = sapply(ox, function(i) {
        aa = strsplit(modified_sequence[i], '')[[1]]
        sapply(strsplit(peptides$modification[i], ';'), function(mods) {
          sapply(mods, function(s) {
            if (grepl('^M[0-9]+\\(Oxidation\\)$', s)) {
              pos = as.integer(sub('^M([0-9]+)\\(Oxidation\\)$', '\\1', s))
              if (aa[pos] == 'M') {
                aa[pos] <<- 'M(ox)'
              }
            }
          })
        })
        paste0(aa, collapse = '')
      })
    }
    
    peptides = data.frame(
      modified_sequence = modified_sequence,
      collision_energy = ce,
      precursor_charge = charge,
      stringsAsFactors = FALSE
    )
    peptides = unique(peptides)
    peptides = subset(peptides, nchar(modified_sequence) <= 30)
    
    write.csv(
      peptides,
      sub('\\.peptide\\.csv$', paste0('_charge', charge, '_ce', ce, '.prosit_peptide.csv'), file), 
      row.names = FALSE
    )
  })
})