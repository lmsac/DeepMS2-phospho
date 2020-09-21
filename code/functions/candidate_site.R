get.candidate.sites = function(sequence, mod.name = 'Oxidation', aa = 'M', max.count = 3) {
  positions = stringr::str_locate_all(sequence, paste0('[', paste0(aa, collapse = ''), ']'))[[1]]
  if (nrow(positions) == 0) {
    return(NULL)
  }
  
  lapply(1:max.count, function(modCount) {
    if (nrow(positions) < modCount) {
      return(NULL)
    }
    indexes = combn(nrow(positions), modCount)
    modifications = apply(indexes, 2, function(idx) {
      paste0(sapply(idx, function(j) {
        paste0(substring(sequence, positions[j, 1], positions[j, 1]), positions[j, 1], '(', mod.name,')')
      }), collapse = ';')
    })
  })
}