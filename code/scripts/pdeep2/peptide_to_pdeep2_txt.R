library(readr)

lapply(list.files(pattern = '\\.peptide\\.csv$'), function(f) {
  peptides = read_csv(f)
  
  lapply(2:3, function(charge) {
    peptidetxt = data.frame(
      peptide = peptides$sequence,
      modification = paste(
        sapply(strsplit(peptides$modification, ';'), function(s) {
          paste(sapply(s, function(x) {
            aa = sub('^([A-Z])([0-9]+)\\(([A-Za-z0-9\\-_]+)\\)$', '\\1', x)
            pos = sub('^([A-Z])([0-9]+)\\(([A-Za-z0-9\\-_]+)\\)$', '\\2', x)
            name = sub('^([A-Z])([0-9]+)\\(([A-Za-z0-9\\-_]+)\\)$', '\\3', x)
            paste0(pos, ',', name, '[', aa, ']')
          }), collapse = ';')
        })
      ),
      charge = charge,
      stringsAsFactors = FALSE
    )
    
    write.table(peptidetxt, sub('\\.peptide\\.csv$', paste0('_charge', charge, '.peptide.txt'), f), sep = '\t', row.names = FALSE, quote = FALSE)
  })
})
