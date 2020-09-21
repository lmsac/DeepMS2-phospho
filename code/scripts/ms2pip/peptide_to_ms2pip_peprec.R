library(readr)

lapply(list.files(pattern = '\\.peptide\\.csv$'), function(f) {
  peptides = read_csv(f)
  
  lapply(2:3, function(charge) {
    peprec = data.frame(
      spec_id = paste0(
        peptides$sequence, '_', charge, '_',
        peptides$modification
      ),
      modifications = sub('\\|$', '', paste(
        sapply(strsplit(peptides$modification, ';'), function(s) {
          paste(sapply(s, function(x) {
            aa = sub('^([A-Z])([0-9]+)\\(([A-Za-z0-9\\-_]+)\\)$', '\\1', x)
            pos = sub('^([A-Z])([0-9]+)\\(([A-Za-z0-9\\-_]+)\\)$', '\\2', x)
            name = sub('^([A-Z])([0-9]+)\\(([A-Za-z0-9\\-_]+)\\)$', '\\3', x)
            if (name == 'Phospho') {
              name = paste0(name, aa)
            }
            paste0(pos, '|', name)
          }), collapse = '|')
        }),
        sapply(strsplit(peptides$sequence, ''), function(x) {
          paste(sapply(which(x == 'C'), function(i) paste0(i, '|Carbamidomethyl')), collapse = '|')
        }),
        sep = '|'
      )),
      peptide = peptides$sequence,
      charge = charge,
      stringsAsFactors = FALSE
    )
    
    write.table(peprec, sub('\\.peptide\\.csv$', paste0('_charge', charge, '.peprec.txt'), f), row.names = FALSE, quote = FALSE)
  })
})

