library(readr)

peptide.files = list.files(pattern = '[STY]phospho1\\.peptide\\.csv$')
peptide.lists = lapply(peptide.files, read_csv)
names(peptide.lists) = peptide.files


lapply(list.files(pattern = '\\.prosit\\.csv$'), function(assay.file) {
  assay = read_csv(assay.file)
  
  print(table(assay$FragmentType))
  print(table(assay$FragmentCharge))
  print(table(assay$FragmentLossType))
  
  
  precursor.indexes = which(sapply(1:nrow(assay), function(i) {
    if (i == 1) {
      return(TRUE)
    }
    
    if (assay$StrippedPeptide[i] != assay$StrippedPeptide[i - 1] ||
        assay$ModifiedPeptide[i] != assay$ModifiedPeptide[i - 1] ||
        assay$LabeledPeptide[i] != assay$LabeledPeptide[i - 1] ||
        assay$PrecursorCharge[i] != assay$PrecursorCharge[i - 1]) {
      TRUE
    }
    else {
      FALSE
    }
  }))
  
  ions = lapply(1:length(precursor.indexes), function(i) {
    if (i < length(precursor.indexes)) {
      rowIndexes = precursor.indexes[i]:precursor.indexes[i + 1]
    }
    else {
      rowIndexes = precursor.indexes[i]:nrow(assay)
    }
    
    sequence = assay$StrippedPeptide[rowIndexes[1]]
    charge = assay$PrecursorCharge[rowIndexes[1]]
    
    ions = list(
      b1 = rep(0, nchar(sequence) - 1),
      #bn1 = rep(0, nchar(sequence) - 1),
      #bo1 = rep(0, nchar(sequence) - 1),
      b2 = rep(0, nchar(sequence) - 1),
      #bn2 = rep(0, nchar(sequence) - 1),
      #bo2 = rep(0, nchar(sequence) - 1),
      y1 = rep(0, nchar(sequence) - 1),
      #yn1 = rep(0, nchar(sequence) - 1),
      #yo1 = rep(0, nchar(sequence) - 1),
      y2 = rep(0, nchar(sequence) - 1)#,
      #yn2 = rep(0, nchar(sequence) - 1),
      #yo2 = rep(0, nchar(sequence) - 1)
    )
    
    lapply(rowIndexes, function(row) {
      type = paste0(
        assay$FragmentType[row],
        ifelse(
          assay$FragmentLossType[row] == 'noloss', '',
          #ifelse(
          #  assay$FragmentLossType[row] == 'NH3', 'n',
          #  ifelse(
          #    assay$FragmentLossType[row] == 'H2O', 'o', 
          #    '?'
          #  )
          #)
          '?'
        ),
        assay$FragmentCharge[row]
      )
      num = as.integer(assay$FragmentNumber[row])
      intensity = as.numeric(assay$RelativeIntensity[row])
      if (type %in% names(ions) && num > 0 && num <= nchar(sequence) - 1 && !is.na(intensity)) {
        ions[[type]][num] <<- intensity
      }
    })
    
    result = list(
      peptide = sequence,
      modifiedPeptide = assay$ModifiedPeptide[rowIndexes[1]],
      charge = charge,
      ions = ions
    )
    result
  })
  
  
  peptide.file.indexes = grep(
    sub('_charge[0-9]+_ce[0-9]+\\.prosit\\.csv$', '', assay.file),
    peptide.files
  )
  lapply(peptide.file.indexes, function(i) {
    peptide.file = peptide.files[i]
    peptides = peptide.lists[[i]]
    
    indexes = local({
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
                  aa[pos] <<- 'M[Oxidation (M)]'
                }
              }
            })
          })
          paste0(aa, collapse = '')
        })
      }
      
      modified_sequence = gsub('C', 'C[Carbamidomethyl (C)]', modified_sequence)
      modified_sequence = paste0('_', modified_sequence, '_')
      
      match(modified_sequence, sapply(ions, function(x) x$modifiedPeptide))
    })
    
    ions = lapply(1:length(indexes), function(i) {
      if (!is.na(indexes[i])) {
        x = ions[[indexes[i]]]
        x$modification = peptides$modification[i]
        x$modifiedPeptide = NULL
        x
      }
      else {
        NULL
      }
    })
    
    ions = ions[!sapply(ions, is.null)]
    
    writeLines(
      rjson::toJSON(ions, indent = 2), 
      paste0(
        sub('\\.peptide\\.csv$', '', peptide.file),
        stringr::str_extract(assay.file, '_charge[0-9]+_ce[0-9]+'),
        '.prosit.ions.json'
      )
    )
  })
})
