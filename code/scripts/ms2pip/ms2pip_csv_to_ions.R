library(readr)
lapply(list.files(pattern = '\\.ms2pip\\.csv$'), function(ms2pip.file) {
  ms2pip = read_csv(ms2pip.file)
  
  psm.indexes = which(sapply(1:nrow(ms2pip), function(i) {
    if (i == 1) {
      return(TRUE)
    }
    
    if (ms2pip$spec_id[i] != ms2pip$spec_id[i - 1]) {
      TRUE
    }
    else {
      FALSE
    }
  }))
  
  
  extracted.ions = lapply(1:length(psm.indexes), function(i) {
    if (i < length(psm.indexes)) {
      rowIndexes = psm.indexes[i]:psm.indexes[i + 1]
    }
    else {
      rowIndexes = psm.indexes[i]:nrow(ms2pip)
    }
    
    #sequence = sub('_[0-9]+$', '', ms2pip$spec_id[rowIndexes[1]])
    sequence = sub('^([A-Z]+)_[0-9]+.*$', '\\1', ms2pip$spec_id[rowIndexes[1]])
    
    # modification = strsplit(sub('^[A-Z]+_[0-9]+(.*)$', '\\1', ms2pip$spec_id[rowIndexes[1]]), '_')[[1]]
    # modification = modification[nchar(modification) > 0]
    # modification = ifelse(
    #   grepl('^[STY][0-9]+$', modification), paste0(modification, '(Phospho)'),
    #   ifelse(
    #     grepl('^[M][0-9]+$', modification), paste0(modification, '(Oxidation)'),
    #     ifelse(
    #       grepl('^[C][0-9]+$', modification), paste0(modification, '(Carbamidomethyl)'),
    #       modification
    #     )
    #   )
    # )
    # if (length(modification) > 0) {
    #   modification = paste(modification, collapse = ';')
    # }
    modification = sub('^_', '', sub('^[A-Z]+_[0-9]+(.*)$', '\\1', ms2pip$spec_id[rowIndexes[1]]))
    
    charge = ms2pip$charge[rowIndexes[1]]
    
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
      type = ifelse(
        ms2pip$ion[row] == 'B', 'b1',
        ifelse(
          ms2pip$ion[row] == 'Y', 'y1',
          ifelse(
            ms2pip$ion[row] == 'B2', 'b2',
            ifelse(
              ms2pip$ion[row] == 'Y2', 'y2', NA
            )
          )
        )
      )
      num = as.integer(ms2pip$ionnumber[row])
      intensity = as.numeric(ms2pip$prediction[row])
      if (type %in% names(ions) && num > 0 && num <= nchar(sequence) - 1 && !is.na(intensity)) {
        ions[[type]][num] <<- intensity
      }
    })
    
    list(
      specId = ms2pip$spec_id[rowIndexes[1]],
      peptide = sequence,
      modification = modification,
      charge = charge,
      ions = ions
    )
  })
  
  writeLines(
    rjson::toJSON(extracted.ions),
    paste0(sub('\\.ms2pip\\.csv', '.ms2pip.ions.json', ms2pip.file))
  )
})
