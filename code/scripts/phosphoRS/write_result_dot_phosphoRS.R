library(readr)

psm = read_csv(list.files(pattern = 'psm\\.csv$')[1])

phosphoRS_result = read_csv(list.files(pattern = '\\.phosphoRS\\.csv$')[1], trim_ws = TRUE)

dot_result = read_csv(list.files(pattern = '\\.result\\.csv$')[1], trim_ws = TRUE)


phosphoRS_result = local({
  result = phosphoRS_result
  result$Spectrum.Name = sub(' File:.*$', '', result$Spectrum.Name)
  result
})

phosphoRS_result = local({
  result = data.frame(
    file = sub('\\.[0-9]+\\.[0-9]+\\.[0-9]*$', '', phosphoRS_result$Spectrum.Name),
    scan = as.integer(sub('^.*\\.[0-9]+\\.([0-9]+)\\.[0-9]*$', '\\1', phosphoRS_result$Spectrum.Name)),
    charge = phosphoRS_result$Spectrum.PrecursorCharge,
    peptide = sapply(1:nrow(phosphoRS_result), function(i) {
      if (!is.na(i) && 
          !is.na(phosphoRS_result$Peptide.Sequence[i]) &&
          !is.na(phosphoRS_result$Isoform.Sites[i])) {
        residues = strsplit(phosphoRS_result$Peptide.Sequence[i], '')[[1]]
        sites = phosphoRS_result$Isoform.Sites[i]
        lapply(strsplit(sites, ' ')[[1]], function(s) {
          index = as.integer(sub('^[A-Z]', '', s))
          residues[index] <<- paste0(residues[index], '*')
        })
        paste0(residues, collapse = '')
      }
      else {
        NA
      }
    }),
    phosphoRS_score = phosphoRS_result$Isoform.Score,
    phosphoRS_probability = phosphoRS_result$Isoform.Probability,
    stringsAsFactors = FALSE
  )
  result = subset(result, !is.na(peptide))
  result
})


dot_result = local({
  result = data.frame(
    file = sub('\\.[0-9]+\\.[0-9]+$', '', dot_result$spectrum),
    scan = as.integer(sub('^.*\\.[0-9]+\\.([0-9]+)$', '\\1', dot_result$spectrum)),
    charge = dot_result$charge,
    peptide = gsub('@', '', dot_result$peptide),
    modified_peptide = dot_result$peptide,
    dot_score = dot_result$score,
    delta_score = dot_result$deltaScore,
    stringsAsFactors = FALSE
  )
  result
})

dot_phosphoRS_result = local({
  result = merge(
    dot_result, phosphoRS_result,
    by = c('file', 'scan', 'charge', 'peptide'),
    suffixes = c('_dot', '_phosphoRS'),
    all = FALSE
  )
})

dot_phosphoRS_result = local({
  result = dot_phosphoRS_result
  
  result = result[order(
    result$file, 
    result$scan
  ), ]
  
  result$relative_dot_score = NA
  rowIndexes = which(!duplicated(paste(result$file, result$scan)))
  rowIndexes = cbind(rowIndexes, c(rowIndexes[-1] - 1, nrow(result)))
  lapply(1:nrow(rowIndexes), function(i) {
    result$relative_dot_score[rowIndexes[i, 1]:rowIndexes[i, 2]] <<- result$dot_score[rowIndexes[i, 1]:rowIndexes[i, 2]] / max(result$dot_score[rowIndexes[i, 1]:rowIndexes[i, 2]])
  })
  
  result$rank = NA
  result$score = result$relative_dot_score * result$phosphoRS_score
  result$probability = 10 ^ (result$score / 10)
  
  result = result[order(
    result$file, 
    result$scan, 
    -result$probability
  ), ]
  
  result$ambiguous = FALSE
  rowIndexes = which(!duplicated(paste(result$file, result$scan)))
  rowIndexes = cbind(rowIndexes, c(rowIndexes[-1] - 1, nrow(result)))
  lapply(1:nrow(rowIndexes), function(i) {
    probability = result$probability[rowIndexes[i, 1]:rowIndexes[i, 2]] / sum(result$probability[rowIndexes[i, 1]:rowIndexes[i, 2]])
    result$probability[rowIndexes[i, 1]:rowIndexes[i, 2]] <<- probability
    
    result$rank[rowIndexes[i, 1]:rowIndexes[i, 2]] <<- 
      1:length(rowIndexes[i, 1]:rowIndexes[i, 2])
    
    result$ambiguous[rowIndexes[i, 1]:rowIndexes[i, 2]] <<- probability == c(probability[-1], 0)
  })
  
  result
})

dot_phosphoRS_result = local({
  result = dot_phosphoRS_result
  result = subset(result, rank == 1)
  result
})

dot_phosphoRS_result = local({
  index = match(
    sprintf('%s.%d.%d', psm$file, psm$scan, psm$scan),
    sprintf('%s.%d.%d', dot_phosphoRS_result$file, dot_phosphoRS_result$scan, dot_phosphoRS_result$scan)
  )
  result = data.frame(
    file = psm$file,
    scan = psm$scan,
    charge = psm$charge,
    original = psm$peptide,
    peptide = dot_phosphoRS_result$modified_peptide[index],
    score = dot_phosphoRS_result$score[index],
    probability = dot_phosphoRS_result$probability[index],
    dot_score = dot_phosphoRS_result$dot_score[index],
    phosphoRS_score = dot_phosphoRS_result$phosphoRS_score[index],
    ambiguous = dot_phosphoRS_result$ambiguous[index],
    stringsAsFactors = FALSE
  )
  result$site = sapply(result$peptide, function(x) {
    if (is.na(x)) {
      NA
    }
    else {
      paste0(sub('(Phospho)', '', stringr::str_extract_all(
        get.sequence.modification(x)$modification, 
        '([STY][0-9]+)\\(Phospho\\)'
      )[[1]], fixed = TRUE), collapse = ';')
    }
  })
  result$site[result$ambiguous] = 'ambiguous'
  result$ambiguous = NULL
  result
})

write.csv(dot_phosphoRS_result, 'dot_phosphoRS_result.csv', row.names = FALSE)
