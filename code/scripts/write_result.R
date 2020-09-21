library(readr)

psm = read_csv(list.files(pattern = 'psm\\.csv$')[1])

deltadot_result = read_csv(list.files(pattern = '\\.result\\.csv$')[1], trim_ws = TRUE)

deltadot_result = local({
  result = deltadot_result
  result = subset(result, rank == 1)
  result
})

deltadot_result = local({
  index = match(
    sprintf('%s.%d.%d', psm$file, psm$scan, psm$scan),
    deltadot_result$spectrum
  )
  result = data.frame(
    file = psm$file,
    scan = psm$scan,
    charge = psm$charge,
    original = psm$peptide,
    peptide = deltadot_result$peptide[index],
    score = deltadot_result$score[index],
    delta_score = deltadot_result$deltaScore[index],
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
  result$site[result$delta_score == 0] = 'ambiguous'
  result
})

write.csv(deltadot_result, sub('\\.result\\.csv$', '_result.csv', file), row.names = FALSE)
