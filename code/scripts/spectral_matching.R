psmFile = list.files(pattern = 'psm.csv')[1] #'{PATH_PSM_FILE}'
referencePath = 'msp' # {PATH_TO_MSP}'
spectraPath = 'mgf' # {PATH_TO_MGF}'
resultFile = 'spectral_matching.result.csv'

fragmentTolerance = 50


library(readr)

psm = read_csv(psmFile)
psm = subset(psm, !grepl('\\*.+\\*', peptide))
psm = subset(psm, charge == 2 | charge == 3)
psm = subset(psm, grepl('[STY].*[STY]', peptide))
psm = psm[order(psm$file), ]

wd = getwd()

setwd(referencePath)
reference.spectra = do.call(c, lapply(list.files(pattern = '\\.msp$', recursive = T), function(msp) {
  message(msp)
  read.msp(msp)
}))

reference.spectra = reference.spectra[!sapply(reference.spectra, function(x) is.na(max(x$ions[, 2])))]

precursorCharge = sapply(reference.spectra, function(x) as.integer(x$charge))
precursorMZ = sapply(reference.spectra, function(x) as.numeric(x$comment$Parent))
strippedSequence = sapply(reference.spectra, function(x) x$peptide)

modifiedSequence = sapply(reference.spectra, function(x) {
  if (length(x$mods$pos) == 0) {
    x$peptide
  }
  else {
    aa = strsplit(x$peptide, '')[[1]]
    sapply(1:length(x$mods$pos), function(i) {
      pos = x$mods$pos[i]
      name = x$mods$name[i]
      
      if (name == 'Phospho' || name %in% c('PhosphoS', 'PhosphoT', 'PhosphoY')) {
        aa[pos] <<- paste0(aa[pos], '*')
      }
      else if (name == 'Oxidation') {
        aa[pos] <<- paste0(aa[pos], '@')
      }
    })
    paste0(aa, collapse = '')
  }
})
modificationCount = list()
invisible(lapply(c('*', '@'), function(symbol) {
  modificationCount[[symbol]] <<- sapply(strsplit(modifiedSequence, ''), function(s) {
    sum(s == symbol)
  })
}))



setwd(spectraPath)
spectra.file = ''
spectra = NULL
spectrumTitle = NULL
result = do.call(rbind, lapply(1:nrow(psm), function(i) {
  file = paste0(psm$file[i], '.mgf')
  if (file != spectra.file || is.null(spectra)) {
    spectra.file <<- file
    message(paste0(Sys.time(), ': loading ', spectra.file))
    spectra <<- read.mgf(spectra.file)
    message(paste0(Sys.time(), ': ', spectra.file, ' loaded'))
    spectrumTitle <<- sapply(spectra, function(x) sub('^(.*\\.[0-9]+\\.[0-9]+)\\.[0-9]* .*', '\\1', x$title))
  }
  
  queryTitle = sprintf('%s.%d.%d', psm$file[i], psm$scan[i], psm$scan[i])
  spectrumIndex = match(queryTitle, spectrumTitle)
  if (is.na(spectrumIndex)) {
    message(paste0(Sys.time(), ': ', queryTitle, ' spectrum not found '))
    return(NULL)
  }
  
  charge = psm$charge[i]
  querySequence = gsub('[*@]', '', psm$peptide[i])
  queryModificationCount = sapply(c('*', '@'), function(symbol) {
    sum(strsplit(psm$peptide[i], '')[[1]] == symbol)
  })
  
  referenceIndexes = which(
    precursorCharge == charge &
      strippedSequence == querySequence &
      do.call('&', lapply(c('*', '@'), function(symbol) {
        modificationCount[[symbol]] == queryModificationCount[symbol]
      }))
  )
  if (length(referenceIndexes) == 0) {
    # message(paste0(
    #   Sys.time(), ': ', queryTitle, ' peptide not found in spectra library ', querySequence, '|',
    #   paste0(paste0(names(queryModificationCount), queryModificationCount), collapse = ''), '|', charge
    # ))
    return(NULL)
  }
  
  actualMZ = spectra[[spectrumIndex]]$pepMass[1]
  actualRT = spectra[[spectrumIndex]]$rtInSeconds / 60
  querySpectrum = spectra[[spectrumIndex]]$ions
  
  peakTable = align.peaklists(
    c(list(querySpectrum), lapply(referenceIndexes, function(j) {
      reference.spectra[[j]]$ions
    })),
    tolerance = fragmentTolerance
  )
  peakIntensityTable = subset(
    peakTable$intensity,
    !is.na(rowMeans(
      subset(peakTable$mz, select = -1),
      na.rm = TRUE
    ))
  )
  # peakIntensityTable = peakTable$intensity
  
  scores = sapply(1:length(referenceIndexes), function(j) {
    dot.product(peakIntensityTable[, 1], peakIntensityTable[, j + 1])
  }) 
  
  # scores = sapply(referenceIndexes, function(j) {
  #   referenceSpectrum = reference.spectra[[j]]$ions
  # 
  #   peakTable = assign.peaks(
  #     querySpectrum, referenceSpectrum[, 1],
  #     tolerance = fragmentTolerance, tolerance.unit = 'ppm',
  #     assignment = 'mostintense'
  #   )
  # 
  #   dot.product(referenceSpectrum[, 2], peakTable[, 3])
  # })
  
  rank = order(scores, decreasing = TRUE)
  scores = scores[rank]
  deltaScores = scores - c(scores[-1], NA)
  indexes = referenceIndexes[rank]
  
  data.frame(
    index = rep(i, length(rank)),
    spectrum = rep(queryTitle, length(rank)),
    charge = rep(charge, length(rank)),
    actualRT = rep(actualRT, length(rank)),
    actualMZ = rep(actualMZ, length(rank)),
    rank = 1:length(rank),
    peptide = sapply(indexes, function(j) modifiedSequence[j]),
    score = scores,
    deltaScore = deltaScores,
    theoreticalMZ = precursorMZ[indexes],
    stringsAsFactors = FALSE
  )
}))


setwd(wd)
rm(wd)
write.csv(result, resultFile, row.names = FALSE)
message(paste0(Sys.time(), 'result saved'))
