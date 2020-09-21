library(XML)
library(readr)
psm.file = list.files(pattern = '\\.psm\\.csv$')[1]
psm = read_csv(psm.file)

activationTypes = 'HCD'
massTolerance = 0.02
scoreNeutralLoss = FALSE

spectra.file = ''
spectra = NULL
spectrumTitle = NULL

spectraNode = newXMLNode('Spectra', lapply(1:nrow(psm), function(i) {
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
  
  name = as.character(spectra[[spectrumIndex]]$title)
  charge = as.integer(psm$charge[i])
  sequence = gsub('[*@]', '', psm$peptide[i])
  modification = paste0('0.', gsub('[A-Z]', '0', gsub('[A-Z]([*@])', '\\1', psm$peptide[i])), '.0')
  
  peaks = paste0(sapply(1:nrow(spectra[[spectrumIndex]]$ions), function(i) {
    paste0(spectra[[spectrumIndex]]$ions[i, 1:2], collapse = ':')
  }), collapse = ',')
  
  newXMLNode(
    'Spectrum', attrs = list(ID = i, Name = name, PrecursorCharge = charge, ActivationTypes = activationTypes),
    newXMLNode('Peaks', newXMLTextNode(peaks)),
    newXMLNode('IdentifiedPhosphorPeptides', newXMLNode(
      'Peptide', attrs = list(ID = i, Sequence = sequence, ModificationInfo = modification)
    ))
  )
}))

phosphoRSInput = newXMLNode(
  'phosphoRSInput',
  newXMLNode('MassTolerance', attrs = list(Value = massTolerance)),
  newXMLNode('Phosphorylation', attrs = list(Symbol = '*')),
  newXMLNode('ScoreNeutralLoss', attrs = list(Value = tolower(as.character(scoreNeutralLoss)))),
  spectraNode,
  newXMLNode(
    'ModificationInfos', 
    newXMLNode('ModificationInfo', attrs = list(
      Symbol = '*', 
      Value = '1:Phospho:Phospho:79.966331:PhosphoLoss:97.976896:STY'
    )),
    newXMLNode('ModificationInfo', attrs = list(
      Symbol = '@', 
      Value = '2:Oxidation:Oxidation:15.994919:null:0:M'
    ))
  )
)

saveXML(phosphoRSInput, file = 'phosphoRSInput.xml')
