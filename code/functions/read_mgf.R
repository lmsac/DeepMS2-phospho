read.mgf = function(file) {
  lines = readLines(con = file)
  
  begin.indexs = grep(pattern = 'BEGIN IONS', x = lines)
  end.indexs = grep(pattern = 'END IONS', x = lines)
  
  lapply(1:length(begin.indexs), function(i) {
    internal.lines = lines[begin.indexs[i]:end.indexs[i]]
    
    fields.indexs = grep(pattern = '=', internal.lines)
    fields = strsplit(internal.lines[fields.indexs], '=')
    fields.name = sapply(fields, function(s) s[1])
    fields = substring(internal.lines[fields.indexs], nchar(fields.name) + 2)
    names(fields) = fields.name
    
    title = fields['TITLE']
    charge = fields['CHARGE']
    pepMass = as.double(strsplit(fields['PEPMASS'], ' ')[[1]])
    rtInSeconds = as.double(strsplit(fields['RTINSECONDS'], ' ')[[1]])
    
    ion.indexs = -c(1, length(internal.lines), grep(pattern = '=', internal.lines))
    records = strsplit(internal.lines[ion.indexs], split = ' |\t')
    if(length(records) > 0)
      ions = do.call(rbind, lapply(records, function(record) {
        as.numeric(record[1:2])
      }))
    else
      ions = matrix(ncol = 2, nrow = 0)
    colnames(ions) = c('mz', 'intensity')
    
    list(ions = ions, title = title, charge = charge, pepMass = pepMass, rtInSeconds = rtInSeconds, fields = fields)
  })
}