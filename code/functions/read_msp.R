read.msp = function(file) {
  lines = readLines(file)
  begin.indexs = grep(pattern = '^Name: ', x = lines)
  end.indexs = c(begin.indexs[-1] - 1, length(lines))
  
  lapply(1:length(begin.indexs), function(i) {
    internal.lines = lines[begin.indexs[i]:end.indexs[i]]
    internal.lines = internal.lines[nchar(internal.lines) > 0]
    
    name.index = grep(pattern = '^Name: ', internal.lines)
    name = substring(internal.lines[name.index], first = 7)
    peptide = strsplit(name, '/')[[1]][1]
    charge = as.numeric(strsplit(name, '/')[[1]][2])
    
    mw.index = grep(pattern = '^MW: ', internal.lines)
    mw = as.numeric(substring(internal.lines[mw.index], first = 4))
    
    comment.index = grep(pattern = '^Comment: ', internal.lines)
    comment = local({
      line = substring(internal.lines[comment.index], first = 10)
      space.indexs = gregexpr(pattern = ' ', line)[[1]]
      quote.indexs = gregexpr(pattern = '"', line)[[1]]
      if (length(quote.indexs) == 1 && quote.indexs == -1) {
        split.indexs = space.indexs
      }
      else {
        quote.start.indexs = quote.indexs[1:length(quote.indexs) %% 2 == 1]
        quote.end.indexs = quote.indexs[1:length(quote.indexs) %% 2 == 0]
        if (length(quote.start.indexs) > length(quote.end.indexs))
          quote.end.indexs = c(quote.end.indexs, nchar(line) + 1)
        split.indexs = space.indexs[sapply(space.indexs, function(space.index) {
          !any(space.index > quote.start.indexs &
                 space.index < quote.end.indexs)
        })]
      }
      start.indexs = c(1, split.indexs + 1)
      end.indexs = c(split.indexs - 1, nchar(line))
      
      fields = strsplit(substring(line, start.indexs, end.indexs), '=')
      result = lapply(fields, function(field)
        field[2])
      names(result) = lapply(fields, function(field)
        field[1])
      result
    })
    
    num.index = grep(pattern = '^Num peaks: ', internal.lines)
    num.peaks = as.numeric(substring(internal.lines[num.index], first = 12))
    
    ion.indexs = (length(internal.lines) - num.peaks + 1):length(internal.lines)
    records = strsplit(internal.lines[ion.indexs], split = '\t')
    ions = do.call(rbind, lapply(records, function(record) {
      as.numeric(record[1:2])
    }))
    colnames(ions) = c('mz', 'intensity')
    ion.labels = sapply(records, function(record) {
      record[3]
    })
    
    mods = local({
      pos.name = strsplit(strsplit(comment$Mods, '/')[[1]][-1], ',')
      pos = sapply(pos.name, function(x) {
        as.numeric(x[1]) + 1
      })
      name = sapply(pos.name, function(x) {
        x[3]
      })
      list(pos = pos, name = name)
    })
    
    list(
      peptide = peptide,
      charge = charge,
      mods = mods,
      ions = ions,
      ion.labels = ion.labels,
      mw = mw,
      comment = comment
    )
  })
}
