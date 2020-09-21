normalize.intensity = function(data, max.intensity = 100) {
  normalized.intensity = data[, 2] / max(data[, 2]) * max.intensity
  normalized.intensity[normalized.intensity < 0] = 0
  new.data = data
  new.data[, 2] = normalized.intensity
  new.data
}