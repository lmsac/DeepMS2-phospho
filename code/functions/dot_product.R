dot.product = function(x, y) {
  x = sqrt(x)
  y = sqrt(y)
  prod = sum(x * y)
  if (prod == 0)
    0
  else
    prod / sqrt(sum(x * x) * sum(y * y))
}