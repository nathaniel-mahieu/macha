localMinima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}


localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-Inf, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

lassign = function(..., pos = 1, envir = as.environment(pos)) {
  as = list(...)

  for (i in seq_along(as)) {
    assign(names(as[i]), as[[i]], envir = envir)
  }
}
