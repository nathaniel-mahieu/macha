fastrangedt = function(DT, range, key) {
  if (key(DT) != key) stop("data.table key not set. setkey() for fast subsetting.")
  if (class(DT[[key]]) != class(range)) stop("Class mismatch between range and key column.")
  if (diff(range) <= 0) return(DT[numeric()])

  range[2] = range[2]+10^(-10) #Small error introduced. Necessary to make upper bound inclusive.

  ind <- DT[.(range), which=TRUE, roll=-Inf, rollends=c(F,F)]

  isna = is.na(ind)
  ind[isna] = c(1, nrow(DT)+1)[isna]

  if (diff(ind) == 0) return(DT[numeric()])
  DT[(ind[1]):(ind[2]-1)]
}
