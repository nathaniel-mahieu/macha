collect_errors = function(l, names = NULL, .rbind = T) {

  if (!is.null(names)) names(l) = names

  cat(names(l))

  ise = sapply(l, function(x) {
    inherits(x, "simpleError") | inherits(x, "try-error")
  })

  list = l[!ise]
  if (.rbind) list = do.call(rbind, list)


  list(
    list = list,
    errors = l[ise]
  )

}

