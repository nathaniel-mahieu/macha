collect_errors = function(l, names = NULL, .rbind = T) {

  if (!is.null(names)) names(l) = names

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


sink_output = function(x, file) {
  f = file(file)

  write(paste("File:", file),file=file,append=TRUE)
  write(paste("Time Ended:", Sys.time()),file=file,append=TRUE)
  write(paste("Output:\n\n```"),file=file,append=TRUE)

  sink(file, append=T)
  x
  sink()

  write(paste("```"),file=file,append=TRUE)
  close(f)

  invisible()
  }
