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

  write(paste("File:", file), file=f,append=TRUE)
  write(paste("Time Ended:", Sys.time()), file=f,append=TRUE)
  write(paste("Output:\n\n```"), file=f,append=TRUE)

  sink(f, append=T)
  print(x)
  sink()

  write(paste("```"), file=f,append=TRUE)
  close(f)

  invisible()
  }
