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

x = 1:10
file="foo.txt"

sink_output = function(x, file) {
  f = file(file, open="a")

  write(paste("File:", file), file=f)
  write(paste("Time Ended:", Sys.time()), file=f)
  write(paste("Output:\n```"), file=f)

  sink(f)
  print(x)
  sink()

  write(paste("```\n\n\n"), file=f)
  close(f)

  invisible()
  }
