
nmacha = function(...) {
  Nmacha = list()
  Nmacha$m = list(...)

  Nmacha$m.c = lapply(seq_along(Nmacha$m), function(j) {
    Nmacha$m[[j]]$c[,m:=j]
  }) %>% do.call(what=rbind)
  Nmacha$m.c[,m.c.:=seq_along(c)]

  Nmacha$m.r = lapply(seq_along(Nmacha$m), function(j) {
    Nmacha$m[[j]]$r[,m:=j]
  }) %>% do.call(what=rbind)

  class(Nmacha) <- c("mach", class(Nmacha))

  Nmacha
  }

