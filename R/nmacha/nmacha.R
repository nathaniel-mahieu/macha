
nmacha = function(...) {
  Nmacha = list()
  Nmacha$m = list(...)

  Nmacha$m.c = lapply(seq_along(Nmacha$m), function(j) {
    Nmacha$m[[j]]$c[,m:=j]
  }) %>% do.call(what=rbind)

  Nmacha$m.r = lapply(seq_along(Nmacha$m), function(j) {
    Nmacha$m[[j]]$r[,m:=j]
  }) %>% do.call(what=rbind)

  Nmacha
  }

