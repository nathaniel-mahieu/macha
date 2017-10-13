
getgroup = function(Nmacha, .g) {
  gl = Nmacha$m.c[Nmacha$m.c_g,,on="m.c."][g==.g]

  roi = lapply(split(gl, gl[,paste(m,r)]), function (g) {
    getroi(Nmacha$m[[g$m[1]]], g$r[1])[,m:=g$m[1]]
  }) %>% do.call(what=rbind)

  list(m.c = gl, roi=roi)
}

