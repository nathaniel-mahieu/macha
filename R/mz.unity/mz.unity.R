
breakzgroup = function(r.n, z, rz) {
  if (length(r.n) < 2) return(rep(1, length(r.n)))

  df = which(arr.ind = T, outer(z, z, "-") == outer(r.n, r.n, "-") * rz) %>% as.data.frame
  mem = graph.data.frame(df) %>% clusters %>% '[['("membership")

  v = vector(length = length(r.n))
  v[as.numeric(names(mem))] = mem
  v[is.na(v)] = seq(from = max(v,na.rm=T), to = sum(is.na(v)) +max(v,na.rm=T))

  v
}

breakigroup = function(r.n, i, abmin, abmax) {
  if (length(r.n) < 2) return(rep(1, length(r.n)))


  w = outer(i, i, "/")  %>% { . < abmax & . > abmin & outer(r.n, r.n, "-") == 1 }

  df = which(arr.ind = T, w)
  mem = graph.data.frame(df) %>% clusters %>% '[['("membership")

  v = vector(length = length(r.n))
  v[as.numeric(names(mem))] = mem
  v[is.na(v)] = seq(from = max(v,na.rm=T), to = sum(is.na(v)) +max(v,na.rm=T))

  v
}

unityme = function(Nmacha, rules, zs = 1:6, ppmwid = 5, rtwid = 0.75, factor = 1) {
  library(mz.unity)

  scales = c(ppmwid * 700 / 1E6, rtwid)


  MZ = lapply(zs, function(.z) copy(Nmacha$cc[,z:=.z])) %>% do.call(what=rbind)
  MZ[, m :=mz*abs(z)][,mz:=NULL]

  # For each rule
  knots = list()
  for (.r in seq_len(nrow(rules))) {
    cat(sep="","\rWorking on rule \"", rules[.r, name.string], "\" (", round(.r/nrow(rules)*100,1), "%)        ")

    rulea = copy(rules[.r])

    MZ[, r := m %% rulea$m]
    MZ[, r.n := m %/% rulea$m]

    # Iterative splitting
    gcols = paste0("g", 1:4)
    MZ[,(gcols) := "1"]

    ngs = 1
    repeat {
      MZ[,g2 := paste(g2, breakgroup(r.n, 1)), by=gcols]
      MZ[,g3 := paste(g3, breakzgroup(r.n, z, rulea$z)), by=gcols]
      MZ[,g4 := paste(g4, breakigroup(r.n, i, rulea$abmin, rulea$abmax)), by=gcols]
      MZ[,g1 := paste(g1, clustgroup(cbind(r, rtpeak), scales, factor)), by=gcols]

      l = length(unique(MZ[,do.call(paste, .SD), .SDcols = gcols])); if (l == ngs) break else ngs = l
    }
    MZ[,rg := as.integer(factor(paste(do.call(paste, .SD)))), .SDcols = gcols][,(gcols) :=NULL]

    numg = MZ[,.(.N, nu = length(unique(r.n))),by="rg"]

    knots = c(knots, list(copy(MZ[numg[N>1 & nu > 1],,on="rg"][,.(cc, z, rg, rule = rep(rulea$rule, length(cc)))])))
  }
  Knots = do.call(what=rbind, knots)[,rg := as.numeric(factor(paste(rg, rule)))]



  ### Determine Charge
  uzs = unique(Knots$z)

  zcounts = matrix(nrow = max(Knots$cc), ncol = length(uzs))
  levels = unique(Knots$rule)

  Knots[,ninrg := .N, by="rg"]
  counts = Knots[, sum(ninrg), by=c("cc","z")]
  zcounts[counts[,.(cc,as.numeric(factor(z, levels = uzs)))] %>% as.matrix] = counts$V1
  zcounts = zcounts / rep(seq_len(ncol(zcounts)), each = nrow(zcounts))

  w = which(zcounts == rowMaxs(zcounts,na.rm=T),arr.ind=T)
  zdef = vector(length = max(Knots$cc))
  zdef[w[,1]] = uzs[w[,2]]

  Knots = Knots[zdef[cc] == z][,n:=.N,by="rg"][n>1][,n:=NULL]

  Nmacha$unity.rules = rules
  Nmacha$unity = Knots[,.(cc, z, rg, rule)]

  Nmacha
}
