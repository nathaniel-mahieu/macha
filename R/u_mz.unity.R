
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
clustgroup = function(mat, scales, factor) {
  if (nrow(mat) < 2) return(rep(1, nrow(mat)))
  
  cutree(hclust(dist(mat/rep(scales, each = nrow(mat))), method="single"), h = factor)
}

breakgroup = function(x, breaksize = 1) {
  o = order(x)
  ass = c(0, which(diff(x[o]) > breaksize), length(x))
  
  as.integer(cut(seq_along(x), breaks = ass))[order(o)]
}

unityme = function(Nmacha, rules, zs = 1:6, ppmwid = 5, rtwid = 0.75, factor = 1) {
  scales = c(ppmwid * 700 / 1E6, rtwid)

  MZ = lapply(zs, function(.z) copy(Nmacha$cc[,z:=.z])) %>% do.call(what=rbind)
  MZ[, m :=mz*abs(z)][,mz:=NULL]

  # For each rule
  knots = foreach (r = seq_len(nrow(rules)), .packages = c("data.table", "macha"), .options.redis=list(chunkSize=50)) %dopar% {
    .r = r; rm(r)
    cat(sep="","\rWorking on rule \"", rules[.r, name.string], "\" (", round(.r/nrow(rules)*100,1), "%)        ")

    rulea = copy(rules[.r])

    MZ[, r := m %% rulea$m]
    MZ[, r.n := m %/% rulea$m]

    # Iterative splitting
    gcols = paste0("g", 1:5)
    MZ[,(gcols) := "1"]
    
    
    g5 = MZ[,dbscan(cbind(r, rt)/rep(scales, each=length(rt)), eps = 1, minPts = 2)$cluster]
    g5[g5 == 0] = seq_len(sum(g5==0)) + max(g5)
    g5x = g5
    MZ[,g5 := as.character(g5x)]

    MZ[,N:=.N,by=gcols]
    checkme = which(MZ$N>1)
    
    ngs = 1
    repeat {
      MZ[checkme,g1 := paste(g1, clustgroup(cbind(r, rt), scales, factor)), by=gcols]
      MZ[,N:=.N,by=gcols]
      checkme = which(MZ$N>1)

      MZ[checkme,g2 := paste(g2, breakgroup(r.n, 1)), by=gcols]
      MZ[checkme,g3 := paste(g3, breakzgroup(r.n, z, rulea$z)), by=gcols]      
      MZ[checkme,g4 := paste(g4, breakigroup(r.n, i, rulea$abmin, rulea$abmax)), by=gcols]
      
      l = length(unique(MZ[,do.call(paste, .SD), .SDcols = gcols])); if (l == ngs) break else ngs = l
    }
    MZ[,rg := as.integer(factor(paste(do.call(paste, .SD)))), .SDcols = gcols][,(gcols) :=NULL]

    numg = MZ[,.(.N, nu = length(unique(r.n))),by="rg"]

    copy(MZ[numg[N>1 & nu > 1],,on="rg"][,.(cc, z, rg, rule = rep(rulea$rule, length(cc)))])
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
