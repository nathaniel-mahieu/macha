
baseline.test = function(Macha, .rs, ppmwin, lambda1s, lambda2s) {

  cat("\nStarted baselining.")

  featsd = macha$k[,.(mz, s, i)]
  roisd = macha$k[macha$k_r[r%in%.rs],,on="k"]

  # Aggregate ROIs which are close in mass but possibly separate in retention time.
  cat("\nAggregating mass channels.")
  massdt = roisd[,mean(mz),by=r]

  o = order(massdt$V1)
  m = massdt$V1[o]
  jumps = diff(m) / m[-1] * 1E6

  gs = c(cumsum(c(0,jumps) > ppmwin - 0.5))
  toolarge = sapply(split(m, gs), function(x) { diff(range(x))/mean(x) * 1E6 }) %>% '>'(ppmwin) %>% which

  maxg = max(gs)
  for (gn in toolarge - 1) {
    wtl = which(gs == gn)
    gs[wtl] = maxg + 1 + seq_along(wtl)
    maxg = maxg + length(wtl)
  }

  massdt[,mchan:=as.numeric(factor(gs))[o]]
  roisd[massdt, mchan := mchan, on="r"]

  mchans = roisd[,.(mz = mean(mz)),by=mchan]
  ranges = cbind(mchans$mz-mchans$mz*(ppmwin+0.5)/1E6, mchans$mz + mchans$mz*(ppmwin+0.5)/1E6)


  cat("\nExtracting mass channels.\n")  # Put all ROIs into a matrix and baseline together (~20x faster than individual baselining.)

  setkey(macha$s, "s")
  setkey(featsd, "mz")
  mat = matrix(ncol = nrow(macha$s), nrow = nrow(ranges), dimnames = list(mchans$mchan, macha$s$s))

  nms = nrow(ranges)
  for (j in seq_len(nrow(ranges))) {
    if (j %% 50 == 0) { cat(paste0("\r", j, " of ", nms, " mass channels analyzed. (", round(j/nms*100), "%)              ")) }

    range = ranges[j,]

    #Stack exchange soltuion here: http://stackoverflow.com/questions/40665673/does-data-table-implement-fast-range-subsetting-based-on-binary-search-what-is
    #microbenchmark::microbenchmark(
    #  featsd[{ind <- featsd[.(c(range[1], range[2])), which=TRUE, roll=TRUE]; (ind[1]+1):ind[2]}],
    #  featsd[mz > range[1] & mz < range[2]]
    #  )

    m.inrange = featsd[{ind <- featsd[.(c(range[1], range[2])), which=TRUE, roll=TRUE]; ind[1][is.na(ind[1])] = 1; ind[2][is.na(ind[2])] = nrow(featsd); (ind[1]+1):ind[2]}]

    mat[j,] =  m.inrange[macha$s, sum(i, na.rm=T), roll=T, on="s", mult="all", by=.EACHI]$V1
  }

  cat("\nCalculating baselines.")


  vars = expand.grid(lambda1=lambda1s, lambda2=lambda2s)

  bl.l =lapply(seq_len(nrow(vars)), function(i) {
    bl = baseline::baseline(mat, lambda1 = vars[i,1], lambda2 = vars[i,2])@baseline

    m= melt(bl) %>% as.matrix
    m = cbind(m, matrix(unlist(vars[i,]), nrow = nrow(m), ncol=2, byrow = T))

    colnames(m) = c("roi", "scan", "ii", colnames(vars))
    m %>% as.matrix
  }) %>% do.call(what="rbind") %>% as.data.frame

  df = melt(mat)
  colnames(df) = c("roi", "scan", "ii")

  bl.l$scan = as.numeric(colnames(mat))[bl.l$scan]

  bl.l$lambda1 = factor(bl.l$lambda1)
  bl.l$lambda2 = factor(bl.l$lambda2)

  for (.r in unique(bl.l$roi)){
    g =  ggplot() + geom_line(data=subset(df,roi == .r), mapping=aes(x = scan, y = ii)) + geom_line(data = subset(bl.l,roi == .r), aes(x=scan, y = ii), colour="red") + facet_grid(lambda1 ~ lambda2) + theme_nate()
    print(g)
  }


}
