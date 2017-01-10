warpgroup.nmacha = function(Nmacha, ugs = NULL, maxdriftrt = 5, maxdriftppm = 1, fracobs = 0.5, fraccontrib = 0.5, refit.var = c(2, 0.5, .1), do.plot = F) {
  cat("\nStarting aggressive refinement with warpgroup.\n")
  library(warpgroup)


  if (is.null(ugs)) ugs = unique(Nmacha$m.c$g)

  glen = length(ugs)
  cc.l = foreach(g.g = ugs, j = icount(), .errorhandling = 'pass') %do% {
    if (j %% 1 == 0) { cat(paste0("\r", j, " of ", glen, " mass channels analyzed. (", round(j/glen*100), "%)              ")) }

    .g = g.g

    Group = getgroup(Nmacha, .g)
    if (do.plot) plot.group(Nmacha, .g)

    #Search through each macha and find ROIs which could possibly contain peaks in this group (mainly concerned about missed peaks here)
    rtrange = c(min(Group$cs$rtmin.g)-maxdriftrt, max(Group$cs$rtmax.g) + maxdriftrt)
    maxdriftmz = maxdriftppm * mean(Group$cs$mz.g)/ 1E6
    mzrange = c(min(Group$cs$mz.g)-maxdriftmz, max(Group$cs$mz.g) + maxdriftmz)

    putativerois = Nmacha$m.r[minmz < mzrange[2] & maxmz > mzrange[1] & minrt < rtrange[2] & maxrt > rtrange[1]]

    # Research each macha to warn of mass peaks which did not qualify to be an ROI but are in our interesting range.
    if (F) {
      missingrois = lapply(seq_along(Nmacha$m), function(j) {
        M = Nmacha$m[[j]]

        mzrange.u = cormz(mzrange, j, Nmacha$gmz, uncorrect = T)
        rtrange.u = corrt(rtrange, j, Nmacha$grt, uncorrect = T)

        M$r[M$k,,nomatch = NA, on="k"][M$s,,on="s"][mz < mzrange.u[2] & mz > mzrange.u[1] & rt < rtrange.u[2] & rt > rtrange.u[1]]$r %>% unique
      })
      testmissing = sapply(missingrois, function(x) which(is.na(x))) %>% unlist
      if (any(testmissing)) warning(length(testmissing), " peaks were detected and could be the analyte but were not included in ROIs")
    }

    #Build the inputs for warpgroup
    rate = mean(Group$rs[,diff(rt),by=c("m", "r")]$V1)*0.9
    rtsout = seq(rtrange[1], rtrange[2], by = rate)

    eic.mat = lapply(seq_len(nrow(putativerois)), function(j) {
      R = putativerois[j,]

      rid = R$r
      roi = Group$rs[r == rid & m == R$m]
      if (nrow(roi) < 1) roi = getroi(Nmacha$m[[R$m]], rid)

      y = approx(corrt(roi$rt, R$m, Nmacha$grt), roi$ii, xout = rtsout)$y
      #y = approx(roi$rt, roi$ii, xout = rtsout)$y

      y[is.na(y)] = min(y, na.rm=T)
      y
    }) %>% do.call(what = cbind)

    ps = matrix(c(approx(rtsout, seq_along(rtsout), xout=unlist(Group$cs[,.(rtpeak.g, rtmin.g, rtmax.g)]))$y,unlist(Group$cs[,m])), nrow = nrow(Group$cs), dimnames=list(NULL, c("sc", "scmin", "scmax", "sample")))
    ps[,1] = round(ps[,1]); ps[,2] = floor(ps[,2]); ps[,3] = ceiling(ps[,3])

    if (ncol(eic.mat) != length(Nmacha$m)) ps[,"sample"] = as.numeric(factor(ps[,"sample"], levels = unique(putativerois$m)))
    # Call Warpgroup
    wgs = warpgroup(ps, eic.mat, sc.aligned.lim = 3, pct.pad = 0.1, min.peaks = 1)
    if (ncol(eic.mat) != length(Nmacha$m)) wgs = lapply(wgs, function(x) { x[,"sample"] = putativerois$m[x[,"sample"]]; x } )
#1
    ## Remove warpgroups spawned by fewer groups than fraccontrib
    contribs = sapply(wgs, function(x) { sum(!is.na(x[,"n"])) })
    wgs = wgs[contribs >= fraccontrib * length(Nmacha$m)]

    if (length(wgs) == 0) return(NULL)

    wgs = lapply(seq_along(wgs), function(j) {
      x = wgs[[j]]
      x = x[!duplicated(x[,"sample"]),,drop=F]
      cbind(x, m = putativerois$m[x[,"sample"]], r = putativerois$r[x[,"sample"]], c=Group$cs$c[x[,"n"]], wg = j)
    })


    ## Find a file which has the most seeds
    seed.origins = lapply(wgs, function(x) x[!is.na(x[,"n"]),"sample"])
    heuristicfiles = which(table(unlist(lapply(seed.origins, unique))) %>% { . == max(.)} ) %>% names %>% as.numeric
    heuristicfile = sample(heuristicfiles, 1)

    if (F) { # TODO: Add better seed selection code?

      pseeds.h = sapply(wgs, function(x) {
        x[x[,"sample"] == heuristicfile,,drop=F][1,,drop=F]
      }) %>% aperm
      colnames(pseeds.h) = colnames(wgs[[1]])

      overlapallowed = 2

      outer(pseeds.h[,"scmin"], pseeds.h[,"sc"], "-") > overlapallowed & outer(pseeds.h[,"scmax"], pseeds.h[,"sc"], "-") < overlapallowed
    }
#2


    # Define fitting constraints (so we wind up with reasonably similar peaks)
    .cs = sapply(wgs, function (x) x[x[,"sample"] %in% heuristicfiles,"n"])

    seed.constraints = apply(.cs, 2, function(x) {
      Group$cs[x][,.(location, scale, shape, factor)] %>% colMeans
      }) %>% aperm


    seed.constraints.var = rep(c(refit.var,0), each = nrow(seed.constraints))

    # Build input for fitting and Fit
    dt = wgs %>% do.call(what = rbind) %>% data.table

#3
    # For each file
    cc.l = lapply(seq_along(putativerois$r), function(.roin) {

      # BUild input for fitting
      seedsdt = dt[sample==putativerois$m[.roin]]
      seeds = seedsdt[,.(scmin, scmax, sc)] %>% as.matrix

      eic = cbind(rt = rtsout, ii = eic.mat[,.roin])

#3.6
      # Set file specific constraints
      local.const = seed.constraints
      local.const[,"location"] = rtsout[seeds[,"sc"]]

      const.lower = local.const - seed.constraints.var
      const.upper = local.const + seed.constraints.var
      const.lower[,4] = 0
      const.upper[,4] = 20
      #3.5

      # Refit
      components = fitseeds(eic, seeds = seeds, unrelated.dist = 30, const.upper = const.upper, const.lower = const.lower, do.plot = do.plot)

      #Store results as consensus compounds .cc
      # cc
      # newcs, mz, r, c, m, g, wg, cc

      mz = c(apply(components, 2, function(x) {
        range = c(x["rtmin"], x["rtmax"])
        range = corrt(range, putativerois$m[.roin], Nmacha$grt, uncorrect = T)

        roi = Group$rs[r == putativerois$r[.roin]]
        if (nrow(roi) < 1) roi = getroi(Nmacha$m[[putativerois$m[.roin]]], putativerois$r[.roin])

        proi = roi[{ind <- roi[.(c(range[1], range[2])), which=TRUE, roll=TRUE, on="rt"]; ind[1][is.na(ind[1])] = 1; ind[2][is.na(ind[2])] = nrow(roi); (ind[1]+1):ind[2]}]

        sum(matrixStats::rowProds(proi[,.(mz, i)] %>% as.matrix),na.rm=T)/sum(proi$i,na.rm=T)
      }))

      rbind(components %>% as.matrix, mz, r = putativerois$r[.roin], m = putativerois$m[.roin], c = Group$cs[seedsdt$n]$c, g = .g, wg = seedsdt$wg) %>% aperm
    }) %>% do.call(what = rbind)

  }
  cat("\nAggregating results.")

  ise = sapply(cc.l, function(x) "error" %in% class(x))
  if (any(ise)) {Nmacha$cc_error = cc.l[ise]; names(Nmacha$cc_error) = as.character(ugs[ise]); warning(paste0(sum(ise), " warpgroup refinement errors. (Stored in $cc_error.)")) } else { Nmacha = Nmacha[!names(Nmacha) == "cc_error"] }

  ism = sapply(cc.l, function(x) "matrix" %in% class(x))
  Nmacha$m.c_cc = data.table(do.call(rbind, cc.l[ism]))
  if (nrow(Nmacha$m.c_cc) > 0) {
    Nmacha$m.c_cc[,cc:=as.numeric(factor(paste(g,wg)))]
    Nmacha$m.c_cc[,g:=NULL]

    Nmacha$cc = Nmacha$m.c_cc[,.(mz = mean(mz, na.rm=T), rt = mean(rtpeak), rtmin = mean(rtmin), rtmax = mean(rtmax), i = mean(intpeak), cv=sd(intpeak)/mean(intpeak), n = length(mz), detected = sum(!is.na(c))), by="cc"]
  }

  isn = sapply(cc.l, is.null)

  if (any(!(ise | ism | isn))) warning(paste0(sum(!(ise | ism | isn)), " rois discarded. This is a bug."))

  cat("\nFinished warpgrouping.")

  Nmacha
}
