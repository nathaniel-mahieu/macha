warpgroup.nmacha.iter = function(Nmacha, ugs = NULL, maxdriftrt = 5, maxdriftppm = 1) {

  if (any(sapply(Nmacha$m, function(x) is.null(x$roi_cache)))) {
    warning("Set keys on individual machas for faster processing 'makeroicache()'.")
  }

  if (!all(c("m", "r") %in% key(Nmacha$m.r))) {
    warning("Set keys on $m.r to 'm' and 'r' for faster processing.")
  }

  do.plot=F
  if (is.null(ugs)) ugs = unique(Nmacha$m.c$g)
  ugs.iter = iter(ugs)

  nextEl = function() {
    g.g = nextElem(ugs.iter)
    .g = g.g

    Group = getgroup(Nmacha, .g)
    if (do.plot) plot.group(Nmacha, .g)

    Group$cs[is.na(rtmin.g),rtmin.g:=0]


    #Search through each macha and find ROIs which could possibly contain peaks in this group (mainly concerned about missed peaks here)
    rtrange = c(min(Group$cs$rtmin.g)-maxdriftrt, max(Group$cs$rtmax.g) + maxdriftrt)
    maxdriftmz = maxdriftppm * mean(Group$cs$mz.g)/ 1E6
    mzrange = c(min(Group$cs$mz.g)-maxdriftmz, max(Group$cs$mz.g) + maxdriftmz)



    putativerois = Nmacha$m.r[minmz < mzrange[2] & maxmz > mzrange[1] & minrt < rtrange[2] & maxrt > rtrange[1]]

    # Research each macha to warn of mass peaks which did not qualify to be an ROI but are in our interesting range.
    if (F) {
      missingrois = lapply(seq_along(Nmacha$m), function(j) {
        M = Nmacha$m[[j]]

        mzrange.u = cormz(mzrange, Nmacha$gmz[[j]], uncorrect = T)
        rtrange.u = corrt(rtrange, Nmacha$grt[[j]], uncorrect = T)

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
      if (nrow(roi) < 1) {
        roi = getroi(Nmacha$m[[R$m]], rid)[,m:=R$m]

        Group$rs <<- rbind(Group$rs, roi)
        }

      y = approx(corrt(roi$rt, Nmacha$grt[[R$m[1]]]), roi$ii, xout = rtsout)$y
      #y = approx(roi$rt, roi$ii, xout = rtsout)$y

      y[is.na(y)] = min(y, na.rm=T)
      y
    }) %>% do.call(what = cbind)

    bb.mat = lapply(seq_len(nrow(putativerois)), function(j) {
      R = putativerois[j,]

      rid = R$r
      roi = Group$rs[r == rid & m == R$m]
      if (nrow(roi) < 1) roi = getroi(Nmacha$m[[R$m]], rid)

      y = approx(corrt(roi$rt, Nmacha$grt[[R$m[1]]]), roi$bb, xout = rtsout)$y
      #y = approx(roi$rt, roi$ii, xout = rtsout)$y

      y[is.na(y)] = min(y, na.rm=T)
      y
    }) %>% do.call(what = cbind)

    ps = matrix(c(approx(rtsout, seq_along(rtsout), xout=unlist(Group$cs[,.(rtpeak.g, rtmin.g, rtmax.g)]))$y,unlist(Group$cs[,m])), nrow = nrow(Group$cs), dimnames=list(NULL, c("sc", "scmin", "scmax", "sample")))
    ps[,1] = round(ps[,1]); ps[,2] = floor(ps[,2]); ps[,3] = ceiling(ps[,3])

    colnames(eic.mat) = putativerois$m
    ps[,"sample"] = which(!duplicated(colnames(eic.mat)))[as.numeric(factor(ps[,"sample"], levels = unique(colnames(eic.mat))))] # This currently assigns a peak in the first roi supplied per file - shoudl assign to the detected ROI


    list(eic.mat = eic.mat, bb.mat = bb.mat, ps = ps, putativerois=copy(putativerois), Group = copy(Group), rtsout = rtsout, .g = .g, m.n = length(Nmacha$m))
  }


  obj <- list(nextElem=nextEl)
  class(obj) <- c('ixts', 'abstractiter', 'iter')
  obj
}

#' Collectively analyze the peaks detected across samples to determine consensus components
#'
#' \code{warpgroup.nmacha} Combines information from all independent peak detection rounds to determine a set of consensus components, greatly increasing the consistency of feature detection and integration across samples.
#'
#' After warpgrouping of the mass traces raw data is refit with seeded features.  The fits are constrained such that fit peaks are similar, but allow for variance from run to run. The allowed variance should reflect the observed variance in peak shape across runs.
#'
#' @param Nmacha An Nmacha object containing peaks grouped between samples.
#' @param ugs Integer. Vector of group numbers to warpgroup allowing a subset of the data to analyzed. If NULL all are warpgrouped.
#' @param	warpgroup.nmacha_data_l List. Pre-computed data to be distributed parallely.  Allows caching of an expensive step and faster recovery in the case of errors.
#' @param sc.aligned.lim Integer. The maximum distance (in scans) two points can remain while still being considered the same point. Applies to chromatographic peak start, end, and peak.  Is used to determine weather to consolidate or split up detected peaks.
#' @param pct.pad Numeric. The EIC is padded by this factor of the original length.  Padding decays from end value to zero over this length.
#' @param min.peaks Integer. The minimum number of peaks across samples in a group to apply warpgrouping. Saves computation time. Don't set higher than fracobs/fraccontrib
#' @param maxdriftrt Nuemric. The maximum observed retention time drift you wish to correct. EICs are expanded by this amount.  Drastically increases processing time.
#' @param fraccontrib Numeric. The minimum number of originally detected peaks remaining in a warpgroup in order to retain it.
#' @param refit.var Numeric. Vector with three values. Limits on the parameters of the refitted peaks. The fit peaks will be constrained to the mean of the contributing peaks parameters plus or minus refit.var. If peaks change throughout a run these values should eb larger.  The contraint ranges refer to the following parameters c(location, scale, shape)
#' @param do.plot Boolean. Plots a bunch of steps.  Slow if T.
#'
#' @return List (an Nmacha object) containing additional list named cc (consensus components), m.c_cc (mapping between m.c and cc).
#'
#' @examples
#' \dontrun{
#' warpgroup.nmacha(Nmacha, ugs = ugs, warpgroup.nmacha_data_l = warpgroup.nmacha_data_l, sc.aligned.lim = 9, pct.pad = 0.1, min.peaks = min.peaks, maxdriftrt = 1, maxdriftppm = 1, fracobs = 0.3, fraccontrib = 0.3, refit.var = c(2, 0.5, 1), do.plot = F)
#' }
#'
#' @export
#'
warpgroup.nmacha = function(Nmacha, ugs = NULL, warpgroup.nmacha_data_l = NULL, sc.aligned.lim = 9, pct.pad = 0.1, min.peaks = 1, maxdriftrt = 5, maxdriftppm = 1, fracobs = NULL, fraccontrib = 0.5, refit.var = c(2, 0.5, .1), do.plot = F) {
  cat("\nStarting aggressive refinement with warpgroup.\n")

  library(warpgroup)

  if (is.null(ugs)) ugs = unique(Nmacha$m.c$g)

  if (is.null(warpgroup.nmacha_data_l)) {
    warpgroup.nmacha_data_l = warpgroup.nmacha.iter(Nmacha, ugs, maxdriftppm, maxdriftrt)
    }

  # Call Warpgroup
  glen = length(ugs)

  grt = Nmacha$grt

  cc.l = foreach (params = warpgroup.nmacha_data_l, j = icount(), .errorhandling = "pass", .packages=c("macha", "warpgroup"), .noexport=c("Nmacha"), .options.redis=list(chunkSize=50)) %dopar% {
    #if (j %% 1 == 0) { cat(paste0("\r", j, " of ", glen, " mass channels analyzed. (", round(j/glen*100), "%)              ")) }

    ps = params[["ps"]]; eic.mat = params[["eic.mat"]]; bb.mat = params[["bb.mat"]]; putativerois = params[["putativerois"]]; Group = params[["Group"]]; rtsout = params[["rtsout"]]; .g = params[[".g"]]; grt = grt

    wgs = warpgroup(ps, eic.mat, sc.aligned.lim = sc.aligned.lim, pct.pad = pct.pad, min.peaks = min.peaks)

    wgs = lapply(wgs, function(x) { x[,"sample"] = as.numeric(colnames(eic.mat))[x[,"sample"]]; x } )

    #1
    ## Remove warpgroups spawned by fewer groups than fraccontrib
    contribs = sapply(wgs, function(x) { sum(!is.na(x[,"n"])) })
    wgs = wgs[contribs >= fraccontrib * length(params[['m.n']])]

    if (length(wgs) == 0) return(NULL)

    wgs = lapply(seq_along(wgs), function(j) {
      x = wgs[[j]]
      x = x[!duplicated(x[,"sample"]),,drop=F]
      cbind(x, m = putativerois$m[x[,"sample"]], r = putativerois$r[x[,"sample"]], c=Group$cs$c[x[,"n"]], wg = j)
    })


    ## Find a file which has the most seeds
    seed.origins = lapply(wgs, function(x) x[!is.na(x[,"n"]),"sample"])
    heuristicfiles = which(table(unlist(lapply(seed.origins, unique))) %>% { . == max(.)} ) %>% names %>% as.numeric
    resample <- function(x, ...) x[sample.int(length(x), ...)]
    heuristicfile = resample(heuristicfiles, 1)

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
    .cs = sapply(wgs, function (x) x[x[,"sample"] %in% heuristicfile,"n"])
    .cs = matrix(.cs, ncol = length(wgs))


    seed.constraints = apply(.cs, 2, function(x) {
      Group$cs[x][,.(location, scale, shape, factor)] %>% colMeans
    }) %>% aperm


    wgs = wgs[!is.na(seed.constraints[,1])]
    seed.constraints = seed.constraints[!is.na(seed.constraints[,1]),,drop=F]

    seed.constraints.var = rep(c(refit.var,0), each = nrow(seed.constraints))

    # Build input for fitting and Fit
    dt = wgs %>% do.call(what = rbind) %>% data.table

    #3
    # For each file
    wgpeaks = lapply(seq_along(putativerois$r), function(.roin) {

      # BUild input for fitting
      seedsdt = dt[sample==putativerois$m[.roin]]
      seeds = seedsdt[,.(scmin, scmax, sc)] %>% as.matrix

      bl = range(seeds) %>% { mean(bb.mat[.[1]:.[2],.roin]) }

      eic = cbind(rt = rtsout, ii = eic.mat[,.roin] - bl)

      seeds[seeds > nrow(eic)] = nrow(eic)
      seeds[seeds < 1] = 1

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

      components["baseline",] = components["baseline",] + bl

      mz = c(apply(components, 2, function(x) {
        range = c(x["rtmin"], x["rtmax"])
        range = corrt(range, grt[[putativerois$m[.roin]]], uncorrect = T)

        roi = Group$rs[r == putativerois$r[.roin]]
        #if (nrow(roi) < 1) roi = getroi(Nmacha$m[[putativerois$m[.roin]]], putativerois$r[.roin])

        proi = roi[{ind <- roi[.(c(range[1], range[2])), which=TRUE, roll=TRUE, on="rt"]; ind[1][is.na(ind[1])] = 1; ind[2][is.na(ind[2])] = nrow(roi); (ind[1]+1):ind[2]}]

        sum(matrixStats::rowProds(proi[,.(mz, i)] %>% as.matrix),na.rm=T)/sum(proi$i,na.rm=T)
      }))

      rbind(components %>% as.matrix, mz, r = putativerois$r[.roin], m = putativerois$m[.roin], c = Group$cs[seedsdt$n]$c, g = .g, wg = seedsdt$wg) %>% aperm
    }) %>% do.call(what = rbind) %>% data.table

    pickbest = function(x) {
      colns = colnames(x)
      if (!any(duplicated(x$m))) return(x)

      poss = which(x[,intpeak == max(intpeak),by="m"]$V1 + x[,!is.na(c)] + x[,!is.na(mz)] > 1)

      meanint = mean(x[poss,intpeak])

      x = x[,.SD[which.min(abs(intpeak - meanint))],by="m"]
      setcolorder(x, colns)
      x
      }

    wgpeaks[,pickbest(.SD),by="wg"] %>% as.matrix
  } #####


  #Store results as consensus compounds .cc
  # cc
  # newcs, mz, r, c, m, g, wg, cc
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
