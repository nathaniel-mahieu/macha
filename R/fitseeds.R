#' Find a set of parameterized curves which describe the ROI well
#'
#' \code{fitandreseed} Iteritively fits a sum of Skew-Normal distibutions to the ROI until quality parameters are met.
#'
#' Used to provide a parameterized rrepresentation of the ROI
#'
#' @param roi data.table. ROI
#' @param seeds matrix. Columns "start", "end", "peak". Rows are individual initialization points for ROI fitting. Values correspond to the index within the supplied ROI.
#' @param unrelated.dist Integer. The distance two seeds must be to fit them individually (much faster).
#' @param min.peakwidth numeric. The minimum width a component can be to be retained. In seconds.
#' @param sn.adjust numeric. A multiplier for the signal to noise limit of retained components. Calculated from the raw data and the supplied baseline.  Used to determine the number of observations as well.
#' @param min.sharpness numeric. A lower limit on peak.height/peak.width.  Units Intensity/Second.
#' @param min.fractionobs numeric. A lower limit on the fraction of expected observations present in the peak region.
#' @param do.plot boolean. Plot a bunch of random stuff.
#'
#' @return List of matrices. Matrices include the parameterization for each component. Eg. location, dispersion, skew, scaling factor, baseline. Additional 95th percentile values are also returned as rtmin and rtmax.
#'
#' @export
#'

# These should all be in RT
fitandreseed = function(roi, seeds, unrelated.dist=30, min.peakwidth = 3, sn.adjust = 1, min.sharpness = 1, min.fracobs = 0.4, do.plot = F) {
  #Make EIC Matrix
  eic = as.matrix(roi)
  setkey(roi, "rt")

  meanrtd = mean(diff(roi$rt))
  min.peakwidth.s = min.peakwidth / meanrtd

  if (do.plot) {
    plot(eic[,c("rt","i")], type="l", ylim=c(0,max(eic[,"i"],na.rm=T)))
    lines(eic[,c("rt","b")])
  }

  #Split if they can be fitted individually (for speed)
  o = order(seeds[,3])
  peaksplits = which(diff(seeds[o,3]) > unrelated.dist)
  peaksectiongroups = split(o, cut(seq_along(o), c(-Inf, peaksplits, Inf)))



  seeds.l = lapply(unname(peaksectiongroups), function(x) { seeds[x,,drop=F] })

  component.mat.l = lapply(seeds.l, function(seeds) {
    iteration = 0
    repeat {
      iteration = iteration + 1
      if (iteration > 8) stop("Did not converge fitting.")
      if (nrow(seeds) < 1) { return(NULL) }

      #Fit peaks to our seeds
      component.mat = fitseeds(eic, seeds, do.plot=do.plot, unrelated.dist = unrelated.dist)

      #Learn about fits
      qc = apply(component.mat, 2, function(x) {
        range = x[c("rtmin", "rtmax")]
        pregion = roi[{ind <- roi[.(c(range[1], range[2])), which=TRUE, roll=TRUE, on="rt"]; ind[1][is.na(ind[1])] = 1; ind[2][is.na(ind[2])] = nrow(roi); (ind[1]+1):ind[2]}]

        fitted = curvemany(x[1:5], pregion$rt)
        fracexpl = sum((pregion$ii - fitted)^2)^0.5 / sum(pregion$ii - x[5])

        obs = sum(pregion$i > pregion$bb * sn.adjust,na.rm=T)
        wid = unname(range[2] - range[1])
        sharpness = unname(x[["intpeak"]]/wid)

        c(wid = wid, obs = obs, fracobs = obs/nrow(pregion), fracexpl = fracexpl, sharpness = sharpness)
        })

      # Remove terrible fits
      keeps = qc["wid",] > min.peakwidth.s/6 & component.mat["factor",] > 0
      component.mat = component.mat[,keeps,drop=F]
      if (sum(keeps) < 1) { return(NULL) }

      #Find local maxima in fitted chromatogram (Lets leverage all this work we did.)
      fitted = curvemany(c(component.mat[1:4,], 0), eic[,"rt"])
      #plot(v[min(seeds):max(seeds)], type="l"); points(fitted, col="red")
      ploci = localMaxima(fitted)
      ploci = ploci[ploci != 1 & ploci != nrow(eic)]
      fitinv = -fitted + max(fitted)
      vloci = localMaxima(fitinv)

      #Find if components have corresponding fitted local maxima
      dists = apply(abs(outer(component.mat["rtpeak",], eic[ploci,"rt"], "-")),2, which.min)

      # If we have a stable set of peaks refine them
      stable.tf = length(unique(dists)) == ncol(component.mat) & all(keeps)
      if (stable.tf) {
        #If we have suitably refined peaks return them
        passrefinement = qc["fracobs",keeps] > min.fracobs & qc["sharpness",keeps] > min.sharpness & qc["wid",keeps] > min.peakwidth
        if (all(passrefinement)) return(component.mat)
        if (all(!passrefinement)) return(NULL)

        component.mat = component.mat[,passrefinement, drop=F]

        fitted = curvemany(c(component.mat[1:4,], 0), eic[,"rt"])
        #plot(v[min(seeds):max(seeds)], type="l"); points(fitted, col="red")
        ploci = localMaxima(fitted)
        ploci = ploci[ploci != 1 & ploci != nrow(eic)]
        fitinv = -fitted + max(fitted)
        vloci = localMaxima(fitinv)
      }

      #Otherwise find new peak regions to fit based on our fitted distribution
      vloci = c(0,vloci[order(vloci)],Inf)
      peaksections = cut(ploci, breaks = vloci) %>% as.numeric
      seeds = cbind(start=vloci[peaksections], end = vloci[peaksections+1], peak = ploci)
    }
    if (iteration > 5) warning(paste0("Fitting took ", iteration, " iterations."))

    component.mat
  })

  component.mat.l
}


fitseeds = function(eic, seeds, unrelated.dist = 0, parscale = c(0.1, 0.01, 0.01, 10, 100), const.lower=NULL, const.upper=NULL, usebaseline = T, do.plot=F) {



  #if (is.null(const.lower) | is.null(const.upper)) { stop("Only one set of bounds supplied") }

  # Find initial parameters for fitting the groups of peaks by fitting each individual chromatogram section
  parhat.mat = matrix(ncol = nrow(seeds), nrow = 5)

  for (i in seq_len(nrow(seeds))) {
    inds = seeds[i,1]:seeds[i,2]

    y = eic[inds,"ii"] / max(eic[inds,"ii"])
    x = eic[inds,"rt"]

    par = c(eic[seeds[i,3],"rt"], 1, 0, 1, 0)
    #Starting curve
    #plot(x, y, type="l")
    #points(x, par[4] * dsn(x, dp = par[1:3]) + par[5],col="red")



    if (!is.null(const.lower)) {
      const.lowerx = c(const.lower[i,], 0)
      const.upperx = c(const.upper[i,], 20)

      parhat = tryCatch(
        optim(par, fitmany, x=x, y=y, method = "L-BFGS-B", lower = const.lowerx, upper = const.upperx, control = list(parscale = parscale*10))$par,
        error = function(e) {
          warning("Fitting BFGS failed in stage 1. Fell back to slow method. Unbounded!!");
          optim(par, fitmany, x=x, y=y, method = "Nelder-Mead", control = list(parscale = parscale*10))$par
        })
    } else {
      parhat =  tryCatch(
        optim(par, fitmany, x=x, y=y, method = "BFGS", control = list(parscale = parscale*10))$par,
        error = function(e) {
          warning("Fitting BFGS failed in stage 1. Fell back to slow method.");
          optim(par, fitmany, x=x, y=y, method = "Nelder-Mead", control = list(parscale = parscale*10))$par
        })
    }

    #Fitted single curve

    if (do.plot) {
      plot(x,y, type="l")
      yhat = parhat[4] * dsn(x, dp = parhat[1:3]) + parhat[5]
      points(x, yhat)
    }

    parhat.mat[,i] = parhat
  }



  #Make individual peaks scaled relative to one another
  scale.factor = eic[seeds[,3],"ii"]/max(eic[seeds[,3],"ii"])
  parhat.mat[4,] = parhat.mat[4,] * scale.factor
  parhat.mat[5,] = parhat.mat[5,] * scale.factor


  grow.range = range(c(seeds-unrelated.dist*0.4, seeds + unrelated.dist*0.4)) %>% { .[. < 1] = 1; . } %>% { .[. > nrow(eic)] = nrow(eic); . }
  inds = grow.range[1]:grow.range[2]

  scale.factor = 1/max(eic[seeds[,3],"ii"])

  y = eic[inds,"ii"] * scale.factor
  x = eic[inds,"rt"]


  if (do.plot) {
    plot(y, type="l")
    foreach (par = iter(parhat.mat, by="col")) %do% {
      lines(curvemany(c(par[1:4],0),x), col="red")
    }
    lines(curvemany(c(c(parhat.mat[1:4,]),.01),x), col="blue")
  }

  #inst/development/optimize_optim.R

  if (!is.null(const.lower)) {
    const.lower = c(const.lower %>% aperm, 0)
    const.upper = c(const.upper %>% aperm, 0.0001)

    opt = tryCatch(
      optim(c(parhat.mat[1:4,], 0), fitmany, x=x, y=y, method = "L-BFGS-B", lower = const.lower, upper = const.upper, control = list(parscale = c(rep(parscale[1:4], ncol(parhat.mat)), parscale[5]))),
      error = function(e) {
        warning("Fitting BFGS failed in stage 2. Fell back to slow method. Unbounded!!");
        optim(c(parhat.mat[1:4,], 0), fitmany, x=x, y=y, method = "Nelder-Mead", control = list(parscale = parscale))
      })
  } else if (usebaseline) {
    b = mean(eic[inds,"bb"],na.rm=T) * scale.factor

    opt = tryCatch(
      optim(c(parhat.mat[1:4,]), fitmanyb, b=b, x=x, y=y, method = "BFGS", control = list(parscale = c(rep(parscale[1:4], ncol(parhat.mat))))),
      error = function(e) {
        warning("Fitting BFGS failed in stage 2. Fell back to slow method.");
        optim(c(parhat.mat[1:4,]), fitmanyb, b=b, x=x, y=y, method = "Nelder-Mead", control = list(parscale = parscale))
      })
    opt$par = c(opt$par, b)
  } else {
    opt = tryCatch(
      optim(c(parhat.mat[1:4,], 0), fitmany, x=x, y=y, method = "BFGS", control = list(parscale = c(rep(parscale[1:4], ncol(parhat.mat))))),
      error = function(e) {
        warning("Fitting BFGS failed in stage 2. Fell back to slow method.");
        optim(c(parhat.mat[1:4,], 0), fitmany, x=x, y=y, method = "Nelder-Mead", control = list(parscale = parscale))
      })
  }



  components = matrix(opt$par[-length(opt$par)],nrow = 4)
  components[4,] = components[4,] / scale.factor
  baseline = opt$par[length(opt$par)] / scale.factor

  #Fitted Multiple Curves, Unscaled

  if (do.plot) {
    plot(eic[inds,c("rt","ii")], type="l")
    foreach (par = iter(matrix(components, nrow = 4), by="col")) %do% {
      lines(x, curvemany(c(par,0),x), col="red")
    }
    lines(x, curvemany(c(components,baseline),x), col="blue")
  }

  # Record some information about the fitted peaks. This could be recaptured elewhere as well.
  peaks = apply(components, 2, function(c) {
    v4 = dsn(x, dp = c[1:3])
    wm = which.max(v4)

    se = c(0,0)
    try({se = round(qsn(c(0.02, 0.98), dp = c[1:3]), 2)}, silent = T)

    c(v4[wm]*c[4]+baseline, x[wm], se)
  })

  component.mat = rbind(components, baseline, peaks)
  rownames(component.mat) = c("location", "scale", "shape", "factor", "baseline", "intpeak", "rtpeak", "rtmin", "rtmax")


  if (do.plot) {
    plot(eic[,c("rt", "ii")], type="l", main="plot")
    lines(x,curvemany(c(component.mat[1:4,], 0),x), col="red")
    for (i in seq_len(ncol(component.mat))) {
      lines(x,curvemany(c(component.mat[1:4,i], 0),x), col="blue")
    }
  }

  component.mat
  }


fitmany = function(parvec, x, y) {
  parmat = matrix(parvec[-length(parvec)], nrow = 4)

  #1 xi: location
  #2 omega: scale
  #3 alpha: slant
  #4 height
  #5 baseline


  #yhat = apply(parmat, 2, function(par) {
  #  par[4] * dsn(x, xi = par[1], omega=par[2], alpha=par[3])
  #})
  #sum((y-rowSums(yhat) - parvec[length(parvec)])^2)

  yhat = parmat[4,] * suppressWarnings( dsn(rep(x, each = ncol(parmat)), xi = parmat[1,], omega = parmat[2,], alpha = parmat[3,]) )
  sum((y-colSums(matrix(yhat, nrow = ncol(parmat))) - parvec[length(parvec)])^2)
}

fitmanyb = function(parvec, x, y, b) {
  parmat = matrix(parvec, nrow = 4)

  yhat = parmat[4,] * suppressWarnings( dsn(rep(x, each = ncol(parmat)), xi = parmat[1,], omega = parmat[2,], alpha = parmat[3,]) )
  sum((y-colSums(matrix(yhat, nrow = ncol(parmat))) - b)^2)
}

curvemany = function(parvec, x) {
  parmat = matrix(parvec[-length(parvec)], nrow = 4)

  #yhat = apply(parmat, 2, function(par) {
  #  par[4] * dsn(x, xi = par[1], omega=par[2], alpha=par[3])
  #})

  #rowSums(yhat) + parvec[length(parvec)]

  yhat = parmat[4,] * suppressWarnings( dsn(rep(x, each = ncol(parmat)), xi = parmat[1,], omega = parmat[2,], alpha = parmat[3,]) )
  colSums(matrix(yhat, nrow = ncol(parmat))) + parvec[length(parvec)]
}

