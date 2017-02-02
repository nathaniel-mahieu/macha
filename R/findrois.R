#' ROI detection based on density estimation.
#'
#' \code{findrois} takes mass spectra and assigns peaks to ROIs
#'
#' See: \code{\link[centroidgroups]}.  In addition this function assures that only one point per scan is in a ROI.
#'
findrois = function(features, minlength = 5, ppm = 1, rtwid = 10) {
  cat("Finding ROIs\n")

  rois = features$k[features$s[,.(s,rt)],.(mz, rt, i, k, k), on="s"] %>% { .[complete.cases(.),] } %>% as.matrix
  colnames(rois)[4] = "r"

  gl = centroidgroups(rois, ppm = ppm, rtwid = rtwid, minlength = minlength)

  cat("\nFiltering ROIs")
  rois[,4] = rep(seq_along(gl),sapply(gl, length))[order(unlist(gl))]
  rois = rois[rois[,4] %in% which(sapply(gl, length) >= minlength),]

  cat("\nMerging Centroids: This currently throws away all but the highest mass peak per retention time.")
  lenspre = sapply(split(seq_along(rois[,4]), rois[,4]), length)

  rois = mergecentroids(rois)

  lenspost = sapply(split(seq_along(rois[,4]), rois[,4]), length)
  cat(" - ROIs merged:", sum(lenspre - lenspost > 0))

  rois[,4] = rois[,4] %>% factor %>% as.numeric
  rois = rois[order(rois[,1], rois[,2]),]

  cat("\nComplete: Found", max(rois[,4]), "ROIs")
  rois

  features$k_r = data.table(rois[,c("k", "r"),drop=F])
  features$r = features$k[features$s,,on="s"][features$k_r,,on="k"][, .(minrt = min(rt), maxrt = max(rt), meanmz = mean(sum(mz*i)/sum(i)), maxmz = max(mz), minmz = min(mz)),by=r]

  features
}

kernelsplitmass = function(features, ppm=4) {
  masses = features[,1]
  intensi = features[,3]

  curr.bw = max(masses) * ppm / 1E6
  n = c(2^floor(log2(length(masses))-2), 1.5*diff(range(masses))/curr.bw)
  n = if (n[1] < 100) { n[2] } else { n[1] }

  if (n < 2) { return(rep(1,length(masses))) }

  all.mass.den<-density(masses, weights=intensi/sum(intensi), bw=curr.bw, n=n)
  #all.mass.den<-density(masses, bw=curr.bw*2, n=n*10)
  lm = localMinima(all.mass.den$y)
  breaks = approx(all.mass.den$x, xout = lm)$y
  #plot(masses); abline(h = breaks)

  cut(masses, breaks = c(-Inf, breaks, Inf)) %>% as.numeric %>% rank(ties.method = "min")
}

kernelsplitrt = function(features, rtwid = 10) {
  rts = features[,2]
  intensi = features[,3]

  curr.bw = rtwid
  n = c(2^floor(log2(length(rts))-2), 1.5*diff(range(rts))/curr.bw)
  n = if (n[1] < 100) { n[2] } else { n[1] }

  if (n < 2) { return(rep(1,length(rts))) }

  all.mass.den<-density(rts, weights=intensi/sum(intensi), bw=curr.bw*2, n=n)
  lm = localMinima(all.mass.den$y)
  breaks = approx(all.mass.den$x, xout = lm)$y

  cut(rts, breaks = c(-Inf, breaks, Inf)) %>% as.numeric %>% rank(ties.method = "min")
}

#' Speed minded function for iteratively splitting a group of points by two density estimates.
#'
#' \code{centroidgroups} takes mass spectra and assigns peaks to ROIs.
#'
#'
#'
#' @param features Matrix. Rows are observed mass peaks. Columns are ordered to contain "mz", "rt", "i".
#' @param ppm Numeric. Bandwidth of kernel for estimating density in the mz dimension.
#' @param rtwid Numeric. Bandwidth of kernel for estimating density in the rt dimension.
#' @param minlength Integer. If an ROI is shorter than this it is discarded.
#'
#'
#' @return Integer vector containing assignment of each point to a group.
#'
centroidgroups = function(features, ppm, rtwid, minlength) {
  warning("centroid groups is deprecated. use dengroup.ppm().")

  dengroup.ppm(features, ppm, rtwid, minlength)
  }
dengroup.ppm = function (features, ppm, rtwid, minlength) {
  l = nrow(features)

  groups = vector(mode = "list", length = l)
  action = vector(mode = 'numeric', length = l)
  unchanged = vector(mode = 'numeric', length = l)

  groups[[1]] = seq(l)
  action[1] = 0
  unchanged[1] = 0
  tg = 0
  i = 1
  groupmaxindex = 1
  repeat {
    if (i %% 100 == 0) { cat("\rProgress indicator:", groupmaxindex - i, "groups remaining to analyze. (This will grow before starting to decrease.)            ") }

    tg = groups[[i]]
    tg <<- tg

    if (length(tg) >= minlength - 1) {

      if (action[i] == 0) {
        ng = kernelsplitmass(features[tg,,drop=F], ppm)
        action[i] = 1
      } else {
        ng = kernelsplitrt(features[tg,,drop=F], rtwid)
        action[i] = 0
      }

      if (max(ng) > 1) {
        newgroups = split(tg, ng) %>% unname
        ngroups = length(newgroups)

        indices_to_replace = c(i,(groupmaxindex+1):(groupmaxindex+ngroups-1))

        groups[indices_to_replace] = unname(newgroups)
        action[indices_to_replace] = action[i]
        unchanged[indices_to_replace] = 0
        groupmaxindex = groupmaxindex + ngroups - 1

      } else if (unchanged[i] == 3) {
        i = i+1

      } else {
        unchanged[i] = unchanged[i] + 1
      }
    } else { i = i+1 }

    if ((i == groupmaxindex & i > 1) | (i == groupmaxindex & unchanged[i] > 2)) {
      cat("\rProgress indicator: 0 groups remaining to analyze. (This will grow before starting to decrease.)            ")
      return(groups[1:groupmaxindex])
    }
  }
}

filtercentroids = function(features, roi.foldrise = 2) { # Mostly crap - fitting needs to distinguish
  groups = features[,4]
  glt = split(seq_along(groups), groups)

  characteristics = sapply(glt[1:5000], function(x) {
    ints = features[x,3]

    sg.window = 3
    rollav = filter(ints, 1/rep(sg.window, sg.window)) %>% { .[is.na(.)] = ints[is.na(.)]; . }
    foo = ints - rollav
    sd = TTR::runSD(foo, sg.window) %>% { .[is.na(.)] = 0; . }

    (rollav/sd) %>% plot

    rollav %>% plot
    sd %>% plot


    medint = median(ints,na.rm=T)
    q25th = quantile(ints, 0.25)[[1]]

    c(
      medint = medint,
      q25th = q25th,
      foldrise = diff(range(ints))/q25th,
      foldrise2 = diff(range(ints)) / mean(sd)
    )
  }) %>% aperm

  hist(characteristics[,"foldrise2"],xlim = c(0, 30), breaks = 5000)

 features[characteristics[,"foldrise"],]
}

mergecentroids = function(features) { # Keeps highest intesnity
  groups = features[,4]
  glt = split(seq_along(groups), groups)

  gl2 = lapply(glt, function(x) {
    fs = features[x,,drop=F]
    keep = !matrixStats::colAnys(outer(fs[,2], fs[,2], "==") & outer(fs[,3], fs[,3], ">"))

    x[keep]
  })

  features = features[unlist(gl2),]
  features[,4] = rep(seq_along(gl2),sapply(gl2, length))

  features
}
