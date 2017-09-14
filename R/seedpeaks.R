#' Find initial guesses for peaks and valleys
#'
#' \code{seedpeaks} Uses "scale space peak picking" to quickly find robust local maxima and minima
#'
#' Used to provide an initialization for further peak modeling and component recovery.
#'
#' @param roi Matrix. Rows are observed mass peaks. Columns named "mz", "rt", "s", "b", "i".  Column "b" is the baseline at that point.
#' @param S Integer vector. The smoothing scales to apply when searching for peaks.
#' @param seed.maxdist Integer. The maximum distance in scans the local maximum of the smoothed trace can travel and still positively reinforce a peak.
#' @param seed.minwidth Numeric. Minimum width in scans of a putative region to be returned as a seed.
#' @param seed.maxdensity Float. Maximum density in which to seed peaks in peaks/scan. Eg 1/5.
#' @param seed.sn.perpeak Numeric vector. Indices of the vector correspond to the number of observations which must be made of the mass in the peak region defined by seed.sn.range. The value at each index defines the minimum intensity necessary for an observation to count in value * baseline. Eg. c(Inf, 10) means a seed with only 1 observation will not be returned whereas a seed with 2 observations at 10*baseline will pass. 
#' @param seed.sn.range Integer. Range in scans around a peak loci in which to look for observations qualifying as in seed.sn.perpeak.
#' @param do.plot Boolean. Plot a bunch of random things.
#'
#'
#' @return Matrix. Columns "start", "end", "peak", correspond to the indices of the peak and its flanking valleys.
#'
#' @seealso \code{\link[sspp]}
#'
#' @export
#'

#These should all be in RT
seedpeaks = function(roi, S=3:8, seed.maxdensity = 1/5, seed.maxdist=4, seed.sn.perpeak = c(Inf, 50, 5, 3, 2), seed.sn.range = 2, seed.minwidth = 3, do.plot=F) {
  # roi = data.table with: factored rt, i, b, d

  nullresult = matrix(nrow = 0, ncol = 3, dimnames=list(NULL, c("start", "end", "peak")))

  if (do.plot) plot.roi(roi)
  eic = as.matrix(roi)

  v = eic[,"ii"]
  orig.region = 1:length(v) + max(S)
  v.pad = c(rep(0,max(S)), v, rep(0, max(S)))


  #Check Inputs
  S = S[order(S)]
  if (length(v) < min(S)) { return(nullresult) }
  if (length(v.pad) < max(S)) { S = seq(S[1], length(v), by=1) }

  #Find stable local maxima and minima with scale space peak picking
  C = sspp(v.pad, S = S, maxdist = seed.maxdist, do.plot = do.plot)[orig.region]
  vinv = -v.pad + max(v.pad)
  Cv = sspp(vinv, S = S, maxdist = seed.maxdist, do.plot = do.plot)[orig.region]

  if (do.plot) {
    scores = cbind(which(C>1), C[C>1])
    scores[,2] = scores[,2] * max(v)/max(scores[,2])

    par(mfcol=c(1,1))
    plot(v, type="l")
    lines(scores, type="h", col = "red")

    scores = cbind(which(Cv>1), Cv[Cv>1])
    scores[,2] = scores[,2] * max(v)/max(scores[,2])
    lines(scores, type="h", col = "blue")
    }

  # Find EIC regions to seed peaks
  ploci = which(C > 0)
  vloci = which(Cv > 0)

  if (length(ploci) == 0) return(nullresult)

  #Filter peaks for SN
  plocisn = matrix(rep(seq(-seed.sn.range, seed.sn.range), each = length(ploci)) + ploci, ncol=1+2*seed.sn.range) + seed.sn.range
  multbl = c(rep(0,seed.sn.range), eic[,"i"]/(eic[,"bb"]), rep(0,seed.sn.range))
  multbl[is.na(multbl)] = 0
  plocisn[] = multbl[plocisn]

  keep.sn = matrixStats::rowAnys(sapply(seq_along(seed.sn.perpeak), function (i) {
    rowSums(plocisn > seed.sn.perpeak[i]) >= i
  }) %>% { if (is.null(dim(.))) { dim(.) = c(nrow(plocisn), length(seed.sn.perpeak)) }; . })

  if (sum(keep.sn) == 0) return(nullresult)

  vloci = c(0,vloci[order(vloci)],Inf)
  peaksections = cut(ploci, breaks = vloci) %>% as.numeric

  peaksections.mat = cbind(start=vloci[peaksections], end = vloci[peaksections+1], peak = ploci)
  peaksections.mat[peaksections.mat > length(v)] = length(v)

  #peaksections.mat = cbind(peaksections.mat, score = C[ploci])  #Do we want to return scores? Can we provide any information to help guide the fitting?

  #Throw away seeds which do not meet our quality metrics
  peaksections.mat = peaksections.mat[keep.sn,,drop=F]
  denregions = cut(peaksections.mat[,"peak"], unique(c(seq(0, length(v), by = 1/seed.maxdensity), length(v)))) %>% as.numeric

  peaksections.mat = lapply(split(seq_len(nrow(peaksections.mat)), denregions), function(inds) {
    if (length(inds) < 2) return(peaksections.mat[inds,,drop=F])

    x = peaksections.mat[inds,,drop=F]

    cbind(start = min(x), end = max(x), peak = floor(mean(x[,"peak"])))
    }) %>% do.call(what=rbind)

  keep.wid = matrixStats::rowDiffs(peaksections.mat[,1:2,drop=F]) > seed.minwidth
  peaksections.mat = peaksections.mat[keep.wid,,drop=F]

  peaksections.mat
}


