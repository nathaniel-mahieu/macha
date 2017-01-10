#' Find initial guesses for peaks and valleys
#'
#' \code{seedpeaks} Uses "scale space peak picking" to quickly find robust local maxima and minima
#'
#' Used to provide an initialization for further peak modeling and component recovery.
#'
#' @param roi Matrix. Rows are observed mass peaks. Columns named "mz", "rt", "s", "b", "i".  Column "b" is the baseline at that point.
#' @param S Integer vector. The smoothing scales to apply when searching for peaks.
#' @param seed.maxdist Integer. The maximum distance in scans the local maximum of the smoothed trace can travel and still positively reinforce a peak.
#' @param sn.adjust Numeric. Multiplied by the signal to noise limit.  The limit is variable but a value of 1 here corresponds to a SN range of 1.5 to 2.5.  See below.
#' @param min.seedwidth Numeric. Minimum width in scans of a putative region to be returned as a seed.
#' @param total.scans Integer. The maximum possible length of a mass channel. Used to estimate how much noise is in a channel and thereby scale the SN cutoff.
#' @param noisy.chan.n Integer. The number of scans a mass-channel must be observed in to be assigned the maximum SN limit.  Defaults to 0.8*total.scans
#' @param do.plot Boolean. Plot a bunch of random things.
#'
#'
#' @return Matrix. Columns "start", "end", "peak", correspond to the indices of the peak and its flanking valleys.
#'
#' @seealso \code{\link[sspp]}
#'
#' @export
#'

seedpeaks = function(roi, total.scans, S=3:8, seed.maxdist=2.5, noisy.chan.n=NULL, sn.adjust = 1, min.seedwidth = 2, do.plot=F) {
  # roi = data.table with: factored rt, i, b, d


  #meanrtd = mean(diff(roi$rt))
  #minwidth = round(minwidth/meanrtd)
  #maxdist = round(maxdist/meanrtd)

  #Make EIC Matrix
  eic = as.matrix(roi)
  if (do.plot) {
    plot(eic[,c("rt","i")], type="l", ylim=c(0,max(eic[,"i"])))
    lines(eic[,c("rt","b")])
  }

  v = eic[,"ii"]

  #Check Inputs
  S = S[order(S)]
  if (length(v) < min(S)) { return(matrix(nrow = 0, ncol = 3)) }
  if (length(v) < max(S)) { S = seq(S[1], length(v), by=1) }
  if (is.null(noisy.chan.n)) { noisy.chan.n = 0.8 * total.scans}


  #Find stable local maxima and minima with scale space peak picking
  C = sspp(v, S = S, maxdist = seed.maxdist, do.plot = F)
  vinv = -v + max(v)
  Cv = sspp(vinv, S = S, maxdist = seed.maxdist, do.plot = F)

  #Set the S/N value needed for a peak to be seeded.  This is scaled to the density in the mass channel.
  d = max(roi[,"d"], na.rm=T)
  cutoff = if (d > noisy.chan.n) { 2.5 * sn.adjust } else { (1.5 + d/total.scans) * sn.adjust }
  C[v/eic[,"b"] < cutoff] = 0

  plot(eic[,"i"]/eic[,"bb"])
  plot(eic[,"ii"]/eic[,"bb"])



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


  if (F) { #Functionality to iteratively insert peaks based on scores. Possibly could insert peaks until sufficiently fit and then terminate. Unused.
    starting.valleys = length(v)/6

    iterations = 10

    peakranges = seq(max(C), 0, -max(C)/iterations)
    peakgroups = rev(split(seq_along(C), cut(C, peakranges)))

    vranges = seq(max(Cv), 0, -max(Cv)/iterations)
    vgroups = rev(split(seq_along(Cv), cut(Cv, vranges)))


    iter = 10
    vsections = unlist(vgroups[1:iter]); vsections = vsections[order(vsections)]
    peakpoints = unlist(peakgroups[1:iter])
    peaksections = cut(peakpoints, breaks = vsections) %>% as.numeric

  }

  # Find EIC regions to seed peaks
  ploci = which(C > 0)
  vloci = which(Cv > 0)

  vloci = c(0,vloci[order(vloci)],Inf)
  peaksections = cut(ploci, breaks = vloci) %>% as.numeric

  peaksections.mat = cbind(start=vloci[peaksections], end = vloci[peaksections+1], peak = ploci)
  peaksections.mat[peaksections.mat > length(v)] = length(v)

  peaksections.mat = cbind(peaksections.mat, score = C[ploci])

  #Throw away seeds which are too narrow
  keep = which(matrixStats::rowDiffs(peaksections.mat[,1:2,drop=F]) > min.seedwidth)
  peaksections.mat = peaksections.mat[keep,,drop=F]

  peaksections.mat
}


