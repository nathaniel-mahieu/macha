getroi.iter = function(macha, .roil) {
  it <- iter(.roil)

  nextEl = function() {
    getroi(macha, roi.id=nextElem(it))
  }

  obj <- list(nextElem=nextEl)
  class(obj) <- c('ixts', 'abstractiter', 'iter')
  obj
  }

#' Analyze ROIs for retention-time-domain peaks satisfying parameters
#'
#' \code{findcomponents} Sequentially analyzes detected ROIs for peaks. This proceeds by several steps.
#'
#' 1. Peaks are seeded by a fast, scale-space peak detection algorithm
#' 2. Skew-normal distributions are initialized and quickly fit to the seed regions.3
#' 3. Distributions within \code{unrelated.dist} of each other are fit simultaneously to the raw data.
#' 4. Fit distributions are checked by the supplied parameters and peaks not satisfying those parameters are removed.  The distribution is then refit.
#' 5. This repeats until the fit distributions stabilize.
#'
#' @param S The smoothing scales to apply when searching for peaks. Dependent on the amount of signal variation.  Suggested values: smooth peaks, 1:3; rough peaks, 3:(1/2 peakwidth in scans)
#' @param seed.maxdensity The maximum density at which to seed peaks in units of peaks/scan. Suggested value: 1/(One third the peak width in scans)
#' @param	seed.maxdist Integer. The maximum distance in scans the local maximum of the smoothed trace can travel and still positively reinforce a peak. Suggested value: One quarter the peak width in scans
#' @param	seed.minwidth Numeric. Minimum width in scans of a putative region to be returned as a seed.
#' @param	seed.sn.perpeak Required signal to noise (intensity/baseline) for a seed to be retained.  Each index of the vector corresponds to the number of peaks needed above the specified level.  c(Inf, 50, 10) requires a single point to be 50 times the baseline but a peak with three observed scans to be 10x the baseline.
#' @param	seed.sn.range Integer. The number of scans around a seed to look for peaks satisfying the signal to noise limits
#'
#'
#' @param unrelated.dist Integer. The distance two seeds must be apart to fit them individually (much faster). Suggested value: large peakwidth * 3
#' @param min.peakwidth numeric. The minimum width a component can be to be retained. In seconds.
#' @param min.sharpness numeric. A lower limit on peak.height/peak.width.  Units Intensity/Second.
#' @param min.fractionobs numeric. A lower limit on the fraction of expected observations present in the peak region.
#'
#' @return List (a Macha object) containing additional list named c containing detected components and the parameters of the skew normal distribution describing that component.
#'
#' @examples
#' \dontrun{
#' #Can run this locally in parallel to speed processing.  Avoids lots of network IO for each job and caches progress in case of a bug or interruption.
#' .roil = foreach(r = unique(macha$r$r), .packages="macha") %dopar% {
#'   nextElem(getroi.iter(macha, r))
#'   }
#' macha = findcomponents(
#'   macha, .roil = .roil,
#'   S = 3:7, seed.maxdensity=1/7, seed.maxdist=4, seed.sn.perpeak =c(Inf, 10, 7, 3, 2.5, 2), seed.sn.range = 3, seed.sn.adjust = 1, seed.minwidth = 4,
#'   unrelated.dist = 40, min.peakwidth = 3, sn.adjust.comp = 1, min.sharpness = 6E3, min.fracobs = .4, do.plot = F
#' )
#' }
#'
#' @export
#'
findcomponents = function(
  macha, .roil = NULL,
  S = 3:8, seed.maxdensity = 1/10, seed.maxdist=4, seed.sn.perpeak = c(Inf, 10, 5, 3, 2), seed.sn.range = 2, seed.sn.adjust = 1, seed.minwidth = 3,
  unrelated.dist = 30, min.peakwidth, sn.adjust.comp = 1, min.sharpness, min.fracobs = 0.4,
  do.plot=F
  ) {
  cat("\nStarting component detection.\n")

  if (is.null(.roil)) .roil =  unique(macha$r$r)
  if (is.numeric(.roil)) {
    roi.job.input = getroi.iter(macha, .roil)
    groups_run = .roil
  } else {
    roi.job.input = .roil
    groups_run = seq_along(.roil)
    }
  ngs = length(.roil)

  cat("\rProgress indicator:", 1, "of", ngs, "(ROIs analyzed)            ")

  components = foreach (roi = roi.job.input, i = icount(), .errorhandling = "pass", .packages=c("macha"), .options.redis=list(chunkSize=50)) %dopar% {
    if (i %% 5 == 0) { cat(paste0("\r", i, " of ", ngs, " ROIs analyzed. (", round(i/ngs*100), "%)              ")) }

    if (do.plot) { plot.roi(roi) }

    seeds = seedpeaks(roi, S = S, seed.maxdensity = seed.maxdensity, seed.maxdist=seed.maxdist, seed.sn.perpeak = seed.sn.perpeak, seed.sn.range = seed.sn.range, seed.sn.adjust = seed.sn.adjust, seed.minwidth = seed.minwidth, do.plot = do.plot)

    if (nrow(seeds) > 0) {
      components = fitandreseed(roi, seeds, unrelated.dist = unrelated.dist, min.peakwidth = min.peakwidth, sn.adjust = sn.adjust.comp, min.sharpness= min.sharpness, min.fracobs = min.fracobs, do.plot=do.plot)

      components = do.call(cbind, components)
      if (is.null(components)) {return(NULL)}

      mz = c(apply(components, 2, function(x) {
        range = c(x["rtmin"], x["rtmax"])
        proi = roi[{ind <- roi[.(c(range[1], range[2])), which=TRUE, roll=TRUE, on="rt"]; ind[1][is.na(ind[1])] = 1; ind[2][is.na(ind[2])] = nrow(roi); (ind[1]+1):ind[2]}]

        sum(matrixStats::rowProds(proi[,.(mz, i)] %>% as.matrix),na.rm=T)/sum(proi$i,na.rm=T)
      }))

      rbind(components, mz, r = roi$r[[1]]) %>% aperm
    }
  }

  ise = sapply(components, function(x) "error" %in% class(x))
  if (any(ise)) {
    macha$c_error = components[ise]; names(macha$c_error) = as.character(groups_run[ise]); warning(paste0(sum(ise), " component detection errors. (Stored in $c.error.)")) } else { macha = macha[!names(macha) == "c_error"] }

  ism = sapply(components, function(x) "matrix" %in% class(x))
  macha$c = data.table(do.call(rbind, components[ism]))
  if (nrow(macha$c) > 0) macha$c[,c:=seq_len(nrow(macha$c))]

  isn = sapply(components, is.null)

  if (any(!(ise | ism | isn))) warning(paste0(sum(!(ise | ism | isn)), " rois discarded. This is a bug."))

  cat("\nFinished component detection.\n\n")
  macha

}
