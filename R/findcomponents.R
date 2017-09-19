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
#' @param unrelated.dist Integer. The distance two seeds must be apart to fit them individually (much faster). Suggested value: large peakwidth * 3
#' @param min.peakwidth numeric. The minimum width a component can be to be retained. In seconds.
#' @param min.sharpness numeric. A lower limit on peak.height/peak.width.  Units Intensity/Second. (Not very diagnostic.)
#' @param min.fractionobs numeric. A lower limit on the fraction of expected observations present in the peak region.
#'
#' @return List (a Macha object) containing additional list named $c containing detected components and the parameters of the skew normal distribution describing that component. (and $c_error)
#'
#'
#' @export
#'
findcomponents = function(
  macha,
  unrelated.dist = 40, min.peakwidth = 3, min.sharpness=2, min.fracobs = 0.4
  ) {
  cat("\nStarting component detection.\n")
  if (!"seeds" %in% names(macha)) { stop("Findcomponents requires peak seeds in macha$seeds.") }

  seed.l = split(macha$seeds, by="r")
  r. = unique(macha$seeds[,.(r)])
  roi.l = split(macha$roi_cache[r.,,on="r",nomatch=0], by="r")
  
  output = foreach (
    roi = roi.l[names(seed.l)], seeds.dt = seed.l,
    .packages="macha", .options.redis=list(chunkSize=50),
    .errorhandling = 'pass', .final = function(x) collect_errors(x, names = names(seed.l))
  ) %dopar% {
    
    seeds = seeds.dt[,1:3] %>% as.matrix
    components = fitandreseed(roi, seeds, unrelated.dist = unrelated.dist, min.peakwidth = min.peakwidth, min.sharpness= min.sharpness, min.fracobs = min.fracobs, do.plot=F) %>% do.call(what=cbind)
    
    if (!is.matrix(components)) return(NULL)
    
    mz = c(apply(components, 2, function(x) {
      range = c(x["rtmin"], x["rtmax"])
      proi = roi[{ind <- roi[.(c(range[1], range[2])), which=TRUE, roll=TRUE, on="rt"]; ind[1][is.na(ind[1])] = 1; ind[2][is.na(ind[2])] = nrow(roi); (ind[1]+1):ind[2]}]
      
      sum(matrixStats::rowProds(proi[,.(mz, i)] %>% as.matrix),na.rm=T)/sum(proi$i,na.rm=T)
    }))
    
    rbind(components, mz, r = roi$r[[1]]) %>% aperm %>% data.table
  }
  
  macha$c = output$list
  macha$c_error = output$error
  
  macha$c[,c:=seq_along(macha$c$mz)]
  
  cat("\nFinished component detection.\n\n")
  macha
}
