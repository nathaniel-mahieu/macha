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
  macha,
  unrelated.dist = 30, min.peakwidth, sn.adjust.comp = 1, min.sharpness, min.fracobs = 0.4,
  do.plot=F
  ) {
  cat("\nStarting component detection.\n")
  if (!"seeds" %in% names(macha)) { stop("Findcomponents requires peak seeds in macha$seeds.") }
  
  
  groups_run = unique(macha$seeds$r)
  components = foreach (r.r = groups_run, .errorhandling = "pass", .packages=c("macha"), .options.redis=list(chunkSize=50)) %dopar% {
    cat("\r", r.r)
    
    roi = getroi(macha, r.r)
    seeds = macha$seeds[r==r.r,1:3] %>% as.matrix
    
    components = fitandreseed(roi, seeds, unrelated.dist = unrelated.dist, min.peakwidth = min.peakwidth, sn.adjust = sn.adjust.comp, min.sharpness= min.sharpness, min.fracobs = min.fracobs, do.plot=do.plot)
    
    components = do.call(cbind, components)
    if (is.null(components)) { return(NULL) }
    
    mz = c(apply(components, 2, function(x) {
      range = c(x["rtmin"], x["rtmax"])
      proi = roi[{ind <- roi[.(c(range[1], range[2])), which=TRUE, roll=TRUE, on="rt"]; ind[1][is.na(ind[1])] = 1; ind[2][is.na(ind[2])] = nrow(roi); (ind[1]+1):ind[2]}]
      
      sum(matrixStats::rowProds(proi[,.(mz, i)] %>% as.matrix),na.rm=T)/sum(proi$i,na.rm=T)
    }))
    
    rbind(components, mz, r = roi$r[[1]]) %>% aperm
  }
  
  ise = sapply(components, function(x) "error" %in% class(x))
  if (any(ise)) {
    macha$c_error = components[ise]; names(macha$c_error) = as.character(groups_run[ise]); warning(paste0(sum(ise), " component detection errors. (Stored in $c.error.)")) } else { macha = macha[!names(macha) == "c_error"] 
    }
  
  ism = sapply(components, function(x) "matrix" %in% class(x))
  macha$c = data.table(do.call(rbind, components[ism]))
  if (nrow(macha$c) > 0) macha$c[,c:=seq_len(nrow(macha$c))]
  
  isn = sapply(components, is.null)
  
  if (any(!(ise | ism | isn))) warning(paste0(sum(!(ise | ism | isn)), " rois discarded. This is a bug."))
  
  cat("\nFinished component detection.\n\n")
  macha
}
