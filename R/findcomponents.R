getroi.iter = function(macha, .roil) {
  it <- iter(.roil)

  nextEl = function() {
    getroi(macha, roi.id=nextElem(it))
  }

  obj <- list(nextElem=nextEl)
  class(obj) <- c('ixts', 'abstractiter', 'iter')
  obj
  }

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
  } else {
    roi.job.input = .roil
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
    macha$c_error = components[ise]; names(macha$c_error) = as.character(.roil[ise]); warning(paste0(sum(ise), " component detection errors. (Stored in $c.error.)")) } else { macha = macha[!names(macha) == "c_error"] }

  ism = sapply(components, function(x) "matrix" %in% class(x))
  macha$c = data.table(do.call(rbind, components[ism]))
  if (nrow(macha$c) > 0) macha$c[,c:=seq_len(nrow(macha$c))]

  isn = sapply(components, is.null)

  if (any(!(ise | ism | isn))) warning(paste0(sum(!(ise | ism | isn)), " rois discarded. This is a bug."))

  cat("\nFinished component detection.\n\n")
  macha

}
