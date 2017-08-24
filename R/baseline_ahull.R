#' Estimate the baseline for each m/z channel
#'
#' \code{macha.baseline.ahull} estimates the baseline using the convex alpha-hull of the mass trace
#'
#' Developed to robustly estimate baselines in a variety of situations. Primarily: 1. cases where the baseline is not always observed, as in QE data; 2. Cases where a single ROI does not capture a large enough region to estimate the baseline.
#' To achieve this, ROIs of similar mass are grouped based on \code{ppmwin}. For each ROI group all the signal throughout the dataset is aggregated into an EIC. This EIC roll-joined, such that 0 values are filled with the last observed intensity. From this EIC a baseline is estimated.
#' Estimation is performed using th alpha-hull of the trace to identify points to be considered as baseline - points within the variance of the interpolated alpha-hull are included in a spline fit to the baseline.
#'
#' @param macha Macha list containing list items k and s.
#' @param ppmwin Numeric. The ppm window within which to sum to denote the baseline.
#' @param pw.scan Integer. A generous estimate of the peak width at baseline in scans. (Excludes peaks from baselline estimation).
#' @param a Numeric. Generally 2, higher numbers lowers sensitivity to peaks.
#'
#'
#' @return List (a Macha object) containing additional list named k_b containing the calculated baseline intensity at every mass peak. (Note that this refers to peak IDs much like a relational database.)
#'
#' @seealso \code{\link[alphahull::ashape]}
#'
#' @examples
#' \dontrun{
#' macha.baseline.ahull(macha, pw.scan = 20, a = 2, ppmwin = 3)
#' }
#'
#' @export
#'

macha.baseline.ahull = function(macha, pw.scan, a, ppmwin=3, plot.summary=F) {
  cat("\nStarted baselining.")

  featsd = macha$k[,.(mz, s, i)]
  roisd = macha$k[macha$k_r,,on="k"]

  # Aggregate ROIs which are close in mass but possibly separate in retention time.
  cat("\nAggregating mass channels.")
  massdt = roisd[,.(meanmz = mean(mz)),by=r]
  setkey(massdt, meanmz)

  jumps = diff(massdt$meanmz) / massdt$meanmz[-1] * 1E6
  massdt[,gs:=cumsum(c(0,jumps) > ppmwin - 0.5) + 1]

  granges = massdt[,.(grange = diff(range(meanmz))/mean(meanmz)*1E6),by=gs]

  for (gn in granges[grange > ppmwin]$gs) {
    massdt[gs == gn, gs := max(massdt$gs) + 1 + seq_along(gs)]
  }

  massdt[,mchan:=as.numeric(factor(gs))]
  roisd[massdt, mchan := mchan, on="r"]

  mchans = roisd[,.(mz = mean(mz)),by=mchan][,':='(mzmin = mz - mz*(ppmwin+0.5)/1E6, mzmax = mz + mz*(ppmwin+0.5)/1E6)]
  setkey(mchans, mchan)

  dupdt = roisd[,.(dups = sum(duplicated(s)), N=.N),by=mchan]
  cat("\nMass channels with duplicates:", sum(dupdt$dups>0), "- Fraction of mass channels with duplicates:", round(sum(dupdt$dups>0)/nrow(dupdt),2), "- Fraction of duplicated scans in mass channels with duplicates:", round(mean(dupdt[dups>0, dups/N]),2))
  cat("\nExtracting mass channels.")  # Put all ROIs into a matrix and baseline together (~20x faster than individual baselining.)

  setkey(macha$s, "s")
  setkey(featsd, "mz")
  mat = matrix(ncol = max(macha$s$s), nrow = nrow(mchans), dimnames = list(mchans$mchan, macha$s$s))

  cat("\nCalculating baselines.")

  nms = nrow(mchans)


  cat("\nProcessing backend used for foreach(baselines):", getDoParName())
  start  = Sys.time()
  bldt = foreach (j=seq_len(nrow(mchans)), .packages = "macha") %dopar% {
    cat(paste0("\r", j, " of ", nms, " mass channels analyzed. (", round(j/nms*100), "%)              "))

    range = mchans[j,.(mzmin, mzmax)] %>% as.numeric

    m.inrange = fastrangedt(featsd, range, "mz")
    s.inrange = fastrangedt(macha$s, range(m.inrange$s), "s")

    trace = m.inrange[s.inrange, sum(i, na.rm=T), roll=T, on="s", mult="all", by=.EACHI]$V1

    bl = baseline.ahull(x=s.inrange$s, y=trace, a=a, x.var=pw.scan, smooth.n = 5)

    cbind(s = s.inrange$s, b = bl, mchan = mchans[j, mchan], d = 0)
  } %>% do.call(what=rbind) %>% data.table
  cat("\nFinished baselining.", round((Sys.time() - start),1), "minutes.")

  #bldt = data.table(b = c(mat), d = rowSums(mat>0, na.rm=T), mchan=as.integer(rownames(mat)), s = rep(as.integer(colnames(mat)), each=nrow(mat)))
  macha$k_b = roisd[bldt,.(k,b,d),nomatch=0, on=.(mchan,s)]
  macha$k_b[b<0, b:=0]

  macha
}

baseline.ahull = function(x, y, a, x.var, y.var=NULL, smooth.n=5, do.plot = F) {

  if (length(y) < 5) { warning("Not enough points to calculate baseline.  Returning supplied points."); return(y) }

  # Remove isolated, low points
  k = min(ceiling(length(y)*0.25), 20)
  which(y == -zoo::rollmax(-y, k, fill=NA)) %>% { y[.] = rowSums(cbind(y[.-1], y[.+1], na.rm=T))/2 }

  # Smooth
  y.sm = filter_all_the_way_down(y, smooth.n)

  if (is.null(y.var)) y.var = quantile(abs(y-y.sm) %>% { .[.>1] }, 0.8)[[1]]

  # Find alpha-hull and return only bottom of hull
  abp = ahull_bottom(x=x/x.var, y=y.sm/(y.var*10), a=a, do.plot=F)

  # Points to include in baseline estimation: within our variance of the bottom of the hull
  reasonbl = which(abs(y - approx(x[abp], y[abp], x)$y) < y.var/2)


  df = diff(range(x))/x.var/5 # spline df: number of possible features/4
  bb = predict(smooth.spline(x[reasonbl], y[reasonbl], df = df, spar = 0.4)$fit, x)$y

  if (do.plot) {
    roi = data.table(rt = x, ii = y, bb = bb)
    print(ggplot(roi) + geom_line(aes(x=rt, y=ii)) + geom_point(data = roi[abp], aes(x=rt, y=ii), colour="red")+ geom_point(data = roi[reasonbl], aes(x=rt, y=ii), colour="blue", alpha = 0.2) +geom_line(aes(x=rt, y=bb), colour="red"))
    }

  bb
  }

filter_all_the_way_down = function(x, n) {
  if ((n-1) %% 2 != 0) {stop("n must be odd")}

  x2 = rep(NA, length(x))
  for (nn in seq(n, 1, by=-2)) {
    nas = is.na(x2)
    x2[nas] = filter(x,rep(1/nn,nn))[nas]
  }
  x2
}

ahull_bottom = function(x,y, a=3, do.plot=F) { # Almost all of this is dedicated to finding the bottom of the hull - faster way?

  ah = alphahull::ashape(x, y, a)
  dt = data.table(ah$edges)
  if (do.plot) plot(ah)

  #Find x-overlapping edges
  dt[dt$x2 - dt$x1 < 0, ':='(x1=x2, x2=x1, y1=y2, y2=y1, ind1=ind2, ind2=ind1)]
  overlap = ( sign(outer(dt$x1, dt$x2, "-")) != sign(outer(dt$x2, dt$x1, "-")) )
  diag(overlap)= F

  #Find adjacent edges
  adjacent = (outer(dt$ind1, dt$ind2, "==") | outer(dt$ind2, dt$ind1, "==") | outer(dt$ind2, dt$ind2, "==") | outer(dt$ind1, dt$ind1, "=="))

  #Find y-overlapping edges
  dt[dt$y2 - dt$y1 < 0, ':='(x1=x2, x2=x1, y1=y2, y2=y1, ind1=ind2, ind2=ind1)]
  strictlytaller =  sign(outer(dt$y1, dt$y2, "-")) == sign(outer(dt$y2, dt$y1, "-")) & outer(dt$y1, dt$y2, '>=')
  possiblytaller =  sign(outer(dt$y1, dt$y2, "-")) != sign(outer(dt$y2, dt$y1, "-"))

  #Throw away lower overlaps
  too_do = which(overlap & possiblytaller & !adjacent, arr.ind=T)

  #speedups
  dt2 = as.matrix(dt)
  #i.y1 = 4
  #i.y2 = 6
  #i.x1 = 3
  #i.x2 = 5

  nokeep.posstall = numeric()
  if(nrow(too_do) > 0) {
    nokeep.posstall = apply(too_do, 1, function(pair) {
      r1 = dt2[pair[1]]
      r2 = dt2[pair[1]]

      xs = c(r1[3], r1[5], r2[3], r2[5])
      seg = xs[order(xs)[2:3]]



      h1 = sum(r1[4] + (r1[6]-r1[4])/(r1[5]-r1[3]) * (seg - r1[3]))
      h2 = sum(r2[4] + (r2[6]-r2[4])/(r2[5]-r2[3]) * (seg - r2[3]))

      pair[which.max(c(h1, h2))]
    })
  }

  nokeep = c(
    which(overlap & strictlytaller, arr.ind=T)[,1],
    nokeep.posstall
  ) %>% unique

  if (length(nokeep)==0) return(dt[, c("ind1", "ind2")] %>% unlist %>% unique)
  dt[-nokeep, c("ind1", "ind2")] %>% unlist %>% unique
}
