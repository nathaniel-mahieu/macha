#' Estimate the baseline for each m/z channel
#'
#' \code{macha.baseline.ahull} estimates the baseline using the convex alpha-hull of the mass trace
#'
#' Developed to robustly estimate baselines in a variety of situations. Primarily: 1. cases where the baseline is not always observed, as in QE data; 2. Cases where a single ROI does not capture a large enough region to estimate the baseline.
#' To achieve this, ROIs of similar mass are grouped based on \code{ppmwin}. For each ROI group all the signal throughout the dataset is aggregated into an EIC. This EIC roll-joined, such that 0 values are filled with the last observed intensity. From this EIC a baseline is estimated.
#' Estimation is performed using th alpha-hull of the trace to identify points to be considered as baseline - points within the variance of the interpolated alpha-hull are included in a spline fit to the baseline.
#'
#' @param macha macha.
#' @param pw.scan Numeric. Scales the x-dimension making the baseline less sensitive to rises on the scale of pw.scan.  For example, when multiple peaks elute the baseline may disappear (rise) for 3 or three peak widths.  Recommended: 2-3 * peakwidth near baseline in scans.
#' @param a Numeric. Scales both the X and Y dimenstions as above.  Recommended: 2.5
#' @param long.n Integer. Length of trace before extra filters for noisy channels is applies.  Filtering excludes points which could be due to ion suppression.
#'
#' @return data.table with columns bb and v.
#'
#' @seealso \code{\link{alphahull::ashape}}
#'
#'
#' @export
#'


macha.baseline.ahull = function(macha, pw.scan, a=2.5, long.n = 150) {
  cat("\nStarted baselining.")
  cat("\nProcessing backend used for foreach(baselines):", getDoParName())
  cat("\nFilling mass trace gaps.")

  traces = macha$k[macha$k_r,,on="k",nomatch=0][macha$r_mchan,,on="r",nomatch=0][,.(k,s,i,mchan)]
  setkey(traces, "mchan","s")
  traces = traces[,.SD[macha$s[,.(s)],,roll=T, rollends=F, on="s", nomatch=0],by="mchan"]
  traces[duplicated(k), k:= NA]

  setkey(traces, "mchan","s")
  trace.l = split(traces, by="mchan")


  start  = Sys.time()
  nms = length(trace.l)
  output = foreach (
    trace=trace.l, i = icount(), .packages = "macha",
    .options.redis=list(chunkSize=50),
    .errorhandling = 'pass', .final = function(x) collect_errors(x, names = names(trace.l))
    ) %dopar% {
    cat(paste0("\r", i, " of ", nms, " mass channels analyzed. (Fraction: ", round(i/nms, 4), ")              "))

    bl = baseline.ahull(x=trace$s, y=trace$i, a=a, x.var=pw.scan, smooth.n = 5, long.n = 150, do.plot=T)

    data.table(s = trace$s, b = bl[,"bb"], v = bl[,"v"], mchan = trace$mchan, d = rep(0,length(trace$s)))
  }
  cat("\nFinished baselining.", round(((Sys.time() - start)/60),1), "minutes.")

  bldt = output$list
  macha$k_b_error = output$error

  #bldt = data.table(b = c(mat), d = rowSums(mat>0, na.rm=T), mchan=as.integer(rownames(mat)), s = rep(as.integer(colnames(mat)), each=nrow(mat)))
  macha$k_b = traces[bldt,.(k,b,d,v),nomatch=0, on=.(mchan,s)][!is.na(k)]
  macha$k_b[b<0, b:=0]

  macha
}

variance_est = function(y.d) {
  n = 31

  if (length(y.d) < (n+3)){
    n = floor((length(y.d)-2)/2)
    }


  y.d.q = zoo::rollapply(abs(y.d) %>% { .[.<1] = NA; . }, n, quantile, 0.8, fill=NA, na.rm=T)
  fnna = which(!is.na(y.d.q))[1]
  lnna = which(!is.na(y.d.q)) %>% tail(n=1)
  y.d.q[1:(fnna-1)] = y.d.q[fnna]
  y.d.q[(lnna+1):length(y.d)] = y.d.q[lnna]
  y.d.q[is.na(y.d.q)] = mean(y.d.q, na.rm=T)

  ydq.f = filter(y.d.q, rep(1/n, n))
  fnna = which(!is.na(ydq.f))[1]
  lnna = which(!is.na(ydq.f)) %>% tail(n=1)
  ydq.f[1:(fnna-1)] = ydq.f[fnna]
  ydq.f[(lnna+1):length(y.d)] = ydq.f[lnna]
  ydq.f %>% as.numeric
  }

baseline.ahull.bak = function(x, y, a, x.var, y.var=NULL, smooth.n=11, do.plot = F) {

  if (length(y) < 5) { warning("Not enough points to calculate baseline.  Returning supplied points."); return(y) }
  if (smooth.n > length(y)) smooth.n = length(y)
  smooth.n = smooth.n %>% { . + (. %% 2-1) }

  # Remove isolated, low points
  k = min(ceiling(length(y)*0.25), 20)
  y = which(y == -zoo::rollmax(-y, k, fill=NA)) %>% { y[.] = rowSums(cbind(y[.-1], y[.+1], na.rm=T))/2; y }

  # Smooth
  y.sm = filter_all_the_way_down(y, smooth.n)
  y.d = y - y.sm

  y.v = variance_est(y.d)

  if (!is.null(y.var)) y.v = rep(y.v, length(y.sm))

  # Find alpha-hull and return only bottom of hull
  #abp = ahull_bottom(x=x/x.var, y=y.sm/(y.v*10), a=a, do.plot=T)
  abp = ahull_bottom(x=x/x.var, y=y.sm/(mean(y.v)*10), a=a, do.plot=T)

  # Points to include in baseline estimation: within our variance of the bottom of the hull
  reasonbl = which(abs(y - approx(x[abp], y[abp], x)$y) < y.v/2)


  npts = length(reasonbl)
  if (npts > 4) {
    df = diff(range(x))/x.var/5 # spline df: number of possible features/4
    bb = predict(smooth.spline(x[reasonbl], y[reasonbl], df = df, spar = 0.4)$fit, x)$y
  } else if (npts == 1) {
    bb = rep(y[reasonbl],length(x))
  } else if (npts == 0) {
    bb = rep(NA, length(x))
  } else {
    bb = approx(x[reasonbl], y[reasonbl], xout = x)$y
    }

  if (do.plot) {
    roi = data.table(rt = x, ii = y, bb = bb)
    print(ggplot(roi) + geom_line(aes(x=rt, y=ii)) + geom_point(data = roi[abp], aes(x=rt, y=ii), colour="red")+ geom_point(data = roi[reasonbl], aes(x=rt, y=ii), colour="blue", alpha = 0.2) +geom_line(aes(x=rt, y=bb), colour="red"))
    }

  cbind(bb=bb, v = y.v)
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

  overlap_taller = function(r1, r2) {
    xs = c(r1[3], r1[5], r2[3], r2[5])
    seg = xs[order(xs)[2:3]]

    h1 = sum(r1[4] + (r1[6]-r1[4])/(r1[5]-r1[3]) * (seg - r1[3]))
    h2 = sum(r2[4] + (r2[6]-r2[4])/(r2[5]-r2[3]) * (seg - r2[3]))

    which.max(c(h1, h2))
  }

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
      r1 = dt2[pair[1],]
      r2 = dt2[pair[2],]

      pair[overlap_taller(r1, r2)]
    })
  }

  nokeep = c(
    which(overlap & strictlytaller, arr.ind=T)[,1],
    nokeep.posstall
  ) %>% unique

  if (do.plot) {
    plot(ah)

    points(dt[nokeep,.(x1, y1)],col="blue")
    points(dt[nokeep,.(x2, y2)],col="blue")
    points(dt[-nokeep,.(x1, y1)],col="red")
    points(dt[-nokeep,.(x2, y2)],col="red")
  }


  if (length(nokeep)==0) {
    remaining = dt
  } else {
    remaining = dt[-nokeep]
  }

  remaining[,':='(end = F, start=F, remove = F)]

  ends = c(remaining[remaining[,x2] %>% which.min,ind2], remaining[remaining[,x1] %>% which.max,ind1])
  remaining[ind1 == ends[2], end := T]
  remaining[ind2 == ends[1], start := T]

  if (sum(remaining$start) > 1) {
    inty = overlap_taller(remaining[start==T,][1] %>% unlist, remaining[start==T,][2] %>% unlist)
    remaining[start==T][inty][,remove:=T]
    }
  if (sum(remaining$end) > 1) {
    inty = overlap_taller(remaining[end==T,][1]%>% unlist, remaining[end==T,][2] %>% unlist)
    remaining[end==T][inty][,remove:=T]
  }

  if (do.plot) {
    plot(ah)

    points(remaining[,.(x1, y1)],col="red")
    points(remaining[,.(x2, y2)],col="red")
  }

  if (nrow(remaining)==0) return(dt[, c("ind1", "ind2")] %>% unlist %>% unique)
  remaining[, c("ind1", "ind2")] %>% unlist %>% unique
}
