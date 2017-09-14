#' Estimate the baseline for each m/z channel
#'
#' \code{macha.baseline.ahull} estimates the baseline using the convex alpha-hull of the mass trace
#'
#' Developed to robustly estimate baselines in a variety of situations. Primarily: 1. cases where the baseline is not always observed, as in QE data; 2. Cases where a single ROI does not capture a large enough region to estimate the baseline.
#' To achieve this, ROIs of similar mass are grouped based on \code{ppmwin}. For each ROI group all the signal throughout the dataset is aggregated into an EIC. This EIC roll-joined, such that 0 values are filled with the last observed intensity. From this EIC a baseline is estimated.
#' Estimation is performed using th alpha-hull of the trace to identify points to be considered as baseline - points within the variance of the interpolated alpha-hull are included in a spline fit to the baseline.
#'
#' @param x Numeric. Integer scan
#' @param y Numeric. Intensity
#' @param a Numeric. Sets the distance between points for them to be considered different groups.  Practically, sets how sensitive the baseline is to peaks. Higher is less sensitive/confused by peaks.
#' @param x.var Numeric. Scales the x-dimension. Set to approximately a peak width * 4. Again, sets the sensitivity to peaks.
#' @param y.var Numeric. Scales the y-dimension. Set to approximately the signal noise intensity. Leave NULL for point-wise calculation (slower, better).
#'
#' @param smooth.n Integer, odd. Smoothing window applied before baselining.
#' @param long.n Integer. Length of trace before extra filters for noisy channels is applies.  Filtering excludes points which could be due to ion suppression.
#'
#'
#' @return data.table with columns bb and v.
#'
#' @seealso \code{\link[alphahull::ashape]}
#'
#'
#' @export
#'

baseline.ahull = function(x, y, a, x.var, y.var=NULL, smooth.n=11, long.n = Inf, do.plot = F) {

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
  abp = ahull_bottom(x=x/x.var, y=y.sm/(y.v*10), a=a, do.plot=do.plot)
  #abp = ahull_bottom(x=x/x.var, y=y.sm/(mean(y.v)*10), a=a, do.plot=T)

  disq=c()
  if (length(y) >= long.n) {
    y.smhard = filter_all_the_way_down(y, 101)

    wrm = which(y[abp] < y.smhard[abp]*0.85)
    wrm = wrm[!wrm %in% c(1, length(abp))]
    disq = abp[wrm]
    abp = abp[!abp %in% disq]
    }

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
    roi[,c("hullbot", "included", "removed"):=list(F,F,F)]
    roi[abp,hullbot:=T]
    roi[reasonbl,included:=T]
    roi[disq, removed:=T]

    print(ggplot(roi) +
            geom_line(aes(x=rt, y=ii)) +
            geom_point(data = roi[hullbot==T], aes(x=rt, y=ii), colour="red")+
            geom_point(data = roi[included==T], aes(x=rt, y=ii), colour="blue", alpha = 0.2) +
            geom_line(aes(x=rt, y=bb), colour="red") +
            geom_point(data = roi[removed==T], aes(x=rt, y=ii), colour="green", alpha = 0.2))
    }

  cbind(bb=bb, v = y.v)
  }

