plot.roi = function(roi) {
  .pardefault <- par(no.readonly = T)
  par(mfrow = c(2,2), mar=c(5,4,1,1)+0.1)

  plot(roi[,.(rt, i)], type="p", main = roi$r[1])
  lines(roi[,.(rt, ii)], type="l", col="grey")
  lines(roi[,.(rt, bb)], type="l", col="red")
  plot(roi[,.(rt, b)], type="p")
  lines(roi[,.(rt, bb)], type="l", col="grey")
  lines(roi[,.(rt, i)], type="l", col="blue")

  mmz = mean(roi$mz, na.rm=T)

  plot(roi[,.(rt, mz)], type="p")
  abline(h = mmz + mmz/1E6, col="red")
  abline(h = mmz - mmz/1E6, col="red")
  plot(roi[,.(log10.i = i%>%log10, mz)], type="p")
  abline(h = mmz + mmz/1E6, col="red")
  abline(h = mmz - mmz/1E6, col="red")

  par(.pardefault)
  }


getroi = function(macha, roi.id) {
  if (is.null(macha[["roi_cache"]])) {
    warning("Set roi_cache for faster processing 'makeroicache()'.")

    macha$roi_cache = makeroicache(macha)
  }

  roi = macha$roi_cache[r == roi.id]
  setkey(roi, "s")

  roi
  }


#' Compile ROIs
#'
#' \code{makeroicache} Compiles ROIs rolling missing values and aggregating baselines and other calculated information.
#' 
#' @param macha macha.
#' 
#' @return data.table
#'
#' @export
#'
makeroicache = function(macha) {
  macha$roi_cache = macha$k[macha$k_r,,on="k"]
  
  if (!is.null(macha$k_b)){
    macha$roi_cache = macha$roi_cache[macha$k_b,,on="k", nomatch=0]
  }
  
  setkey(macha$roi_cache, r, s)
  rolledmacha = macha$roi_cache[,.SD[macha$s,,roll=T, rollends=F, on="s", nomatch=0],by="r"]
  
  macha$roi_cache = macha$roi_cache[rolledmacha[,.(s,rt,r,polarity,ii = i, bb = b, vv = v)],,on=c("r", "s")]
  setkey(macha$roi_cache, "r", "s")
  macha$roi_cache[,.(s, k, r, rt, i, mz, b, d, v, ii, bb, vv, polarity)]
  }
