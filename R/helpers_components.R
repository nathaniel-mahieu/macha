plot.components = function(macha, roin) {

  peaks = macha$c[r==roin]
  roi = getroi(macha, roin)

  plot(roi[,.(rt,i)], type="l", main=paste0("ROI Number: ", roin, ". ", nrow(peaks), " peaks."))
  for (i in seq_len(nrow(peaks))) {
    lines(roi$rt, curvemany(unlist(peaks[i,.(location, scale, shape, factor, baseline)]), roi$rt), col = "blue")
  }
  lines(roi$rt, roi$bb, col="red")

  span = diff(range(roi$i, na.rm=T))
  text(peaks$rtpeak, peaks$intpeak, paste0("Scale: ", round(peaks$scale,1)), pos = 3, offset = 1, col="darkgrey")
  text(peaks$rtpeak, peaks$intpeak, paste0("Shape: ", round(peaks$shape,1)), pos = 3, offset = 0, col="darkgrey")
  text(peaks$rtpeak, peaks$intpeak, paste0("SN: ", round(peaks$intpeak/peaks$baseline,1)), pos = 3, offset = -1, col="darkgrey")
}
