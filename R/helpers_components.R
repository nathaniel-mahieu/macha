plot.components = function(macha, roin) {

  peaks = macha$c[r==roin]
  roi = getroi(macha, roin)

  plot(roi[,.(rt,i)], type="l", main=paste0("ROI Number: ", roin, ". ", nrow(peaks), " peaks."))
  for (i in seq_len(nrow(peaks))) {
    lines(roi$rt, curvemany(unlist(peaks[i,.(location, scale, shape, factor, baseline)]), roi$rt), col = "blue")
  }
  lines(roi$rt, roi$bb, col="red")

  span = diff(range(roi$i, na.rm=T))
  try({
    text(peaks$rtpeak, peaks$intpeak, paste0("Scale: ", round(peaks$scale,1)), pos = 3, offset = 1, col="darkgrey")
    text(peaks$rtpeak, peaks$intpeak, paste0("Shape: ", round(peaks$shape,1)), pos = 3, offset = 0, col="darkgrey")
    text(peaks$rtpeak, peaks$intpeak, paste0("SN: ", round(peaks$intpeak/peaks$baseline,1)), pos = 3, offset = -1, col="darkgrey")
  })
}

plot.components2 = function(macha, .r = NULL) {

  roi = getroi(macha, .r)
  peaks = macha$c[r==.r]
  peaks[,label:=numeric()]


  peakcurves = data.table(c=numeric(), rt = numeric(), ii = numeric())
  if (nrow(peaks) > 0) {
    peakcurves = lapply (peaks$c, function(.c) {
      rt = roi[peaks[c==.c, rtmin]-10 < rt & peaks[c==.c, rtmax]+10 > rt]$rt
      data.table(rt = rt, ii = curvemany(unlist(peaks[c==.c,.(location, scale, shape, factor, baseline)]), rt), c = .c)
    }) %>% do.call(rbind, .)

    csamp = sample(peaks$c) %>% as.character
    peakcurves[,c:=factor(c,csamp)]
    peaks[,c:=factor(c,csamp)]

    peaks[,label := paste0("Scale: ",round(scale,1),"\nShape: ",shape%>%round(1), "\nSN: ", round(intpeak/baseline,1))]
    }


  ggplot(roi[,.(rt,ii,bb)], aes(x=rt)) + geom_line(aes(y=ii), alpha = 0.5) + geom_line(aes(y=bb), colour="red", alpha = 0.5, size = 0.5) +
    geom_line(data=peakcurves, aes(x=rt, y=ii, colour = c), size = 0.8) +
    theme(legend.position="none") +
    ggtitle(paste0("ROI Number: ", .r, ". ", nrow(peaks), " peaks.")) +
    geom_label_repel(data = peaks, aes(location, intpeak, colour=c, label=label), nudge_y = 0.25*max(roi$ii,na.rm=T), box.padding = unit(1, "lines"))


  }

filterkeepends = function(x, filt) {
  foo = filter(x, filt)
  foo[is.na(foo)] = x[is.na(foo)]
  foo
}
