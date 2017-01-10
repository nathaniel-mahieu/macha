plot.roi = function(roi) {
  .pardefault <- par(no.readonly = T)
  par(mfrow = c(2,2), mar=c(5,4,1,1)+0.1)

  plot(roi[,.(rt, i)], type="p")
  lines(roi[,.(rt, ii)], type="l", col="grey")
  lines(roi[,.(rt, bb)], type="l", col="red")
  plot(roi[,.(rt, b)], type="p")
  lines(roi[,.(rt, bb)], type="l", col="grey")

  mmz = mean(roi$mz, na.rm=T)

  plot(roi[,.(rt, mz)], type="p")
  abline(h = mmz + mmz/1E6, col="red")
  abline(h = mmz - mmz/1E6, col="red")
  plot(roi[,.(i, mz)], type="p")
  abline(h = mmz + mmz/1E6, col="red")
  abline(h = mmz - mmz/1E6, col="red")

  par(.pardefault)
  }

getroi = function(macha, roi.id) {
  roi = macha$k[macha$k_r[r==roi.id],,on="k"][macha$k_b,,on="k", nomatch=0]

  region = macha$s[s %in% min(roi$s):max(roi$s)]

  roi = roi[region,,nomatch=0,on="s"][,rt:=NULL]
  foo = roi[region[,.(s,rt)], .(i,b,s, rt), roll=T, on="s"][,ii:=i][,bb:=b][,b:=NULL][,i:=NULL]
  roi = roi[foo,,on="s"][,r:=roi.id]
  setkey(roi, "s")

  roi
  }
