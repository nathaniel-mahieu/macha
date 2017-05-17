

grtgmzcor = function(Nmacha, shaperng = 2, fracobs=0.6) {

  tmpsize = length(Nmacha$m)
  gcounts = Nmacha$m.c[,length(unique(m)) == length(m) & length(m) >= tmpsize*fracobs & diff(range(shape)) < shaperng, by = g]

  ftg = Nmacha$m.c[gcounts[V1==T],,on="g"]


  # Correct RTs
  rtmat = matrix(NA, nrow = length(unique(ftg$g)), ncol = length(Nmacha$m))
  assmat = cbind(as.numeric(factor(ftg$g)), ftg$m)

  rtmat[assmat] = ftg$rtpeak

  maxrt = max(Nmacha$m.c$rtmax) + 1000

  rtmat = rbind(rep(-100,length(Nmacha$m)), rep(0,length(Nmacha$m)), rtmat, rep(maxrt+100, length(Nmacha$m)))
  meanrts = rowMeans(rtmat, na.rm=T)

  grt = list()
  for (coln in seq_len(ncol(rtmat))) {
    df = data.frame(rt = rtmat[,coln], grt = meanrts)
    lo  = loess(grt ~ rt, df, degree=2, span = 0.1) #grt from rt

    #plot(rtd ~ rt, data = df); lines(lo$x[o], lo$fitted[o], col = "red")

    rtouts = seq(-30, max(Nmacha$m[[coln]]$s$rt)+30, by = mean(diff(Nmacha$m[[coln]]$s$rt))*3)

    grts = predict(lo, rtouts)

    grt = c(grt, list(data.table(rt = rtouts, grt = grts)))
    }


  # Correct mzs
  rtmat = matrix(NA, nrow = length(unique(ftg$g)), ncol = length(Nmacha$m))
  assmat = cbind(as.numeric(factor(ftg$g)), ftg$m)

  rtmat[assmat] = ftg$mz


  maxmz = max(Nmacha$m.c$mz) + 1000
  rtmat = rbind(rep(-100,length(Nmacha$m)),rep(0,length(Nmacha$m)),rtmat,  rep(maxmz+100, length(Nmacha$m)))
  meanrts =  rowMeans(rtmat, na.rm=T)
  #rtdmat[2,] = rtdmat[nrow(rtdmat),]

  gmz = list()
  for (coln in seq_len(ncol(rtmat))) {
    df = data.frame(mz = rtmat[,coln], gmz = meanrts)
    lo  = loess(gmz ~ mz, df, degree=2, span = 0.1) # predicts = rt - mean(rt); rt - prediction = mean(rt)

    #plot(rtd ~ rt, data = df); lines(lo$x[o], lo$fitted[o], col = "red")

    rtouts = seq(-2, max(Nmacha$m[[coln]]$k$mz)+2, by = mean(diff(Nmacha$m[[coln]]$k$mz))*3)

    gmzs = predict(lo, rtouts)

    gmz = c(gmz, list(data.table(mz = rtouts, gmz = gmzs)))
  }

  Nmacha$gcs = ftg[,.(c,m,g)]
  Nmacha$grt = grt
  Nmacha$gmz = gmz
  Nmacha
}



corrt = function(rt, grt, uncorrect = F) {
  if (uncorrect) approx(grt[,.(grt, rt)], xout = rt)$y
  else approx(grt[,.(rt, grt)], xout = rt)$y
}

cormz = function(rt, grt, uncorrect = F) {
    if (uncorrect) approx(grt[,.(gmz, mz)], xout = rt)$y
    else approx(grt[,.(mz, gmz)], xout = rt)$y
}


