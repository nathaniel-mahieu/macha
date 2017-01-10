

grtgmzcor = function(Nmacha, shaperng = 2, fracobs=0.6) {

  tmpsize = length(Nmacha$m)
  gcounts = Nmacha$m.c[,length(unique(m)) == length(m) & length(m) >= tmpsize*fracobs & diff(range(shape)) < shaperng, by = g]

  ftg = Nmacha$m.c[gcounts[V1==T],,on="g"]


  # Correct RTs
  rtmat = matrix(NA, nrow = length(unique(ftg$g)), ncol = length(Nmacha$m))
  assmat = cbind(as.numeric(factor(ftg$g)), ftg$m)

  rtmat[assmat] = ftg$rtpeak

  maxrt = max(Nmacha$m.c$rtmax) + 1000

  rtmat = rbind(rep(0,length(Nmacha$m)), rep(maxrt, length(Nmacha$m)), rtmat)
  rtdmat =  rtmat - rowMeans(rtmat)
  rtdmat[2,] = rtdmat[nrow(rtdmat),]

  df = data.frame(rt = c(rtmat), rtd = c(rtdmat), g = rep(seq_len(dim(rtdmat)[1]), dim(rtdmat)[2]), m = rep(seq_len(dim(rtdmat)[2]), each = dim(rtdmat)[1]))
  Nmacha$grt = lapply(seq_len(length(Nmacha$m)), function (j) {
    lo  = loess(rtd ~ rt, df, subset = m == j | m==0, degree=2, span = 0.44)

    rts = c(0, Nmacha$m[[j]]$s$rt, maxrt)
    grts = predict(lo, rts)

    data.table(rt = rts, grt = rts+grts)
  })


  Nmacha$gcs = ftg[,.(c,m,g)]

  # Correct mzs
  rtmat = matrix(NA, nrow = length(unique(ftg$g)), ncol = length(Nmacha$m))
  assmat = cbind(as.numeric(factor(ftg$g)), ftg$m)

  rtmat[assmat] = ftg$mz


  maxmz = max(Nmacha$m.c$mz) + 1000
  rtmat = rbind(rep(0,length(Nmacha$m)), rep(maxmz, length(Nmacha$m)), rtmat)
  rtdmat =  rtmat - rowMeans(rtmat)
  rtdmat[2,] = rtdmat[nrow(rtdmat),]

  df = data.frame(rt = c(rtmat), rtd = c(rtdmat), g = rep(seq_len(dim(rtdmat)[1]), dim(rtdmat)[2]), m = rep(seq_len(dim(rtdmat)[2]), each = dim(rtdmat)[1]))
  Nmacha$gmz = lapply(seq_len(length(Nmacha$m)), function (j) {
    lo  = loess(rtd ~ rt, df, subset = m == j | m==0, degree=2, span = 0.44)

    mzr = range(Nmacha$m[[j]]$k$mz)
    mzs = c(0, seq(mzr[1], mzr[2], by = 0.1), maxmz)
    grts = predict(lo, mzs)

    data.table(mz = mzs, gmz = mzs+grts)
  })


  Nmacha
}




corrt = function(rt, .m, grt, uncorrect = F) {
  rtout = vector(mode="numeric", length = length(rt))

  idl = split(seq_along(rt), .m)
  .ml = as.numeric(names(idl))

  for (j in seq_along(idl)) {
    if (uncorrect) rtout[idl[[j]]] = approx(grt[[.ml[j]]][,.(rt, grt)], xout = rt[idl[[j]]])$y
    else rtout[idl[[j]]] = approx(grt[[.ml[j]]][,.(grt, rt)], xout = rt[idl[[j]]])$y
  }

  rtout
}

cormz = function(rt, .m, grt, uncorrect = F) {
  rtout = vector(mode="numeric", length = length(rt))

  idl = split(seq_along(rt), .m)
  .ml = as.numeric(names(idl))

  for (j in seq_along(idl)) {
    if (uncorrect) rtout[idl[[j]]] = approx(grt[[.ml[j]]][,.(mz, gmz)], xout = rt[idl[[j]]])$y
    else rtout[idl[[j]]] = approx(grt[[.ml[j]]][,.(gmz, mz)], xout = rt[idl[[j]]])$y
  }

  rtout
}

