#' Use well behaved peak groups to correct drift in both retention time and mass across samples
#'
#' \code{grtmzcor} Determines global correction curves (one correction per scan) for both mz and rt and builds a mapping between the global retention time (rt.g) and observed retention time (rt)
#' 
#' A loess curve is fit to the difference between the observed value and the mean of all values against the retention time.
#' 
#' @param Nmacha An Nmacha object containing peaks grouped between samples.
#' @param shaperng Peaks within a group must exhibit a range of shape parameters less than this value to be included in the fit.
#' @param	fracobs Peak groups must contain at least this fraction of peaks to be included in the fit.
#' 
#' @return List (an Nmacha object) containing additional list named grt (mapping between grt and rt for each file), gmz (mapping between gmz and mz for each file), gc. 
#' 
#' @examples
#' \dontrun{
#' grtgmzcor(Nmacha, shaperng = 2, fracobs = .4)
#' }
#'
#' @export
#'

grtgmzcor = function(Nmacha, shaperng = 2, fracobs=0.6) {

  needed.obs = length(Nmacha$m) * fracobs
  good.groups = Nmacha$m.c[Nmacha$m.c_g,,on="m.c."][,length(unique(m)) == length(m) & length(m) >= needed.obs & diff(range(shape)) < shaperng, by = g]

  
  ftg = Nmacha$m.c[Nmacha$m.c_g,,on="m.c."][good.groups[V1==T],,on="g"]


  # Correct RTs
  rtmat = matrix(NA, nrow = length(unique(ftg$g)), ncol = length(Nmacha$m))
  assmat = cbind(as.numeric(factor(ftg$g)), ftg$m)
  rtmat[assmat] = ftg$location

  rtmat = rbind(rep(-200,length(Nmacha$m)), rep(-100,length(Nmacha$m)), rtmat, rep(max(Nmacha$m.c$rtmax)+100, length(Nmacha$m)), rep(max(Nmacha$m.c$rtmax)+200, length(Nmacha$m)))
  meanrts = rowMeans(rtmat, na.rm=T)

  Nmacha$grt = lapply (seq_len(ncol(rtmat)), function(coln) {
    df = data.frame(rt = rtmat[,coln], grt = meanrts)
    df = df[!is.na(df$rt),]
    smooth.spline(df$grt, df$rt, df = 2, spar = .6)
    })


  # Correct mzs
  rtmat = matrix(NA, nrow = length(unique(ftg$g)), ncol = length(Nmacha$m))
  assmat = cbind(as.numeric(factor(ftg$g)), ftg$m)
  rtmat[assmat] = ftg$mz


  rtmat = rbind(rep(-10,length(Nmacha$m)),rep(-5,length(Nmacha$m)), rtmat,  rep(max(Nmacha$m.c$mz)+5, length(Nmacha$m)), rep(max(Nmacha$m.c$mz)+10, length(Nmacha$m)))
  meanrts =  rowMeans(rtmat, na.rm=T)

  Nmacha$gmz = lapply (seq_len(ncol(rtmat)), function(coln) {
    df = data.frame(mz = rtmat[,coln], gmz = meanrts)
    df = df[!is.na(df$mz),]
    smooth.spline(df$gmz, df$mz, df = 2, spar = 0.6)
  })


  Nmacha
}



cor.global = function(rts, grt, uncorrect = F) {
  if (!uncorrect) predict(grt, rts)$y
  else approx(grt$data$y, grt$data$x, xout = rts)$y
}

