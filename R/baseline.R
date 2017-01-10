#' Estimate the baseline for each m/z channel
#'
#' \code{baseline} estimates the baseline noise level at each point supplied in \code{roi}
#'
#' Developed to robustly estimate baselines in a variety of situations. Primarily: 1. cases where the baseline is not always observed, as in QE data; 2. Cases where a single ROI does not capture a large enough region to estimate the baseline.
#' To achieve this, ROIs of similar mass are grouped based on \code{ppmwin}. For each ROI group all the signal throughout the dataset is aggregated into an EIC. This EIC roll-joined, such that 0 values are filled with the last observed intensity. From this EIC a baseline is estimated.
#'
#' @param rois Matrix. Rows are observed mass peaks. Columns named "mz", "rt", "g", "i".  Column "g" denotes which ROI the mass peak belongs to.
#' @param features Matrix. Rows are observed mass peaks. Columns named "mz", "rt", "i". This should include all noise for baseline estimation.
#' @param rts Character. The retention times of each scan. Used to convert "rt" above into scan number.
#' @param ppmwin Numeric. The ppm window within which to sum to denote the baseline.
#' @param lambda1,lambda2 Integer. See \code{\link[baseline::baseline.irls]}
#'
#'
#' @return Matrix. Column "b" is the baseline intensity. Column "d" is the fraction of the mass channel had signal in it. Each row corresponds to the same row in the supplied \code{roi}.
#'
#' @seealso \code{\link[baseline::baseline.irls]}
#'
#' @export
#'

baseline = function(macha, ppmwin = 3, lambda1 = 3, lambda2 = 6) {
  cat("\nStarted baselining.")

  featsd = macha$k[,.(mz, s, i)]
  roisd = macha$k[macha$k_r,,on="k"]

  # Aggregate ROIs which are close in mass but possibly separate in retention time.
  cat("\nAggregating mass channels.")
  massdt = roisd[,mean(mz),by=r]

  o = order(massdt$V1)
  m = massdt$V1[o]
  jumps = diff(m) / m[-1] * 1E6

  gs = c(cumsum(c(0,jumps) > ppmwin - 0.5))
  toolarge = sapply(split(m, gs), function(x) { diff(range(x))/mean(x) * 1E6 }) %>% '>'(ppmwin) %>% which

  maxg = max(gs)
  for (gn in toolarge - 1) {
    wtl = which(gs == gn)
    gs[wtl] = maxg + 1 + seq_along(wtl)
    maxg = maxg + length(wtl)
  }

  massdt[,mchan:=as.numeric(factor(gs))[o]]
  roisd[massdt, mchan := mchan, on="r"]

  mchans = roisd[,.(mz = mean(mz)),by=mchan]
  ranges = cbind(mchans$mz-mchans$mz*(ppmwin+0.5)/1E6, mchans$mz + mchans$mz*(ppmwin+0.5)/1E6)


  cat("\nExtracting mass channels.\n")  # Put all ROIs into a matrix and baseline together (~20x faster than individual baselining.)

  setkey(macha$s, "s")
  setkey(featsd, "mz")
  mat = matrix(ncol = nrow(macha$s), nrow = nrow(ranges), dimnames = list(mchans$mchan, macha$s$s))

  nms = nrow(ranges)
  for (j in seq_len(nrow(ranges))) {
    if (j %% 50 == 0) { cat(paste0("\r", j, " of ", nms, " mass channels analyzed. (", round(j/nms*100), "%)              ")) }

    range = ranges[j,]

    #Stack exchange soltuion here: http://stackoverflow.com/questions/40665673/does-data-table-implement-fast-range-subsetting-based-on-binary-search-what-is
    #microbenchmark::microbenchmark(
    #  featsd[{ind <- featsd[.(c(range[1], range[2])), which=TRUE, roll=TRUE]; (ind[1]+1):ind[2]}],
    #  featsd[mz > range[1] & mz < range[2]]
    #  )

    m.inrange = featsd[{ind <- featsd[.(c(range[1], range[2])), which=TRUE, roll=TRUE]; ind[1][is.na(ind[1])] = 1; ind[2][is.na(ind[2])] = nrow(featsd); (ind[1]+1):ind[2]}]

    mat[j,] =  m.inrange[macha$s, sum(i, na.rm=T), roll=T, on="s", mult="all", by=.EACHI]$V1
  }

  cat("\nCalculating baselines.")

  bl = baseline::baseline(mat, lambda1 = lambda1, lambda2 = lambda2)@baseline


  cat("\nFinished baselining.")

  bldt = data.table(b = c(bl), d = rowSums(mat>0), mchan=as.integer(rownames(mat)), s = rep(as.integer(colnames(mat)), each=nrow(mat)))
  macha$k_b = roisd[bldt,.(k,b,d),nomatch=0, on=.(mchan,s)]

  macha
}
