#' Estimate the baseline for each m/z channel
#'
#' \code{baseline} estimates the baseline noise level at each point supplied in \code{roi}
#'
#' Developed to robustly estimate baselines in a variety of situations. Primarily: 1. cases where the baseline is not always observed, as in QE data; 2. Cases where a single ROI does not capture a large enough region to estimate the baseline.
#' To achieve this, ROIs of similar mass are grouped based on \code{ppmwin}. For each ROI group all the signal throughout the dataset is aggregated into an EIC. This EIC roll-joined, such that 0 values are filled with the last observed intensity. From this EIC a baseline is estimated.
#'
#' @param macha Macha list containing list items k and s.
#' @param ppmwin Numeric. The ppm window within which to sum to denote the baseline.
#' @param lambda1,lambda2 Integer. See \code{\link[baseline::baseline.irls]}
#'
#'
#' @return List (a Macha object) containing additional list named k_b containing the calculated baseline intensity at every mass peak. (Note that this refers to peak IDs much like a relational database.)
#' 
#' @seealso \code{\link[baseline::baseline.irls]}
#'
#' @examples
#' \dontrun{
#' baseline(macha, ppmwin = 3, lambda1 = 6, lambda2 = 7)
#' }
#'
#' @export
#'

baseline = function(macha, ppmwin = 3, lambda1 = 3, lambda2 = 6, plot.summary=F) {
  cat("\nStarted baselining.")

  featsd = macha$k[,.(mz, s, i)]
  roisd = macha$k[macha$k_r,,on="k"]

  # Aggregate ROIs which are close in mass but possibly separate in retention time.
  cat("\nAggregating mass channels.")
  massdt = roisd[,.(meanmz = mean(mz)),by=r]
  setkey(massdt, meanmz)
  
  jumps = diff(massdt$meanmz) / massdt$meanmz[-1] * 1E6
  massdt[,gs:=cumsum(c(0,jumps) > ppmwin - 0.5) + 1]
  
  granges = massdt[,.(grange = diff(range(meanmz))/mean(meanmz)*1E6),by=gs]
  
  for (gn in granges[grange > ppmwin]$gs) {
    massdt[gs == gn, gs := max(massdt$gs) + 1 + seq_along(gs)]
  }

  massdt[,mchan:=as.numeric(factor(gs))]
  roisd[massdt, mchan := mchan, on="r"]

  mchans = roisd[,.(mz = mean(mz)),by=mchan][,':='(mzmin = mz - mz*(ppmwin+0.5)/1E6, mzmax = mz + mz*(ppmwin+0.5)/1E6)]
  setkey(mchans, mchan)

  cat("\nExtracting mass channels.\n")  # Put all ROIs into a matrix and baseline together (~20x faster than individual baselining.)

  setkey(macha$s, "s")
  setkey(featsd, "mz")
  mat = matrix(ncol = nrow(macha$s), nrow = nrow(mchans), dimnames = list(mchans$mchan, macha$s$s))

  nms = nrow(mchans)
  for (j in seq_len(nrow(mchans))) {
    if (j %% 50 == 0) { cat(paste0("\r", j, " of ", nms, " mass channels analyzed. (", round(j/nms*100), "%)              ")) }

    range = mchans[j,.(mzmin, mzmax)] %>% as.numeric

    #Stack exchange soltuion here: http://stackoverflow.com/questions/40665673/does-data-table-implement-fast-range-subsetting-based-on-binary-search-what-is
    #microbenchmark::microbenchmark(
    #  featsd[{ind <- featsd[.(c(range[1], range[2])), which=TRUE, roll=TRUE]; (ind[1]+1):ind[2]}],
    #  featsd[mz > range[1] & mz < range[2]]
    #  )

    m.inrange = fastrangedt(featsd, range, "mz")

    mat[j,] =  m.inrange[macha$s, sum(i, na.rm=T), roll=T, on="s", mult="all", by=.EACHI]$V1
  }

  cat("\nCalculating baselines.")

  bl = baseline::baseline(mat, lambda1 = lambda1, lambda2 = lambda2)@baseline

  cat("\nFinished baselining.")

  bldt = data.table(b = c(bl), d = rowSums(mat>0), mchan=as.integer(rownames(mat)), s = rep(as.integer(colnames(mat)), each=nrow(mat)))
  macha$k_b = roisd[bldt,.(k,b,d),nomatch=0, on=.(mchan,s)]

  try({if (is.character(plot.summary)) {
    pdf(file="test_baseline.pdf", width = 12, height=7)
    
    dt = macha$k[macha$k_r,,on="k"][macha$k_b,,on="k"]
    
    rlens = dt[,.N,by=r]
    
    roisd[,.(dups = sum(duplicated(s))),by=mchan]$dups %>% hist(main="Multiple mass peaks per scan in each mass channel for baselining.")
    roisd[,.(dups = sum(duplicated(s))),by=mchan][dups > 0]$dups %>% hist(main="Multiple mass peaks per scan in each mass channel for baselining.")
    
    for (.r in sample(rlens[N > quantile(N, 0.75)]$r, 20)) {
      { dt[r == .r] %>% ggplot() + geom_line(aes(x = s, y = i)) + geom_line(aes(x = s, y = b), colour = "red") + ggtitle("Example Long Baseline", paste0("ROI ID: ", i)) } %>% print
    }
    for (.r in sample(rlens[N < quantile(N, 0.75)]$r, 20)) {
      { dt[r == .r] %>% ggplot() + geom_line(aes(x = s, y = i)) + geom_line(aes(x = s, y = b), colour = "red") + ggtitle("Example Short Baseline", paste0("ROI ID: ", i)) } %>% print
    }
    dev.off()
  }})
  
  macha
}
