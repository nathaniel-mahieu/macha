rawdata = function(file, splits = NULL) {
  cat("\n\nLoading and merging raw data from", file, "\n")

  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("Package 'mzR' needed for this function to work. Please install it.", call. = FALSE)
    }

  library(mzR)

  ramp = mzR:::rampOpen(file)
  data = mzR:::rampRawData(ramp)
  mzR:::rampClose(ramp)

  if (!is.null(splits)) {
    nsplits = nrow(splits)

    seglengths = diff(c(data$scanindex, length(data$mz)))
    seg = rep(rep(seq(nsplits), ceiling(length(data$rt)/nsplits))[1:length(seglengths)], seglengths)
    keepmat = outer(splits[,1], data$mz, "<") & outer(splits[,2], data$mz, ">") & outer(seq(nsplits), seg, "==")
    keepvec = matrixStats::colAnys(keepmat)

    newindices = cumsum(keepvec)[data$scanindex+1]-1
    newindices = newindices[seq(from = 1, to = length(data$scanindex), by = nsplits)]

    newis = seq(from = 1, to = length(data$rt), by = nsplits)

    data$mz = data$mz[keepvec]
    data$intensity = data$intensity[keepvec]
    data$scanindex = as.integer(newindices)
    data$rt = data$rt[newis]
    data$acquisitionNum = data$acquisitionNum[newis]
    data$polarity = data$polarity[newis]
    }

  list(
    s = data.table(data.frame(s = data$acquisitionNum, rt = data$rt, polarity = data$polarity)),
    k = data.table(data.frame(mz = data$mz, i = data$intensity, s = rep(data$acquisitionNum, diff(c(data$scanindex, length(data$mz)))), k = seq_along(data$mz))),
    metadata = list(file=file, samplenumber=NA, sampletype = factor(levels = c("pooledqc", "sample", "blankinjection", "blankextraction", "blankchromatography")), rundate = NA, description = NA)
    )
  }
