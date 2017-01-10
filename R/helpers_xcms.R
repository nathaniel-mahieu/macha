machatoxcmsraw = function(macha) {
  if (any(c("c", "r", "b") %in% names(macha))) warning("machatoxcmsraw() currently only converts the raw data. ROIs, Baselines, and Components are not saved in the retruned object.")

  features = macha$p[macha$s,,nomatch=0,on="s"] %>% as.matrix %>% .[order(.[,"rt"], .[,"mz"]),]

  xr <- new("xcmsRaw")

  xr@env$mz = features[,"mz"]
  xr@env$intensity = features[,"i"]

  xr@scantime = unique(features[,"rt"])
  xr@acquisitionNum = seq_along(xr@scantime)

  scans = as.numeric(as.factor(features[,"rt"])) %>% diff %>% '>'(0) %>% which
  xr@scanindex = c(scans, length(xr@env$mz))

  xr@tic = 0
  xr@env$profile =  matrix(0)

  xr
}
