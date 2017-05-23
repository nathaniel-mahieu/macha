#devtools::install_github("nathaniel-mahieu/macha", auth_token = "1fb2171d23be84346ef42d8654440acde687301f") #macha-machine token

library(macha)
library(data.table)

setwd("X:/Nate/2016/05172017_cphB/data_positive")
main_dir = getwd()

folders = c("all_data")
if (!exists("main_dir")) main_dir = getwd()

for (folder in folders) {
  cat("Working in folder", main_dir, "\n")
  cat("Starting folder", folder, "\n")

  setwd(file.path(main_dir, folder))
  files = list.files("./", full.names = T, pattern = "\\.mzXML$")

  library(doRedis)
  registerDoRedis("nate", "128.252.109.58")

  #########
  #### Machas
  #########


  #### Local

  library(doParallel)
  cl <- makeCluster(6)
  registerDoParallel(cl)

  cat("Working on Baselines", "\n")

  # ROIs, Baselines
  foreach (file = files, .packages = "macha", .errorhandling = "pass") %dopar% {
    cat(file)

    macha = rawdata(file, rbind(c(0, 400), c(400, 3000)))

    macha = findrois(macha, minlength = 15, ppm = 2, rtwid = 7)

    macha = baseline(macha, ppmwin = 3, lambda1 = 6, lambda2 = 7)

    saveRDS(macha, file = file.path("./", paste0(basename(file), ".macha.rds")))
  }

  cat("Caching findComponents", "\n")

  # Cache findComponents
  foreach (file = files, .errorhandling = "pass") %do% {
    cat("  ", file, "\n")

    macha = readRDS(file.path("./", paste0(basename(file), ".macha.rds")))

    macha = makeroicache(macha)

    .roil = foreach(r = unique(macha$r$r), .packages="macha") %dopar% {
      nextElem(getroi.iter(macha, r))
    }

    saveRDS(.roil, file = file.path("./", paste0(basename(file), ".macha_findcomponentsinput.rds")))
  }
  
  stopCluster(cl)
  #### /Local

  #### Remote

  library(doRedis)
  registerDoRedis("nate", "128.252.109.58")

  cat("Running findComponents", "\n")

  # findComponents
  foreach (file = files, .errorhandling = "pass") %do% {
    cat("  ", file, "\n")

    macha = readRDS(file.path("./", paste0(basename(file), ".macha.rds")))

    .roil = readRDS(file.path("./", paste0(basename(file), ".macha_findcomponentsinput.rds")))

    macha = findcomponents(
      macha, .roil = .roil,
      S = 3:7, seed.maxdensity=1/7, seed.maxdist=4, seed.sn.perpeak =c(Inf, 10, 7, 3, 2.5, 2), seed.sn.range = 3, seed.sn.adjust = 1, seed.minwidth = 4,
      unrelated.dist = 40, min.peakwidth = 3, sn.adjust.comp = 1, min.sharpness = 6E3, min.fracobs = .4, do.plot = F
    )

    saveRDS(macha, file = file.path("./", paste0(basename(file), ".macha.rds")))

    pdf(file = file.path("./", paste0(basename(file), ".macha.pdf")), width = 10, height = 10)
    try({
      for (i in 1:10) { plot.components(macha, sample(macha$r$r, 1)) }
      for (i in 1:10) { plot.components(macha, sample(which(macha$c$r %>% table >1) %>% names %>% as.numeric, 1)) }
    })
    dev.off()

  }

  #### /Remote



  #########
  #### Nmacha
  #########

  #### Local
  files = list.files("./", full.names = T, pattern = "\\.mzXML\\.macha\\.rds$")
  filename_replicates = paste0(strsplit(basename(files[[1]]), "_")[[1]][[1]], "_", basename(folder))

  cat("Working on Nmacha", "\n")


  Nmacha = do.call(what=nmacha, lapply(files, readRDS))

  gl = dengroup.ppm(Nmacha$m.c[,.(mz, rtpeak, intpeak)] %>% as.matrix, ppm = 2, rtwid = 1, minlength = 32)
  Nmacha$m.c[,g:=rep(seq_along(gl),sapply(gl, length))[order(unlist(gl))]]

  Nmacha = grtgmzcor(Nmacha, shaperng = 1, fracobs = .7)

  # Cache corrected retention times
  Nmacha$m.c[,rtpeak.g := corrt(rtpeak, Nmacha$grt[[m[1]]]), by="m"]
  Nmacha$m.c[,rtmin.g := corrt(rtmin, Nmacha$grt[[m[1]]]), by="m"]
  Nmacha$m.c[,rtmax.g := corrt(rtmax, Nmacha$grt[[m[1]]]), by="m"]
  Nmacha$m.c[,mz.g := cormz(mz, Nmacha$gmz[[m[1]]]), by="m"]

  Nmacha$m.r[,minrt.g := corrt(minrt, Nmacha$grt[[m[1]]]), by="m"]
  Nmacha$m.r[,maxrt.g := corrt(maxrt, Nmacha$grt[[m[1]]]), by="m"]
  Nmacha$m.r[,minmz.g := cormz(minmz, Nmacha$gmz[[m[1]]]), by="m"]
  Nmacha$m.r[,maxmz.g := cormz(maxmz, Nmacha$gmz[[m[1]]]), by="m"]
  Nmacha$m.r[,meanmz.g := cormz(meanmz, Nmacha$gmz[[m[1]]]), by="m"]

  # Regroup with corrected mass and retention time
  gl = dengroup.ppm(Nmacha$m.c[,.(mz.g, rtpeak.g, intpeak)] %>% as.matrix, ppm = 1, rtwid = 1, minlength = 1)
  Nmacha$m.c[,g:=rep(seq_along(gl),sapply(gl, length))[order(unlist(gl))]]

  Nmacha$m.c = Nmacha$m.c[!is.na(mz.g + rtpeak.g + rtmin.g + rtmax.g)]

  pdf(file = file.path(paste0(filename_replicates, "_grt.nmacha.pdf")), width = 10, height = 10)
  try({ plotgrt(Nmacha) }); gc()
  dev.off()
  
  #pdf(file = file.path(paste0(filename_replicates, "_gmz.nmacha.pdf")), width = 10, height = 10)
  #try({ plotgmz(Nmacha) }); gc()
  #dev.off()
  
  pdf(file = file.path(paste0(filename_replicates, "_groups.nmacha.pdf")), width = 10, height = 10)
  try({
    for (i in sample(Nmacha$m.c$g %>% table %>% '>'(ceiling(length(files)/3)) %>% which %>% names %>% as.numeric, 10)) { plot.group(Nmacha,i) }
    for (i in sample(Nmacha$m.c$g %>% table %>% '>'(ceiling(length(files))) %>% which %>% names %>% as.numeric, 10)) { plot.group(Nmacha,i) }
  }); gc()
  dev.off()


  library(doParallel)
  cl <- makeCluster(6)
  registerDoParallel(cl)

  maxdriftppm = 1
  maxdriftrt = 1
  min.peaks = 21
  
  ugs = unique(Nmacha$m.c[,.N,by="g"][N>min.peaks, g])

  Nmacha$m = lapply(Nmacha$m, makeroicache)
  setkey(Nmacha$m.r, "m", "r")

  cat("Caching warpgroup Nmacha", "\n")


  # Cache warpgroup
  warpgroup.nmacha_data_l = foreach(ug = ugs, .packages = "macha") %dopar% {
    f = nextElem(warpgroup.nmacha.iter(Nmacha, ug, maxdriftppm, maxdriftrt))
    gc()
    f
  }
  saveRDS(warpgroup.nmacha_data_l, paste0(filename_replicates, ".nmacha_data_l.rds"))

  stopCluster(cl)
  #### /Local

  #### Remote

  library(doRedis)
  registerDoRedis("nate", "128.252.109.58")

  cat("Warpgrouping Nmacha", "\n")
  # Warpgroup
  Nmacha = warpgroup.nmacha(Nmacha, ugs = ugs, warpgroup.nmacha_data_l = warpgroup.nmacha_data_l, sc.aligned.lim = 4, pct.pad = 0.1, min.peaks = min.peaks, maxdriftrt = 1, maxdriftppm = 1, fraccontrib = 0.6, refit.var = c(1.5, 0.25, 0.5), do.plot = F)
  saveRDS(Nmacha, paste0(filename_replicates, ".nmacha.rds"))

  pdf(file = file.path("./", paste0(filename_replicates, "_wgroup1.nmacha.warpgroup.pdf")), width = 10, height = 10)
  try({ for (i in sample(which(table(Nmacha$m.c$g) > ceiling(length(files)*1)))) { plot.wgroup(Nmacha, i) } }); gc()
  dev.off()
  
  pdf(file = file.path("./", paste0(filename_replicates, "_wgroup3.nmacha.warpgroup.pdf")), width = 10, height = 10)
  try({ for (i in sample(which(table(Nmacha$m.c$g) > ceiling(length(files)*3)))) { plot.wgroup(Nmacha, i) } }); gc()
  dev.off()
  
  rm(Nmacha)

  cat("Done")
  #### /Remote

  }

