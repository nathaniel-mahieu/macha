warpcombine = function(Nmacha, rt.padding = 10) {
  if (is.null(Nmacha$trace_cache)) stop("Please populate Nmacha$trace_cache first.")

  cat("Starting. ", 0);   lt = Sys.time()
  m.c_mchan = Nmacha$m.c[,.SD[Nmacha$m[[m[1]]]$r_mchan,,on="r",nomatch=0],by="m"]
  #Nmacha$m.c_g[,.N,by="g"][N==6]

  comps = putative = list()
  srate = min(sapply(Nmacha$m, function(x) mean(diff(x$s$rt))))*0.8
  rtouts.g = sapply(Nmacha$m, function(x) range(x$s$rt)) %>% {seq(-srate*2, max(.[,2]) + srate*2, srate)}

  #m.c_g.bak = Nmacha$m.c_g
  #Nmacha$m.c_g = Nmacha$m.c_g[g %in% sample(unique(g), 100)]

  ms = seq_along(Nmacha$m)

  cat("Splitting Up Jobs. ", lt - Sys.time());   lt = Sys.time()
  mchan_m_g = m.c_mchan[Nmacha$m.c_g,,on="m.c."]
  cs.l = split(mchan_m_g, by="g")

  trace_cache.l = split(Nmacha$trace_cache[mchan_m_g[,.(mchan, m, g)] %>% { .[!duplicated(.)] },,on=c("mchan","m")],by="g")
  trace_cache.l.o = match(names(cs.l), names(trace_cache.l))

  cat("Aggregating trace and group information. ", lt - Sys.time());   lt = Sys.time()
  trace_ranges = Nmacha$trace_cache[, .(rtmin = min(rt, na.rm=T), rtmax = max(rt, na.rm=T), mzmin = min(mz, na.rm=T), mzmax = max(mz, na.rm=T)),by=c("m", "mchan")]

  g_ranges = mchan_m_g[,.(rtmin = min(rtmin, na.rm=T), rtmax = max(rtmax, na.rm=T), mzmin = min(mz, na.rm=T), mzmax = max(mz, na.rm=T)),by=c("g", "m")]

  setkey(g_ranges, m, rtmin, rtmax, mzmin, mzmax)
  setkey(trace_ranges, m, rtmin, rtmax, mzmin, mzmax)

  cat("Searching for missing traces. ", lt - Sys.time());  lt = Sys.time()
  group_traces = trace_ranges[g_ranges,
     .(m, g, mchan),
     on = .(rtmin <= rtmax, rtmax >= rtmin, mzmin <= mzmax, mzmax >= mzmin),
     nomatch = F, allow.cartesian=T, mult="all"
     ]
  group_traces = group_traces %>% { .[!duplicated(.)] }

  cat("Splitting up original traces. ", lt - Sys.time());  lt = Sys.time()
  raw_traces.l = split(Nmacha$trace_cache[group_traces,,on=c("mchan", "m"), allow.cartesian = T],by="g")
  raw_traces.l.o = match(names(cs.l), names(raw_traces.l))

  cat("Starting foreach loop. ", lt - Sys.time());  lt = Sys.time()

  ug. = Nmacha$m.c_g$g %>% unique; lug. = length(ug.)
  #g. = "231"; lassign(cs = cs.l[[g.]], trace_cache = trace_cache.l[[g.]], raw_traces = raw_traces.l[[g.]])
  output = foreach (
    cs = cs.l, trace_cache = trace_cache.l[trace_cache.l.o], raw_traces = raw_traces.l[raw_traces.l.o], i = icount(),
    .packages = c("macha", "dtw"), .options.redis=list(chunkSize=10), .noexport = c("Nmacha", "m.c_mchan"),
    .errorhandling = 'pass', .final = function(x) collect_errors(x, names = names(cs.l), .rbind=F)
  ) %dopar% {
    cat("\rWarpcombine: ", round(i/lug.,2), "       ")

    trange = cs[,.(rtmin,rtmax)] %>% range
    trange.lim = c( trange[1]-rt.padding, trange[2]+rt.padding)
    mrange = cs[,(mz/1E6*1) %>% { c(mz - ., mz+.) }] %>% range

    #if (!cs$m %>% unique %>% length == length(Nmacha$m)) warning("Missing peaks from some files.")
    if (cs[,.(m, mchan)] %>% unique %>% nrow > cs$m %>% unique %>% length) warning("Multiple mass channel from some files")


    rois = by(cs, cs$m, function(rows) {
      rows = unique(rows[,.(m,mchan)])

      r = trace_cache[rows,,on=.(mchan, m)]
      r = r[rt > trange.lim[1] & rt < trange.lim[2]]

      # Merge multiple ROIs
      ds = duplicated(r$s)
      if (any(ds)) {
        mergek = function(x) x[,.(rt = mean(rt), i = sum(i), mz = mean(mz), b = sum(b), polarity=mean(polarity))]

        rm = r[,mergek(.SD),by=s]
        r = rbind(r[!s %in% rm$s], rm, fill=T)
      }

      r
    }) %>% do.call(what=rbind)


    # Find ROIs in files where no peak was detected.
    missing.m = which(!(ms %in% cs$m)) %>% unique


    rois = lapply (missing.m, function(m.) { #cat(m., " ")
      newroi = raw_traces[m == m.]
      newroi = newroi[rt > trange.lim[1] & rt < trange.lim[2]]

      if (any(duplicated(newroi$s))) { stop("Ambiguous ROI choice.") }

      newroi
    }) %>% { c(., list(rois)) } %>% do.call(what=rbind)

    # Ask: We have ROIs but do those encompass all the necessary region?
    # Only necessary if peaks were not found for a file.
    rs = rois[,{tmp = range(rt); list(min=tmp[1],max=tmp[2])},by="m"]
    look.for.missing = cs[,.(rtmin,rtmax)] %>% range %>% {rs$min < .[1] + 2 & rs$max > .[2] - 2}
    #if (any(!look.for.missing)) warning("Missing some portion of the expected regions")

    if (length(unique(rois$m)) < 2) stop("Trace was found in only one file.")

    # Ask: We have ROIs but do we have ROIs for all files?
    if (unique(rois$m) %>% length < length(ms)) { warning(paste(sep=": ", "Missing ROI for some files", cs$g[[1]], unique(rois$m))) }

    # Ask: We have ROIs but did those ROIs miss any peaks?

    wfi = which(findInterval(rtouts.g, range(cs[,.(rtmin-5, rtmax+5)]))==1) %>% { c(.[1]-1, ., .[length(.)] + 1) }
    rtouts = rtouts.g[wfi] %>% {.[!is.na(.)]}

    # Align EICs
    #rois.bak = copy(rois)

    #rois = copy(rois.bak)
    eic.l = split(rois[,.(rt = rt, intensity = i, scan = s, m = m)], by="m")
    eic.m = sapply(simplify="array", eic.l, function (x) {
      cbind(rt = rtouts, intensity = approx(x$rt, x$i, xout = rtouts)$y %>% {./max(., na.rm=T); .[is.na(.)] = min(.,na.rm=T); . })
    })


    for (i in seq_len(dim(eic.m)[3])[-2]) { #cat(m.)
      m. = names(eic.l)[i] %>% as.numeric
      dt = dtw(
        c(rep(0, 20), eic.m[,"intensity",i], rep(0,20)),
        c(rep(0, 20), eic.m[,"intensity",2], rep(0, 20)),
        open.end = F, open.begin = F, window.type = "none", keep=F, step.pattern = asymmetricP2)

      #plot(dt$index1, dt$index2)
      i1 = dt$index1 %>% { .[20:(length(.)-20)] - 20 } %>% { .[.<1] = 1; .[.>length(eic.m[,"rt",i])] = length(eic.m[,"rt",i]); . }
      i2 = dt$index2 %>% { .[20:(length(.)-20)] - 20 } %>% { .[.<1] = 1; .[.>length(eic.m[,"rt",i])] = length(eic.m[,"rt",i]); . }
      sms = smooth.spline(eic.m[i1,"rt",i], eic.m[i2,"rt",1], spar = .5)
      #plot(eic.m[m.,i1,"rt"], eic.m[1,i2,"rt"]); lines(sms)

      #plot(predict(s,eic.m[1,,"rt"])$y,eic.m[1,,"intensity"], type="l", col = "black")
      #lines(eic.m[2,,"rt"],eic.m[2,,"intensity"], type="l", col = "red")

      cs[m==m., ':='(location=predict(sms, location)$y, rtpeak=predict(sms,rtpeak)$y, rtmin=predict(sms,rtmin)$y, rtmax=predict(sms,rtmax)$y)]
      rois[m==m., rt:= predict(sms,rt)$y]
    }
    #ggplot(rois) + geom_line(aes(x = rt, y = i, color = factor(m)))

    tf = function(rt, i, rto=rtouts) approx(x=rt, y=i, xout=rtouts, rule=1)$y
    roicor = rois[,.(rt = rtouts, i = as.integer(tf(rt, i)), mz = tf(rt, mz), b = as.integer(tf(rt, b))),by="m"]

    #grid.arrange(
    #  ggplot(rois.bak) + geom_line(aes(x=rt, y = i, colour = factor(m))) + facet_wrap(~m),
    #  ggplot(rois) + geom_line(aes(x=rt, y = i, colour = factor(m))) + facet_wrap(~m)
    #)


    #rtsd = roicor[,.(mean = mean(i, na.rm=T), sd = sd(i, na.rm=T)),by="rt"]
    #grid.arrange(
    #  ggplot(rtsd) + geom_line(aes(x=rt, y = mean)) + geom_ribbon(aes(x = rt, ymin = mean - sd, ymax = mean + sd), alpha = 0.2),
    #  ggplot(roicor) + geom_line(aes(x=rt, y = i, colour = factor(m)))
    #  )


    list(composite_groups = roicor[,g:=cs$g[1]], putative_peaks = cs[,.(m, m.c., g, location, scale, shape, factor, baseline, mz, rtmin, rtmax)])
  }

  cat("\n\n")
  Nmacha$warpcombine_error = output$error
  Nmacha$composite_groups = data.table::rbindlist(lapply(output$list, '[[', 'composite_groups'))
  Nmacha$putative_peaks = data.table::rbindlist(lapply(output$list, '[[', 'putative_peaks'))

  return(Nmacha)
}

# Parallelize!
warpcombine_peaks = function(Nmacha, refit_constraints_range = c(1, 0.5, 0.1, Inf), peak_group_bw = 1, min_peaks = 2) {
  # Choose Peaks
    ug. = Nmacha$putative_peaks$g %>% unique

    putative_peaks.l = split(Nmacha$putative_peaks, by="g")
    composite_groups.l = split(Nmacha$composite_groups, by="g")
    composite_groups.l.o = match(names(putative_peaks.l), names(composite_groups.l))

    l.groups = length(putative_peaks.l)

    output = foreach (
      g. = as.numeric(names(putative_peaks.l)), pps = putative_peaks.l, cgs = composite_groups.l[composite_groups.l.o], i = icount(),
      .packages = "macha", .options.redis=list(chunkSize=10), .noexport = c("Nmacha", "putative_peaks.l", "composite_groups.l"),
      .errorhandling = 'pass', .final = function(x) collect_errors(x, names = names(putative_peaks.l))
    ) %dopar% {
      cat("\rWarpcombine peaks: ", round(i/l.groups,2), "      ")

      setkey(cgs, rt)
      cgs.l = split(cgs, by="m")

      for (r in seq_len(nrow(pps))) { #Find actual peak RTs
        x = pps[r]
        rt = cgs.l[[x$m]][abs(rt - x$location)<1][i == max(i,na.rm=T)]$rt
        pps[r,location := rt]
      }



      #k = pps[order(location)][,{ sum(diff(location)>4)+1 }]
      #d = pps[,.(location, mz)]

      #d_clust <- Mclust(as.matrix(d), G=(k-1):(k+2))
      #k <- dim(d_clust$z)[2]
      #plot(d_clust)


      d = density(pps$location, bw = peak_group_bw)
      lms = d$y %>% { .[.<0.001] = 0; . } %>% localMaxima
      k = length(lms)

      #multid = sapply(seq(0.1, 3, 0.2), function(bw) {
      #  density(pps$location, bw = bw)$y %>% { .[.<0.001] = 0; . } %>% localMaxima %>% length
      #  })

      if (k == nrow(pps)) {
        km = list(cluster = seq_len(nrow(pps)), centers = cbind(location=pps$location))
      } else {
        #km = kmeans(pps[,.(location, mz)], k)
        km = kmeans(pps[,.(location, mz)], centers = cbind(d$x[lms], mean(pps$mz,na.rm=T)))
      }
      pps[,k:=km$cluster]



      if (F) {
        plot(d)
        grid.arrange(
          ggplot(pps) + geom_point(aes(x = location, y=mz, colour = factor(m)), alpha = 0.5, size = 3),
          ggplot(pps) + geom_segment(aes(x = rtmin, xend = rtmax, y=mz, yend=mz, colour = factor(m)), alpha = 0.3, size = 2),
          ggplot(cgs) + geom_line(aes(x=rt, y = i, colour = factor(m))),
          ggplot(pps) + geom_point(aes(x = location, y=mz, colour = k %>% factor), size = 3, alpha = 0.5),
          ncol=2
        )
      }


      # Refit Peaks
      #components = fitseeds(eic, seeds = seeds, unrelated.dist = 30, const.upper = const.upper, const.lower = const.lower, do.plot = do.plot)

      keep.center = table(km$cluster) %>% { .[order(as.numeric(names(.)))] }  %>% { . >= min_peaks }
      pps = pps[k %in% which(keep.center)]

      seeds = km$centers[keep.center,"location"]


      okm = order(seeds)
      rt.ts = cgs$rt %>% unique %>% {.[order(.)]}

      seeds = diff(seeds) %>% { cbind(seeds, seeds-c(max(./2),./2), seeds+c(./2,max(./2))) }
      seeds[] = seeds %>% { approx(rt.ts, seq_along(rt.ts), xout = .)$y } %>% round
      seeds[is.na(seeds[,2]),2] = 0
      seeds[is.na(seeds[,3]),3] = length(rt.ts)
      seeds = as.matrix(seeds)

      setkey(cgs, m, rt)

      seed.constraints = pps[,.(loc = mean(location), scl = mean(scale), shp = mean(shape), fct = mean(factor), bsl = mean(baseline)),by="k"]
      seed.constraints = pps[,.(loc = mean(location), scl = 1, shp = 0, fct = mean(factor), bsl = mean(baseline)),by="k"]
      #seed.constraints$shape = mean(shape)

      refit_constraints_range = c(1, 0.25, 0.25, Inf)
      seed.constraints.var.m = rep(refit_constraints_range, each = nrow(seed.constraints))
      consts = seed.constraints[,.(loc, scl, shp, fct)] %>% as.matrix %>% { list(.-seed.constraints.var.m, . + seed.constraints.var.m) }

      consts[[1]][,"fct"] = 0
      #consts[[1]][,"loc"] = max(0, consts[[1]][,"loc"])
      #consts[[2]][,"loc"] = min(consts[[2]][,"loc"], max(rt.ts, na.rm=T))


      cgs[is.na(i), i:= as.integer(min(cgs$i,na.rm=T)/2)]

      components = lapply(split(cgs, by ="m"), function(x) {
        comps = fitseeds(x[,.(rt=rt, i = i, ii = i, b = b, bb = b)] %>% as.matrix, seeds = seeds, unrelated.dist = 30, const.lower = consts[[1]], const.upper = consts[[2]], do.plot = F)

        mzs = sapply(seq_len(ncol(comps)), function(c) {
          sum(curvemany(c(comps[1:3,c],1,0), x$rt) %>% { ./sum(.,na.rm=T) } * x$mz, na.rm=T)
          })

        comps %>%  rbind(., m = rep(x$m[[1]], ncol(.)), mz = mzs, ccclust = g., cc = okm)

      }) %>% do.call(what = cbind) %>% aperm %>% data.table

      components
    }
  Nmacha$cc = output$list
  Nmacha$cc_error = output$errors

  Nmacha
}
