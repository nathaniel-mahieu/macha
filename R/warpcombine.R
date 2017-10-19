warpcombine = function(Nmacha, rt.padding = 10) {
  cat("Building trace cache.\n")
  trace_cache = make_trace_cache(Nmacha)

  m.c_mchan = Nmacha$m.c[,.SD[Nmacha$m[[m[1]]]$r_mchan,,on="r",nomatch=0],by="m"]
  #Nmacha$m.c_g[,.N,by="g"][N==6]

  comps = putative = list()
  srate = min(sapply(Nmacha$m, function(x) mean(diff(x$s$rt))))*0.8
  rtouts.g = sapply(Nmacha$m, function(x) range(x$s$rt)) %>% {seq(-srate*2, max(.[,2]) + srate*2, srate)}

  #m.c_g.bak = Nmacha$m.c_g
  #Nmacha$m.c_g = Nmacha$m.c_g[g %in% sample(unique(g), 100)]

  ug. = Nmacha$m.c_g$g %>% unique; lug. = length(ug.)
  output = foreach (
    g. = ug., i = icount(),
    .packages = "macha", .options.redis=list(chunkSize=50),
    .errorhandling = 'pass', .final = function(x) collect_errors(x, names = ug., .rbind=F)
  ) %dopar% {
    cat("\rWarpcombine: ", round(i/lug.,2), "       ")

    cs = m.c_mchan[Nmacha$m.c_g[g==g.],,on="m.c."]

    trange = cs[,.(rtmin,rtmax)] %>% range
    trange.lim = c( trange[1]-rt.padding, trange[2]+rt.padding)
    mrange = cs[,(mz/1E6*1) %>% { c(mz - ., mz+.) }] %>% range

    #if (!cs$m %>% unique %>% length == length(Nmacha$m)) warning("Missing peaks from some files.")
    if (cs[,.(m, mchan)] %>% unique %>% nrow > cs$m %>% unique %>% length) warning("Multiple mass channel from some files")


    rois = by(cs, cs$m, function(rows) {
      rows = unique(rows[,.(m,mchan)])

      r = trace_cache[[rows$m[1]]][rows[,.(mchan)], on="mchan"]
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
    missing.m = which(!(seq_along(Nmacha$m) %in% cs$m)) %>% unique


    rois = lapply (missing.m, function(m.) { #cat(m., " ")
      rs = Nmacha$m[[m.]]$r[,.(r = r, r1i = findInterval(minrt, trange), r2i = findInterval(maxrt, trange), m1i = findInterval(maxmz, mrange), m2i =  findInterval(minmz, mrange))][(r1i != r2i | r1i == 1 & r2i == 1) & (m1i != m2i | m1i == 1 & m2i == 1)]
      mchans = Nmacha$m[[m.]]$r_mchan[r %in% rs$r]

      newroi = trace_cache[[m.]][mchan %in% mchans$mchan][rt > trange.lim[1] & rt < trange.lim[2]]
      if (any(duplicated(newroi$s))) { stop("Ambiguous ROI choice.") }

      newroi
    }) %>% { c(., list(rois)) } %>% do.call(what=rbind)

    # Ask: We have ROIs but do those encompass all the necessary region?
    # Only necessary if peaks were not found for a file.
    rs = rois[,{tmp = range(rt); list(min=tmp[1],max=tmp[2])},by="m"]
    look.for.missing = cs[,.(rtmin,rtmax)] %>% range %>% {rs$min < .[1] + 2 & rs$max > .[2] - 2}
    #if (any(!look.for.missing)) warning("Missing some portion of the expected regions")

    # Ask: We have ROIs but do we have ROIs for all files?
    if (unique(rois$m) %>% length < length(Nmacha$m)) { warning(paste(sep=": ", "Missing ROI for some files", g., m.)) }

    # Ask: We have ROIs but did those ROIs miss any peaks?

    wfi = which(findInterval(rtouts.g, range(cs[,.(rtmin, rtmax)]))==1) %>% { c(.[1]-1, ., .[length(.)] + 1) }
    rtouts = rtouts.g[wfi]

    # Align EICs
    rois.bak = copy(rois)

    rois = copy(rois.bak)
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
  Nmacha$composite_groups = lapply(output$list, '[[', 'composite_groups') %>% do.call(what=rbind)
  Nmacha$putative_peaks = lapply(output$list, '[[', 'putative_peaks') %>% do.call(what=rbind)

  return(Nmacha)
}


warpcombine_peaks = function(Nmacha, refit_constraints_range = c(1, 0.5, 0.1, Inf)) {
  # Choose Peaks

    ug. = Nmacha$putative_peaks$g %>% unique
    ug. = sample(ug., 20);
    lug. = length(ug.)
    #ug. = 1665
    output = foreach (
      g. = ug., i = icount(),
      .packages = "macha", .options.redis=list(chunkSize=50),
      .errorhandling = 'pass', .final = function(x) collect_errors(x, names = ug.)
    ) %dopar% {
      cat("\rWarpcombine peaks: ", round(i/lug.,2), "      ")

      #g. =  Nmacha$putative_peaks[,.N,by="g"][N>5]$g %>% sample(1)
      #g. = sample(Nmacha$putative_peaks$g, 1)
      pps = Nmacha$putative_peaks[g==g.]
      cgs = Nmacha$composite_groups[g==g.]

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


      d = density(pps$location, bw = 0.3)
      lms = d$y %>% { .[.<0.001] = 0; . } %>% localMaxima
      k = length(lms)

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
      seeds = km$centers[,"location"]
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
        fitseeds(x[,.(rt=rt, i = i, ii = i, b = b, bb = b)] %>% as.matrix, seeds = seeds, unrelated.dist = 30, const.lower = consts[[1]], const.upper = consts[[2]], do.plot = F) %>%
          rbind(., m = rep(x$m[[1]], ncol(.)), ccclust = g., cc = okm)
      }) %>% do.call(what = cbind) %>% aperm %>% data.table

      components
    }
  Nmacha$cc = output$list
  Nmacha$cc_error = output$errors

  Nmacha
}
