#' Split up peaks based on scan and ppm distance.
#'
#' \code{rectroi} Splits peaks based on mass and scan spacing. Conceptually, rectangles are drawn around groups of peaks and if rectangles do not overlap peaks are split into distinct groups.
#'
#' Functional as a fast, first pass ROI grouping.
#'
#' @param k data.table. Three columns: mz, s, k and (optional) g. s is scan, k is a unique integer ID for each peak, g is the initial group assignment.
#' @param ppm numeric. Mass distance in ppm to consider groups of peaks distinct. (If using a ROI refinement like \code{\link[webtrace]} a generous value is acceptable.) Suggested: 2*"instrument ppm"
#' @param scan integer. Scan distance to consider groups of peaks distinct. Consider the number of unobserved values to infer were missing. Suggested: 4.
#'
#' @return Integer vector. Each integer is a group assignment. In the same order as supplied k.
#'
#' @seealso \code{\link[webtrace]}
#'
#' @export
#'

rectroi = function(k, ppm, scan, plot.summary = F) {
  #k: mz, s, k
  if (any(duplicated(k$k))) { stop("Peak IDs not unique (column 'k').") }

  k = copy(k)
  korder = k$k

  if ("g" %in% colnames(k)) {
    message("Starting from supplied groups in k$g")
  } else {
    k[,g:=1]
  }

  keeprunning = T

  prevgn = 0
  currgn = 1
  while (prevgn != currgn) { #OPTIONAL: Speed up by caching completed groups
    prevgn = currgn
    cat("\rNumber of Groups:", currgn, "       ")

    setkey(k, g, mz)
    x = k$mz

    ass = c(0, which(diff(x)/x[-1] * 1E6 > ppm), length(x))
    inds = as.integer(cut(seq_along(x), breaks = ass))
    k[,g:= as.numeric(factor(paste(g, inds)))]


    setkey(k, g, s)
    x = k$s

    ass = c(0, which(diff(x) > scan), length(x))
    inds = as.integer(cut(seq_along(x), breaks = ass))
    k[,g:= as.numeric(factor(paste(g, inds)))]

    currgn = max(k$g)
    cat("\rGroup Number:", currgn, "       ")
    }

  ####
  #Plots
  try({
  if (is.character(plot.summary)) {
    ppmrange = k[,.(n = .N, mz = mean(mz), s = mean(s), ppm = diff(range(mz))/mean(mz)*1E6, s = diff(range(s)), dups = sum(duplicated(s)), dupfrac = 2*sum(duplicated(s)) / length(s)),by="g"]

    pdf(paste0(plot.summary,"_rectroi.pdf"), height=7, width = 12) #PDF
    
    hist(ppmrange$ppm, breaks = 200, main="ppm range of groups")
    hist(ppmrange$s, breaks = 200, main = "scan range of groups")
    hist(ppmrange$dupfrac %>% log10, main="fraction of peaks whcih are duplicated per group")

    ggplot(ppmrange[n>5 & dupfrac > 0.1]) + geom_point(aes(x = s, y = ppm)) + facet_wrap(~cut(dupfrac,4)) + ggtitle("Groups broken out by fraction of duplicated peak and plotted as ppm-span vs scan-span of group.")

    for (i in 1:5) {
      kg = k[g==ppmrange[dupfrac > .30 & n > 30]$g %>% sample(1)]
      mzmin = min(kg$mz) - 0.01
      mzmax = max(kg$mz) + 0.01
      smin = min(kg$s) - 40
      smax = max(kg$s) + 40
      { ggplot(k[mz > mzmin & mz < mzmax & s > smin & s < smax]) + geom_point(aes(y = mz, x = s, colour = duplicated(s))) + geom_text(aes(x = median(s), y = median(mz), label = paste(unique(g), collapse = ", ")), check_overlap = T) + ggtitle("Group with > 30 observations and > 30% duplicated.") + coord_cartesian(xlim = c(smin, smax), ylim = c(mzmin, mzmax)) + theme(legend.position = "none") +
          geom_segment(aes(yend =mzmin, y = mzmin + mzmin * ppm/1E6, xend = smin+5+0.5*scan, x = smin+5+0.5*scan)) + 
          geom_segment(aes(yend =mzmin+ mzmin * ppm/1E6*0.5, y = mzmin + mzmin * ppm/1E6*0.5, xend = smin+5, x = smin+5+scan)) 
      } %>%
        print
    }

    for (i in 1:5) {
      kg = k[g==ppmrange[n > 30]$g %>% sample(1)]
      mzmin = min(kg$mz) - 0.01
      mzmax = max(kg$mz) + 0.01
      smin = min(kg$s) - 40
      smax = max(kg$s) + 40
      { ggplot(k[mz > mzmin & mz < mzmax & s > smin & s < smax]) + geom_point(aes(y = mz, x = s, colour = factor(g))) + geom_text(aes(x = median(s), y = median(mz), label = g), check_overlap = T) + 
          ggtitle("Random group with > 30 observations.") + coord_cartesian(xlim = c(smin, smax), ylim = c(mzmin, mzmax)) + theme(legend.position = "none") + 
          geom_segment(aes(yend =mzmin, y = mzmin + mzmin * ppm/1E6, xend = smin+5+0.5*scan, x = smin+5+0.5*scan)) + 
          geom_segment(aes(yend =mzmin+ mzmin * ppm/1E6*0.5, y = mzmin + mzmin * ppm/1E6*0.5, xend = smin+5, x = smin+5+scan)) 
        } %>% print
    }
    
    mzs = ppmrange[s>100,.(floor(mz/2), floor(s/50))] %>% table %>% { which(. == max(.), arr.ind=T) }
    


    ks = k[mz > mzs[1]*2 - 0.5*5 & mz < mzs[1]*2 +0.5*5 & s > mzs[2]*50 - 100 & s < mzs[2]*50 + 100]
    rectdt = k[ks[,.SD[.N > 1],by=g][!duplicated(g),.(g)], on="g"][,.(ymin = min(mz), ymax = max(mz), xmin = min(s), xmax = max(s)),by=g]
    
    { ks %>%
      ggplot() + geom_point(aes(y = mz, x = s), size = 0.1) + geom_rect(data=rectdt, fill = NA, aes(colour = factor(g, levels = sample(g)), xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), size = 0.1) +
      coord_cartesian(xlim = c(min(ks$s), max(ks$s)), ylim = c(min(ks$mz), max(ks$mz))) + theme(legend.position="none") + 
      ggtitle("Overview of detected rect groups") } %>% print
   
    
    ks = k[mz > mzs[1]*2 - 2*0.5 & mz < mzs[1]*2 +2*0.5 & s > mzs[2]*50 - 100 & s < mzs[2]*50 + 100]
    rectdt = k[ks[,.SD[.N > 1],by=g][!duplicated(g),.(g)], on="g"][,.(ymin = min(mz), ymax = max(mz), xmin = min(s), xmax = max(s)),by=g]
    
    { ks %>%
      ggplot() + geom_point(aes(y = mz, x = s), size = 0.1) + geom_rect(data=rectdt, size = 0.1, fill = NA, aes(colour = factor(g, levels = sample(g)), xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
      coord_cartesian(xlim = c(min(ks$s), max(ks$s)), ylim = c(min(ks$mz), max(ks$mz))) + theme(legend.position="none") + 
      ggtitle("Overview of detected rect groups") } %>% print
    
    
     ks = k[mz > mzs[1]*2 - 0.5 & mz < mzs[1]*2 +0.5 & s > mzs[2]*50 - 100 & s < mzs[2]*50 + 100]
    rectdt = k[ks[,.SD[.N > 1],by=g][!duplicated(g),.(g)], on="g"][,.(ymin = min(mz), ymax = max(mz), xmin = min(s), xmax = max(s)),by=g]
    
    { ks %>%
      ggplot() + geom_point(aes(y = mz, x = s), size = 0.1) + geom_rect(data=rectdt, size = 0.1, fill = NA, aes(colour = factor(g, levels = sample(g)), xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
      coord_cartesian(xlim = c(min(ks$s), max(ks$s)), ylim = c(min(ks$mz), max(ks$mz))) + theme(legend.position="none")  + 
      ggtitle("Overview of detected rect groups") } %>% print
    
    dev.off() #/PDF
  }})
  #/Plots
  ####

  k[,.(k,g)]
  }



#' Split up peaks based on weighted graph of nearby mass observations
#'
#' \code{webtrace} Splits peaks based on mass and scan spacing. A graph is formed by linking masses of nearest m/z values within a scan window.  The mincut of the graph which segregates conflicting ROIs (two mass peaks in a single scan) is determined.  The graph is iteratively cut until no conflicts remain.
#'
#' Functional as a slower, robust ROI detection sensitive to shifts in mass. Necessary for high resolution and high-variance mass data.
#'
#' @param k data.table. Three columns: mz, s, k, g. s is scan, k is a unique integer ID for each peak, g is the initial group assignment.
#' @param scan integer. Scan window within to look for similar mass peaks. Suggested: 10.
#' @param neighbors integer. Within each window, for each mass, how many nearest neighbors to link. Suggested: 3
#'
#' @return Integer vector. Each integer is a group assignment. In the same order as supplied k.
#'
#' @seealso \code{\link[rectroi]}
#'
#' @export
#'

webtrace = function(k, scan, neighbors = 3, plot.summary=F) {
  cat("\rWebtrace group:", k$g[[1]])
  k = copy(k)
  korder = k$k
  setkey(k, s, mz)

  if (scan < 2) { stop("Parameter scan too small.") }

  #Group to each peak
  matches = lapply(unique(k$s), function(s.i) {
    subk = copy(k[abs(s-s.i) <= scan])

    d.mz = abs(outer(subk$mz, subk$mz, "-"))
    d.s = abs(outer(subk$s, subk$s, "-")) > 0
    d.centeronly = abs(outer(rep(1, length(subk$s)), subk$s, "*") - s.i) <= scan-1

    d.mz[!d.s] = Inf
    diag(d.mz) = Inf

    matches = lapply(seq_len(neighbors), function(x) {
      d.tf = d.mz == rowMins(d.mz) & d.mz < Inf & d.centeronly
      mats = which(d.tf, arr.ind = T)
      d.mz[mats] <<- Inf

      mats[] = as.character(subk$k[c(mats)])
      mats
    }) %>% do.call(what = rbind)
  }) %>% do.call(what = rbind)


  colnames(matches) = c("v1", "v2")
  matchdt = data.table(matches)

  matchdt[,v1.s := k[match(v1,k), s]]
  matchdt[,v1.mz := k[match(v1,k), mz]]
  matchdt[,v2.s := k[match(v2,k), s]]
  matchdt[,v2.mz := k[match(v2,k), mz]]
  matchdt[,weight := abs((v1.mz - v2.mz) * 1E3)]
  suppressWarnings(matchdt[,weight := max(weight)-weight + 1])

  # Trim links with high ppm
  #matchdt = matchdt[matchdt$weight > quantile(matchdt$weight, .01)]
  matchdt = matchdt[weight > min(weight)]
  

  #Try Min Cut Approaches
  g = igraph::graph.data.frame(matchdt[,.(as.character(v1), as.character(v2), weight)], directed = F, vertices = as.character(k$k))

  while(T) {
    # More resolution necessary?
    o.s = outer(k$s, k$s, "==")
    o.s[upper.tri(o.s, T)] = F
    storesolv_gs = storesolv = which(o.s, arr.ind = T)
    storesolv_gs[] = storesolv[] = k$k[storesolv] %>% as.character
    
    comps = igraph::membership(igraph::components(g))
    storesolv_gs[] = comps[match(storesolv_gs, names(comps))] %>% as.character
    storesolv = storesolv[storesolv_gs[,1] == storesolv_gs[,2],,drop=F]
    
    cat("\rWebtrace group:", k$g[[1]], "-", nrow(storesolv), "conflicts remaining.      ")
    
    if (nrow(storesolv) < 1) {break}

    mincuts = lapply(seq_len(nrow(storesolv)), function(row) {
      igraph::max_flow(g, storesolv[row,1], storesolv[row,2])
    })

    mincutt = lapply(mincuts, function(x) {
      x$cut %>% sapply(as.numeric)
    }) %>% unlist %>% table
    
    delnum = max(3, floor(length(mincutt)/5), sum(mincutt > .5*nrow(storesolv)))

    g = igraph::delete.edges(g, head(as.numeric(names(mincutt[order(mincutt, decreasing = T)])), n = delnum))
    }

  ####
  #Plots
  try({if (is.character(plot.summary)) {
    pdf(paste0(plot.summary, "_webtrace.pdf"), width = 12, height = 7)
    
    ids = igraph::membership(igraph::components(g))[match(as.numeric(names(igraph::membership(igraph::components(g)))), k$k)] %>% unname
    k[,g := ids]

    #{ plot(g, vertex.label = NA, vertex.size = 3, edge.arrow.size = 0, main="Graph of peak relations") } %>% print
    { ggplot(k) + geom_point(aes(y = mz, x = s, colour = duplicated(s))) + geom_text(aes(x = median(s), y = median(mz), label = g), check_overlap = T) + ggtitle("Multiple masses in a single scan.") } %>% print
    { ggplot(k) + geom_point(aes(y = mz, x = s, colour = g %>% factor)) + geom_text(aes(x = median(s), y = median(mz), label = g), check_overlap = T) + ggtitle("Peak Group Assignments") + theme(legend.position="none") } %>% print
    { ggplot(k) + geom_point(aes(y = mz, x = s, colour = g %>% factor)) + geom_text(aes(x = median(s), y = median(mz), label = g), check_overlap = T) +
        geom_segment(data = matchdt, aes(x = v1.s, xend = v2.s, y = v1.mz, yend = v2.mz), position = "jitter") + ggtitle("Peak assignments with web pre trimming") + theme(legend.position="none") } %>% print

    dev.off()
    }})
  #/Plots
  #####
  mem = igraph::membership(components(g))
  return(data.table(k = as.numeric(names(mem)), g = unname(mem)))
  }




#' Mass trace detection first applying rectroi and then webroi to resolve conflicts
#'
#' \code{rectwebtrace} Splits peaks based on mass and scan spacing.
#'
#' Functional as a slower, robust ROI detection sensitive to shifts in mass. Necessary for high resolution and high-variance mass data.
#'
#' See other functions for parameter definitions.
#'
#' @param ppm.rect numeric.
#' @param scan.rect integer.
#' @param scan.web integer.
#' @param neighbors.web integer.
#' @return Integer vector. Each integer is a group assignment. In the same order as supplied k.
#'
#' @seealso \code{\link[rectroi]} \code{\link[webtrace]}
#'
#' @export
#'

rectwebtrace = function (k, ppm.rect = 8, scan.rect = 4, scan.web = 20, neighbors.web = 3, plot.summary=F) {
  k = copy(k)

  cat("Running rectroi.\n")
  foo = rectroi(k, ppm.rect, scan.rect, plot.summary = plot.summary)
  k[foo, g:=g, on="k"]
  cat("\nCompleted rectroi.\n")

  k[,webg := g]
  dups = k[,sum(duplicated(s)),by="g"][V1 > 0]

  cat("Running webtrace.\n")
  if (is.character(plot.summary)) { print.these = sample(seq_len(nrow(dups)), 5) } else { print.these = 0 }
  for (i in seq_len(nrow(dups))) {
    dg = dups[i]$g
    ksub = k[g == dg]

    cat("                           Resolving conflicting group", i, "of", nrow(dups), "with webtrace.", nrow(ksub), "peaks in group.           ")
    
    if (i %in% print.these) {plot.summary2 =paste0(plot.summary, "_", dg) } else {plot.summary2 = F }
    
    foo = webtrace(ksub, scan.web, neighbors.web, plot.summary = plot.summary2)

    #p = profvis({
    #  webtrace(ksub, scan.web, neighbors.web)
    #  })
    
    k[foo, webg := i.g, on="k"]
    i = i+1
  }
  cat("\nCompleted webtrace.")
  
  k[,r:= as.numeric(factor(paste(g, webg)))]


  k[,.(k,r)]
  }

