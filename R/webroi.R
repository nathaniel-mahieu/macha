intfactor = function(...) {
  intl <- list(...)

  powers = sapply(intl, function(x) {
    c(small=floor(log10(min(x))), big=ceiling(log10(max(x))), len = length(x))
    })

  if (length(unique(powers["len",])) != 1) stop("Ints must be same length to factor.")
  #if (length(intl) < 2) stop("Must supply at least two vectors")

  overlap = ( sign(outer(powers["small",], powers["big",], "-")) != sign(outer(powers["big",], powers["small",], "-")) )
  overlap[upper.tri(overlap, T)] = F

  overlaps = which(overlap, arr.ind=T)

  maxp = max(powers[1:2,])
  for (r in seq_len(nrow(overlaps))) {
    intl[[overlaps[r,2]]] = intl[[overlaps[r,2]]] + 10^(maxp+1)
    maxp = ceiling(log10(max(intl[[overlaps[r,2]]])))
    }

  x = intl[[1]]
  for (intv in intl[-1]) {
    x = x + intv
    }

  o = order(x)
  ass = which(diff(x[o])>0)
  rep(seq_len(length(ass)+1),c(ass, length(x)) - c(0,ass))[o]
  }

#' Split up peaks based on scan and ppm distance.
#'
#' \code{rectroi} Splits peaks based on mass and scan spacing. Conceptually, rectangles are drawn around groups of peaks and if rectangles do not overlap peaks are split into distinct groups.
#'
#' Functional as a fast, first pass ROI grouping.
#'
#' @param k data.table. Three columns: mz, s, k and (optional) g. s is scan, k is a unique integer ID for each peak, g is the initial group assignment.
#' @param ppm numeric. Mass distance in ppm to consider groups of peaks distinct. (If using a ROI refinement like \code{\link[webtrace]} a generous value is acceptable.) Suggested: 2*"instrument ppm"
#' @param scan integer. Scan distance to consider groups of peaks distinct. Consider the number of unobserved values to infer were missing. Suggested: 6.
#'
#' @return Integer vector. Each integer is a group assignment. In the same order as supplied k.
#'
#' @seealso \code{\link[webtrace]}
#'
#' @export
#'

rectroi = function(k, ppm, scan) {
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
    #inds = as.integer(cut(seq_along(x), breaks = ass))
    inds = rep(seq_len(length(ass)+1),c(ass, length(x)) - c(0,ass))
    #k[,g:= as.numeric(factor(paste(g, inds)))]
    k[,g:= as.numeric(intfactor(g, inds))]


    setkey(k, g, s)
    x = k$s

    ass = c(0, which(diff(x) > scan), length(x))
    #inds = as.integer(cut(seq_along(x), breaks = ass))
    inds = rep(seq_len(length(ass)+1),c(ass, length(x)) - c(0,ass))
    #k[,g:= as.numeric(factor(paste(g, inds)))]
    k[,g:= as.numeric(intfactor(g, inds))]

    currgn = max(k$g)
    cat("\rGroup Number:", currgn, "       ")
    }

  k[,.(k,r=g)]
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

webtrace = function(k, scan, neighbors = 3) {
  cat("\rWebtrace group:", k$r[[1]], "starting.                                              ")
  k = copy(k)
  korder = k$k
  setkey(k, s, mz)

  if (scan < 2) { stop("Parameter scan too small.") }

  #Group to each peak
  matches = lapply(unique(k$s), function(s.i) {
    subk = k[abs(s-s.i) <= scan]

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

    cat("\rWebtrace group:", k$r[[1]], "-", nrow(storesolv), "conflicts remaining.      ")

    # Non-ideal time savings
    storesolv.n = nrow(storesolv)
    if (storesolv.n < 1) {
      break
    } else if (storesolv.n < 60) {
      resolve_subset = seq_len(storesolv.n)
    } else {
      resolve_subset = sample.int(storesolv.n, 60)
    }

    mincuts = lapply(resolve_subset, function(row) {
      igraph::max_flow(g, storesolv[row,1], storesolv[row,2])
    })

    mincutt = lapply(mincuts, function(x) {
      x$cut %>% sapply(as.numeric)
    }) %>% unlist %>% table

    delnum = max(5, floor(length(mincutt)/5), sum(mincutt > .5*nrow(storesolv)))

    g = igraph::delete.edges(g, head(as.numeric(names(mincutt[order(mincutt, decreasing = T)])), n = delnum))
    }

  mem = igraph::membership(components(g))
  return(data.table(k = as.numeric(names(mem)), r = unname(mem)))
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

rectwebtrace = function (k, ppm.rect = 8, scan.rect = 4, scan.web = 20, neighbors.web = 3) {
  k = copy(k)

  cat("Running rectroi.\n")
  foo = rectroi(k, ppm.rect, scan.rect)
  k[foo, g:=g, on="k"]
  cat("\nCompleted rectroi.\n")

  k[,webg := g]
  dups = k[,sum(duplicated(s)),by="g"][V1 > 0]

  cat("Running webtrace.\n")
  if (is.character(plot.summary)) { print.these = sample(seq_len(nrow(dups)), 5) } else { print.these = 0 }
  for (i in seq_len(nrow(dups))) {
    dg = dups[i]$g
    ksub = k[g == dg]

    cat("                      Resolving conflicting group", i, "of", nrow(dups), "with webtrace.", nrow(ksub), "peaks in group.           ")

    foo = webtrace(ksub, scan.web, neighbors.web)

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

