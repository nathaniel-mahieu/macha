dengroup = function (features, splits, minlength=1) {
  if (ncol(features) != length(splits)) stop("Incorrect number of splits supplied to dengroup.")

  ns = length(splits)

  l = nrow(features)
  groups = vector(mode = "list", length = l)
  action = vector(mode = 'numeric', length = l)

  groups[[1]] = seq_len(l)
  action[1] = 0
  i = 1
  groupmaxindex = 1

  repeat {
    #if (i %% 100 == 0) { cat("\rProgress indicator:", groupmaxindex - i, "groups remaining to analyze. (This will grow before starting to decrease.)            ") }

    if (i > groupmaxindex) {
     # cat("\rProgress indicator: 0 groups remaining to analyze. (This will grow before starting to decrease.)            \n")
      gl = groups[1:groupmaxindex]
      return(rep(seq_along(gl),sapply(gl, length))[order(unlist(gl))])
    }


    tg = groups[[i]]

    if (length(tg) >= minlength - 1) {

      act = action[i] %% ns + 1
      action[i] = action[i] + 1

      ng = kernelsplit(features[tg, act], splits[act])

      if (max(ng) > 1) {
        newgroups = split(tg, ng) %>% unname
        ngroups = length(newgroups)

        indices_to_replace = c(i,(groupmaxindex+1):(groupmaxindex+ngroups-1))

        groups[indices_to_replace] = newgroups
        action[indices_to_replace] = 0
        groupmaxindex = groupmaxindex + ngroups - 1

      } else if (action[i] > ns) { i = i+1 }

    } else { i = i+1 }
  }
}


kernelsplit = function(v, split) {
  curr.bw = split

  n = c(2^floor(log2(length(v))-2), 1.5*diff(range(v))/curr.bw)
  n = if (n[1] < 100) { n[2] } else { n[1] }
  if (n < 2) return(rep(1,length(v)))

  all.mass.den = density(v, bw=curr.bw, n=n)
  lm = localMinima(all.mass.den$y)
  breaks = approx(all.mass.den$x, xout = lm)$y

  cut(v, breaks = c(-Inf, breaks, Inf)) %>% as.numeric %>% rank(ties.method = "min")
  }


breakgroup = function(x, breaksize = 1) {
  o = order(x)
  ass = c(0, which(diff(x[o]) > breaksize), length(x))

  as.integer(cut(seq_along(x), breaks = ass))[order(o)]
}
