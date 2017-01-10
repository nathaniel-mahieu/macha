
sspp = function(v, S = 1:20, maxdist = 5, do.plot = T) {

  dC = C = vector(mode = "numeric", length=length(v))
  O = NULL
  N = length(v)
  v.original = v.prev = v

  #dC.m Might ofer additional diagnostic information.  Could exclude peaks that moved around.  Probably just unnecesarily complex
  #dC.m = matrix(ncol = length(dC), nrow = length(S))

  for (i in seq_along(S)) {
    s = S[i]
    #s = round(max(c(1,i*N / (10 * S))))

    W = signal::hamming(s)
    W = W/sum(W)

    v.prev = v

    v = as.numeric(stats::filter(v,W))
    v[is.na(v)] = v.prev[is.na(v)]

    P = localMaxima(v)

    if (length(P) == 0) break

    if (i == 1) {
      dC[P] = v[P]
      O = P
    } else {
      dC[] = 0

      dist = abs(outer(O, P, "-")) < maxdist

      ints = matrix(v.original[O], ncol = ncol(dist), nrow = nrow(dist))
      ints[!dist] = NA
      maxints = matrix(matrixStats::colMaxs(ints, na.rm=T), ncol = ncol(ints), nrow = nrow(ints), byrow = T)

      neighbors = which(ints == maxints, arr.ind =T)

      dC[O[neighbors[,1]]] = v[P[neighbors[,2]]] * s^1.4
    }
    C = C + dC
    #dC.m[i,] = dC
    #O = union(O, P)
  }

  if (do.plot) {
    scores = cbind(which(C>1), C[C>1])
    scores[,2] = scores[,2] * max(v.original)/max(scores[,2])

    par(mfcol=c(1,1))
    plot(v.original, type="l")
    lines(scores, type="h", col = "red")
    lines(v, type="l", col="grey")
    points(y=v.original[O], x=O, col="red")
    }

  C
}
