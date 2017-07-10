library(macha)
macha = readRDS("Q:/Datashare/05172017_cphB/data_positive/all_data/05172017_ExpB_C8_B18_21.mzXML.macha.rds")

foo = copy(macha$k)
foo[,g:=1]

###


setkey(foo, g, mz)
x = foo$mz
ass = c(0, which(diff(x)/x[-1] * 1E6 > 8), length(x))
inds = as.integer(cut(seq_along(x), breaks = ass))

foo[,g:= as.numeric(factor(paste(g, inds)))]


setkey(foo, g, s)
x = foo$s
ass = c(0, which(diff(x)/2 > 4), length(x))
inds = as.integer(cut(seq_along(x), breaks = ass))

foo[,g:= as.numeric(factor(paste(g, inds)))]

length(unique(foo$g))


##

ppmrange = foo[,.(n = .N, ppm = diff(range(mz))/mean(mz)*1E6, s = diff(range(s)), dups = sum(duplicated(s)), dupfrac = 2*sum(duplicated(s)) / length(s)),by="g"]

hist(ppmrange$ppm, breaks = 200)
hist(ppmrange$s, breaks = 200)
hist(ppmrange$dupfrac %>% log10)

plot(ppmrange[n > 5, .(ppm, s)])

library(ggplot2)
ggplot(ppmrange[n>5]) + geom_point(aes(x = s, y = ppm, colour = dups))
ggplot(ppmrange[n>5]) + geom_point(aes(x = s, y = ppm, colour = dupfrac)) + facet_wrap(~cut(dupfrac,4))
ggplot(ppmrange[n>5 & dupfrac < 0.1]) + geom_point(aes(x = s, y = ppm, colour = dups)) + facet_wrap(~cut(dupfrac,4))
ggplot(ppmrange[n > 2]) + geom_point(aes(x = n, y = dupfrac))
ppmrange[dupfrac > 0]$dupfrac %>% hist(breaks = 50)

ggplot(foo[g==ppmrange[dupfrac > .30 & n > 30]$g %>% sample(1)]) + geom_point(aes(y = mz, x = s, colour = duplicated(s))) + geom_text(aes(x = median(s), y = median(mz), label = g), check_overlap = T)

#534857
#10537
#17912

bar = copy(foo[g == 17912])
bar
ggplot(bar) + geom_point(aes(y = mz, x = s, colour = duplicated(s))) + geom_text(aes(x = median(s), y = median(mz), label = g), check_overlap = T)

setkey(bar, s, mz)
subbar = copy(bar[abs(s-3201) < 10])


#Group to peak densities?
den = subbar$mz %>% density(adjust = 0.2)
den %>% plot
ps = den$x[which(diff(sign(diff(den$y)))==-2)+1]

diffs = abs(outer(subbar$mz, ps, "-"))
mats = which(rowMins(diffs) == diffs, arr.ind=T)

if (length(na.omit(mats[,1])) != length(subbar$mz) | anyDuplicated(mats[,1])) { stop("Clustering edge case.") }

cbind(subbar$k, mats[order(mats[,1]), 2])


#Group to each peak
matches = lapply(unique(bar$s), function(s.i) {
  subbar = copy(bar[abs(s-s.i) <= 10*2])
  
  nrow(subbar) %>% print
  
  d.mz = abs(outer(subbar$mz, subbar$mz, "-"))
  d.s = abs(outer(subbar$s, subbar$s, "-")) > 0
  d.centeronly = abs(outer(rep(1, length(subbar$s)), subbar$s, "*") - s.i) <= 8*2

  d.mz[!d.s] = Inf
  diag(d.mz) = Inf
  
  n = 3
  matches = lapply(seq_len(n), function(x) {
    d.tf = d.mz == rowMins(d.mz) & d.mz < Inf & d.centeronly
    mats = which(d.tf, arr.ind = T)
    d.mz[mats] <<- Inf
    
    mats[] = as.character(subbar$k[c(mats)])
    mats
    }) %>% do.call(what = rbind)
  }) %>% do.call(what = rbind)


colnames(matches) = c("v1", "v2")
matchdt = data.table(matches)

matchdt[,v1.s := bar[match(v1,as.numeric(k)), s]]
matchdt[,v1.mz := bar[match(v1,as.numeric(k)), mz]]
matchdt[,v2.s := bar[match(v2,as.numeric(k)), s]]
matchdt[,v2.mz := bar[match(v2,as.numeric(k)), mz]]
matchdt[,weight := abs((v1.mz - v2.mz) * 1E3)]
matchdt[,weight := max(weight)-weight + 1]

#With ppm weighted trimming
matchdt[,v1.mz - v2.mz] %>% abs %>% hist (breaks = 500)
matchdt[,v1.mz - v2.mz] %>% abs %>% quantile(.98) %>% {abline(v=.)}

matchdt = matchdt[matchdt$weight > quantile(matchdt$weight, .02)]


#Try Min Cut Approaches
library(igraph)
g = graph.data.frame(matchdt[,.(as.character(v1), as.character(v2), weight)], directed = F, vertices = as.character(bar$k))


o.s = outer(bar$s, bar$s, "==")
o.s[upper.tri(o.s, T)] = F
storesolv = which(o.s, arr.ind = T)
storesolv[] = bar$k[storesolv] %>% as.character


mincuts = lapply(seq_len(nrow(storesolv)), function(row) {
  max_flow(g, storesolv[row,1], storesolv[row,2])
  })

mincutt = lapply(mincuts, function(x) {
  x$cut %>% sapply(as.numeric)
  }) %>% unlist %>% table

mincutt[order(mincutt, decreasing = T)]

g = delete.edges(g, head(as.numeric(names(mincutt[order(mincutt, decreasing = T)])), n = 3))





# Visualize
graph =g
bar[,g2 := membership(components(graph))[match(k, as.numeric(names(membership(components(graph)))))]]

ggplot(bar) + geom_point(aes(y = mz, x = s, colour = g2 %>% factor)) + geom_text(aes(x = median(s), y = median(mz), label = g), check_overlap = T)

plot(g, vertex.label = NA, vertex.size = 3, edge.arrow.size = 0)

ggplot(bar) + geom_point(aes(y = mz, x = s, colour = g2 %>% factor)) + geom_text(aes(x = median(s), y = median(mz), label = g), check_overlap = T) +
  geom_segment(data = matchdt, aes(x = v1.s, xend = v2.s, y = v1.mz, yend = v2.mz), position = "jitter")

ggplot(bar) + geom_point(aes(y = mz, x = s, colour = g2 %>% factor)) + geom_text(aes(x = median(s), y = median(mz), label = g), check_overlap = T) +
  geom_segment(data = matchdt, aes(x = v1.s, xend = v2.s, y = v1.mz, yend = v2.mz), position = "jitter") + xlim(3250, 3300)
