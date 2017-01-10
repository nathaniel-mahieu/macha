
getgroup = function(Nmacha, .g) {
  gl = Nmacha$m.c[g==.g]

  ugl = split(gl, gl[,paste(m,r)])



  #if (F) { # Group into single DT?
  group = list(cs = gl)

  group$rs = lapply(split(gl, gl[,paste(m,r)]), function (g) {
    getroi(Nmacha$m[[g$m[1]]], g$r[1])[,m:=g$m[1]]
  }) %>% do.call(what=rbind)
  #}

  group$r = Nmacha$m.r[group$cs[,.(m,r)],,on=c("m", "r"), nomatch=0]

  group
}


plot.group = function(Nmacha, .g) {

  group = getgroup(Nmacha, .g)

  library(ggplot2)
  library(gridExtra)


  df.rs = data.frame(group$rs)
  df.cs = data.frame(group$cs)
  df.ccs = lapply(seq_len(nrow(group$cs)), function(j) {

    roi = group$rs[group$cs[j,.(m, r)],,on=c("m","r")]

    cbind(rt = roi$rt, m=group$cs$m[j], c = group$cs$c[j], r = group$cs$r[j], i = curvemany(unlist(group$cs[j,.(location, scale, shape, factor, baseline)]), roi$rt))

  }) %>% do.call(what=rbind) %>% data.frame

  df.ccs$uid = paste(df.ccs$m, df.ccs$c)
  df.cs$uid = paste(df.cs$m, df.cs$c)

  df.cs$rtmean = sapply(seq_len(nrow(df.cs)), function(j) getmean.sn(df.cs[j,"location"], df.cs[j,"scale"], df.cs[j,"shape"]))

  xlim = range(group$cs[,.(rtmin, rtmax)]); xlim[1] = xlim[1] - 3; xlim[2] = xlim[2]+3

  all.o = ggplot() +
    geom_line(data = df.rs, aes(x=rt, y = i, group = m, colour = (m))) +
    geom_line(data=df.ccs, aes(x = rt, y = i, group = uid), colour="red") +
    theme_nate()  + theme(legend.position = "none") + scale_y_continuous(labels=fancy_scientific) + xlim(xlim)


  cs.o = ggplot() +
    geom_line(data=df.ccs, aes(x = rt, y = i, group = uid, colour = m)) +
    theme_nate()  + theme(legend.position = "none") + scale_y_continuous(labels=fancy_scientific) + xlim(xlim)
  is.o =   ggplot() +
    geom_line(data = df.rs, aes(x=rt, y = i, group = m, colour = (m))) +
    theme_nate()  + theme(legend.position = "none") + scale_y_continuous(labels=fancy_scientific) + xlim(xlim)

  peakpos = ggplot(df.cs) + geom_vline(aes(xintercept = rtpeak, colour=m), size = 1) + theme_nate() + theme(legend.position = "none") + scale_y_continuous(labels=fancy_scientific)


  all.w = ggplot() +
    geom_line(data = df.rs, aes(x=rt, y = i, group = m, colour = (m))) +
    geom_line(data=df.ccs, aes(x = rt, y = i, group = uid), colour="red") +
    theme_nate() + facet_wrap(~m, nrow = round(length(unique(df.cs$m))/3)) + theme(legend.position = "none") + scale_y_continuous(labels=fancy_scientific) + xlim(xlim)

  grid.arrange(arrangeGrob(all.o, cs.o, is.o, peakpos, ncol=1), all.w, ncol=2, nrow =1, top = paste0("Group: ", .g), layout_matrix=cbind(c(1), 2,2))

}
