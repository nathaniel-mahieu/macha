plotgrt = function(Nmacha) {

  df = Nmacha$grt %>% do.call(what=rbind)
  df[,m:=factor(rep(seq_along(Nmacha$grt), sapply(Nmacha$grt, nrow)))]

  df = df %>% as.data.frame
  df = subset(df, !(rt == 0 | rt == max(rt)))

  curves = ggplot(df) + geom_line(aes(x = rt, y = rt - grt, group = m, colour = m)) + theme_nate() + guides(col = guide_legend(nrow = 1))



  Gcs = Nmacha$m.c[Nmacha$gcs[,.(c,m)],,on=.(c,m)]
  tgs = sapply(split(Gcs$g, cut(Gcs$rtpeak, 4)), sample, size=1)

  preps = lapply(tgs, function(.g) {

    Group = getgroup(Nmacha, .g)
    Group$rs[,rt.g := corrt(rt, Nmacha$grt[[m[1]]]),by="m"]

    xlim = range(Group$cs[,.(rtmin, rtmax)])

    df.rois = Group$rs %>% as.data.frame

    df.rois$m = factor(df.rois$m)

    list(
      ggplot() +
        geom_line(data = df.rois, aes(x=rt, y = i, group = m, colour = (m))) +
        theme_nate()  + theme(legend.position = "none") + scale_y_continuous(labels=fancy_scientific) + xlim(xlim),
      ggplot() +
        geom_line(data = df.rois, aes(x=rt.g, y = i, group = m, colour = (m))) +
        theme_nate()  + theme(legend.position = "none") + scale_y_continuous(labels=fancy_scientific) + xlim(xlim)
    )

    })


  g = do.call(arrangeGrob, c(sapply(preps, '[', 1), list(nrow = 1)))
  g2 = do.call(arrangeGrob, c(sapply(preps, '[', 2), list(nrow = 1)))
  grid.arrange(curves, g, g2, nrow = 3, top="RT correction curves and randomly selected peak groups.", layout_matrix = cbind(c(1,1,2,3)))

}

plotgmz = function(Nmacha) {

  df = Nmacha$gmz %>% do.call(what=rbind)
  df[,m:=factor(rep(seq_along(Nmacha$gmz), sapply(Nmacha$gmz, nrow)))]

  df = df %>% as.data.frame
  df = subset(df, !(mz == 0 | mz == max(mz)))
  curves = ggplot(df) + geom_line(aes(x = mz, y = mz - gmz, group = m, colour = m)) + theme_nate() + guides(col = guide_legend(nrow = 1)) +
    geom_hline(data =data.frame(y = c(200/1E6, -200/1E6)), mapping = aes(yintercept=y), colour="red", alpha=0.5) + geom_hline(data =data.frame(y = c(200/1E6 * 0.5, -200/1E6 * 0.5)), mapping = aes(yintercept=y), colour="yellow", alpha=0.5)


  Gcs = Nmacha$m.c[Nmacha$gcs[,.(c,m)],,on=.(c,m)]
  Gcs[,mz.g:=cormz(mz, Nmacha$gmz[[m[1]]],by="m")]


  tgs = sapply(split(Gcs$g, cut(Gcs$rtpeak, 4)), sample, size=1)

  preps = lapply(tgs, function(.g) {

    Group = getgroup(Nmacha, .g)
    Group$cs[,mz.g := cormz(mz, Nmacha$gmz[[m[1]]], by = "m")]

    xlim = range(Group$cs[,.(mz, mz.g)])

    df.mzs = Group$cs %>% as.data.frame
    df.mzs$m = factor(df.mzs$m)

    list(
      ggplot() +
        geom_vline(data = df.mzs, aes(xintercept=mz, colour = (m))) +
        theme_nate()  + theme(legend.position = "none") + scale_y_continuous(labels=fancy_scientific) + xlim(xlim),
      ggplot() +
        geom_vline(data = df.mzs, aes(xintercept=mz.g, colour = (m))) +
        theme_nate()  + theme(legend.position = "none") + scale_y_continuous(labels=fancy_scientific) + xlim(xlim)
    )

  })

  g = do.call(arrangeGrob, c(sapply(preps, '[', 1), list(nrow = 1)))
  g2 = do.call(arrangeGrob, c(sapply(preps, '[', 2), list(nrow = 1)))
  grid.arrange(curves, g, g2, nrow = 3, top="Mz correction curves and randomly selected peak groups.", layout_matrix = cbind(c(1,1,2,3)))

}
