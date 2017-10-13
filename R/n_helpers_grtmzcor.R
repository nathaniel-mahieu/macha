plotgmz = function(Nmacha) {

  df = Nmacha$gmz %>% do.call(what=rbind)
  df[,m:=factor(rep(seq_along(Nmacha$gmz), sapply(Nmacha$gmz, nrow)))]

  df = df %>% as.data.frame
  df = subset(df, !(mz == 0 | mz == max(mz)))
  curves = ggplot(df) + geom_line(aes(x = mz, y = mz - gmz, group = m, colour = m)) + theme_nate() + guides(col = guide_legend(nrow = 1)) +
    geom_hline(data =data.frame(y = c(200/1E6, -200/1E6)), mapping = aes(yintercept=y), colour="red", alpha=0.5) + geom_hline(data =data.frame(y = c(200/1E6 * 0.5, -200/1E6 * 0.5)), mapping = aes(yintercept=y), colour="yellow", alpha=0.5)


  Gcs = Nmacha$m.c[Nmacha$gcs[,.(c,m)],,on=.(c,m)]
  Gcs[,mz.g:=cormz(mz, Nmacha$gmz[[m[1]]]),by="m"]


  tgs = sapply(split(Gcs$g, cut(Gcs$rtpeak, 4)), sample, size=1)

  preps = lapply(tgs, function(.g) {

    Group = getgroup(Nmacha, .g)
    Group$cs[,mz.g := cormz(mz, Nmacha$gmz[[m[1]]]), by = "m"]

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


explore_grouping = function(Nmacha, scans = seq(1, 20, 4), ppms = seq(1, 10, 2)) { # Unfinished function to help choose grouping parameters.

  grid = expand.grid(scans, ppms)
  dt = Nmacha$m.c[,.(k = m.c., mz = mz, s = location)]

  dt.g = apply(grid, 1, function(x) {
    cbind(g = rectroi(dt, scan = x[1], ppm = x[2])[,.(m.c. = k, g = r)]$g %>% table, ppm = x[1], scan = x[2])
  }) %>% do.call(what=rbind) %>% data.table

  dt.g.h = dt.g[,.N, by="g,ppm,scan"]
  grid.arrange(
    ggplot(dt.g) + geom_histogram(aes(x = g), stat="count") + facet_grid(ppm ~ scan) + xlim(0, length(Nmacha$m)*2),
    ggplot(dt.g) + geom_histogram(aes(x = g), stat="count") + facet_grid(ppm ~ scan) + scale_y_log10(),
    ncol = 2
  )
}
