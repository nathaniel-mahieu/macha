
plot.wgroup = function (Nmacha, .g_den) {

  .ccs = Nmacha$m.c_cc[Nmacha$m.c[g==.g_den,.(m, c)],cc,on=c("m", "c")]

  ps = Nmacha$m.c_cc[cc %in% .ccs]

  Group = getgroup(Nmacha, .g_den)

  Group$rs[,rt.g := corrt(rt, Nmacha$grt[[m[1]]]),by="m"]

  df.rs = data.frame(Group$rs)
  df.rs$m = factor(df.rs$m)

  xlim = range(ps[,.(rtmin, rtmax)])

  df.ccs = lapply(seq_len(nrow(ps)), function(j) {

    roi = Group$rs[ps[j,.(m, r)],,on=c("m","r")]
    if (nrow(roi) == 0) roi = getroi(Nmacha$m[[ps$m[[j]]]], ps$r[[j]])

    cbind(rt = roi$rt.g, m=ps$m[j], wg=ps$wg[j], cc = ps$cc[j], r = ps$r[j], i = curvemany(unlist(ps[j,.(location, scale, shape, factor, baseline)]), roi$rt.g))

  }) %>% do.call(what=rbind) %>% data.frame

  df.ccs$m = factor(df.ccs$m)

  eachfile=ggplot(df.rs) + geom_line(aes(x = rt.g, y = i, group = m), colour = "grey") + geom_line(data=df.ccs, aes(x = rt, y = i, colour = factor(wg), group = cc)) + theme_nate() + facet_wrap(~m, ncol = round(length(unique(df.rs$m))/3)) + xlim(xlim)+ theme(legend.position = "none")
  overlaycs = ggplot(df.ccs) + geom_line(aes(x = rt, y = i, colour = factor(m), group = factor(paste(cc, m))))  + theme_nate() + xlim(xlim)+ theme(legend.position = "none")
  overlay = ggplot(df.rs) + geom_line(aes(x = rt.g, y = i, colour = m, group = m)) + theme_nate() + xlim(xlim)+ theme(legend.position = "none")
  wrapcs = ggplot(df.ccs[!is.na(df.ccs[,"rt"]),,drop=F]) + geom_line(aes(x = rt, y = i, colour = factor(m), group = m)) + facet_wrap(~cc, ncol = 3, scales = "free_y") + theme_nate() + xlim(xlim)+ theme(legend.position = "none")

  grid.arrange(arrangeGrob(overlay, overlaycs, wrapcs, ncol=1), eachfile, ncol = 2, top = paste0("Group: ",  .g_den,". Warpgroups: ", paste(collapse=", ", unique(ps$wg))))
}
