function() {
  kq = macha$k[,.(q05 = quantile(i, .05), q10 = quantile(i, .1), q30 = quantile(i, .3), q50 = quantile(i, .5), q70 = quantile(i, .7), q90 = quantile(i, .9), q95 = quantile(i, .95)),by="s"]
  ggplot(kq %>% melt(id.vars="s")) + geom_line(aes(x = s, y=value %>% log10, colour= variable)) + ggtitle("Quantiles of detected peaks across the run.")
  
  s. = seq(from = min(macha$s$s), to = max(macha$s$s), by = floor(diff(range(macha$s$s))/11))
  ggplot(macha$k[s %in% s.]) + geom_density(aes(x = i %>% log10)) + facet_wrap(~factor(s)) + ggtitle("Density of log10(intensity) at selected scans") 
  }
