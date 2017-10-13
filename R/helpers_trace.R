make_trace_cache = function(Nmacha) {
  trace_cache = lapply(seq_along(Nmacha$m), function(j) {
    macha = Nmacha$m[[j]]

    traces = macha$k[macha$k_r,,on="k",nomatch=0][macha$r_mchan,,on="r",nomatch=0][macha$k_b,,on="k", nomatch=0][,.(k,s,i,r,mz,b,mchan)]
    setkey(traces, "mchan","s")

    mrng = traces[,{tmp = range(s); .(mins = min(s), maxs = max(s))},by="mchan"]
    mchanscan = mrng[,.(s = seq(mins, maxs)),by="mchan"]
    setkey(mchanscan, mchan, s)

    traces = traces[mchanscan,,roll=T, rollends=F, on=c("mchan","s"), nomatch=0]
    traces[duplicated(k), ':='(k=NA, mz=NA,i=NA)]
    traces[,m:=j]
    traces = traces[macha$s,,on="s"]


    setkey(traces, r, m, mchan)
    traces
  })

  trace_cache
  }
