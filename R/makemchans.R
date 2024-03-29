
#' Reconstruct ROIs split due to missing observations
#'
#' \code{makemchans} Find's ROIs which are likely the same mass channel.
#'
#' Combining ROIs which were previously split due to missing values or conflicts in the webtrace step allows more complete information to be supplied to the baselining step.
#'
#' @param macha macha. Three columns: mz, s, k and (optional) g. s is scan, k is a unique integer ID for each peak, g is the initial group assignment.
#' @param ppmwin numeric. Mass distance in ppm to disqualify ROI aggregation. Recommended: instrument_ppm
#' @param rtwin numeric. RT distance in seconds to disqualify ROI aggregation. Recommended: 1/5 of run length
#'
#' @return macha. macha object with additional item $r_mchan
#'
#'
#' @export
#'
makemchans = function(macha, ppmwin = 3, rtwin = 30) {
  cat("\nAggregating mass channels and reconstructing traces.\n")
  
  mr = macha$k[macha$s,,on="s"][macha$k_r,,on="k"][, .(minrt = min(rt), maxrt = max(rt), meanmz = mean(sum(mz*i)/sum(i)), maxmz = max(mz), minmz = min(mz)),by=r] 
  setkey(mr, "meanmz")
  mr = mr %>% as.matrix
  mr.nrow = nrow(mr)
  massspans = mr[,"meanmz"] %>% { c(0,diff(.))/.*1E6 } %>% cumsum
  
  mchans = mr[,"r"]
  rs = mr[,"r"]
  maxmchan = max(mchans)
  laststartmass = lastmchan = i = 0
  startmass = 1
  
  #profvis::profvis({
  while(startmass <= length(massspans) + 1) {
    
    cons = { massspans - massspans[startmass] } %>% { . >= 0 & . < ppmwin } %>% which
    mchans.tf = mchans %in% mchans[cons]
    ts = mr[mchans.tf,,drop=F]
    if (nrow(ts) < 2) { startmass = startmass + 1; next }
    
    i = i+1
    #if (i == 51039) break;
    cat("\rFraction Complete:", round(max(cons)/mr.nrow, 3), i, "         ")
    
    o1 = outer(ts[,"maxrt"], ts[,"minrt"], "-")
    o2 = outer(ts[,"minrt"], ts[,"maxrt"], "-")
    overlap = ( sign(o1) != sign(o2) ) #& ( o1 != 0 & o2 != 0 )
    dist = outer(ts[,"maxrt"], ts[,"minrt"], "-")
    
    samemchan = outer(mchans[mchans.tf], mchans[mchans.tf], "==") %>% { .[upper.tri(., T)] = NA; . } %>% which(arr.ind=T)
    for (r in seq_len(nrow(samemchan))) {
      overlap[samemchan[r,1],] = overlap[samemchan[r,1],] | overlap[,samemchan[r,2]]
      overlap[samemchan[r,2],] = overlap[samemchan[r,2],] | overlap[,samemchan[r,1]]
      overlap[,samemchan[r,1]] = overlap[,samemchan[r,1]] | overlap[samemchan[r,2],]
      overlap[,samemchan[r,2]] = overlap[,samemchan[r,2]] | overlap[samemchan[r,1],]
    }
    overlap = overlap %>% { .[upper.tri(., T)] = NA; . }
    
    pmerg = which(!overlap & (dist < rtwin | aperm(dist) < rtwin), arr.ind=T)
    
    dist.v = cbind(dist[pmerg], -dist[pmerg[,c(2,1),drop=F]]) %>% abs %>% rowMins
    
    #update = mchans %in% ts[pmerg[which.min(dist),],"mchan"]
    update = mchans[mchans.tf] %in% mchans[mchans.tf][pmerg[which.min(dist.v),]]
    #ts[update,]
    
    if (!any(update)) { 
      startmass = startmass + 1 
    } else {
      maxmchan = maxmchan + 1
      mchans[mchans.tf][update] = rep(maxmchan, sum(update))
      }
  }
  #})
  macha$r_mchan = data.table(r=rs, mchan = mchans)
  macha
}
